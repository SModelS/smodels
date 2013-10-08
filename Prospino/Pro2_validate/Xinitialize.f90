! ===========================================================================================================
module xx_initialize
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  use xx_in_out
  implicit none 
  private :: TEST_SETTINGS_GLOBAL, TEST_SETTINGS_ITERATION
  public  :: PROSPINO_CHECK_FS, PROSPINO_CHECK_HIGGS
  public  :: INIT_GLOBAL, INIT_ITERATION, INIT_VEGAS, INIT_IFAST, INIT_SUSY
  public  :: INIT_ALPHAS, FILL_COMMONS_SQ_GL
contains
! ------------------------------
  subroutine PROSPINO_CHECK_HIGGS(final_state_in)

    character(len=2), intent(in) :: final_state_in

    if ( (final_state_in=='hh').or.(final_state_in=='ht') ) then
       print*, " private code, not part of the official Prospino2 distribution "
       print*, "  "
       print*, " please contact tilman.plehn@cern.ch "
       print*, "  "
       call HARD_STOP
    end if

  end subroutine PROSPINO_CHECK_HIGGS

! ------------------------------
  subroutine PROSPINO_CHECK_FS(final_state_in,ipart1_in,ipart2_in,lfinal)

    character(len=2), intent(in) :: final_state_in
    integer, intent(in)          :: ipart1_in,ipart2_in
    logical, intent(out)         :: lfinal

    lfinal = .true.
    if ( (final_state_in/='nn').and.(ipart2_in/=1) ) lfinal = .false. ! only chi-chi needs both numbers
    if ( (final_state_in/='ll').and.(ipart1_in>8)  ) lfinal = .false. ! only sleptons have numbers>8
    if ( (final_state_in=='gg').and.(ipart1_in>1)  ) lfinal = .false. ! squark and gluino do not care about ipart
    if ( (final_state_in=='sb').and.(ipart1_in>1)  ) lfinal = .false.
    if ( (final_state_in=='sg').and.(ipart1_in>1)  ) lfinal = .false.
    if ( (final_state_in=='ss').and.(ipart1_in>1)  ) lfinal = .false.
    if ( (final_state_in=='tb').and.(ipart1_in>2)  ) lfinal = .false.
    if ( (final_state_in=='bb').and.(ipart1_in>2)  ) lfinal = .false.
    if ( (final_state_in=='lq').and.(ipart1_in>1)  ) lfinal = .false.
    if ( (final_state_in=='le').and.(ipart1_in>1)  ) lfinal = .false.
    if ( (final_state_in=='sb').and.(ipart1_in>1)  ) lfinal = .false.
    if ( (final_state_in=='xx').and.(ipart1_in>2)  ) lfinal = .false.

    if ( .not. lfinal ) print*, " PROSPINO_CHECK_FS: final state not valid "
       
  end subroutine PROSPINO_CHECK_FS

! ------------------------------
  subroutine INIT_GLOBAL(ipart1_in,ipart2_in,lvalid_global)

    integer, intent(in)  :: ipart1_in,ipart2_in
    logical, intent(out) :: lvalid_global

    ii_done(-1:30) = 0                                  ! initialize array as nothing done yet
    ii_done(-2)    = 1
    dim(-1:30)     = 1                                  ! dimension of integration 

    select case (final_state) 
    case ('ng')
       dim(-1:5)  = (/ 3,3,3,3,4,4,5 /)
    case ('ns')
       dim(-1:7)  = (/ 3,3,3,3,4,4,5,4,5 /)
    case ('nn')
       dim(-1:11) = (/ 3,3,3,3,6,4,4,6,6,6,6,6,6 /)
    case ('ll')
       dim(-1:6)  = (/ 3,3,3,3,6,4,6,4 /)
    case ('tb','bb')
       dim(-1:4)  = (/ 3,3,3,3,4,4 /)
    case ('lq')
       dim(-1:4)  = (/ 3,3,3,3,4,4 /)
    case ('le')
       dim(-1:10) = (/ 3,3,3,3,4,4,5,8,4,4,5,5 /)
    case ('hh')
       dim(-1:8)  = (/ 3,3,3,3,3,4,7,5,4,4 /)
    case ('ht')
       dim(-1:9)  = (/ 3,3,3,3,3,4,7,7,4,5,5 /)
    case ('gg','sb','sg','ss','xx')                     ! dummy setting
       dim(-1:11) = (/ 1,1,1,1,1,1,1,1,1,1,1,1,1 /)
    case default
       print*, " INIT_GLOBAL: final_state not set ",final_state
    end select

!----------------------------------------------------------------------------
    iscaling = 0  ! calculate cross sections[0] or scaling functions[1]     !
    eta = 1.e-3                                                             !
!----------------------------------------------------------------------------
    if (iscaling==1) icoll = 4   ! flag for the luminosities
    
!----------------------------------------------------------------------------
    imx = 0       ! negative masses[0] hard wired                           !
!----------------------------------------------------------------------------
    
!----------------------------------------------------------------------------
    isca = 0      ! scale of process average mass[0] , cm energy[1]         !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
    cut     = 1.e-8   ! cut on integration variables [<1.e-8]               !
    eps_sli = 1.e-5   ! cut for phase space slicing [1.e-5 -> 0.0]          !
    eps_sub = 1.e-4   ! for for phase space subtraction [1.e-4 -> 1.e-6]    !
    ewi     = 1.e-2   ! width of intermediate squark [1.e-2 -> 1.e-5]       !
!----------------------------------------------------------------------------
    if (ewi<1.e-3) ewi = 1.e-3            ! minimum value for substitution tan(z) 

    if (final_state=='ng') then 
       if (ipart1_in>=5) then
          idub = 1                        ! set idub for (q,qbar'), needed for ng
       else if (ipart1_in<=4) then 
          idub = 0                        ! set idub for (q,qbar), needed for ng
       end if
    else if (final_state=='ns') then 
       idub = 0
    else if (final_state=='nn') then 
       if ( ((ipart1_in<5).and.(ipart2_in<5)).or.(ipart1_in>4).and.(ipart2_in>4)) then
          idub = 0                        ! set idub=1 for (q,qbar')
       else 
          idub = 1 
       end if
    else if (final_state=='ll') then 
       if ((ipart1_in==4).or.(ipart1_in==5).or.((ipart1_in>=10).and.(ipart1_in<=13)) ) then
          idub = 1
       else 
          idub = 0
       end if
    else if (final_state=='tb') then 
       idub = 0
    else if (final_state=='bb') then 
       idub = 0
    else if (final_state=='lq') then 
       idub = 0
    else if (final_state=='le') then 
       idub = 0
    end if

    if (icoll>3) then                     ! may happen for scaling functions 
       sc = 1.0
    end if
    
    call TEST_SETTINGS_GLOBAL(ipart1_in,ipart2_in,lvalid_global)  ! check everything set in main program
    
  end subroutine INIT_GLOBAL
  
! ------------------------------
! note that routine is only internal to this module, called in INIT_GLOBAL
  subroutine TEST_SETTINGS_GLOBAL(ipart1_in,ipart2_in,lvalid_global)
    
    integer, intent(in)  :: ipart1_in,ipart2_in
    logical, intent(out) :: lvalid_global
    
    lvalid_global = .true.      ! go through unless something happens
    
    if ((final_state/='ng').and.(final_state/='ns').and.                                 & ! test final_state
        (final_state/='nn').and.(final_state/='ll').and.                                 &
        (final_state/='gg').and.(final_state/='sb').and.                                 &
        (final_state/='sg').and.(final_state/='ss').and.                                 & 
        (final_state/='tb').and.(final_state/='bb').and.                                 &
        (final_state/='lq').and.(final_state/='le').and.                                 &
        (final_state/='hh').and.(final_state/='ht').and.                                 &
        (final_state/='xx')                                                       ) then 
       print*, " TEST_SETTINGS_GLOBAL: final_state not correct ",final_state
       lvalid_global = .false.
    end if
       
!tp    if (( ipart1_in < 1).or.(ipart2_in < 1) ) then                                         ! start with a general check
!tp       print*, " TEST_SETTINGS_GLOBAL: ipart1_in or ipart2_in not in range >0 ",ipart1_in,ipart2_in
!tp       lvalid_global = .false.
!tp    end if
       
    if ( ((final_state=='nn').or.(final_state=='ns').or.(final_state=='ng')).and.        & ! test ipart1_in
         ((ipart1_in<1).or.(ipart1_in>8))                                              ) then 
       print*, " TEST_SETTINGS_GLOBAL: ipart1_in not in range 1...8 ",ipart1_in
       lvalid_global = .false.
    else if ( (final_state=='ll').and.                                                   &
         ((ipart1_in<0).or.(ipart1_in>14))                                             ) then 
       print*, " TEST_SETTINGS_GLOBAL: ipart1_in not in range 1...8 ",ipart1_in
       lvalid_global = .false.
    end if
       
    if (  (final_state=='nn').and.                                                       & ! test ipart2_in
         ((ipart2_in<1).or.(ipart2_in>8))                                              ) then 
       print*, " TEST_SETTINGS_GLOBAL: ipart2_in not in range 1...8 ",ipart2_in
       lvalid_global = .false.
    end if
       
    if ((final_state=='nn').and.                           &                               ! charge if final state cannot be 2
        ( ((ipart1_in==5).and.(ipart2_in==5)).or.          &
          ((ipart1_in==6).and.(ipart2_in==6)).or.          &
          ((ipart1_in==7).and.(ipart2_in==7)).or.          &
          ((ipart1_in==8).and.(ipart2_in==8)).or.          &
          ((ipart1_in==5).and.(ipart2_in==6)).or.          &
          ((ipart1_in==6).and.(ipart2_in==5)).or.          &
          ((ipart1_in==7).and.(ipart2_in==8)).or.          &
          ((ipart1_in==8).and.(ipart2_in==7))     ) ) then 
       print*, " TEST_SETTINGS_GLOBAL: final state not possible ",final_state,ipart1_in,ipart2_in
       print*, " TEST_SETTINGS_GLOBAL: final state not possible ",final_state,ipart1_in,ipart2_in
       lvalid_global = .false.
    end if

    if ( (icoll<0).or.(icoll>3) ) then                                                     ! test icoll 
       print*, " TEST_SETTINGS_GLOBAL: icoll not in range 0...1 ",icoll 
       lvalid_global = .false.
    end if
       
  end subroutine TEST_SETTINGS_GLOBAL

! ------------------------------
  subroutine INIT_ITERATION(lvalid_iteration)

    logical, intent(out) :: lvalid_iteration

    n_faulty = 0                                                                           ! bad events counter initialized 

    lvalid_iteration = .true.                                                              ! go through unless something happens

    if (ii<=0) then                                                                        ! initialize the parton densities, either LO or NLO
       call INIT_PDF(0)
    else 
       call INIT_PDF(1)
    end if

    call TEST_SETTINGS_ITERATION(lvalid_iteration)                                         ! check everything set in main program, for now only ii

  end subroutine INIT_ITERATION

! ------------------------------
! note that routine is only internal to this module, called in INIT_GLOBAL
  subroutine TEST_SETTINGS_ITERATION(lvalid_iteration)
    
    logical, intent(out) :: lvalid_iteration
    
    lvalid_iteration = .true.      ! go through unless something happens
    
    if (ii<-1) then                                                                         ! lower limit to ii universal 
       print*, " TEST_SETTINGS_ITERATION: ii not correct ",ii
       lvalid_iteration = .false.
    end if

    select case (final_state)
    case('ng')
       if (ii>5) then 
          lvalid_iteration = .false.
       end if
    case('ns')
       if (ii>1) then 
          print*, " Sorry, the NLO corrections to this process are not yet checked and implemented "
          lvalid_iteration = .false.
       end if
    case('nn')
       if (ii>11) then 
          lvalid_iteration = .false.
       end if
    case('ll')
       if ((ii==3).or.(ii==5).or.(ii>6)) then 
          lvalid_iteration = .false.
       end if
    case('tb','bb')
       if (ii>4) then 
          lvalid_iteration = .false.
       end if
    case('lq')
       if (ii>4) then 
          lvalid_iteration = .false.
       end if
    case('le')
       if (ii>10) then 
          lvalid_iteration = .false.
       end if
    case('hh')
       if (ii>8) then 
          lvalid_iteration = .false.
       end if
    case('ht')
       if (ii>9) then 
          lvalid_iteration = .false.
       end if
    end select
    
  end subroutine TEST_SETTINGS_ITERATION
       
! ------------------------------
  subroutine INIT_VEGAS(i_lofast,ivegas)
    
    integer,              intent(in)  :: i_lofast
    integer,dimension(4), intent(out) :: ivegas

    select case (final_state)                       ! the different ivegas input values to get <1% error for nlo contributions 
    case('ng')

       select case (ii)
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/ 5000  , 4 , 30000  , 4 /)
       case(2)                                      ! virtual
          ivegas(1:4) = (/ 2000  , 4 , 10000  , 4 /)
       case(3)                                      ! real emission
          ivegas(1:4) = (/ 10000 , 8 , 100000 , 4 /)
       case(4)                                      ! crossed matrix elements
          ivegas(1:4) = (/ 10000 , 8 , 200000 , 4 /)
       case(5)                                      ! crossed os subtraction 
          ivegas(1:4) = (/ 10000 , 8 , 200000 , 6 /)
       case default       
          ivegas(1:4) = (/ 100   , 4 ,  1000  , 2 /)
       end select

    case('ns')

       select case (ii)
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/ 5000   , 4 ,  30000  , 4 /)
       case(2)                                      ! virtual
          ivegas(1:4) = (/ 100    , 4 ,  10000  , 4 /)
       case(3)                                      ! real qg matrix element
          ivegas(1:4) = (/ 10000  , 8 , 100000  , 4 /)
       case(4,6)                                    ! crossed matrix elements
          ivegas(1:4) = (/ 10000  , 8 , 200000  , 4 /)
       case(5,7)                                    ! crossed os subtractions
          ivegas(1:4) = (/ 20000  , 8 , 200000  , 6 /)
       case default       
          ivegas(1:4) = (/ 100    , 4 ,  1000   , 2 /)
       end select

    case('nn')

       select case (ii)
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/ 5000  , 4 , 30000  , 4 /)
       case(2)                                      ! virtual
          ivegas(1:4) = (/ 500   , 4 , 10000  , 6 /)
       case(3)                                      ! real matrix element
          ivegas(1:4) = (/ 20000 , 4 , 100000 , 6 /)
       case(4)                                      ! real dipole subtraction 
          ivegas(1:4) = (/ 10000 , 4 , 200000 , 6 /)
       case(5)                                      ! crossed dipole subtraction 
          ivegas(1:4) = (/ 10000 , 4 , 100000 , 6 /)
       case(6,9)                                    ! crossed matrix element
          if ( ( (ipart1<5).and.((ipart2==6).or.(ipart2==8)) ).or.       & ! N-C2 all fucked up because of dominant o-s subtraction 
               ( (ipart2<5).and.((ipart1==6).or.(ipart1==8)) )    ) then 
             ivegas(1:4) = (/ 80000 , 4 ,1500000 ,10 /)
          else 
             ivegas(1:4) = (/ 80000 , 4 , 500000 , 6 /)
          end if
       case(7,8,10,11)                              ! crossed o-s subtraction 
          ivegas(1:4) = (/ 80000 , 4 , 500000 , 6 /)
       case default       
          ivegas(1:4) = (/ 100   , 4 ,  1000  , 2 /)
       end select
       
    case('ll')

       select case (ii)
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/ 5000  , 4 , 30000  , 4 /)
       case(2)                                      ! virtual
          ivegas(1:4) = (/ 500   , 4 , 10000  , 6 /)
       case(3,5,7)                                  ! real emission matrix element
          ivegas(1:4) = (/ 20000 , 4 , 100000 , 6 /)
       case(4,6,8)                                  ! real dipole subtraction 
          ivegas(1:4) = (/ 10000 , 4 , 300000 , 6 /)
       case default       
          ivegas(1:4) = (/ 100   , 4 ,  1000  , 2 /)
       end select
       
    case('tb','bb')

       select case (ii)
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/ 5000  , 4 , 30000  , 4 /)
       case(2)                                      ! virtual
          ivegas(1:4) = (/  500  , 4 , 5000  , 6 /)
       case(3,4)                                    ! real emission matrix element
          ivegas(1:4) = (/ 10000 , 4 , 50000  , 6 /)
       case default       
          ivegas(1:4) = (/ 100   , 4 ,  1000  , 2 /)
       end select

    case('lq')

       select case (ii)
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/ 5000  , 4 , 30000  , 4 /)
       case(2)                                      ! virtual
          ivegas(1:4) = (/  500  , 4 ,  5000  , 6 /)
       case(3,4)                                    ! real emission matrix element
          ivegas(1:4) = (/ 10000 , 4 , 50000  , 6 /)
       case default       
          ivegas(1:4) = (/ 100   , 4 ,  1000  , 2 /)
       end select

    case('le')

       select case (ii)
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/   5000, 4 , 30000 , 4 /)
       case(2)                                      ! virtual
          ivegas(1:4) = (/    500, 4 , 15000 , 6 /)
       case(3)                                      ! qg matrix element 
          ivegas(1:4) = (/   3000 ,4 , 80000 , 6 /)
       case(4)                                      ! gg matrix element 
          ivegas(1:4) = (/   3000 ,4 , 30000 , 6 /)
       case(5)                                      ! gg on-shell term
          ivegas(1:4) = (/  30000 ,4 ,100000 , 6 /)
       case(6)                                      ! gg os subtraction
          ivegas(1:4) = (/  80000 ,6 ,100000 ,12 /)
       case(7)                                      ! qq matrix element 
          ivegas(1:4) = (/  10000 ,4 , 40000 , 6 /)
       case(8)                                      ! qb matrix element 
          ivegas(1:4) = (/  10000 ,4 ,100000 , 6 /)
       case(9)                                      ! qb os subtraction 
          ivegas(1:4) = (/  10000 ,4 , 40000 , 6 /)
       case(10)                                     ! qb os subtraction
          ivegas(1:4) = (/  10000 ,4 ,100000 , 6 /)
       case default       
          ivegas(1:4) = (/    100 ,4 ,  1000 , 2 /)
       end select

    case('hh')

       select case (ii)
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/  5000 , 4 , 30000 , 4 /)
       case(2,3)                                    ! virtual and virtual-susy
          ivegas(1:4) = (/ 10000 , 8 , 20000 , 4 /)
       case(4)                                      ! real emission
          ivegas(1:4) = (/ 10000 , 8 , 80000 , 4 /)
       case(5,7)                                    ! crossed matrix element
          ivegas(1:4) = (/ 10000 , 8 , 50000 , 4 /)
       case(6,8)                                    ! o-s subtraction
          ivegas(1:4) = (/ 10000 , 8 , 50000 , 4 /)
       case default       
          ivegas(1:4) = (/ 100   , 4 ,  1000 , 2 /)
       end select

    case('ht')
       
       select case (ii)  
       case(-1,0,1)                                 ! born
          ivegas(1:4) = (/  5000 , 4 ,  30000 , 4 /)
       case(2,3)                                    ! virtual and virtual-susy
          ivegas(1:4) = (/ 2000  , 8 ,  20000 , 4 /)
       case(4)                                      ! real qg matrix element
          ivegas(1:4) = (/ 20000 , 8 , 100000 , 4 /)
       case(5)                                      ! gg matrix element 
          ivegas(1:4) = (/ 20000 , 8 ,  50000 , 4 /)
       case(6,7)                                    ! qb, qq matrix element 
          ivegas(1:4) = (/ 20000 , 8 , 200000 , 4 /)
       case(8,9)                                    ! qb, gg on-shell subtraction
          ivegas(1:4) = (/ 20000 , 8 , 100000 , 4 /)
       case default       
          ivegas(1:4) = (/ 100   , 4 , 1000   , 2 /)
       end select

    case default
       print *, " INIT_VEGAS: old channels "
    end select

    if ( (i_lofast == 1).and.( ii <= 1 ) ) then                           ! to speed up at expense of precision
       print *, " INIT_VEGAS: speeding up via i_lofast "
       ivegas(1:4) = (/  500  , 4 ,  2000  , 4 /)
    end if

  end subroutine INIT_VEGAS
       
! ------------------------------
  subroutine INIT_IFAST(isq_gl,inext_sq,ifast)
    
    integer,              intent(in)  :: isq_gl,inext_sq
    integer,              intent(out) :: ifast

    ifast = 1                                                             ! make fast mode the default

    if ( isq_gl == 0 ) then                                               ! new channels
       if ((ii == -1 ).and.(inext_sq == 0 )) then                         ! ii=-1 is only appearance of squark loops 
          print *, " INIT_IFAST: first term summing squarks ",isq_gl
          ifast = 0
       end if
       
       if ((ifast == 1).and.(ii_done(ii-1) /= 1)) then                    ! check if channel before has been computed
          print*, " INIT_IFAST: recent ii not used         ",ii,ii_done(ii-1)
          ifast = 0
       end if
       
       if ((ii > -1).and.(ifast == 1).and.(dim(ii) /= dim(ii-1))) then   ! check if dimensionality is the same
          print*, " INIT_IFAST: new number of dimensions   ",ii,dim(ii),dim(ii-1)
          ifast = 0  
       end if
    else if ( isq_gl == 1 ) then
       if (inext_sq == 0 ) then                                           ! ii=-1 is only appearance of squark loops 
          print *, " INIT_IFAST: first term summing squarks ",isq_gl
          ifast = 0
       end if
    else 
       print*, " INIT_IFAST: problem with isq_gl ",isq_gl
       call HARD_STOP
    end if

  end subroutine INIT_IFAST

! -----------------------------------------------------------------------!
! note that INIT_ALPHAS is run inside, in case sugra routine gets called !
! GET_SPECTRUM is the interface to some SUSY spectrum code               !
!                                                                        !
! conventions for weak-scale MSSM parameters used inside Prospino2:      !
!                                                                        !
!       lowmass(0)  mu                                                   !
!       lowmass(1)  m_1                                                  !
!       lowmass(2)  m_2                                                  !
!       lowmass(3)  m_3                                                  !
!                                                                        !
!       lowmass(4)  gluino mass                                          !
!       lowmass(5)  \                                                    !
!       lowmass(6)   \                                                   !
!       lowmass(7)   /  neutralino masses [with sign]                    !
!       lowmass(8)  /                                                    !
!       lowmass(9)  \                                                    !
!       lowmass(10) /   chargino masses                                  !
!                                                                        !
!       lowmass(15) degenerate squark mass (8)                           !
!       lowmass(16) degenerate squark mass (10)                          !
!                                                                        !
!       lowmass(21) a_b                                                  !
!       lowmass(24) a_t                                                  !
!                                                                        !
!       lowmass(30) selectron_l mass                                     !
!       lowmass(31) selectron_r mass                                     !
!       lowmass(32) selectron-neutrino mass                              !
!       lowmass(33) stau_1 mass                                          !
!       lowmass(34) stau_2 mass                                          !
!       lowmass(35) stau-neutrino mass                                   !
!       lowmass(36) a_tau                                                !
!                                                                        !
!       lowmass(40) cp odd higgs mass                                    !
!       lowmass(41) light cp even higgs mass                             !
!       lowmass(42) heavy cp even higgs mass                             !
!       lowmass(43) charged higgs mass                                   !
!       lowmass(44) sin(alpha)                                           !
!       lowmass(45) cos(alpha)                                           !
!                                                                        !
!       like cteq: u,d,s,c,b,t first L then R                            !
!       lowmass(51) sup_L                                                !
!       lowmass(52) sdown_L                                              !
!       lowmass(53) sstrange_L                                           !
!       lowmass(54) scharm_L                                             !
!       lowmass(55) sbottom_1                                            !
!       lowmass(56) stop_1 (used as scalar LQ mass as well)              !
!       lowmass(57) sup_R                                                !
!       lowmass(58) sdown_R                                              !
!       lowmass(59) sstrange_R                                           !
!       lowmass(60) scharm_R                                             !
!       lowmass(61) sbottom_2                                            !
!       lowmass(62) stop_2                                               !
!                                                                        !
!       lowmass(80) unification scale                                    !
!       lowmass(81) unified coupling alpha(m_x)                          !
!                                                                        !
!       lowmass(91) trilinear higgs coupling lambda(1)                   !
!       .......                                                          !
!       lowmass(97) lambda(7)                                            !
!                                                                        !
!       lowmass(99)                                                      !
!                                                                        !
!       bwmix neutralino mixing matrix (bino-wino)                       !
!       pzmix neutralino mixing matrix (photino-zino)                    !
!       uumix chargino mixing matrix u                                   !
!       vvmix chargino mixing matrix v                                   !
!                                                                        !
!-------------------------------------------------------------------------
  subroutine INIT_SUSY(nin1,nin2,ipart1_in,ipart2_in,isq_ng_in,run0,run1,unimass,lowmass)
    
    integer,           intent(in) :: nin1,nin2
    integer,           intent(in) :: ipart1_in,ipart2_in,isq_ng_in
    real(kind=double), intent(in) :: run0,run1
    real(kind=double), dimension(1:20), intent(out) :: unimass
    real(kind=double), dimension(0:99), intent(out) :: lowmass

    integer                              :: n
    real(kind=double), dimension(2,2)    :: uu_in,vv_in
    real(kind=double), dimension(2,2)    :: mst_in,msb_in,msl_in
    real(kind=double), dimension(4,4)    :: bw_in
    real(kind=double)                    :: sin2x
    complex(kind=double), dimension(1:4) :: sig

    unimass(1:20) = 0.0

    call GET_SPECTRUM(nin2,unimass,lowmass,bw_in,uu_in,vv_in,mst_in,msb_in,msl_in)

    uu(1:2,1:2)  = uu_in(1:2,1:2)      ! set the mixing matrices as defined globally in Xvital.f90
    vv(1:2,1:2)  = vv_in(1:2,1:2)
    bw(1:4,1:4)  = bw_in(1:4,1:4)
    mst(1:2,1:2) = mst_in(1:2,1:2)
    msb(1:2,1:2) = msb_in(1:2,1:2)
    msl(1:2,1:2) = msl_in(1:2,1:2)

    tan_b        = unimass(10)
    mu_susy      = lowmass(0)          ! mu parameter
    mg           = lowmass(4)          ! gluino mass 
    ms           = lowmass(15)         ! degenerate squark mass (still needed for old prospino code)
    a_b          = lowmass(21)         ! tri-scalar coupling
    a_t          = lowmass(24)         ! tri-scalar coupling
    mh1          = lowmass(41)         ! light scalar higgs mass 
    mh2          = lowmass(42)         ! heavy scalar higgs mass 
    mch          = lowmass(43)         ! charged higgs mass 
    sin_a        = lowmass(44)
    cos_a        = lowmass(45)

    smass_n(1:6) = lowmass(5:10)       ! neutralino/chargino masses 
    smass_n(7:8) = smass_n(5:6)

    if (isq_ng_in == 1 ) then          ! allow for non-degenerate squark masses
       msq(-6:-1) = lowmass(56:51:-1)  ! the squarkL masses: u,d,s,c,b,t
       msq( 1: 6) = lowmass(57:62:1)   ! the squarkR masses: u,d,s,c,b,t
    else if (isq_ng_in == 0 ) then     ! or if switched off set all masses to the average value
       msq(-6:-1) = ms
       msq( 1: 6) = ms
    end if

!tp    mg = mg * 1.5**run0                  ! examples for the use of run0 and run1
!tp    mg = mg + run1
!tp    msq(+6) = msq(+6) + run0
!tp    msq(-6) = msq(-6) + run0
!tp    print*, " mass shift ",mg

    if (imx==1) then                   ! neutralino masses, switch only works for Born term
       mass_n(1:8) = abs(smass_n(1:8)) ! positive mass array for complex mixing matrix
       sig(1:4)  = (1.0,0.0)           ! set correction factor eta 
       do n=1,4 
          if (smass_n(n)<0) sig(n) = (0.0,1.0)  
          zz(n,1:4) = sig(n) * bw(n,1:4) 
       end do
    else if (imx==0) then
       mass_n(1:8) = smass_n(1:8)
       zz(1:4,1:4) = bw(1:4,1:4)
    end if

    if (final_state=='nn') then        ! this is for the neutralino pairs only
       if (ipart1_in < ipart2_in) then ! ordered by charge: +1 -2
          ipart1 = ipart1_in
          ipart2 = ipart2_in
       else 
          ipart2 = ipart1_in
          ipart1 = ipart2_in
       end if
    else 
       ipart1 = ipart1_in
       ipart2 = ipart2_in              ! where ipart2 is not really needed
    end if
       
    mass_s(1:4) = 0.0
    if (final_state=='ll') then 
       select case (ipart1)
       case(0) 
          mass_s(1)   = ( lowmass(30) + lowmass(31) )/2.D0 
          mass_s(2)   = ( lowmass(30) + lowmass(31) )/2.D0 
       case(1)
          mass_s(1)   = lowmass(30)
          mass_s(2)   = lowmass(30)
       case(2) 
          mass_s(1)   = lowmass(31)
          mass_s(2)   = lowmass(31)
       case (3)
          mass_s(1)   = lowmass(32)
          mass_s(2)   = lowmass(32)
       case(4,5)
          mass_s(1)   = lowmass(30)
          mass_s(2)   = lowmass(32)
       case(6)
          mass_s(1)   = lowmass(33)
          mass_s(2)   = lowmass(33)
          call COMPUTE_SCALAR_ANGLE(msl,sin2x)
          mass_s(3)   = sin2x
       case(7)
          mass_s(1)   = lowmass(34)
          mass_s(2)   = lowmass(34)
          call COMPUTE_SCALAR_ANGLE(msl,sin2x)
          mass_s(3)   = sin2x
       case(8)
          mass_s(1)   = lowmass(33)
          mass_s(2)   = lowmass(34)
          call COMPUTE_SCALAR_ANGLE(msl,sin2x)
          mass_s(3)   = sin2x
       case(9)
          mass_s(1)   = lowmass(35)
          mass_s(2)   = lowmass(35)
       case(10,11)
          mass_s(1)   = lowmass(33)
          mass_s(2)   = lowmass(35)
          call COMPUTE_SCALAR_ANGLE(msl,sin2x)
          mass_s(3)   = sin2x
       case(12,13)
          mass_s(1)   = lowmass(34)
          mass_s(2)   = lowmass(35)
          call COMPUTE_SCALAR_ANGLE(msl,sin2x)
          mass_s(3)   = sin2x
       case(14)                                      ! charged Higgs pairs
          mass_s(1)   = lowmass(43)
          mass_s(2)   = lowmass(43)
       end select
    else if ((final_state=='tb').or.(final_state=='xx')) then 
       mass_s(1) = msq(-6)                           ! different from the stau syntax!!!
       mass_s(2) = msq(+6)
       mass_s(3) = 2*mst(1,1)*mst(1,2)               ! like in the Form code
       mass_s(4) = mst(1,1)**2 - mst(1,2)**2 
       mass_x(1) = msq(-5)                           ! also fix the sbottom sector  
       mass_x(2) = msq(+5)
       mass_x(3) = 2*msb(1,1)*msb(1,2)               ! like in the Form code
       mass_x(4) = msb(1,1)**2 - msb(1,2)**2 
    else if (final_state=='lq') then 
       mass_s(1) = msq(-6)         
       mass_s(2) = msq(-6)
    else if (final_state=='le') then 
       mass_s(1) = msq(-6)         
    else if (final_state=='bb') then 
       mass_s(1) = msq(-5)                           ! different from the stau syntax!!!
       mass_s(2) = msq(+5)
       mass_s(3) = 2*msb(1,1)*msb(1,2)               ! like in the Form code
       mass_s(4) = msb(1,1)**2 - msb(1,2)**2 
       mass_x(1) = msq(-6)                           ! also fix the stop sector  
       mass_x(2) = msq(+6)
       mass_x(3) = 2*mst(1,1)*mst(1,2)               ! like in the Form code
       mass_x(4) = mst(1,1)**2 - mst(1,2)**2 
    end if

    if (final_state=='hh') then                      ! report removed s-channel coupligs
       if ( 2.0*mch < mz  ) print*, " INIT_SUSY: removed coupling for on-shell decay  Z->H+H- "
       if ( 2.0*mch < mh1 ) print*, " INIT_SUSY: removed coupling for on-shell decay h0->H+H- "
       if ( 2.0*mch < mh2 ) print*, " INIT_SUSY: removed coupling for on-shell decay H0->H+H- "
    end if
       
    call DECOUPLE_SPECTRUM(ipart1_in,ipart2_in)
    
  end subroutine INIT_SUSY

! ------------------------------
! routine to take care of critical mass choices 
! solve by overwriting in xx_pass_integ
  subroutine DECOUPLE_SPECTRUM(ipart1_in,ipart2_in)

    integer, intent(in) :: ipart1_in,ipart2_in
    integer             :: i1
    real(kind=double)   :: mn_min,ml_min,ms_min,ms3_min

    if (abs(ms-mg) < 1.0) then
       ms = ms + 1.0
       print*, " DECOUPLE_SPECTRUM: spectrum degenerate, change ms ",ms
    end if

    if ( (final_state=='ht') .and. (abs(mch-mt) < 1.0) ) then               ! only a problem for charged Higgs
       mch = mch + 1.0
       print*, " DECOUPLE_SPECTRUM: spectrum degenerate, change mH ",mch
    end if

    mn_min = min( abs(mass_n(ipart1_in)) ,abs(mass_n(ipart2_in)) )          ! smallest final-state mass for nn channel
    ml_min = min( abs(mass_s(1)), abs(mass_s(2)) )                          ! smallest final-state mass for ll channel
    ms3_min = mass_s(ipart1_in)                                             ! final-state mass for bb and tb channels

    ms_min = 1.e10
    do i1=-6,6
       if (i1==0) cycle
       ms_min = min(msq(i1),ms_min) 
    end do

    if ( (final_state=='gg') .and. ( (ms-mg)/mg > 10) ) then                ! case with virtual-only squarks in Prospino1
       ms = (10+1) * mg 
       print*, " DECOUPLE_SPECTRUM: squark mass decoupled: ms=",ms
    end if

    if ( ((final_state=='ss') .or.                                         &
          (final_state=='sb')     ).and. ( (mg-ms_min)/ms_min > 10) ) then  ! case with virtual-only gluinos in Prospino1
       mg = (10+1) * ms_min
       print*, " DECOUPLE_SPECTRUM: gluino mass decoupled: mg=",mg
    end if

    if (final_state=='bb') then                                             ! case with virtual-only gluinos in Prospino2
!tp       print*, " ratio mg ",(mg-ms3_min)/ms3_min
!tp       print*, " ratio ms ",(ms-ms3_min)/ms3_min
       if ( (mg-ms3_min)/ms3_min > 2)  then
          mg = (2+1) * ms3_min
          print*, " DECOUPLE_SPECTRUM: gluino mass decoupled: mg=",mg
       end if
       if ( (ms-ms3_min)/ms3_min > 20)  then
          ms = (20+1) * ms3_min
          print*, " DECOUPLE_SPECTRUM: squark mass decoupled: mg=",ms
       end if
    end if

    if (final_state=='tb') then                                             ! case with virtual-only gluinos in Prospino2
       if ( (mg-ms3_min)/ms3_min > 5) then
          mg = (5+1) * ms3_min
          print*, " DECOUPLE_SPECTRUM: gluino mass decoupled: mg=",mg
       end if
       if ( (ms-ms3_min)/ms3_min > 20) then
          ms = (20+1) * ms3_min
          print*, " DECOUPLE_SPECTRUM: squark mass decoupled: ms=",ms
       end if
    end if

    if (final_state=='nn') then                                             ! case with mixed Born/virtual squarks in Prospino2
       if ( (ms-mn_min)/mn_min > 10) then                                   ! already there for the t channel propagator, using ms for NLO part
          ms = (10+1) * mn_min
          print*, " DECOUPLE_SPECTRUM: squark mass decoupled: ms=",ms
       end if
       if ( (mg-mn_min)/mn_min > 15) then                                   ! only in the actual loop
          mg = (15+1) * mn_min
          print*, " DECOUPLE_SPECTRUM: gluino mass decoupled: mg=",mg
       end if
    end if

    if (final_state=='ll') then                                             ! case with virtual-only squarks in Prospino2
       if ( (ms-ml_min)/ml_min > 15) then                                   ! only in the actual loop
          ms = (15+1) * ml_min
          print*, " DECOUPLE_SPECTRUM: squark mass decoupled: ms=",ms
       end if
       if ( (mg-ml_min)/ml_min > 15) then                                   ! only in the actual loop
          mg = (15+1) * ml_min
          print*, " DECOUPLE_SPECTRUM: gluino mass decoupled: mg=",mg
       end if
    end if

    if (final_state=='ng') then                                             ! case with mixed Born/virtual squarks in Prospino2
       if ( (ms-mg)/mg > 10) then                                   ! already there for the t channel propagator, using ms for NLO part
          ms = (10+1) * mg
          print*, " DECOUPLE_SPECTRUM: squark mass decoupled: ms=",ms
       end if
    end if

  end subroutine DECOUPLE_SPECTRUM

! ------------------------------
! note that in the old prospino there is routine INIT_ALPHAS_ARG(inlo) for an independent outside call
  subroutine INIT_ALPHAS
    
    integer           :: nq,inlo_alphas               ! note that inlo_alphas is an internal parameter
    real(kind=double) :: lam_dummy,acc,lam_qcd,mcq,mbq,mtq

    if (ii<=0) then 
       inlo_alphas = 0 
    else 
       inlo_alphas = 1 
    end if

    call GET_LAMBDA_QCD(inlo_alphas,lam_dummy)

    lam_qcd = lam_dummy
    acc = 1.e-8
    mcq = 1.5              
    mbq = mb
    mtq = 175000.0                                    ! top quark decoupled from alphas
    nq  = 5 

    call ALSINI(acc,lam_qcd,mcq,mbq,mtq,nq)

    mtq = mt                                          ! top quark with finite mass for running Yukawas
    nq  = 6
    call ALSINI_RUNM(acc,lam_qcd,mcq,mbq,mtq,nq)

  end subroutine INIT_ALPHAS

! ------------------------------
! fill all common blocks used in the old squark/gluino routines
  subroutine FILL_COMMONS_SQ_GL(inlo,scafac)
    
    integer,           intent(in) :: inlo
    real(kind=double), intent(in) :: scafac

    real(kind=double) :: s_1,energy_1,alphas_1,ms_1,mg_1,mt_1
    common/CONST1/       s_1,energy_1,alphas_1,ms_1,mg_1,mt_1                       ! energy and masses filled here, s and alphas filled in hadron**.f

    real(kind=double), dimension(-6:6) :: msq_common    
    integer                            :: isquark1_common,isquark2_common,i_ngtest_common
    common/SQUARKS/                       msq_common,isquark1_common,isquark2_common,i_ngtest_common
    
    real(kind=double) :: s_2,energy_2,alphas_2,mst1_2,mg_2,mt_2,ms_2,mst2_2,sin2t_2
    common/CONST6/       s_2,energy_2,alphas_2,mst1_2,mg_2,mt_2,ms_2,mst2_2,sin2t_2 ! energy and masses filled here, s and alphas filled in hadron**.f

    real(kind=double) :: scale,scafac_1                                             ! scale filled in hadron**.f, all otherse filled here
    integer           ::                icoll_1,iscapt
    common/CONST2/       scale,scafac_1,icoll_1,iscapt

    integer           :: ilo_common, inlo_common, ionlylo_common, idg_common, ing_common
    common/CONST5/       ilo_common, inlo_common, ionlylo_common, idg_common, ing_common  ! first three filled here
    
    integer           :: icharconj
    common/CHARCONJ/     icharconj                                                  ! for now dummy block

    real(kind=double) :: ptmin, ptmax
    common/CUT1/         ptmin, ptmax                                               ! for now dummy block, filled in here

    real(kind=double) :: ymin, ymax
    common/CUT2/         ymin, ymax                                                 ! for now dummy block, filled in here

    integer           :: iflavor, itotal
    common/FLAVOR/       iflavor, itotal

    energy_1 = sqrt( sc )                                                           ! collider energy
    ms_1     = ms                                                                   ! make sure INIT_SUSY is called first
    mg_1     = mg
    mt_1     = mt

    energy_2 = sqrt( sc )                                                           ! collider energy
    ms_2     = ms                                                                   ! make sure INIT_SUSY is called first
    mg_2     = mg
    mt_2     = mt
    mst1_2   = msq(-6)
    mst2_2   = msq( 6)
    call COMPUTE_SCALAR_ANGLE(mst,sin2t_2)                                          ! compute scalar mixing angle from matrix

    i_ngtest_common = i_ngtest                                                      ! on switch (i_ngtest=1) only for testing purpose
    msq_common(-6:6) = msq(-6:6)                                                    ! also needs INIT_SUSY
    isquark1_common  = isquark1                                                     ! needed by PROSPINO1
    if (isquark1_common==9) isquark1_common=1                                       ! overwrite for mass degenerate case

    isquark2_common  = isquark2
    if (isquark2_common==9) isquark2_common=1

    scafac_1 = scafac
    icoll_1  = icoll                                                                ! this is just the Tevatron/LHC flag now, as set in prospino_main.f90
    iscapt   = 0                                                                    ! always use the average final state mass as the scale 

    ilo_common  = 50000                                                                    ! vegas iterations, but no impact here, go to initpdf.f
    inlo_common = 50000

    if (inlo==0) then                                                                ! leading order only flag
       ionlylo_common = 1
    else if (inlo==1) then 
       ionlylo_common = 0 
    else 
       print*, " FILL_COMMONS_SQ_GL: problem with inlo ",inlo
       call HARD_STOP
    end if

    icharconj = 0                                                                   ! for distributions only: with respact to s or b
    
    ptmin = 0.0                                                                     ! for distributions only
    ptmax = energy_1

    ymin =   0.0                                                                    ! for distributions only
    ymax = 100.0

    iflavor = 0                                                                     ! all initial state flavors available
    itotal  = 1                                                                     ! total cross section without cuts, faster version 
    
  end subroutine FILL_COMMONS_SQ_GL

end module xx_initialize



