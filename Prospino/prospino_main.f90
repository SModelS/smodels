program main
  use xx_kinds
  use xx_prospino_subroutine
  implicit none

  integer                              :: inlo,isq_ng_in,icoll_in,i_error_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in
  logical                              :: lfinal
  character(len=2)                     :: final_state_in

!----------------------------------------------------------------------------
  inlo = 1          ! specify LO only[0] or complete NLO (slower)[1]        !
!                   ! results: LO     - leading order, degenerate squarks   !
!                   !          NLO    - NLO, degenerate squarks             !
!                   !          LO_ms  - leading order, free squark masses   !
!                   !          NLO_ms - NLO, free squark masses             !
!                   ! all numerical errors (hopefully) better than 1%       !
!                   ! follow Vergas iteration on screen to check            !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  isq_ng_in = 1     ! specify degenerate [0] or free [1] squark masses      !
                    ! [0] means Prospino2.0 with average squark masses      !
                    ! [0] invalidates isquark_in switch                     !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  icoll_in = 3      ! collider : tevatron[0], lhc14[1], lhc7[2], lhc8[3]    !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  i_error_in = 0    ! with central scale [0] or scale variation [1]         !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  final_state_in = 'll'                                                     !
!                                                                           !
!                   ng     neutralino/chargino + gluino                     !
!                   ns     neutralino/chargino + squark                     !
!                   nn     neutralino/chargino pair combinations            !
!                   ll     slepton pair combinations                        !
!                   sb     squark-antisquark                                !
!                   ss     squark-squark                                    !
!                   tb     stop-antistop                                    !
!                   bb     sbottom-antisbottom                              !
!                   gg     gluino pair                                      !
!                   sg     squark + gluino                                  !
!                   lq     leptoquark pairs (using stop1 mass)              !
!                   le     leptoquark plus lepton (using stop1 mass)        !
!                   hh     charged Higgs pairs (private code only!)         !
!                   ht     charged Higgs with top (private code only!)      !
!                                                                           !
!  squark and antisquark added, but taking into account different sb or ss  !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  ipart1_in = 1                                                             !
  ipart2_in = 1                                                             !
!                                                                           !
!  final_state_in = ng,ns,nn                                                !
!  ipart1_in   = 1,2,3,4  neutralinos                                       !
!                5,6      positive charge charginos                         !
!                7,8      negative charge charginos                         !
!  ipart2_in the same                                                       !
!      chargino+ and chargino- different processes                          !
!                                                                           !
!  final_state_in = ll                                                      !
!  ipart1_in   = 0        sel,sel + ser,ser  (first generation)             !
!                1        sel,sel                                           !
!                2        ser,ser                                           !
!                3        snel,snel                                         !
!                4        sel+,snl                                          !
!                5        sel-,snl                                          !
!                6        stau1,stau1                                       !
!                7        stau2,stau2                                       !
!                8        stau1,stau2                                       !
!                9        sntau,sntau                                       !
!               10        stau1+,sntau                                      !
!               11        stau1-,sntau                                      !
!               12        stau2+,sntau                                      !
!               13        stau2-,sntau                                      !
!               14        H+,H- in Drell-Yan channel                        !
!                                                                           !
!  final_state_in = tb and bb                                               !
!  ipart1_in   = 1        stop1/sbottom1 pairs                              !
!                2        stop2/sbottom2 pairs                              !
!                                                                           !
!  note: otherwise ipart1_in,ipart2_in have to set to one if not used       !
!                                                                           !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  isquark1_in = 0                                                           !
  isquark2_in = 0                                                           !
!                                                                           !
!  for LO with light-squark flavor in the final state                       !
!  isquark1_in     =  -5,-4,-3,-2,-1,+1,+2,+3,+4,+5                         !
!                    (bL cL sL dL uL uR dR sR cR bR) in CteQ ordering       !
!  isquark1_in     = 0 sum over light-flavor squarks throughout             !
!                      (the squark mass in the data files is then averaged) !
!                                                                           !
!  flavors in initial state: only light-flavor partons, no bottoms          !
!                            bottom partons only for Higgs channels         !
!                                                                           !
!  flavors in final state: light-flavor quarks summed over five flavors     !
!                                                                           !
!----------------------------------------------------------------------------
  
  call PROSPINO_OPEN_CLOSE(0)                                                            ! open all input/output files
  
  call PROSPINO_CHECK_HIGGS(final_state_in)                                              ! lock Higgs final states
  call PROSPINO_CHECK_FS(final_state_in,ipart1_in,ipart2_in,lfinal)                      ! check final state 
  if (.not. lfinal ) then
     print*, " final state not correct ",final_state_in,ipart1_in,ipart2_in
     call HARD_STOP                                                                      ! finish if final state bad
  end if

  call PROSPINO(inlo,isq_ng_in,icoll_in,i_error_in,final_state_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in) ! actual prospino call
        
!----------------------------------------------------------------------------
!  input file: prospino.in.leshouches                                       !
!              use block MASS for masses, plus low-energy mixing matrices   !
!                                                                           !
!  output file: prospino.dat   for compact output format                    !
!               prospino.dat2  for long output including subchannels        !
!               prospino.dat3  lo file for masses, flags, etc               !
!----------------------------------------------------------------------------
  call PROSPINO_OPEN_CLOSE(1)                                                            ! close all input/output files 

end program main


