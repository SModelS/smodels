program main
  use xx_kinds
  use xx_prospino_subroutine
  implicit none

  integer                              :: inlo,isq_ng_in,icoll_in,i_error_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in
  integer                              :: ifs
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
  icoll_in = 1      ! specify the collider : tevatron[0], lhc14[1], lhc7[2] !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  i_error_in = 0    ! with central scale [0] or scale variation [1]         !
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!  final_state_in = ng   ! neutralino/chargino + gluino                     !
!                   ns   ! neutralino/chargino + squark                     !
!                   nn   ! neutralino/chargino pair combinations            !
!                   ll   ! slepton pair combinations                        !
!                   sb   ! squark-antisquark                                !
!                   ss   ! squark-squark                                    !
!                   tb   ! stop-antistop                                    !
!                   bb   ! sbottom-antisbottom                              !
!                   gg   ! gluino pair                                      !
!                   sg   ! squark + gluino                                  !
!                   lq   ! leptoquark pairs (using stop1 mass)              !
!                   le   ! leptoquark plus lepton (using stop1 mass)        !
!                   hh   ! charged Higgs pairs (private code only!)         !
!                   ht   ! charged Higgs with top (private code only!)      !
!                   xx   ! stop-antistop, wrapped Propino1                  !
!                                                                           !
!  squark and antisquark added, but taking into account different sb or ss  !
!----------------------------------------------------------------------------
  final_state_in = 'le'

!----------------------------------------------------------------------------
!  final_state_in = ng,ns,nn                                                !
!  ipart1_in   = 1,2,3,4  neutralinos                                       !
!                5,6      positive charge charginos                         !
!                7,8      negative charge charginos                         !
!  ipart2_in the same                                                       !
!      chargino+ and chargino- different processes                          !
!                                                                           !
!  final_state_in = ll                                                      !
!  ipart1_in   = 0        sel,sel + ser,ser                                 !
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
  ipart1_in = 1
  ipart2_in = 1 
  isquark1_in = 0
  isquark2_in = 0
  
  call PROSPINO_OPEN_CLOSE(0)                                                            ! open all input/output files
  
  do i_error_in=0,1,1                                                                    ! loop for test purposes only
  do icoll_in=0,2,2                                                                      ! loop for test purposes only
  do ifs = 1,15,1                                                                        ! loop and if structure for test purposes only
  if (ifs==1) then          
     final_state_in = 'ng'     
  else if (ifs==2) then     
     final_state_in = 'ns'     
  else if (ifs==3) then     
     final_state_in = 'nn'     
  else if (ifs==4) then 
     final_state_in = 'll'
  else if (ifs==5) then 
     final_state_in = 'gg'
  else if (ifs==6) then 
     final_state_in = 'sb'
  else if (ifs==7) then 
     final_state_in = 'sg'
  else if (ifs==8) then 
     final_state_in = 'ss'
  else if (ifs==9) then 
     final_state_in = 'tb'
  else if (ifs==10) then 
     final_state_in = 'bb'
  else if (ifs==11) then 
     final_state_in = 'lq'
  else if (ifs==12) then 
     final_state_in = 'le'
  else if (ifs==13) then 
     final_state_in = 'hh'
  else if (ifs==14) then 
     final_state_in = 'ht'
  else if (ifs==15) then 
     final_state_in = 'xx'
  end if
  
  do ipart1_in = 1,14,1                                                                   ! loop for test purposes only
  do ipart2_in = ipart1_in,6,1                                                                   ! loop for test purposes only

  call PROSPINO_CHECK_FS(final_state_in,ipart1_in,ipart2_in,lfinal)                      ! check final state 
  if (.not. lfinal ) cycle
     
  call PROSPINO(inlo,isq_ng_in,icoll_in,i_error_in,final_state_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in) ! actual prospino call
        
  end do                                                                                 ! close both ipart loops 
  end do
  end do                                                                                 ! close the ifs dummy loop
  end do                                                                                 ! close the icoll_in loop
  end do                                                                                 ! close the i_error_in loop

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


