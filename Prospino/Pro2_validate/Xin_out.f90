! ===========================================================================================================
module xx_in_out
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ   
  implicit none 
  private
  public :: DAT1,DAT2,DAT3
contains
! ------------------------------
  subroutine DAT1(ndat1,run0,run1,ms1_print,ms2_print,fct,rel,kfac,kng)             ! the short output e.g. for paw

    integer,                             intent(in) :: ndat1
    real(kind=double), dimension(-1:22), intent(in) :: fct,rel
    real(kind=double),                   intent(in) :: run0,run1,ms1_print,ms2_print,kfac,kng
    character(len=120), parameter ::                                             &  ! short output format
         form1="(a2,1x,2(i2,1x),1x,5(f6.1,1x),1(f6.3,1x),4(g9.3,1x),(f6.4,1x),2(g9.3,1x))"
    character(len=120), parameter ::                                             &  ! short output format
         form2="(a2,1x,2(i2,1x),1x,5(f6.1,1x),1(f6.3,1x),4(g9.3,1x),(f6.4,1x),1(g9.3,1x))"

    integer           :: i_zero
    real(kind=double) :: r_zero

    i_zero = 0
    r_zero = 0.0

    select case(final_state)
    case ('ng')
       write(ndat1, fmt=form1) final_state,ipart1,i_zero,run0,run1,scafac,mass_n(ipart1),mg,r_zero             &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('ns')
       write(ndat1, fmt=form1) final_state,ipart1,i_zero,run0,run1,scafac,mass_n(ipart1),ms1_print,r_zero      &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('nn')
       write(ndat1, fmt=form1) final_state,ipart1,ipart2,run0,run1,scafac,mass_n(ipart1),mass_n(ipart2),r_zero &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('ll')
       write(ndat1, fmt=form1) final_state,ipart1,i_zero,run0,run1,scafac,mass_s(1:3)                          &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('gg')
       write(ndat1, fmt=form1) final_state,i_zero,i_zero,run0,run1,scafac,mg,mg,r_zero                         &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('sb','ss')
       write(ndat1, fmt=form1) final_state,i_zero,i_zero,run0,run1,scafac,ms1_print,ms2_print,r_zero           &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('sg')
       write(ndat1, fmt=form1) final_state,i_zero,i_zero,run0,run1,scafac,ms1_print,mg,r_zero                  &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('tb','bb')
       if (ipart1.eq.1) then
       write(ndat1, fmt=form1) final_state,ipart1,i_zero,run0,run1,scafac,mass_s(1),mass_s(1),mass_s(3)        &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
       elseif (ipart1.eq.2) then
       write(ndat1, fmt=form1) final_state,ipart1,i_zero,run0,run1,scafac,mass_s(2),mass_s(2),-mass_s(3)       &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
       else
       print*," DAT1: problem with stop/sbottom sector, ipart1 = ",ipart1
       call HARD_STOP
    end if
    case ('lq')
       write(ndat1, fmt=form1) final_state,ipart1,i_zero,run0,run1,scafac,mass_s(1),mass_s(1),r_zero           &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('le')
       write(ndat1, fmt=form1) final_state,ipart1,i_zero,run0,run1,scafac,mass_s(1),r_zero,r_zero              &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('xx')
       write(ndat1, fmt=form1) final_state,ipart1,i_zero,run0,run1,scafac,mass_s(1),mass_s(1),r_zero           &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(21),fct(22)
    case ('hh')
       write(ndat1, fmt=form2) final_state,ipart1,i_zero,run0,run1,scafac,mch,mch,tan_b                        &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(22)
    case ('ht')
       write(ndat1, fmt=form2) final_state,ipart1,i_zero,run0,run1,scafac,mch,mt,tan_b                         &
                              ,fct(0),rel(0),fct(20),rel(20),kfac,fct(22)
    end select

  end subroutine DAT1

! ------------------------------
  subroutine DAT2(ndat2,run0,run1,ms1_print,ms2_print,fct,rel,kfac,kng)   ! the long output

    integer,                             intent(in) :: ndat2
    real(kind=double), dimension(-1:22), intent(in) :: fct,rel
    real(kind=double),                   intent(in) :: run0,run1,ms1_print,ms2_print,kfac,kng
    character(len=500), parameter ::           &                          ! output format
         form1="(a2,1x,2(i2,3x),2(f10.4,1x),/  &
               &,4(f10.4,1x),/                 &       
               &,4(f10.4,1x),/,2(f10.4,1x),/,/ &       
               &,5(g10.4,1x),/,5(g10.4,1x),/,/ &
               &,4(g10.4,1x),/,4(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::            &
         form2="(a2,1x,2(i2,3x),2(f10.4,1x),/   &
               &,4(f10.4,1x),/                  &       
               &,4(f10.4,1x),/,2(f10.4,1x),/,/  &       
               &,5(g10.4,1x),/,5(g10.4,1x),/,/  &
               &,2(g10.4,1x),/,2(g10.4,1x),/,/  &
               &,3(g10.4,1x),/,3(g10.4,1x),/,/  &
               &,3(g10.4,1x),/,3(g10.4,1x),/,/  &
               &,2(g10.4,1x),/,1(g10.4,11x)     &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::           & 
         form3="(a2,1x,1(i2,3x),2(f10.4,1x),/  &
               &,3(f10.4,1x),/,2(f10.4,1x),/,/ &       
               &,5(g10.4,1x),/,5(g10.4,1x),/,/ &
               &,3(g10.4,1x),/,3(g10.4,1x),/,/ &
               &,3(g10.4,1x),/,3(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::           & 
         form4="(a2,1x,2(f10.4,1x),/           &
               &,2(f10.4,1x),/,/               &       
               &,2(g10.4,1x),/,2(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::           & 
         form5="(a2,1x,5(f10.4,1x),/           &
               &,5(g10.4,1x),/,5(g10.4,1x),/,/ &
               &,1(g10.4,1x),/,1(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::           & 
         form6="(a2,1x,4(f10.4,1x),/           &
               &,5(g10.4,1x),/,5(g10.4,1x),/,/ &
               &,1(g10.4,1x),/,1(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::           & 
         form7="(a2,1x,3(f10.4,1x),/           &
               &,5(g10.4,1x),/,5(g10.4,1x),/,/ &
               &,3(g10.4,1x),/,3(g10.4,1x),/,/ &
               &,4(g10.4,1x),/,4(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::           & 
         form8="(a2,1x,2(f10.4,1x),/           &
               &,5(f10.4,1x),/,/               &       
               &,2(g10.4,1x),/,2(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::           & 
         form9="(a2,1x,4(f10.4,1x),/           &
               &,4(g10.4,1x),/,4(g10.4,1x),/,/ &
               &,5(g10.4,1x),/,5(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"
    character(len=500), parameter ::           & 
        form10="(a2,1x,4(f10.4,1x),/           &
               &,4(g10.4,1x),/,4(g10.4,1x),/,/ &
               &,4(g10.4,1x),/,4(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,2(g10.4,1x),/,/ &
               &,2(g10.4,1x),/,1(g10.4,11x)    &
               &,2(f9.6,1x),/,/)"

    select case(final_state)
    case ('ng','ns')
       write(ndat2, fmt=form1) final_state                                  &
                              ,ipart1,ipart2,run0,run1,mass_n,ms1_print,mg  &
                              ,fct(0),fct(1),fct(21),fct(2),fct(3)          &
                              ,rel(0),rel(1),rel(21),rel(2),rel(3)          &
                              ,fct(4),fct(5),fct(6),fct(7)                  &
                              ,rel(4),rel(5),rel(6),rel(7)                  &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    case ('nn')
       write(ndat2, fmt=form2) final_state                                  &
                              ,ipart1,ipart2,run0,run1,mass_n,ms,mg         &
                              ,fct(0),fct(1),fct(21),fct(2),fct(3)          &
                              ,rel(0),rel(1),rel(21),rel(2),rel(3)          &
                              ,fct(4),fct(5)                                &
                              ,rel(4),rel(5)                                &
                              ,fct(6),fct(7),fct(8)                         &
                              ,rel(6),rel(7),rel(8)                         &
                              ,fct(9),fct(10),fct(11)                       &
                              ,rel(9),rel(10),rel(11)                       &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    case ('ll')
       write(ndat2, fmt=form3) final_state                                  &
                              ,ipart1,run0,run1,mass_s(1:3),ms,mg           &
                              ,fct(0),fct(1),fct(21),fct(2),fct(3)          &
                              ,rel(0),rel(1),rel(21),rel(2),rel(3)          &
                              ,fct(4),fct(5),fct(6)                         &
                              ,rel(4),rel(5),rel(6)                         &
                              ,fct(7),fct(8),fct(9)                         &
                              ,rel(7),rel(8),rel(9)                         & 
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    case ('gg')
       write(ndat2, fmt=form4) final_state                                  &
                              ,run0,run1,ms,mg                              &
                              ,fct(0),fct(21)                               &
                              ,rel(0),rel(21)                               &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng

    case ('sg')
       write(ndat2, fmt=form4) final_state                                  &
                              ,run0,run1,ms1_print,mg                       &
                              ,fct(0),fct(21)                               &
                              ,rel(0),rel(21)                               &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng

    case ('sb','ss')
       write(ndat2, fmt=form4) final_state                                  &
                              ,run0,run1,ms1_print,ms2_print                &
                              ,fct(0),fct(21)                               &
                              ,rel(0),rel(21)                               &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng

    case ('tb','bb')
       write(ndat2, fmt=form5) final_state                                  &
                              ,run0,run1,mass_s(1:3)                        &
                              ,fct(0),fct(1),fct(2),fct(21),fct(3)          &
                              ,rel(0),rel(1),rel(2),rel(21),rel(3)          &
                              ,fct(4)                                       &
                              ,rel(4)                                       &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    case ('lq')
       write(ndat2, fmt=form6) final_state                                  &
                              ,run0,run1,mass_s(1:2)                        &
                              ,fct(0),fct(1),fct(2),fct(21),fct(3)          &
                              ,rel(0),rel(1),rel(2),rel(21),rel(3)          & 
                              ,fct(4)                                       &
                              ,rel(4)                                       &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    case ('le')
       write(ndat2, fmt=form7) final_state                                  &
                              ,run0,run1,mass_s(1)                          &
                              ,fct(0),fct(1),fct(2),fct(21),fct(3)          &
                              ,rel(0),rel(1),rel(2),rel(21),rel(3)          &
                              ,fct(4),fct(5),fct(6)                         &
                              ,rel(4),rel(5),rel(6)                         &
                              ,fct(7),fct(8),fct(9),fct(10)                 &
                              ,rel(7),rel(8),rel(9),rel(10)                 &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    case ('xx')
       write(ndat2, fmt=form8) final_state                                  &
                              ,run0,run1,mass_s(1:3),ms,mg                  &
                              ,fct(0),fct(21)                               &
                              ,rel(0),rel(21)                               &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    case ('hh')
       write(ndat2, fmt=form9) final_state                                  &
                              ,run0,run1,mch,tan_b                          &
                              ,fct(0),fct(1),fct(2),fct(21)                 &
                              ,rel(0),rel(1),rel(2),rel(21)                 &
                              ,fct(4),fct(5),fct(6),fct(7),fct(8)           &
                              ,rel(4),rel(5),rel(6),rel(7),rel(8)           &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    case ('ht')
       write(ndat2, fmt=form10)final_state                                  &
                              ,run0,run1,mch,tan_b                          &
                              ,fct(0),fct(1),fct(2),fct(21)                 &
                              ,rel(0),rel(1),rel(2),rel(21)                 &
                              ,fct(4),fct(5),fct(6),fct(7)                  &
                              ,rel(4),rel(5),rel(6),rel(7)                  &
                              ,fct(8),fct(9)                                &
                              ,rel(8),rel(9)                                &
                              ,fct(20),fct(22)                              &
                              ,rel(20),kfac,kng
    end select

  end subroutine DAT2

! ------------------------------
  subroutine DAT3(ndat3,lowmass,unimass)               ! the complete output for documentation 

    integer,                            intent(in)   :: ndat3
    real(kind=double), dimension(0:99), intent(in)   :: lowmass
    real(kind=double), dimension(1:20), intent(in)   :: unimass

    write(ndat3, fmt="(a50)"           )  " ======================================================= "
    write(ndat3, fmt="(a15,a2)"        )  " final state = ",final_state
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a15,f10.2)"     )  " cm energy   = ",sqrt(sc)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,f12.5)"     )  " tan_beta = ",unimass(10)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,f12.5)"     )  " m_gluino = ",lowmass(4)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_neut1  = ",lowmass(5)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_neut2  = ",lowmass(6)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_neut3  = ",lowmass(7)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_neut4  = ",lowmass(8)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_char1  = ",lowmass(9)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_char2  = ",lowmass(10)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sdl    = ",lowmass(11)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sdr    = ",lowmass(12)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sul    = ",lowmass(13)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sur    = ",lowmass(14)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sq8    = ",lowmass(15)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sq10   = ",lowmass(16)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sb1    = ",lowmass(17)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sb2    = ",lowmass(18)
    write(ndat3, fmt="(a12,f12.5)"     )  " A_b      = ",lowmass(21)
    write(ndat3, fmt="(a12,f12.5)"     )  " sin(2b)  = ",lowmass(22)
    write(ndat3, fmt="(a12,f12.5)"     )  " cos(2b)  = ",lowmass(23)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_st1    = ",lowmass(19)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_st2    = ",lowmass(20)
    write(ndat3, fmt="(a12,f12.5)"     )  " A_t      = ",lowmass(24)
    write(ndat3, fmt="(a12,f12.5)"     )  " sin(2t)  = ",lowmass(25)
    write(ndat3, fmt="(a12,f12.5)"     )  " cos(2t)  = ",lowmass(26)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sel    = ",lowmass(30)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_ser    = ",lowmass(31)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sne    = ",lowmass(32)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_stau1  = ",lowmass(33)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_stau2  = ",lowmass(34)
    write(ndat3, fmt="(a12,f12.5)"     )  " m_sntau  = ",lowmass(35)
    write(ndat3, fmt="(a12,f12.5)"     )  " A_l      = ",lowmass(36)
    write(ndat3, fmt="(a12,f12.5)"     )  " sin(2l)  = ",lowmass(37)
    write(ndat3, fmt="(a12,f12.5)"     )  " cos(2l)  = ",lowmass(38)
    write(ndat3, fmt="(a12)     "      ) " ---------- "
    write(ndat3, fmt="(a12,f12.5)"     ) " m_ha     = ",lowmass(40)
    write(ndat3, fmt="(a12,f12.5)"     ) " m_hl     = ",lowmass(41)
    write(ndat3, fmt="(a12,f12.5)"     ) " m_hh     = ",lowmass(42)
    write(ndat3, fmt="(a12,f12.5)"     ) " m_hc     = ",lowmass(43)
    write(ndat3, fmt="(a12,f12.5)"     ) " sin(alp) = ",lowmass(44)
    write(ndat3, fmt="(a12,f12.5)"     ) " cos(alp) = ",lowmass(45)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,8(f9.6,2x))")  " Z_bw     = ",real( zz(1,1:4) ),aimag( zz(1,1:4) )
    write(ndat3, fmt="(a12,8(f9.6,2x))")  "            ",real( zz(2,1:4) ),aimag( zz(2,1:4) )
    write(ndat3, fmt="(a12,8(f9.6,2x))")  "            ",real( zz(3,1:4) ),aimag( zz(3,1:4) )
    write(ndat3, fmt="(a12,8(f9.6,2x))")  "            ",real( zz(4,1:4) ),aimag( zz(4,1:4) )
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,2(f9.6,2x))")  " U        = ",uu(1,1:2)
    write(ndat3, fmt="(a12,2(f9.6,2x))")  "            ",uu(2,1:2)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,2(f9.6,2x))")  " V        = ",vv(1,1:2)
    write(ndat3, fmt="(a12,2(f9.6,2x))")  "            ",vv(2,1:2)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,2(f9.6,2x))")  " F (st)   = ",mst(1,1:2)
    write(ndat3, fmt="(a12,2(f9.6,2x))")  "            ",mst(2,1:2)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,2(f9.6,2x))")  " F (sb)   = ",msb(1,1:2)
    write(ndat3, fmt="(a12,2(f9.6,2x))")  "            ",msb(2,1:2)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a12,2(f9.6,2x))")  " F (sl)   = ",msl(1,1:2)
    write(ndat3, fmt="(a12,2(f9.6,2x))")  "            ",msl(2,1:2)
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a18)"           )  " flags : "
    write(ndat3, fmt="(a12,i2)"        )  " isq_ng   = ",isq_ng
    write(ndat3, fmt="(a12,i2)"        )  " ipart1   = ",ipart1
    write(ndat3, fmt="(a12,i2)"        )  " ipart2   = ",ipart2
    write(ndat3, fmt="(a12,i2)"        )  " isquark1 = ",isquark1
    write(ndat3, fmt="(a12,i2)"        )  " isquark2 = ",isquark2
    write(ndat3, fmt="(a12,i2)"        )  " icoll    = ",icoll
    write(ndat3, fmt="(a12,i2)"        )  " isca     = ",isca
    write(ndat3, fmt="(a12,i2)"        )  " imx      = ",imx
    write(ndat3, fmt="(a12)"           )  " ---------- "
    write(ndat3, fmt="(a18)"           )  " numerical cuts : "
    write(ndat3, fmt="(a17,f14.8)"     )  " all integrals : ",cut
    write(ndat3, fmt="(a17,f14.8)"     )  " slicing       : ",eps_sli
    write(ndat3, fmt="(a17,f14.8)"     )  " subtraction   : ",eps_sub
    write(ndat3, fmt="(a17,f14.8)"     )  " small width   : ",ewi


  end subroutine DAT3

end module xx_in_out



