! ===========================================================================================================
module xx_kinds      ! stolen from thorsten's vamp_bundle.f90
  implicit none
  integer, parameter, public :: &
       double = selected_real_kind (precision (1.0) + 1, range (1.0) + 1)
  character(len=*), public, parameter :: KINDS_RCS_ID = &
       "$Id: kinds.nw,v 1.17 1998/05/07 10:43:48 ohl Exp $"
end module xx_kinds


! ===========================================================================================================
module xx_pass_integ     ! transfer variable values past vegas 
  use xx_kinds 
  implicit none

  integer, dimension(-1:30):: dim                        ! set in INIT_GLOBAL : dimensionality of integrals
  integer, dimension(-2:30):: ii_done                    ! set in IFCT_**_X12 : keep track of partonic channels 
  integer                  :: ipart1, ipart2             ! set in INIT_SUSY   : internal variable for final state particles 
  integer                  :: icoll                      ! set in PROSPINO    : Tevatron vs. LHC
  integer                  :: isq_ng                     ! set in PROSPINO    : use mass degenerate squarks
  integer                  :: i_ngtest                   ! set in PROSPINO    : go into testing mode (not recommended)
  integer                  :: isquark1,isquark2          ! set in PROSPINO    : squark flavor kept for LO results
  integer                  :: ii                         ! set in PROSPINO    : LO-NLO partonic channels
  integer                  :: idub                       ! set in INIT_GLOBAL : doublet structure of incoming quarks 
  integer                  :: isca                       ! set in INIT_GLOBAL : scale choice 
  integer                  :: imx                        ! set in INIT_GLOBAL : complex neutralino mixing vs. negative mass
  integer                  :: iscaling                   ! set in INIT_GLOBAL : scaling functions option [only for test use]
  integer                  :: ifla_le                    ! set in PROSPINO    : incoming flavor for LQ+lepton production
  integer                  :: n_faulty                   ! set in PROSPINO    : counts the number of non-valid calls

  character(len=2)         :: final_state                ! set in PROSPINO

  real(kind=double)        :: cut,eps_sli,ewi,eps_sub    ! set in INIT_GLOBAL

  real(kind=double)                  :: sc               ! set in PROSPINO
  real(kind=double)                  :: mg,ms            ! set in INIT_SUSY
  real(kind=double)                  :: mh1,mh2,mch      ! set in INIT_SUSY
  real(kind=double)                  :: sin_a,cos_a      ! set in INIT_SUSY
  real(kind=double)                  :: mu_susy,tan_b    ! set in INIT_SUSY
  real(kind=double)                  :: a_b,a_t          ! set in INIT_SUSY
  real(kind=double), dimension(-6:6) :: msq              ! set in INIT_SUSY
  real(kind=double), dimension(1:8)  :: smass_n,mass_n   ! set in INIT_SUSY
  real(kind=double), dimension(2,2)  :: uu,vv            ! set in INIT_SUSY
  real(kind=double), dimension(2,2)  :: mst,msb,msl      ! set in INIT_SUSY
  real(kind=double), dimension(4,4)  :: bw               ! set in INIT_SUSY
  real(kind=double), dimension(1:4)  :: mass_s           ! set in INIT_SUSY
  real(kind=double), dimension(1:4)  :: mass_x           ! set in INIT_SUSY
  real(kind=double)                  :: scafac           ! set in PROSPINO 
  real(kind=double)                  :: eta              ! set in INIT_GLOBAL

  complex(kind=double), dimension(4,4) :: zz             ! set in INIT_SUSY
end module xx_pass_integ


! ===========================================================================================================
module xx_public_variables   ! all kind of public sm and mssm parameters, used everywhere   
   use xx_kinds
   implicit none 

   real(kind=double), parameter, public :: gf=1.1663e-5                 ! set all kinds of
   real(kind=double), parameter, public :: pi=3.1415926                 ! public parameters
   real(kind=double), parameter, public :: gevpb=389379660.0
   real(kind=double), parameter, public :: mz=91.187, mw=80.41
   real(kind=double), parameter, public :: mb=4.6,    mt=172.0

   logical, parameter, public           :: ldebug=.false.

end module xx_public_variables



