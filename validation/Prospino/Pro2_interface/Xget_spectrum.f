cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c GET_SPECTRUM(nin2,unimass,lowmass,bw,uu,vv,mst,msb,msl)              c
c             interface calling the SUSY spectrum                      c
c                                                                      c
c   nin2        unit to read from in Les Houches format                c
c                                                                      c
c   output  unimass(10)   listed in Xinitialize.f90                    c 
c           lowmass(0:99) listed in Xinitialize.f90                    c
c           bw            neutralino mixing matrix bino-wino basis     c
c           uu            chargino mixing matrix                       c
c           vv            chargino mixing matrix                       c  
c           mst,msb,msl   stop/sbottom/slepton mixing matrices         c
c                                                                      c
c n.b. make sure INIT_ALPHAS is run first from Xinitialize.f90         c   
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine GET_SPECTRUM(nin2,unimass,lowmass,bw,uu,vv,mst,msb,msl)

      implicit none

      integer nin2,i1,i2
      real*8  unimass(1:20),lowmass(0:99),width(0:99)
      real*8  bw(4,4),pz(4,4),uu(2,2),vv(2,2)
      real*8  msb(2,2),mst(2,2),msl(2,2)
      real*8  mw_dum,mz_dum
      
      do i1=0,99
         lowmass(i1) = 0.D0
      end do

      do i1=1,4
         do i2=1,4
            bw(i1,i2) = 0.D0
            pz(i1,i2) = 0.D0
         end do
      end do

      do i1=1,2
         do i2=1,2
            uu(i1,i2)  = 0.D0
            vv(i1,i2)  = 0.D0
            mst(i1,i2) = 0.D0
            msb(i1,i2) = 0.D0
            msl(i1,i2) = 0.D0
         end do
      end do

      call READ_SLHA_RIP(nin2,unimass,lowmass,width,
     &                      bw,uu,vv,mst,msb,msl,mz_dum,mw_dum)

      return
      end





