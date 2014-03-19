cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c    approximative mSUGRA subroutine a la Drees & Martin               c
c                                                                      c
c                                                                      c
c                   SUGRA                                              c
c   |      |     |         |          |           |        |           c
c XITLA ALSINI RUNM_ORIG B_SPECIAL GAUGINOS DIAGONALIZE SUBH(BORNH)    c
c                                            |          |   |     |    c
c                                    RUNM_ORIG   RUNM_ORIG GFUN ALPHAS c
c                                                          |  |        c
c                                                  RUNM_ORIG  ALPHAS   c
c                                                                      c
c                                                                      c
c    subroutines used : ALSINI in Xalphas_smooth.f                     c
c    [externally]       ALPHAS in Xalphas_smooth.f [initialized]       c
c                       RUNM_ORIG   in Xalphas_smooth.f [initialized]  c
c                                                                      c
c    literature used  : Drees & Martin [madph-95-879]                  c
c                       Collins, Martin, Squires                       c
c                       Kileng's thesis                                c
c                       my thesis                                      c
c                       Spira's program                                c
c                       Carena, Quiros, Wagner                         c
c                                                                      c
c    n.b. :                                                            c
c    gauge coupling unification {12,13} via variable nuni={12,13}      c
c    negative mass squared lead to default value for all masses        c
c    ms,mc,mb,mt,mtau,mz,mw,gf,...  fixed inside the program           c
c                                                                      c
c    23-11-98 : leading order higgs sector included -> BORNH           c
c    08-12-98 : only higgs sector -> HIGGS_MSSM                        c
c    16-12-98 : remove common block MSSM_SUGRA                         c
c    16-03-99 : set on-shell sw for low energy masses                  c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c    input parameters:  ncomphep={0,1}  for CompHEP's weak parameters  c
c                       nlo_h={0,1} for NLO Higgs sector               c
c                                                                      c
c                       unimass(1)  universal scalar mass m_0(1,2)     c
c                       unimass(2)  additional mass term for m_0(1,2)  c
c                       unimass(3)  additional mass term for m_0(3)    c
c                       unimass(4)  additional mass term for m_0(Higgs)c
c                       unimass(5)  universal fermion mass m_1/2       c
c                       unimass(6)  additional mass term for m1        c
c                       unimass(7)  additional mass term for m2        c
c                       unimass(8)  additional mass term for m3        c
c                       unimass(9)  universal mixing parameter a_0     c
c                       unimass(10) tan(beta)                          c
c                       unimass(11) sign of higgs mass parameter mu    c
c                       .............                                  c
c                       unimass(20)                                    c
c                                                                      c
c    output parameters: lowmass(0)  mu                                 c
c                       lowmass(1)  m_1                                c
c                       lowmass(2)  m_2                                c
c                       lowmass(3)  m_3                                c
c                                                                      c
c                       lowmass(4)  gluino mass                        c
c                       lowmass(5)  \                                  c
c                       lowmass(6)   \                                 c
c                       lowmass(7)   /  neutralino masses [with sign]  c
c                       lowmass(8)  /                                  c
c                       lowmass(9)  \                                  c
c                       lowmass(10) /   chargino masses                c
c                                                                      c
c                       lowmass(11) sdown_l mass                       c
c                       lowmass(12) sdown_r mass                       c
c                       lowmass(13) sup_l mass                         c
c                       lowmass(14) sup_r mass                         c
c                       lowmass(15) degenerate squark mass (8)         c
c                       lowmass(16) degenerate squark mass (10)        c
c                       lowmass(17) sbottom-1 mass                     c
c                       lowmass(18) sbottom-2 mass                     c
c                       lowmass(19) stop-1 mass                        c
c                       lowmass(20) stop-2 mass                        c
c                                                                      c
c                       lowmass(21) a_b                                c
c                       lowmass(22) sin(2 theta_b)                     c
c                       lowmass(23) cos(2 theta_b)                     c
c                       lowmass(24) a_t                                c
c                       lowmass(25) sin(2 theta_t)                     c
c                       lowmass(26) cos(2 theta_t)                     c
c                                                                      c
c                       lowmass(30) selectron_l mass                   c
c                       lowmass(31) selectron_r mass                   c
c                       lowmass(32) selectron-neutrino mass            c
c                       lowmass(33) stau_1 mass                        c
c                       lowmass(34) stau_2 mass                        c
c                       lowmass(35) stau-neutrino mass                 c
c                       lowmass(36) a_tau                              c
c                       lowmass(37) sin(2 theta_tau)                   c
c                       lowmass(38) cos(2 theta_tau)                   c
c                                                                      c
c                       lowmass(40) cp odd higgs mass                  c
c                       lowmass(41) light cp even higgs mass           c
c                       lowmass(42) heavy cp even higgs mass           c
c                       lowmass(43) charged higgs mass                 c
c                       lowmass(44) sin(alpha)                         c
c                       lowmass(45) cos(alpha)                         c
c                                                                      c
c                       like cteq: u,d,s,c,b,t first L then R          c
c                       lowmass(51) sup_L                              c
c                       lowmass(52) sdown_L                            c
c                       lowmass(53) sstrange_L                         c
c                       lowmass(54) scharm_L                           c
c                       lowmass(55) sbottom_1                          c
c                       lowmass(56) stop_1                             c
c                       lowmass(57) sup_R                              c
c                       lowmass(58) sdown_R                            c
c                       lowmass(59) sstrange_R                         c
c                       lowmass(60) scharm_R                           c
c                       lowmass(61) sbottom_2                          c
c                       lowmass(62) stop_2                             c
c                                                                      c
c                       lowmass(80) unification scale                  c
c                       lowmass(81) unified coupling alpha(m_x)        c
c                                                                      c
c                       lowmass(91) trilinear higgs coupling lambda(1) c
c                       .......                                        c
c                       lowmass(97) lambda(7)                          c
c                                                                      c
c                       lowmass(99)                                    c
c                                                                      c
c                       bwmix neutralino mixing matrix (bino-wino)     c
c                       pzmix neutralino mixing matrix (photino-zino)  c
c                       uumix chargino mixing matrix u                 c
c                       vvmix chargino mixing matrix v                 c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine SUGRA(ncomphep,nlo_h,unimass,lowmass,bw,pz,uu,vv)

      implicit none 

      integer n,ng,nh,nloop,nq,nuni,ncomphep,nlo_h
      real*8  unimass(20),lowmass(0:99),pz(4,4),bw(4,4),uu(2,2),vv(2,2)
     &       ,XITLA,RUNM_ORIG,B_SPECIAL,lam_qcd,lambda(1:7)
     &       ,m0(3),m12,a0,sgnmu,ms,mc,gmas(3)
     &       ,alsmz,cb,cbb,sb,acc,sw02,sw2,al1,al2,al3
     &       ,b(3),mx,alpx,cof(3),cofq(2),alsq,q
     &       ,qq(2),dd,cgl(2),mer,mel,mnl,mur,mdr,mul,mdl,rmt,rmb
     &       ,xtau,xt,xb,m2er,m2el,m2nl,m2ur,m2dr,m2ul,m2dl
     &       ,mer3,mel3,mnl3,mql3,mur3,mdr3,at,ab,al,mg
     &       ,m2er3,m2el3,m2nl3,m2sn3,m2ql3,m2ur3,m2dr3
     &       ,m2stl,m2str,mlrt,trast,detst,m2st1,m2st2,mst1,mst2,s2t,c2t
     &       ,m2sbl,m2sbr,mlrb,trasb,detsb,m2sb1,m2sb2,msb1,msb2,s2b,c2b
     &       ,m2sel,m2ser,mlre,trase,detse,m2se1,m2se2,mse1,mse2,s2e,c2e
     &       ,msq8,msq10,m1,m2,m3,eps,mh12,mh22,mu,mu2,ma,m2a
     &       ,mneut(4),mchar(2)
     &       ,mhl,mhh,mhc,tanba,mchi,sa,ca
     &       ,pi,mz,mw,sw,swos,gf,alpmz,mb,mt,mtau,tb,default

      external XITLA,RUNM_ORIG,B_SPECIAL

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               all kind of standard model parameters                  c

c               set nuni={12,13} for 1-2 of 1-3 gauge coupling unif.
      nuni = 12

      pi   = 4.D0*atan(1.D0)

      mz   = 91.187D0
      mw   = 80.41D0

      mtau = 1.7771D0
      
      ms = 0.3D0
      mc = 1.5D0
      mb = 5.D0
      mt = 175.D0

c               check the input
      if ((ncomphep.ne.0).and.(ncomphep.ne.1)) then
         print*, ' Xsugra.f: wrong value for ncomphep : ',ncomphep
         stop
      end if 

      if ((nlo_h.ne.0).and.(nlo_h.ne.1)) then
         print*, ' Xsugra.f: wrong value for nlo_h : ',nlo_h
         stop
      end if 

c               weak coupling sector given for 1-3-unification 
      if (nuni.eq.13) then      
         alsmz = 0.118D0
c               everything the on-shell and gf scheme  
         gf    = 1.16639D-5
         sw    = sqrt( 1.D0 - mw**2/mz**2 )
         swos  = sw
         alpmz = gf    * sqrt(2.D0) * sw**2 * mw**2 / pi
      elseif (nuni.eq.12) then 
c               everything in the msbar scheme,
c               except for low energy masses 
         alpmz = 1.D0/127.9D0
         if (ncomphep.eq.0) then 
            sw    = sqrt( 0.2315D0 ) 
            swos  = sqrt( 1.D0 - mw**2/mz**2 )
         else if (ncomphep.eq.1) then
            sw    = 0.4740D0
            swos  = sw
            print*,' Xsugra.f: sw set a la CompHEP '
         end if
         gf    = alpmz / sqrt(2.D0) / sw**2 / mw**2 * pi 
      else 
         print*, ' Xsugra.f: wrong gauge coupling unification ? '
         stop
      end if 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c               translate the input 
      m0(1)   = unimass(1) + unimass(2)
      m0(2)   = unimass(1) + unimass(3)
      m0(3)   = unimass(1) + unimass(4)
      m12     = unimass(5)
      gmas(1) = unimass(6)
      gmas(2) = unimass(7)
      gmas(3) = unimass(8)
      a0      = unimass(9)
      tb      = unimass(10)

c               set the sign of mu 
      if (unimass(11).gt.0.D0) then 
         sgnmu = 1.D0
      else if (unimass(11).lt.0.D0) then 
         sgnmu = -1.D0 
      else 
         print*,' Xsugra.f: wrong unimass(11) : ',unimass(11)
      end if 

c               the accuracy for iterations 
      acc = 1.d-10

c               default mass value for negative mass squared 
      default = 0.D0

c               the beta functions [cms 10.113]
      ng = 3
      nh = 2 
      b(1) = -2.D0 * dble(ng) - 3.D0/10.D0 * dble(nh) 
      b(2) = 6.D0 - 2.D0 * dble(ng) - 1.D0/2.D0 * dble(nh) 
      b(3) = 9.D0 - 2.D0 * dble(ng) 

c               the gut-scale weak mixing angle [cms 7.12]
      sw02 = 3.D0/8.D0

c               the low-scale gauge couplings [cms 7.17]
      al1 = alpmz * (1.D0-sw02)/sw02 /(1.D0-sw**2)

      if (nuni.eq.13) then 
         al3 = alsmz
c               the gut scale [cms 7.14]
         mx = mz * dexp(2.D0*pi/(b(3)-b(1))*(1.D0/al1-1.D0/al3))
c               the unified gauge coupling alpha(mx) [cms 7.14] 
         alpx = al3/(1.D0 + al3*b(3)/(4.D0*pi) * log(mx**2/mz**2))
      else if (nuni.eq.12) then 
         al2 = alpmz/sw**2 
c               the gut scale [cms 7.14]
         mx = mz * dexp(2.D0*pi/(b(2)-b(1))*(1.D0/al1-1.D0/al2))
c               the unified gauge coupling alpha(mx) [cms 7.14] 
         alpx = al2/(1.D0 + al2*b(2)/(4.D0*pi) * log(mx**2/mz**2))
      end if 

c               the different low-scale gauge couplings 
      do n=1,3 
         cof(n) = 1.D0/(1.D0 + alpx*b(n)/(4.D0*pi) * log(mz**2/mx**2))
      end do

c               the alpha_s value for qcd type different masses 
      alsq = cof(3) * alpx

c               the improved lambda_qcd value at q=mz for ALPHAS routine
      nloop   = 2
      q       = mz
      lam_qcd = XITLA(nloop,alsq,q,acc)

c               initialize alpha_s
      if (q.lt.mt) then
         nq = 5
      else
         nq = 6
      end if
      call ALSINI(acc,lam_qcd,mc,mb,mt,nq)

c               the gluino loop parameter [dm 17] for light generations 
      qq(1) = sqrt( m0(1)**2 + 6*m12**2 )
      cofq(1) = 1.D0/(1.D0 + alpx*b(3)/(4.D0*pi) * log(qq(1)**2/mx**2))
      cgl(1)  = 8.D0/9.D0 * (cofq(1)**2-1.D0)
c               same for heavy generations
      qq(2) = sqrt( m0(2)**2 + 6*m12**2 )
      cofq(2) = 1.D0/(1.D0 + alpx*b(3)/(4.D0*pi) * log(qq(2)**2/mx**2))
      cgl(2)  = 8.D0/9.D0 * (cofq(2)**2-1.D0)

c               all the non-mixing scalar masses [dm 16]
      sw2  = swos**2
      cb  = 1.D0 / sqrt(1.D0+tb**2)
      cbb = 2.D0*cb**2 - 1.D0
      sb  = tb * cb

      dd = mz**2 * cbb
      m2er = m0(1)**2 + 0.15D0         * m12**2 -        sw2    *dd
      m2el = m0(1)**2 + 0.52D0         * m12**2 - (0.5D0-sw2)   *dd
      m2nl = m0(1)**2 + 0.52D0         * m12**2 +  0.5D0        *dd
      m2ur = m0(1)**2 + (0.07D0+cgl(1))* m12**2 + 2.D0/3.D0*sw2 *dd
      m2dr = m0(1)**2 + (0.02D0+cgl(1))* m12**2 - 1.D0/3.D0*sw2 *dd
      m2ul = m0(1)**2 + (0.47D0+cgl(1))* m12**2
     &                               + (1.D0/2.D0-2.D0/3.D0*sw2)*dd
      m2dl = m0(1)**2 + (0.47D0+cgl(1))* m12**2
     &                               - (1.D0/2.D0-1.D0/3.D0*sw2)*dd

      if (m2er.gt.0.D0) then
         mer = sqrt(m2er)
      else 
         print*,' Xsugra.f: mer^2 negative '
         mer = default
      end if 

      if (m2el.gt.0.D0) then
         mel = sqrt(m2el)
      else 
         print*,' Xsugra.f: mel^2 negative '
         mel = default
      end if 

      if (m2nl.gt.0.D0) then
         mnl = sqrt(m2nl)
      else 
         print*,' Xsugra.f: mnl^2 negative '
         mnl = default
      end if 

      if (m2ur.gt.0.D0) then
         mur = sqrt(m2ur)
      else 
         print*,' Xsugra.f: mur^2 negative '
         mur = default
      end if 

      if (m2dr.gt.0.D0) then
         mdr = sqrt(m2dr)
      else 
         print*,' Xsugra.f: mdr^2 negative '
         mdr = default
      end if 

      if (m2ul.gt.0.D0) then
         mul = sqrt(m2ul)
      else 
         print*,' Xsugra.f: mul^2 negative '
         mul = default
      end if 

      if (m2dl.gt.0.D0) then
         mdl = sqrt(m2dl)
      else 
         print*,' Xsugra.f: mdl^2 negative '
         mdl = default
      end if 

c               the running top and bottom mass
      rmt = RUNM_ORIG(mt,6) * (1.D0 + 8.D0/3.D0 * alpmz**2)
      rmb = RUNM_ORIG(mt,5) * (1.D0 + 8.D0/3.D0 * alpmz**2)

c               approximate yukawa couplings [dm 25]
      xt = (rmt/150.D0/sb)**2 * 
     &     (  0.9D0*m0(2)**2 + 2.1D0*m12**2 
     &      + (1.D0-(rmt/190.D0/sb)**3) * (0.24D0*a0**2 + a0*m12) )
      xb = (rmb/150.D0/cb)**2 * 
     &     (  0.9D0*m0(2)**2 + 2.1D0*m12**2
     .      + (1.D0-(rmb/190.D0/cb)**3) * (0.24D0*a0**2 + a0*m12) )
      xtau = 1.d-4/cb**2 * ( m0(2)**2 + 0.15D0*m12**2 + 0.33D0*a0**2 )

c               approximate third generation scalar masses [dm 24]
      m2er3 = m0(2)**2 + 0.15D0          * m12**2 - 2.D0/3.D0 * xtau
      m2el3 = m0(2)**2 + 0.52D0          * m12**2 - 1.D0/3.D0 * xtau
      m2nl3 = m2el3 
      m2ql3 = m0(2)**2 + (0.47D0+cgl(2)) * m12**2 - 1.D0/3.D0 * (xb+xt)
      m2ur3 = m0(2)**2 + (0.07D0+cgl(2)) * m12**2 - 2.D0/3.D0 * xt
      m2dr3 = m0(2)**2 + (0.02D0+cgl(2)) * m12**2 - 2.D0/3.D0 * xb

      if (m2er3.gt.0.D0) then
         mer3 = sqrt(m2er3)
      else 
         print*,' Xsugra.f: mer3^2 negative '
         mer3 = default
      end if 

      if (m2el3.gt.0.D0) then
         mel3 = sqrt(m2el3)
      else 
         print*,' Xsugra.f: mel3^2 negative '
         mel3 = default
      end if 

      if (m2nl3.gt.0.D0) then
         mnl3 = sqrt(m2nl3)
      else 
         print*,' Xsugra.f: mnl3^2 negative '
         mnl3 = default
      end if 

      if (m2ql3.gt.0.D0) then
         mql3 = sqrt(m2ql3)
      else 
         print*,' Xsugra.f: mql3^2 negative '
         mql3 = default
      end if 

      if (m2ur3.gt.0.D0) then
         mur3 = sqrt(m2ur3)
      else 
         print*,' Xsugra.f: mur3^2 negative '
         mur3 = default
      end if 

      if (m2dr3.gt.0.D0) then
         mdr3 = sqrt(m2dr3)
      else 
         print*,' Xsugra.f: mdr3^2 negative '
         mdr3 = default
      end if 

c               approximate stop trilinear coupling [dm 26]
      at = (  a0 * (1.D0-(rmt/190.D0/sb)**2)
     &      + m12 * (3.47D0-1.9D0*(rmt/190.D0/sb)**2) )

c               no down-type trilinear coupling 
      ab = 0.D0
      al = 0.D0

c               approximate higgs doublet masses [dm 23]
      mh12 = m0(3)**2 + 0.52D0*m12**2 - xb - xtau/3.D0
      mh22 = m0(3)**2 + 0.52D0*m12**2 - xt

c               approximate higgsino mass parameter [dm 9a]
      mu2 = ( mh22*sb**2 - mh12*cb**2 )/cbb - mz**2/2.D0
      if (mu2.gt.0.D0) then 
         mu = sgnmu * sqrt(mu2)
      else 
         print*,' Xsugra.f: mu^2 negative '
         mu = default
      end if 

c               pseudoscalar higgs mass at tree level [dm 11]
      m2a = mh12 + mh22 + 2.D0*mu2
      if (m2a.gt.0.D0) then 
         ma = sqrt(m2a)
      else 
         print*,' Xsugra.f: ma^2 negative '
         ma = default
      end if 

c               stop mass matrix [dm 12a, thesis 1.8]
      m2stl = mql3**2 + rmt**2 + (0.5D0-2.D0/3.D0*sw2) *dd
      m2str = mur3**2 + rmt**2 +        2.D0/3.D0*sw2  *dd
      mlrt  = -rmt * (at+mu/tb)

      trast = m2stl + m2str 
      detst = m2stl*m2str - mlrt**2 

      m2st1 = ( trast - sqrt( trast**2 - 4.D0*detst ))/2.D0
      m2st2 = ( trast + sqrt( trast**2 - 4.D0*detst ))/2.D0
      s2t   = 2.D0*mlrt / sqrt( trast**2 - 4.D0*detst )
      c2t   = abs(m2stl-m2str) / sqrt( trast**2 - 4.D0*detst ) 
      s2t   = -sign(1.D0,m2stl-m2str) * s2t 

      if (m2st1.gt.0.D0) then 
         mst1 = sqrt(m2st1)
      else 
         print *, 'Xsugra.f: m_stop1^2 negative '
         mst1 = default
      end if 

      if (m2st2.gt.0.D0) then 
         mst2 = sqrt(m2st2)
      else 
         print *, 'Xsugra.f: m_stop2^2 negative '
         mst2 = default
      end if 

c               sbottom mass matrix [dto]
      m2sbl = mql3**2 + rmb**2 + (-0.5D0+1.D0/3.D0*sw2) *dd
      m2sbr = mdr3**2 + rmb**2 -         1.D0/3.D0*sw2  *dd
      mlrb  = -rmb * (ab+mu*tb)

      trasb = m2sbl + m2sbr 
      detsb = m2sbl*m2sbr - mlrb**2 

      m2sb1 = ( trasb - sqrt( trasb**2 - 4.D0*detsb ))/2.D0
      m2sb2 = ( trasb + sqrt( trasb**2 - 4.D0*detsb ))/2.D0
      s2b   = 2.D0*mlrb / sqrt( trasb**2 - 4.D0*detsb )
      c2b   = abs(m2sbl-m2sbr) / sqrt( trasb**2 - 4.D0*detsb ) 
      s2b   = -sign(1.D0,m2sbl-m2sbr) * s2b 

      if (m2sb1.gt.0.D0) then 
         msb1 = sqrt(m2sb1)
      else 
         print *, 'Xsugra.f: m_sbottom1^2 negative '
         msb1 = default
      end if 

      if (m2sb2.gt.0.D0) then 
         msb2 = sqrt(m2sb2)
      else 
         print *, 'Xsugra.f: m_sbottom2^2 negative '
         msb2 = default
      end if 

c               stau mass matrix [dm 12b]
      m2sel = mel3**2 + mtau**2 + (-0.5D0+sw2) *dd
      m2ser = mer3**2 + mtau**2 -         sw2  *dd
      m2sn3 = mnl3**2           +   0.5D0      *dd
      mlre  = -mtau * (al+mu*tb)

      if (m2sn3.gt.0.D0) then 
         mnl3 = sqrt(m2sn3)
      else 
         print *, 'Xsugra.f: m_snu_tau^2 negative '
         mnl3 = default
      end if 

      trase = m2sel + m2ser 
      detse = m2sel*m2ser - mlre**2 

      m2se1 = ( trase - sqrt( trase**2 - 4.D0*detse ))/2.D0
      m2se2 = ( trase + sqrt( trase**2 - 4.D0*detse ))/2.D0
      s2e   = 2.D0*mlre / sqrt( trase**2 - 4.D0*detse )
      c2e   = abs(m2sel-m2ser) / sqrt( trase**2 - 4.D0*detse ) 
      s2e   = -sign(1.D0,m2sel-m2ser) * s2e 

      if (m2se1.gt.0.D0) then 
         mse1 = sqrt(m2se1)
      else 
         print *, 'Xsugra.f: m_stau1^2 negative '
         mse1 = default
      end if 

      if (m2se2.gt.0.D0) then 
         mse2 = sqrt(m2se2)
      else 
         print *, 'Xsugra.f: m_stau2^2 negative '
         mse2 = default
      end if 

c               the mass degenerate squark masses 
      msq8  = ( mur + mdr + mul + mdl )/4.D0
      msq10 = ( 2.D0*(mur + mdr + mul + mdl) + msb1 + msb2 )/10.D0

c               the gaugino masses [dm 18]: non-degenerate 
      m1 = cof(1) * ( m12 + gmas(1) )
      m2 = cof(2) * ( m12 + gmas(2) )
      m3 = cof(3) * ( m12 + gmas(3) )

c               approximate gluino mass [dm 19] 
      eps = 1.d-2
      mg = m3 * ( 1.D0 + alsq/(4.D0*pi) * ( 15.D0 - 18.D0*log(m3/q)
     &           + B_SPECIAL(m3**2,mur,eps,q**2)
     &           + B_SPECIAL(m3**2,mul,eps,q**2)
     &           + B_SPECIAL(m3**2,mdr,eps,q**2)
     &           + B_SPECIAL(m3**2,mdl,eps,q**2)
     &           + B_SPECIAL(m3**2,mur,mc,q**2)
     &           + B_SPECIAL(m3**2,mul,mc,q**2)
     &           + B_SPECIAL(m3**2,mdr,ms,q**2)
     &           + B_SPECIAL(m3**2,mdl,ms,q**2)
     &           + B_SPECIAL(m3**2,mur3,mt,q**2)
     &           + B_SPECIAL(m3**2,mql3,mt,q**2)
     &           + B_SPECIAL(m3**2,mdr3,mb,q**2)
     &           + B_SPECIAL(m3**2,mql3,mb,q**2)
     &    ))

c               all the gaugino masses and mixing matrices 
      call GAUGINO(tb,mw,mz,swos,mu,m1,m2,mchar,mneut,bw,pz,uu,vv)

c               the higgs sector a la carena et al. or to leading order
      if (nlo_h.eq.0) then 
         call BORNH(ma,tb,mw,mz,mhl,mhh,mhc,sa,ca)
         do n=1,7,1
            lambda(n) = 0.D0
         end do
      else if (nlo_h.eq.1) then 
         mchi = mchar(2)
         call SUBH(ma,tb,mw,mz,swos,gf,mt
     &            ,mql3,mur3,mdr3,at,ab,mu,mchi,default
     &            ,mhl,mhh,mhc,sa,ca,tanba,lambda)
      end if 
c               check the diagonalization of the mass matrices 
      call DIAGONALIZE(tb,mw,mz,swos,mt,alpmz,mtau
     &                ,m1,m2,mu,mql3,mur3,mdr3,mel3,mer3
     &                ,at,ab,al,s2t,s2b,s2e,bw,uu,vv)

c               initialize the output array 
      do n=0,99
         lowmass(n) = 0.D0
      end do

c               translate the output 
      lowmass(0)  = mu  
      lowmass(1)  = m1
      lowmass(2)  = m2
      lowmass(3)  = m3

      lowmass(4)  = mg
      lowmass(5)  = mneut(1)
      lowmass(6)  = mneut(2)
      lowmass(7)  = mneut(3)
      lowmass(8)  = mneut(4)
      lowmass(9)  = mchar(1)
      lowmass(10) = mchar(2)

      lowmass(11) = mdl
      lowmass(12) = mdr
      lowmass(13) = mul
      lowmass(14) = mur
      lowmass(15) = msq8
      lowmass(16) = msq10
      lowmass(17) = msb1 
      lowmass(18) = msb2
      lowmass(19) = mst1
      lowmass(20) = mst2

      lowmass(21) = ab
      lowmass(22) = s2b
      lowmass(23) = c2b
      lowmass(24) = at
      lowmass(25) = s2t
      lowmass(26) = c2t
      
      lowmass(30) = mel
      lowmass(31) = mer
      lowmass(32) = mnl
      lowmass(33) = mse1
      lowmass(34) = mse2
      lowmass(35) = mnl3
      lowmass(36) = al
      lowmass(37) = s2e
      lowmass(38) = c2e
      
      lowmass(40) = ma
      lowmass(41) = mhl
      lowmass(42) = mhh
      lowmass(43) = mhc
      lowmass(44) = sa
      lowmass(45) = ca 

      lowmass(51) = mdl
      lowmass(52) = mul
      lowmass(53) = mul
      lowmass(54) = mdl
      lowmass(55) = msb1
      lowmass(56) = mst1
      lowmass(57) = mdr
      lowmass(58) = mur
      lowmass(59) = mur
      lowmass(60) = mdr
      lowmass(61) = msb2
      lowmass(62) = mst2

      lowmass(80) = mx
      lowmass(81) = alpx 

      lowmass(91) = lambda(1)
      lowmass(92) = lambda(2)
      lowmass(93) = lambda(3)
      lowmass(94) = lambda(4)
      lowmass(95) = lambda(5)
      lowmass(96) = lambda(6)
      lowmass(97) = lambda(7)

      return
      end

c ----------------------------------------------------------------------
c          routine calculating all neutralino/chargino
c           mass and mixing matrices; i.e. output bw,pz,uu,vv 
      subroutine GAUGINO(tb,mw,mz,sw,mu,m1,m2,mc,mn,bw,pz,uu,vv)
      
      implicit none 

      integer    i,j,k,l,i1,idummy,n,irem(2),iord(4)
      real*8     mu,m1,m2,mc(2),mn(4)
     &          ,bw(4,4),pz(4,4),uu(2,2),vv(2,2),cm(2,2),dm(2,2)
     &          ,uux(2,2),vvx(2,2)
     &          ,cw,cb,cbb,sb,sbb,eps,xc2,xc3,xc4,xs,xu,x0
     &          ,xmn(4),ymn(4),zx(4,4),yz(4,4)
     &          ,xx0,xx1,sx,tx,c2u,c2v,sigu,sigv,cu,su,cv,sv
     &          ,mz,mw,sw,tb
      complex*16 cxa,cxb,cxc,cxd,cx1,cx2,cx3

c               weak mixing angle from SUGRA routine 
      cw = sqrt(1.D0-sw**2)

c               susy mixing angle 
      cb  = 1.D0 / sqrt(1.D0+tb**2)
      cbb = 2.D0*cb**2 - 1.D0
      sb  = tb * cb
      sbb = 2.D0 * sb * cb 

c               shift into the complex plane 
      eps = -1.d-10

c               neutralino mass matrix similar to b.kileng {-p,q,r}
      xc2 = m1*m2 - mz**2 - mu**2 - 3.D0/8.D0*(m1+m2)**2
      xc3 = - 1.D0/8.D0*(m1+m2)**3
     &      + 1.D0/2.D0*(m1+m2)*(m1*m2-mz**2-mu**2)
     &      + (m1+m2)*mu**2
     &      + (m1*cw**2+m2*sw**2)*mz**2
     &      - mu*mz**2*sbb
      xc4 =   (m1*cw**2+m2*sw**2)*mu*mz**2*sbb
     &      - m1*m2*mu**2
     &      + 1.D0/4.D0*(m1+m2)*( (m1+m2)*mu**2
     &                           +(m1*cw**2+m2*sw**2)*mz**2
     &                           -mu*mz**2*sbb )
     &      + 1.D0/16.D0*(m1+m2)**2*(m1*m2-mz**2-mu**2)
     &      - 3.D0/256.D0*(m1+m2)**4

      xs = - xc3**2 - 2.D0/27.D0*xc2**3 + 8.D0/3.D0*xc2*xc4
      xu = - 1.D0/3.D0*xc2**2 - 4.D0*xc4
      cxd = ( -4.D0*xu**3 - 27.D0*xs**2 ) * dcmplx(1.D0,eps)
      cxc = 1.D0/2.D0 * ( -xs + dcmplx(0.D0,1.D0)*sqrt(cxd/27.D0) )
      cxa = real(cxc**(1.D0/3.D0)) *  dcmplx(1.D0,-eps)
      cxb = 8.D0*cxa - 8.D0/3.D0*xc2* dcmplx(1.D0,-eps)

      x0  = (m1+m2)/4.D0

      cx1 =  cxa/2.D0 - xc2/6.D0 * dcmplx(1.D0,-eps)
      cx2 = -cxa/2.D0 - xc2/3.D0 * dcmplx(1.D0,-eps)
      cx3 = xc3 * dcmplx(1.D0,-eps)/sqrt(cxb)

c               neutralino mass matrix eigenvalues [not ordered]
      xmn(1) = x0 - abs(sqrt(cx1)) + abs(sqrt(cx2+cx3))
      xmn(2) = x0 + abs(sqrt(cx1)) - abs(sqrt(cx2-cx3))
      xmn(3) = x0 - abs(sqrt(cx1)) - abs(sqrt(cx2+cx3))
      xmn(4) = x0 + abs(sqrt(cx1)) + abs(sqrt(cx2-cx3))

c               bw basis mixing matrix yz(4,4) [thesis A.8]
      do i=1,4
         mn(i)  = abs(xmn(i))
         ymn(i) = xmn(i)
         zx(i,2) = -cw/sw * (m1-xmn(i)) / (m2-xmn(i))
         zx(i,3) = ( mu * (m2-xmn(i)) * (m1-xmn(i))
     &              -mz**2*sb*cb * ( (m1-m2)*cw**2+m2-xmn(i)))
     &             /mz/(m2-xmn(i))/sw/(mu*cb+xmn(i)*sb)
         zx(i,4) = (-xmn(i) * (m2-xmn(i)) * (m1-xmn(i))
     &              -mz**2*cb*cb * ( (m1-m2)*cw**2+m2-xmn(i)))
     &             /mz/(m2-xmn(i))/sw/(mu*cb+xmn(i)*sb)
         zx(i,1) = 1.D0/sqrt(1.D0+zx(i,2)**2+zx(i,3)**2+zx(i,4)**2) 
         yz(i,1)=zx(i,1)
         yz(i,2)=zx(i,2)*zx(i,1)
         yz(i,3)=zx(i,3)*zx(i,1)
         yz(i,4)=zx(i,4)*zx(i,1)
      end do

c               ordering the neutralinos
      xx0 = dmin1(mn(1),mn(2),mn(3),mn(4))
      xx1 = dmax1(mn(1),mn(2),mn(3),mn(4))
      idummy = 1
      do i=1,4
         if (mn(i).eq.xx0) then
            iord(1) = i
         elseif (mn(i).eq.xx1) then
            iord(4) = i
         else
            irem(idummy) = i
            idummy = idummy+1
         endif
      end do
      if (mn(irem(1)).le.mn(irem(2))) then
         iord(2) = irem(1)
         iord(3) = irem(2)
      else
         iord(2) = irem(2)
         iord(3) = irem(1)
      end if

c               bw basis mixing matrix; mn(4) eigenvalues 
      do j=1,4
         i = iord(j)
         mn(j) = ymn(i)
         do i1=1,4
            bw(j,i1) = yz(i,i1)
         end do
      end do

c               pz basis mixing matrix 
      do n=1,4 
         pz(n,1) =   cw * bw(n,1) + sw * bw(n,2)
         pz(n,2) = - sw * bw(n,1) + cw * bw(n,2)
         pz(n,3) = bw(n,3)
         pz(n,4) = bw(n,4)
      end do

c               the charginos [kileng 3.116 with sb<->cb]
      sx = (m2**2+mu**2)/2.D0 + mw**2
      tx = (m2**2-mu**2)**2/4.D0 + mw**4*cbb**2  
     &    + mw**2*( m2**2 + mu**2 + 2.D0*mu*m2*sbb )

c               masses ordered [kileng 3.115 with mc* interchanged]
      mc(2) = sqrt(abs( sx + sqrt(tx) ))
      mc(1) = sqrt(abs( sx - sqrt(tx) ))

c               mixing angles [kileng 3.118]
      c2u = ( m2**2 - mu**2 - 2.D0*mw**2*cbb )/(mc(2)**2-mc(1)**2)
      c2v = ( m2**2 - mu**2 + 2.D0*mw**2*cbb )/(mc(2)**2-mc(1)**2)

c               sign of the two angles [kileng 3.122]
      sigu = sign(1.D0,m2*cb+mu*sb)
      sigv = sign(1.D0,m2*sb+mu*cb)
      
      cu  =      sqrt((1.D0+c2u)/2.D0)
      su  = sigu*sqrt((1.D0-c2u)/2.D0)
      cv  =      sqrt((1.D0+c2v)/2.D0)
      sv  = sigv*sqrt((1.D0-c2v)/2.D0)

c               mixing matrices 
      uux(1,1) =  cu 
      uux(1,2) =  su 
      uux(2,1) = -su
      uux(2,2) =  cu

      vvx(1,1) =  cv 
      vvx(1,2) =  sv 
      vvx(2,1) = -sv * sign(1.D0,mu*m2-mw**2*sbb)
      vvx(2,2) =  cv * sign(1.D0,mu*m2-mw**2*sbb)

c               one last sign to be fixed 
      cm(1,1) = m2
      cm(1,2) = sqrt(2.D0)*mw*sb
      cm(2,1) = sqrt(2.D0)*mw*cb
      cm(2,2) = mu
      
      do i=1,2
         do j=1,2
            dm(i,j) = 0.D0
         end do
      end do   

      do i=1,2
         do j=1,2
            do k=1,2
               do l=1,2
                  dm(i,j) = dm(i,j) + uux(i,k)*cm(k,l)*vvx(j,l)
               end do
            end do
         end do
      end do

      if (dm(1,1).lt.0.D0) then 
         do i=1,2
            do j=1,2
               uux(i,j) = -uux(i,j)
            end do
         end do
      end if 

c              order the charginos w.r.t. to their mass 
      if (dm(1,1).gt.dm(2,2)) then 
         do j=1,2
            uu(1,j) = uux(2,j)
            uu(2,j) = uux(1,j)
            vv(1,j) = vvx(2,j)
            vv(2,j) = vvx(1,j)
         end do
      else
         do j=1,2
            uu(1,j) = uux(1,j)
            uu(2,j) = uux(2,j)
            vv(1,j) = vvx(1,j)
            vv(2,j) = vvx(2,j)
         end do
      end if

      return
      end

c ----------------------------------------------------------------------
c        iteration routine to determine improved lambda's
      real*8 function XITLA(nloop,als_fix,q,acc)

      implicit none 

      integer nloop,nf,n
      real*8  als_fix,q,acc,aa,bb,als,xit,pi,a,b,x,lam,diff

c               the over-all prefactor 1/beta0 as usual 
      aa(nf,pi)= 12.D0*pi / (33.D0 - 2.D0*dble(nf))
c               beta1 as usual
      bb(nf,pi)= (153.D0-19.D0*dble(nf))/(33.D0-2.D0*dble(nf))
     &           /(2.D0*pi)

c               alpha_s in terms of x=log(q^2/lambda^2)
      als(nf,x) = aa(nf,pi)/x * ( 1.D0 - bb(nf,pi)*aa(nf,pi)*log(x)/x )

c               the iterated solution
      xit(a,b,x) = a/2.D0 * ( 1D0 + sqrt(1D0-4D0*b*log(x)) )

c               the prefactors in the calculation       
      nf = 5
      pi = 4.D0 * atan(1.D0)
      a = aa(nf,pi)/als_fix
      b = bb(nf,pi)*als_fix

c               the starting value in leading order 
      lam = q * dexp(-aa(nf,pi)/als_fix/2.D0)

c               leading order result 
      if (nloop.eq.1) then 
         XITLA = lam
         return
      end if 

c               the do loop for the iteration 
      do n=1,1000

c               iteration of lambda using initial value 
         x = log(q**2/lam**2)
         lam = q * dexp(-xit(a,b,x)/2.D0)
      
c               improved value for lambda
         x = log(q**2/lam**2) 
         diff = abs( (als(nf,x)-als_fix)/als_fix ) 
         
c               right value for alpha_s
         if (diff.lt.acc) then 
            XITLA = lam
            return
         end if 

      end do

c                in case the accuracy is not reached 
      print*,'Xsugra.f:  iteration for lam_qcd stopped, lambda =',lam
      XITLA = lam
     
      return 
      end

c ----------------------------------------------------------------------
c               spirix' version 
      real*8 function B_SPECIAL(s,m1,m2,mu2)

      implicit none

      real*8     s,m1,m2,mu2,m12,m22,A,B 
      complex*16 zkappa,x1,x2 

      A(m12,s) = m12 * ( 1.D0 - log(m12/s) )

      m12 = m1**2
      m22 = m2**2

c               copied from Xtwopoint_spira.f 
      zkappa = sqrt( dcmplx(s**2 + m12**2 + m22**2
     &                     -2.D0 * (s*m12 + s*m22 + m12*m22) ))

      x1 = dcmplx( (s-m22+m12+zkappa)/(2.D0*s) )
      x2 = dcmplx( (s-m22+m12-zkappa)/(2.D0*s) )

c               B02(s,m1,m2,s)
      if ((real(x1).eq.0.d0).or.(real(x2).eq.0.d0)) then
         B = 2.3D0   ! This is just an average value... -SA
      else
         B = real( 2.D0 + log(s/m22)
     &                + x1 * log(1.D0-1.D0/x1)
     &                + x2 * log(1.D0-1.D0/x2) )
      end if

      B_SPECIAL = ( A(m22,s) - A(m12,s) + (m12-m22-s) * B )/2.D0/s
     &           + log(s/mu2) /2.D0

      return
      end

c ----------------------------------------------------------------------
c               check all the mass matrix diagonalization 
c               all variables are input 
      subroutine DIAGONALIZE(tb,mw,mz,sw,mt,alpmz,mtau
     &                      ,m1,m2,mu,mql,mur,mdr,mel,mer
     &                      ,at,ab,al,s2t,s2b,s2e,bw,uu,vv)

      implicit none 

      integer i,j,k,l
      real*8  m1,m2,s2t,s2b,s2e,mql,bw(4,4),uu(2,2),vv(2,2)
     &       ,cbet,sbet,cw,cthet,sthet,ctheb,stheb,cthel,sthel
     &       ,sw2,rmt,rmb,cbb,at,ab,al,mu,RUNM_ORIG
     &       ,mul,mur,mdl,mdr,mel,mer,c2t,c2b,c2e
     &       ,mstl2,mstr2,msbl2,msbr2,msll2,mslr2
     &       ,xm(4,4),ym(4,4),cm(2,2),dm(2,2)
     &       ,st(2,2),zt(2,2),dt(2,2)
     &       ,sb(2,2),zb(2,2),db(2,2)
     &       ,sl(2,2),zl(2,2),dl(2,2)
     &       ,mz,mw,sw,alpmz,mt,mtau,tb,del_check

      external RUNM_ORIG

c               accuracy for the cheeck
      del_check = 1.D-4

c               the susy mixing angle beta 
      cbet = 1.D0 / sqrt(1.D0+tb**2)
      sbet = tb * cbet
      cbb  = cbet**2 - sbet**2 

c               weak mixing angle 
      cw  = sqrt(1.D0-sw**2)
      sw2 = sw**2

c               all kind of mixing angles 
      c2t = sqrt(1.D0-s2t**2)
      c2b = sqrt(1.D0-s2b**2)
      c2e = sqrt(1.D0-s2e**2)

      sthet = sign(1.D0,s2t) * sqrt((1.D0+c2t)/2.D0)
      cthet =                  sqrt((1.D0-c2t)/2.D0)
      stheb = sign(1.D0,s2b) * sqrt((1.D0+c2b)/2.D0)
      ctheb =                  sqrt((1.D0-c2b)/2.D0)
      sthel = sign(1.D0,s2e) * sqrt((1.D0+c2e)/2.D0)
      cthel =                  sqrt((1.D0-c2e)/2.D0)

c               the doublet masses
      mdl = mql
      mul = mql

c               initialize the neutralino, chargino, sfermion matrices 
      do i=1,4
         do j=1,4            
            xm(i,j) = 0
            ym(i,j) = 0
         end do
      end do
      do i=1,2
         do j=1,2
            dm(i,j) = 0
            st(i,j) = 0
            zt(i,j) = 0
            dt(i,j) = 0
            sb(i,j) = 0
            zb(i,j) = 0
            db(i,j) = 0
            sl(i,j) = 0
            zl(i,j) = 0
            dl(i,j) = 0
         end do
      end do

c               neutralino mass matrix [thesis A.6]
      xm(1,1) =  m1
      xm(1,3) = -mz*sw*cbet
      xm(1,4) =  mz*sw*sbet
      xm(2,2) =  m2
      xm(2,3) =  mz*cw*cbet
      xm(2,4) = -mz*cw*sbet
      xm(3,1) = -mz*sw*cbet
      xm(3,2) =  mz*cw*cbet
      xm(3,4) = -mu
      xm(4,1) =  mz*sw*sbet
      xm(4,2) = -mz*cw*sbet
      xm(4,3) = -mu

c               chargino mass matrix [thesis A.11]
      cm(1,1) = m2
      cm(1,2) = sqrt(2.D0)*mw*sbet
      cm(2,1) = sqrt(2.D0)*mw*cbet
      cm(2,2) = mu

c               diagonalized neutralino matrix
      do i=1,4
         do j=1,4
            do k=1,4
               do l=1,4
                  ym(i,j) = ym(i,j) + bw(i,k)*xm(k,l)*bw(j,l)
               end do
            end do
         end do
      end do

c               diagonalized chargino matrix 
      do i=1,2
         do j=1,2
            do k=1,2
               do l=1,2
                  dm(i,j) = dm(i,j) + uu(i,k)*cm(k,l)*vv(j,l)
               end do
            end do
         end do
      end do

c               the running top and bottom mass
      rmt = RUNM_ORIG(mt,6) * (1.D0 + 8.D0/3.D0 * alpmz**2)
      rmb = RUNM_ORIG(mt,5) * (1.D0 + 8.D0/3.D0 * alpmz**2)

c               stop mass matrix [dm 12a]
      mstl2 = mul**2 + (0.5D0-2.D0/3.D0*sw2) * mz**2*cbb
      mstr2 = mur**2 +        2.D0/3.D0*sw2  * mz**2*cbb
      st(1,1) = mstl2 + rmt**2
      st(1,2) = -rmt * (at+mu/tb)
      st(2,1) = -rmt * (at+mu/tb)
      st(2,2) = mstr2 + rmt**2
c               stop rotation matrix 
      zt(1,1) =  cthet
      zt(1,2) =  sthet
      zt(2,1) = -sthet
      zt(2,2) =  cthet

c               sbottom mass matrix [dto]
      msbl2 = mdl**2 + (-0.5D0+1.D0/3.D0*sw2) * mz**2*cbb
      msbr2 = mdr**2 -         1.D0/3.D0*sw2  * mz**2*cbb
      sb(1,1) = msbl2 + rmb**2
      sb(1,2) = -rmb * (ab+mu*tb)
      sb(2,1) = -rmb * (ab+mu*tb)
      sb(2,2) = msbr2 + rmb**2
c               sbottom rotation matrix 
      zb(1,1) =  ctheb
      zb(1,2) =  stheb
      zb(2,1) = -stheb
      zb(2,2) =  ctheb

c               stau mass matrix [dm 12b]
      msll2 = mel**2 + (-0.5D0+sw2) * mz**2*cbb
      mslr2 = mer**2 -         sw2  * mz**2*cbb
      sl(1,1) = msll2 + mtau**2
      sl(1,2) = -mtau * (al+mu*tb)
      sl(2,1) = -mtau * (al+mu*tb)
      sl(2,2) = mslr2 + mtau**2
c               stau rotation matrix 
      zl(1,1) =  cthel
      zl(1,2) =  sthel
      zl(2,1) = -sthel
      zl(2,2) =  cthel

c               test the rotation for stop,sbottom,stau
      do i=1,2
         do j=1,2
            do k=1,2
               do l=1,2
                  dt(i,j) = dt(i,j) + zt(i,k)*st(k,l)*zt(j,l)
                  db(i,j) = db(i,j) + zb(i,k)*sb(k,l)*zb(j,l)
                  dl(i,j) = dl(i,j) + zl(i,k)*sl(k,l)*zl(j,l)
               end do
            end do
         end do
      end do

c               check if all mass matrices are diagonalized correctly
      do i=1,2
         if (dt(i,i).lt.0.D0) print*,'Xsugra.f:  mass not real: st',i
         if (db(i,i).lt.0.D0) print*,'Xsugra.f:  mass not real: sb',i
         if (dl(i,i).lt.0.D0) print*,'Xsugra.f:  mass not real: sl',i
         if (dm(i,i).lt.0.D0) print*,'Xsugra.f:  mass negative: ch',i
         do j=1,2 
            if (i.ne.j) then 
             if (abs(dt(i,j)).gt.del_check)
     &        print*,'Xsugra.f: st mass matrix non diagonal',i,j,dt(i,j)
             if (abs(db(i,j)).gt.del_check)
     &        print*,'Xsugra.f: sb mass matrix non diagonal',i,j,db(i,j)
             if (abs(dl(i,j)).gt.del_check)
     &        print*,'Xsugra.f: sl mass matrix non diagonal',i,j,dl(i,j)
             if (abs(dm(i,j)).gt.del_check)
     &        print*,'Xsugra.f: ch mass matrix non diagonal',i,j,dm(i,j)
            end if 
         end do
      end do
      
      do i=1,4
         do j=1,4 
            if (i.ne.j) then 
             if (abs(ym(i,j)).gt.del_check)
     &        print*,'Xsugra.f: nn mass matrix non diagonal',i,j,ym(i,j)
            end if 
         end do
      end do

c               all the test output 
c$$$      write(6,*)' testing the diagonalization : '
c$$$      write(6,*)
c$$$      write(6,'(a6,/,4(4(f14.8,1x),/))')'m0:',((xm(i,j),j=1,4),i=1,4)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,4(4(f14.8,1x),/))')' m:',((ym(i,j),j=1,4),i=1,4)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,4(4(f14.8,1x),/))')' z:',((bw(i,j),j=1,4),i=1,4)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.8,1x),/))')'m0:',((cm(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.8,1x),/))')' m:',((dm(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.8,1x),/))')' u:',((uu(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.8,1x),/))')' v:',((vv(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f15.6,1x),/))')'mt0:',((st(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f12.5,1x),/))')'mt:',
c$$$     .                            ((sqrt(abs(dt(i,j))),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.8,1x),/))')'zt:',((zt(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f15.6,1x),/))')'mb0:',((sb(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.5,1x),/))')'mb:',
c$$$     .                            ((sqrt(abs(db(i,j))),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.8,1x),/))')'zb:',((zb(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f15.6,1x),/))')'ml0:',((sl(i,j),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.5,1x),/))')'ml:',
c$$$     .                            ((sqrt(abs(dl(i,j))),j=1,2),i=1,2)
c$$$      write(6,*)
c$$$      write(6,'(a6,/,2(2(f14.8,1x),/))')'zl:',((zl(i,j),j=1,2),i=1,2)
c$$$      write(6,*)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THIS PROGRAM COMPUTES THE RENORMALIZATION GROUP IMPROVED
C     VALUES OF HIGGS MASSES AND COUPLINGS IN THE MSSM.
C
C     INPUT: MA,TANB = TAN(BETA),MQ,MUR,MDR,MTOP,AU,AD,MU,MCHI
C
C     ALL MASSES IN GEV UNITS. MA IS THE CP-ODD HIGGS MASS,
C     MTOP IS THE PHYSICAL TOP MASS, MQ AND MUR/MDR ARE THE SOFT
C     SUPERSYMMETRY BREAKING MASS PARAMETERS OF LEFT HANDED
C     AND RIGHT HANDED STOPS RESPECTIVELY, AU AND AD ARE THE
C     STOP AND SBOTTOM TRILINEAR SOFT BREAKING TERMS,
C     RESPECTIVELY,  AND MU IS THE SUPERSYMMETRIC
C     HIGGS MASS PARAMETER. WE USE THE  CONVENTIONS FROM
C     THE PHYSICS REPORT OF HABER AND KANE: LEFT RIGHT
C     STOP MIXING TERM PROPORTIONAL TO (AU - MU/TANB).
C     MCHI IS THE HEAVIEST CHARGINO MASS. 
C     WE USE AS INPUT TANB DEFINED AT THE SCALE MTOP.

C     OUTPUT: MH,HM,MCH, SA = SIN(ALPHA), CA= COS(ALPHA), TANBA
C     WHERE MHP AND HPM ARE THE LIGHTEST AND HEAVIEST CP-EVEN
C     HIGGS MASSES, MHCH IS THE CHARGED HIGGS MASS AND
C     ALPHA IS THE HIGGS MIXING ANGLE.
C     TANBA IS THE ANGLE TANB AT THE CP-ODD HIGGS MASS SCALE.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Program based on the work by M. Carena, M. Quiros
c       and C.E.M. Wagner, "Effective potential methods and
c       the Higgs mass spectrum in the MSSM", Phys. Lett.
c       B335 (1995) 209. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine SUBH(ma,tb,mw,mz,sw,gf,mt
     &               ,mq,mur,md,au,ad,mu,mchi,default
     &               ,mhl,mhh,mhc,sa,ca,tanba,lambda)

      implicit none
      
      integer i,j
      real*8 ma,mq,mur,md,au,ad,mu,mchi,mhl,mhh,mhc,sa,ca
     &      ,tanba,tanbt,v,sw2,cw2,alpha2,alpha1,alpha3
     &      ,rmb,rmt,tq,tu,td,sinb,cosb,g1,g2,g3,hu,hd
     &      ,tp,tpd,tdp,tdpd,dlambda1,dlambda2,dlambda3,dlambda4
     &      ,lambda(1:7),vh(2,2),m2(2,2),m2p(2,2),mssusy,tchar
     &      ,deltal12,deltal3p4,deltam112,deltam122,deltam222
     &      ,trm2p,detm2p,mhl2,mhh2,mhc2,cos2alpha
     &      ,RUNM_ORIG,ALPHAS
     &      ,pi,mz,mw,sw,gf,mt,tb,default

      external RUNM_ORIG,ALPHAS

      pi = 4.D0 * atan(1.D0)

      v   = 1/sqrt(2.D0*sqrt(2.D0)*gf)
      sw2 = sw**2 
      cw2 = 1.D0-sw2 

      alpha2  = (2.D0*mw/v/sqrt(2.D0))**2/4.D0/pi
      alpha1  = alpha2*sw2/cw2

c               alphas a la spira has to be initialized 
      alpha3  = ALPHAS(mt,2)

c               running quark masses used 
      rmb = RUNM_ORIG(mt,5)
      rmt = RUNM_ORIG(mt,6)

      tq = log((mq**2  + mt**2)/mt**2)
      tu = log((mur**2 + mt**2)/mt**2)
      td = log((md**2  + mt**2)/mt**2)

c               initialize the tan(beta) values at different scales 
      tanba = tb
      tanbt = tb
      sinb  = tb/sqrt(1.D0 + tb**2)
      cosb  = 1.D0/sqrt(1.D0 + tb**2)
      
c               tan(beta) reset for 
      if (ma.gt.mt)
     *     tanba = tb*(1.D0-3.D0/32.D0/pi**2*
     *             (rmt**2/v**2/sinb**2-rmb**2/v**2/cosb**2)*
     *             log(ma**2/mt**2))
      if (ma.le.mt) tanbt = tanba

c               sin/cos(beta) set at scale mt 
      sinb = tanbt/sqrt(1.D0 + tanbt**2)
      cosb = 1.D0/sqrt(1.D0 + tanbt**2)

c               gauge couplings at scale mz
      g1 = sqrt(alpha1*4.D0*pi)
      g2 = sqrt(alpha2*4.D0*pi)
      g3 = sqrt(alpha3*4.D0*pi)

c               yukawa couplings at scale mt 
      hu = rmt/v/sinb
      hd = rmb/v/cosb

      if (mq.gt.mur) tp   = tq - tu
      if (mq.le.mur) tp   = tu - tq
      if (mq.gt.mur) tdp  = tu
      if (mq.le.mur) tdp  = tq
      if (mq.gt.md)  tpd  = tq - td
      if (mq.le.md)  tpd  = td - tq
      if (mq.gt.md)  tdpd = td
      if (mq.le.md)  tdpd = tq

      if (mq.gt.md) then 
         dlambda1 = 6.D0/96.D0/pi**2*hd**2*g1**2*tpd
      else if (mq.le.md) then 
         dlambda1 = 3.D0/32.D0/pi**2*hd**2*(g1**2/3.D0+g2**2)*tpd
      end if 
         
      if(mq.gt.mur) then 
         dlambda2 = 12.D0/96.D0/pi**2*hu**2*g1**2*tp
      else if (mq.le.mur) then 
         dlambda2 =  3.D0/32.D0/pi**2*hu**2*(-g1**2/3.D0+g2**2)*tp
      end if 

      dlambda3 = 0.D0
      dlambda4 = 0.D0

      if (mq.gt.md) then 
         dlambda3 = -1.D0/32.D0/pi**2*hd**2*g1**2*tpd
      else if (mq.le.md) then 
         dlambda3 =  3.D0/64.D0/pi**2*hd**2*(g2**2-g1**2/3.D0)*tpd
      end if 

      if (mq.gt.mur) then 
         dlambda3 =dlambda3-1.D0/16.D0/pi**2*hu**2*g1**2*tp
      else if (mq.le.mur) then 
         dlambda3 =dlambda3+3.D0/64.D0/pi**2*hu**2*(g2**2+g1**2/3.D0)*tp
      end if 

      if (mq.lt.mur) then 
         dlambda4 = -3.D0/32.D0/pi**2*hu**2*g2**2*tp
      end if 

      if (mq.lt.md) then
         dlambda4 = dlambda4 - 3.D0/32.D0/pi**2*g2**2*hd**2*tpd
      end if 

c               all the usual lambdas, as defined in kane's book 
      lambda(1) = ((g1**2+g2**2)/4.D0)
     &           * (1.D0-3.D0*hd**2*(tpd + tdpd)/8.D0/pi**2)
     &        + (3.D0*hd**4.D0/16.D0/pi**2) * tpd 
     &           * (1.D0 + (3.D0*hd**2/2.D0 + hu**2/2.D0 - 8.D0*g3**2)
     &                      * (tpd + 2.D0*tdpd)/16.D0/pi**2            )
     &        + (3.D0*hd**4.D0/8.D0/pi**2) * tdpd
     &           * (1.D0  + (3.D0*hd**2/2.D0 + hu**2/2.D0 - 8.D0*g3**2)
     &                       * tdpd/16.D0/pi**2                        )
     &        + dlambda1 

      lambda(2) = ((g1**2 + g2**2)/4.D0)
     &           *(1.D0-3.D0*hu**2*(tp + tdp)/8.D0/pi**2)
     &        + (3.D0*hu**4.D0/16.D0/pi**2) * tp
     &           *(1.D0 + (3.D0*hu**2/2.D0 + hd**2/2.D0 - 8.D0*g3**2) 
     &                     * (tp + 2.D0*tdp)/16.D0/pi**2               )
     &        + (3.D0*hu**4.D0/8.D0/pi**2) * tdp
     &           * (1.D0 + (3.D0*hu**2/2.D0 + hd**2/2.D0 - 8.D0*g3**2)
     &                      * tdp/16.D0/pi**2                          )
     &        + dlambda2 

      lambda(3) = ((g2**2 - g1**2)/4.D0)
     &           *(1.D0-3.D0*hu**2*(tp + tdp)/16.D0/pi**2 
     &                - 3.D0*hd**2*(tpd + tdpd)/16.D0/pi**2) 
     &        + dlambda3 

      lambda(4) = (- g2**2/2.D0)
     &           *(1.D0-3.D0*(hu**2)*(tp + tdp)/16.D0/pi**2
     &                 -3.D0*(hd**2)*(tpd + tdpd)/16.D0/pi**2) 
     &        + dlambda4

      lambda(5) = 0.D0
      lambda(6) = 0.D0
      lambda(7) = 0.D0

c               the cp even mass matrix 
      m2(1,1) = 2.D0*v**2*(  lambda(1)*cosb**2 
     &                     + 2.D0*lambda(6)*cosb*sinb 
     &                     + lambda(5)*sinb**2) 
     &          + ma**2*sinb**2
      m2(2,2) = 2.D0*v**2*(  lambda(5)*cosb**2 
     &                     + 2.D0*lambda(7)*cosb*sinb 
     &                     + lambda(2)*sinb**2) 
     &          + ma**2*cosb**2
      m2(1,2) = 2.D0*v**2*(  lambda(6)*cosb**2 
     &                     + ( lambda(3) + lambda(4) )*cosb*sinb
     &                     + lambda(7)*sinb**2) 
     &          - ma**2*sinb*cosb
      m2(2,1) = m2(1,2)

c               contribution from light charginos/neutralinos
      mssusy=sqrt(0.5D0*(mq**2+mur**2)+mt**2)
      if (mchi.le.mssusy) then 
         if (mchi.lt.mt) mchi=mt
         
         tchar = log(mssusy**2/mchi**2)

         deltal12  = ( 9.D0/64.D0 /pi**2*g2**4
     &               + 5.D0/192.D0/pi**2*g1**4 )*tchar
         deltal3p4 = ( 3.D0/64.D0 /pi**2*g2**4
     &               + 7.D0/192.D0/pi**2*g1**4
     &               + 4.D0/32    /pi**2*g1**2*g2**2 )*tchar
         deltam112 = 2.D0 * deltal12  *v**2 * cosb**2
         deltam222 = 2.D0 * deltal12  *v**2 * sinb**2
         deltam122 = 2.D0 * deltal3p4 *v**2 * sinb*cosb

         m2(1,1) = m2(1,1)+deltam112
         m2(2,2) = m2(2,2)+deltam222
         m2(1,2) = m2(1,2)+deltam122
         m2(2,1) = m2(2,1)+deltam122
      end if 

      call GFUN(ma,tb,mw,mz,sw,gf,mt,mq,mur,md,au,ad,mu,default,vh)
      
      do i=1,2
         do j=1,2
            m2p(i,j) = m2(i,j) + vh(i,j)
         end do
      end do

c               trace and determinant of the modified mass matrix
      trm2p  = m2p(1,1) + m2p(2,2)
      detm2p = m2p(1,1)*m2p(2,2) - m2p(1,2)*m2p(2,1)

c               light and heavy cp even higgs mass  
      mhl2 = (trm2p - sqrt(trm2p**2 - 4.D0* detm2p))/2.D0
      mhh2 = (trm2p + sqrt(trm2p**2 - 4.D0* detm2p))/2.D0

      if (mhl2.gt.0.D0) then 
         mhl = sqrt(mhl2)
      else 
         print *, 'Xsugra.f: m_higgs_l^2 negative '
         mhl = default
      end if 

      if (mhh2.gt.0.D0) then 
         mhh = sqrt(mhh2)
      else 
         print *, 'Xsugra.f: m_higgs_h^2 negative '
         mhh = default
      end if 

c               charged higgs mass !!!!!!!!!!!!!!!!
      mhc2 = ma**2 + (lambda(5)-lambda(4)) * v**2

      if (mhc2.gt.0.D0) then 
         mhc = sqrt(mhc2)
      else 
         print *, 'Xsugra.f: m_higgs_c^2 negative '
         mhc = default
      end if 

c               the mixing angle alpha 
ctp      sin2alpha = 2.D0*m2p(1,2)/sqrt(trm2p**2-4.D0*detm2p)
      cos2alpha = (m2p(1,1)-m2p(2,2))/sqrt(trm2p**2-4.D0*detm2p)
ctp      if((cos2alpha.gt.0.D0).and.(abs(sin2alpha).lt.1.D0))
ctp     &                       alpha = dasin(sin2alpha)/2.D0
ctp      if((cos2alpha.lt.0.D0).and.(abs(sin2alpha).lt.-1.D0))
ctp     &                       alpha = -pi/2.D0-dasin(sin2alpha)/2.D0
ctp      if(cos2alpha.gt.0.D0) alpha = dasin(sin2alpha)/2.D0
ctp      if(cos2alpha.lt.0.D0) alpha = -pi/2.D0-dasin(sin2alpha)/2.D0

      sa =-sqrt((1.D0-cos2alpha)/2.D0)
      ca = sqrt((1.D0+cos2alpha)/2.D0)

      return
      end

c ----------------------------------------------------------------------
c               non degenerate stop/sbottom effects
c               output: vh
      subroutine GFUN(ma,tanb,mw,mz,sw,gf,mt
     &               ,mq,mur,md,au,ad,mu,default,vh)

      implicit none 

      integer i,j
      real*8 ma,mq,mur,md,au,ad,mu,au1,ad1
     &      ,vh(2,2),vh3t(2,2),vh3b(2,2),al(2,2)
     &      ,x,y,g,mq2,mur2,md2,sinb,cosb,tanba,sinba,cosba
     &      ,v,sw2,cw2,alpha2,alpha1,alpha3,g1,g2,g32
     &      ,rmb,rmt,mst,msusyt,msb,msusyb,tt,tb,ht,htst,hb
     &      ,bt2,bb2,al2,al1,mt4,mt2,mb4,mb2
     &      ,vi,h1i,h2i,h1t,h2t,h1b,h2b,tanbst,sinbt
     &      ,tanbsb,sinbb,cosbb,stop12,stop22,sbot12,sbot22
     &      ,stop1,stop2,sbot1,sbot2,f1t,f1b,f2t,f2b
     &      ,RUNM_ORIG,ALPHAS
     &      ,pi,mz,mw,sw,gf,mt,tanb,default

      external RUNM_ORIG

      g(x,y) = 2.D0 - (x+y)/(x-y)*log(x/y)

      pi = 4.D0 * atan(1.D0)

c               wrong sign of a_{b,t} in squark mass matrix
      au1 = -au 
      ad1 = -ad 

      mq2   = mq**2
      mur2  = mur**2
      md2   = md**2

      sinb  = tanb/sqrt(tanb**2+1.D0)
      cosb  = sinb/tanb

      tanba = tanb
      sinba = tanba/sqrt(tanba**2+1.D0)
      cosba = sinba/tanba        

c               exactly the same as in SUBH
      v  = 1/sqrt(2.D0*sqrt(2.D0)*gf)
      sw2 = sw**2 
      cw2 = 1.D0-sw2

      alpha2  = (2.D0*mw/v/sqrt(2.D0))**2/4.D0/pi
      alpha1  = alpha2*sw2/cw2
      alpha3  = ALPHAS(mt,2)

      g1  = sqrt(alpha1*4.D0*pi)
      g2  = sqrt(alpha2*4.D0*pi)
      g32 = alpha3*4.D0*pi
      
      rmb = RUNM_ORIG(mt,5)
      rmt = RUNM_ORIG(mt,6)

c               WHY NO RUNNING TOP MASS ???
      mst = dmax1(mq,mur)
      msusyt = sqrt(mst**2  + mt**2)

      msb = dmax1(mq,md)
      msusyb = sqrt(msb**2 + rmb**2)

      tt = log(msusyt**2/mt**2)
      tb = log(msusyb**2/mt**2)

c               yukawa couplings using running masses 
      ht   = rmt/v/sinb
      htst = rmt/v
      hb   = rmb/v/cosb

      bt2   = -(8.D0*g32 - 9.D0*ht**2/2.D0 - hb**2/2.D0)/(4.D0*pi)**2
      bb2   = -(8.D0*g32 - 9.D0*hb**2/2.D0 - ht**2/2.D0)/(4.D0*pi)**2
      al2   = 3.D0/8.D0/pi**2 * ht**2
      al1   = 3.D0/8.D0/pi**2 * hb**2

      al(1,1) = al1
      al(1,2) = (al2+al1)/2.D0
      al(2,1) = (al2+al1)/2.D0
      al(2,2) = al2
      
      mt4 = rmt**4 * (1.D0 + 2.D0*bt2*tt - al2*tt)
      mt2 = sqrt(mt4)
      mb4 = rmb**4 * (1.D0 + 2.D0*bb2*tb - al1*tb)
      mb2 = sqrt(mb4)

      if(ma.gt.mt) then
         vi  = v   * (1.D0 + 3.D0/32.D0/pi**2*htst**2*log(mt**2/ma**2))
         h1i = vi * cosba
         h2i = vi * sinba
         h1t = h1i * sqrt(sqrt((1.D0 
     &                  + 3.D0/8.D0/pi**2*hb**2*log(ma**2/msusyt**2))))
         h2t = h2i * sqrt(sqrt((1.D0
     &                  + 3.D0/8.D0/pi**2*ht**2*log(ma**2/msusyt**2))))
         h1b = h1i * sqrt(sqrt((1.D0
     &                  + 3.D0/8.D0/pi**2*hb**2*log(ma**2/msusyb**2))))
         h2b = h2i * sqrt(sqrt((1.D0
     &                  + 3.D0/8.D0/pi**2*ht**2*log(ma**2/msusyb**2))))
      else
         vi  =  v
         h1i = vi*cosb
         h2i = vi*sinb
         h1t = h1i * sqrt(sqrt((1.D0
     &                  + 3.D0/8.D0/pi**2*hb**2*log(mt**2/msusyt**2))))
         h2t = h2i * sqrt(sqrt((1.D0
     &                  + 3.D0/8.D0/pi**2*ht**2*log(mt**2/msusyt**2))))
         h1b = h1i * sqrt(sqrt((1.D0
     &                  + 3.D0/8.D0/pi**2*hb**2*log(mt**2/msusyb**2))))
         h2b = h2i * sqrt(sqrt((1.D0
     &                  + 3.D0/8.D0/pi**2*ht**2*log(mt**2/msusyb**2))))
      end if

      tanbst = h2t/h1t
      sinbt = tanbst/sqrt(1.D0+tanbst**2)

      tanbsb = h2b/h1b
      sinbb = tanbsb/sqrt(1.D0+tanbsb**2)
      cosbb = sinbb/tanbsb

c               note that st2 and sb2 are the light particles  
      stop12 = (mq2 + mur2)*.5D0 + mt2 
     &       + 1.D0/8.D0*(g2**2+g1**2)*(h1t**2-h2t**2)
     &       + sqrt( ( (g2**2-5.D0*g1**2/3.D0)/4.D0 * (h1t**2-h2t**2)
     &                  + mq2 - mur2 )**2*0.25D0 
     &               + mt2*(au1-mu/tanbst)**2 )
      
      stop22 = (mq2 + mur2)*.5D0 + mt2 
     &       + 1.D0/8.D0*(g2**2+g1**2)*(h1t**2-h2t**2) 
     &       - sqrt( ( (g2**2-5.D0*g1**2/3.D0)/4.D0 * (h1t**2-h2t**2)
     &                  + mq2 - mur2 )**2*0.25D0 
     &              + mt2*(au1-mu/tanbst)**2 )

      sbot12 = (mq2 + md2)*.5D0  
     &       - 1.D0/8.D0*(g2**2+g1**2)*(h1b**2-h2b**2)
     &       + sqrt( ( (g1**2/3.D0-g2**2)/4.D0 * (h1b**2-h2b**2) 
     &                  + mq2 - md2)**2*0.25D0 
     &              + mb2*(ad1-mu*tanbsb)**2 )

      sbot22 = (mq2 + md2)*.5D0  
     &       - 1.D0/8.D0*(g2**2+g1**2)*(h1b**2-h2b**2)
     &       - sqrt( ( (g1**2/3.D0-g2**2)/4.D0*(h1b**2-h2b**2)
     &                  + mq2 - md2)**2*0.25D0 
     &              + mb2*(ad1-mu*tanbsb)**2 )

      if ((stop22.lt.0.D0).or.(sbot22.lt.0.D0)) then 
         do i =1,2
            do j = 1,2
               vh(i,j) = -1.d+15
            end do
         end do
         return 
      end if 

      if (stop12.gt.0.D0) then 
         stop1 = sqrt( stop12 )
      else 
         print *,'Xsugra.f: (st1 yukawa contribution)^2 negative '
         stop1 = default 
      end if

      if (stop22.gt.0.D0) then 
         stop2 = sqrt( stop22 ) 
      else 
         print *,'Xsugra.f: (st2 yukawa contribution)^2 negative '
         stop2 = default 
      end if

      if (sbot12.gt.0.D0) then 
         sbot1 = sqrt( sbot12 )
      else 
         print *,'Xsugra.f: (sb1 yukawa contribution)^2 negative '
         sbot1 = default 
      end if

      if (sbot22.gt.0.D0) then 
         sbot2 = sqrt( sbot22 )
      else 
         print *,'Xsugra.f: (sb2 yukawa contribution)^2 negative '
         sbot2 = default 
      end if

c               d-terms
      f1t = (mq2-mur2)/(stop12-stop22) 
     &                * (.5D0-4.D0/3.D0*sw2) * log(stop1/stop2)
     &     + (.5D0-2.D0/3.D0*sw2) * log(stop1*stop2/(mq2+mt2))
     &     +       2.D0/3.D0*sw2  * log(stop1*stop2/(mur2+mt2))
      
      f1b = (mq2-md2)/(sbot12-sbot22)
     &               * (-.5D0+2.D0/3.D0*sw2) * log(sbot1/sbot2)
     &     + (-.5D0+1.D0/3.D0*sw2) * log(sbot1*sbot2/(mq2+mb2))
     &     -        1.D0/3.D0*sw2  * log(sbot1*sbot2/(md2+mb2))

      f2t = sqrt(mt2)*(au1-mu/tanbst)/(stop12-stop22)
     &                * (-.5D0*log(stop12/stop22)
     &     + (4.D0/3.D0*sw2-.5D0) * (mq2-mur2)/(stop12-stop22)
     &                            * g(stop12,stop22))

      f2b = sqrt(mb2)*(ad1-mu*tanbsb)/(sbot12-sbot22)
     &                * (.5D0*log(sbot12/sbot22)
     &     + (-2.D0/3.D0*sw2+.5D0) * (mq2-md2)/(sbot12-sbot22)
     &                             * g(sbot12,sbot22))

      vh3b(1,1) = mb4/(cosbb**2)
     &  * ( log(sbot1**2*sbot2**2/(mq2+mb2)/(md2+mb2)) 
     &     + 2.D0 * (ad1*(ad1-mu*tanbsb)/(sbot1**2-sbot2**2))
     &            * log(sbot1**2/sbot2**2))  
     &     + mb4/(cosbb**2)*(ad1*(ad1-mu*tanbsb)/(sbot1**2-sbot2**2))**2
     &          * g(sbot12,sbot22) 

      vh3t(1,1) = mt4/(sinbt**2)
     &  * (mu*(-au1+mu/tanbst)/(stop1**2-stop2**2))**2*g(stop12,stop22)

      vh3b(1,1) = vh3b(1,1) + mz**2*(2*mb2*f1b-sqrt(mb2)*ad1*f2b)

      vh3t(1,1) = vh3t(1,1) + mz**2*(sqrt(mt2)*mu/tanbst*f2t)  

      vh3t(2,2) = mt4/(sinbt**2)
     &  * ( log(stop1**2*stop2**2/(mq2+mt2)/(mur2+mt2)) 
     &     + 2.D0 * (au1*(au1-mu/tanbst)/(stop1**2-stop2**2))
     &            * log(stop1**2/stop2**2))
     &     + mt4/(sinbt**2)*(au1*(au1-mu/tanbst)/(stop1**2-stop2**2))**2
     &          * g(stop12,stop22) 

      vh3b(2,2) = mb4/(cosbb**2)
     &  * (mu*(-ad1+mu*tanbsb)/(sbot1**2-sbot2**2))**2*g(sbot12,sbot22)

      vh3t(2,2) = vh3t(2,2) + mz**2*(-2*mt2*f1t+sqrt(mt2)*au1*f2t)

      vh3b(2,2) = vh3b(2,2) - mz**2*sqrt(mb2)*mu*tanbsb*f2b  
 
      vh3t(1,2) = -mt4/(sinbt**2)*mu*(au1-mu/tanbst)/(stop1**2-stop2**2)
     &  * ( log(stop1**2/stop2**2)
     &     + au1*(au1 - mu/tanbst)/(stop1**2-stop2**2)*g(stop12,stop22))

      vh3b(1,2) = -mb4/(cosbb**2)*mu*(au1-mu*tanbsb)/(sbot1**2-sbot2**2)
     &  * ( log(sbot1**2/sbot2**2)
     &     + ad1*(ad1 - mu*tanbsb)/(sbot1**2-sbot2**2)*g(sbot12,sbot22))

      vh3t(1,2) = vh3t(1,2)
     &    + mz**2*( mt2/tanbst*f1t-sqrt(mt2)*(au1/tanbst+mu)/2.D0*f2t)

      vh3b(1,2) = vh3b(1,2)
     &    + mz**2*(-mb2*tanbsb*f1b+sqrt(mb2)*(ad1*tanbsb+mu)/2.D0*f2b)

      vh3t(2,1) = vh3t(1,2)
      vh3b(2,1) = vh3b(1,2)

      do i = 1,2
         do j = 1,2
            vh(i,j) = 6.D0/(8.D0 * pi**2 * (h1t**2+h2t**2))
     &                    * vh3t(i,j) * 0.5D0 * (1.D0-al(i,j)*tt/2.D0) 
     &              + 6.D0/(8.D0 * pi**2 * (h1b**2+h2b**2))
     &                    * vh3b(i,j) * 0.5D0 * (1.D0-al(i,j)*tb/2.D0)
         end do 
      end do
       
      return
      end 




c ----------------------------------------------------------------------
c               higgs sector to leading order
c input : ma,tb,mw,mz
      subroutine BORNH(ma,tb,mw,mz,mhl,mhh,mhc,sa,ca)

      implicit none 

      real*8 ma,mhl,mhh,mhc,sa,ca,cb,cbb,delta,caa
     &      ,mw,mz,tb

      cb  = 1.D0 / sqrt(1.D0+tb**2)
      cbb = 2.D0*cb**2 - 1.D0
  
      delta = sqrt( (ma**2+mz**2)**2 - (2.D0*ma*mz*cbb)**2 ) 
      mhl = sqrt( (ma**2 + mz**2 - delta)/2.D0 )
      mhh = sqrt( (ma**2 + mz**2 + delta)/2.D0 )
      mhc = sqrt( ma**2 + mw**2 )

      caa = -cbb * ( (ma**2-mz**2)/(mhh**2-mhl**2) )      
      sa  = -sqrt((1.D0 - caa)/2.D0)
      ca  =  sqrt((1.D0 + caa)/2.D0)

c$$$      caa  = 2.D0*ca**2 - 1.D0
c$$$      saa  = 2.D0*sa*ca
c$$$      sab  = sa*cb + ca*sb
c$$$      cab  = ca*cb - sa*sb
c$$$      cbma = cb*ca + sb*sa
c$$$      sbma = sb*ca - cb*sa    

      return
      end 
