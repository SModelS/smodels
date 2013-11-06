cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  LUMI(inlo,ix,icoll,idub,iq,x1,x2,mu)                                c 
c           for completely exclusive cross sections :                  c
c                                                                      c
c  inlo  = 0,1 for LO or NLO densities                                 c
c                                                                      c
c  icoll = 0 : tevatron (p,pbar collider)                              c
c          1 : lhc 14 tev (p,p collider)                               c
c          2 : lhc  7 tev (p,p collider)                               c
c          3 : lhc  8 tev (p,p collider)                               c
c          4 : LUMI=1.0, only scaling function!!!                      c
c                                                                      c
c  idub = 0 ; iq = +1  (u,ubar) or (u,g) or (g,ubar) + perm.           c
c             iq = -1  (d,dbar) or (d,g) or (g,dbar) + perm.           c
c   i.e. iq the quantum number of the quark involved [unique]          c
c                                                                      c
c  idub = 1 ; iq = +1  (u,dbar) or (u,g) or (g,dbar) + perm.           c
c             iq = -1  (d,ubar) or (d,g) or (g,ubar) + perm.           c
c   i.e. iq the charge of the qb state, gluon replaces one quark       c
c                                                                      c
c  ix = 20 : [ix=21 + ix=22]                                           c
c  ix = 21 : q(x1)*b(x2)            ix = 22 : b(x1)*q(x2)              c
c  ix = 30 : [ix=31 + ix=32]                                           c
c  ix = 31 : q(x1)*g(x2)            ix = 32 : g(x1)*q(x2)              c
c  ix = 40 : [ix=41 + ix=42]                                           c
c  ix = 41 : g(x1)*b(x2)            ix = 42 : b(x1)*g(x2)              c
c  ix = 50 : [ix=51 + ix=52]                                           c
c  ix = 51 : g(x1)*g(x2)/2          ix = 52 : g(x1)*g(x2)/2            c
c  ix = 60 : [ix=61 + ix=62] like 20, but bottoms only                 c
c  ix = 61 : q(x1)*b(x2)            ix = 62 : b(x1)*q(x2)              c
c  ix = 70 : [ix=71 + ix=72] like 30, but bottoms only                 c
c  ix = 71 : q(x1)*g(x2)            ix = 72 : g(x1)*q(x2)              c
c  ix = 80 : [ix=81 + ix=82] like 40, but bottoms only                 c
c  ix = 81 : g(x1)*b(x2)            ix = 82 : b(x1)*g(x2)              c
c  ix = 90 : [ix=91 + ix=92] bottoms only                              c
c  ix = 91 : q(x1)*q(x2)            ix = 92 ; b(x1)*b(x2)              c
c  ix = 100: [ix=101 + ix=102], one bottom                             c
c  ix = 101: q(x1)*b(x2)            ix = 102: b(x1)*q(x2)              c
c  ix = 110: [ix=111 + ix=112], one bottom                             c
c  ix = 111: qB(x1)*q(x2)            ix = 112: q(x1)*qB(x2)            c
c  ix = 120: [ix=121 + ix=122], one bottom                             c
c  ix = 121: bB(x1)*b(x2)            ix = 122: b(x1)*bB(x2)            c
c                                                                      c
c  ix   = 2*    ; idub = 0 ; iq = +1 : (u1,ubar2 + ubar1,u2)           c
c                 idub = 0 ; iq = -1 : (d1,dbar2 + dbar1,d2)           c
c                 idub = 1 ; iq = +1 : (u1,dbar2 + dbar1,u2)           c
c                 idub = 1 ; iq = -1 : (d1,ubar2 + ubar1,d2)           c
c                                                                      c
c  ix   = 3*    ;            iq = +1 : (u1,g2 + g1,u2)                 c
c                            iq = -1 : (d1,g2 + g1,d2)                 c
c                                                                      c
c  ix   = 4*    ; idub = 0 ; iq = +1 : (g1,ubar2 + ubar1,g2)           c
c                 idub = 0 ; iq = -1 : (g1,dbar2 + dbar1,g2)           c
c                 idub = 1 ; iq = +1 : (g1,dbar2 + dbar1,g2)           c
c                 idub = 1 ; iq = -1 : (g1,ubar2 + ubar1,g2)           c
c                                                                      c
c note : iq = 2*t3(iq) = 'charge of the outgoing state'                c
c        cteq4 parton densities with 1,2,3,4,5 as u,d,s,c,b            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function LUMI(inlo,ix,icoll,idub,iq,x1,x2,mu)

      implicit none

c               n = dummy variable in this case      
      integer ix,icoll,idub,iq,inlo 
      real*8  x1,x2,mu,pdf1(-6:6),pdf2(-6:6)
      real*8  u1q,d1q,s1q,c1q,b1q,g1x
      real*8  u1b,d1b,s1b,c1b,b1b
      real*8  u2q,d2q,s2q,c2q,b2q,g2x
      real*8  u2b,d2b,s2b,c2b,b2b

c               for scaling functions only
      if (icoll.eq.4) then
         LUMI = 1.D0
         return 
      end if 

c ---------------------------------------------------------------
c               define the parton densities, pdf array like CTEQ
      call GET_PDF(inlo,x1,mu,pdf1)
      call GET_PDF(inlo,x2,mu,pdf2)

c               assign internal structure for Tevatron (p-pbar)
      if (icoll.eq.0) then 
         u1q = pdf1( 1)
         u1b = pdf1(-1)
         d1q = pdf1( 2)
         d1b = pdf1(-2)
         s1q = pdf1( 3)
         s1b = pdf1(-3)
         c1q = pdf1( 4)
         c1b = pdf1(-4)
         b1q = pdf1( 5)
         b1b = pdf1(-5)
         g1x = pdf1( 0)
         
         u2b = pdf2( 1)
         u2q = pdf2(-1)
         d2b = pdf2( 2)
         d2q = pdf2(-2)
         s2b = pdf2( 3)
         s2q = pdf2(-3)
         c2b = pdf2( 4)
         c2q = pdf2(-4)
         b2b = pdf2( 5)
         b2q = pdf2(-5)
         g2x = pdf2( 0)

c               the same for the LHC
      else if ( (icoll.eq.1).or.(icoll.eq.2).or.(icoll.eq.3) ) then 
         u1q = pdf1( 1)
         u1b = pdf1(-1)
         d1q = pdf1( 2)
         d1b = pdf1(-2)
         s1q = pdf1( 3)
         s1b = pdf1(-3)
         c1q = pdf1( 4)
         c1b = pdf1(-4)
         b1q = pdf1( 5)
         b1b = pdf1(-5)
         g1x = pdf1( 0)

         u2q = pdf2( 1)
         u2b = pdf2(-1)
         d2q = pdf2( 2)
         d2b = pdf2(-2)
         s2q = pdf2( 3)
         s2b = pdf2(-3)
         c2q = pdf2( 4)
         c2b = pdf2(-4)
         b2q = pdf2( 5)
         b2b = pdf2(-5)
         g2x = pdf2( 0)
      end if 

c ---------------------------------------------------------------
c               (q-qbar)+(qbar-q)   for idub=0, charge=0
c               (q-qbar')+(qbar-q') for idub=1, charge=iq 
      if (ix.eq.20) then 

         if (idub.eq.0) then 
            if (iq.eq.1) then 
               LUMI = u1q*u2b + c1q*c2b
     &              + u1b*u2q + c1b*c2q
            else if (iq.eq.-1) then 
               LUMI = d1q*d2b + s1q*s2b + 0*b1q*b2b
     &              + d1b*d2q + s1b*s2q + 0*b1b*b2q
            end if 
         else if (idub.eq.1) then 
            if (iq.eq.1) then 
               LUMI = u1q*d2b + c1q*s2b
     &              + d1b*u2q + s1b*c2q
            else if (iq.eq.-1) then 
               LUMI = d1q*u2b + s1q*c2b
     &              + u1b*d2q + c1b*s2q
            end if 
         end if 

      else if (ix.eq.21) then 

         if (idub.eq.0) then 
            if (iq.eq.1) then 
               LUMI =  u1q*u2b + c1q*c2b
            else if (iq.eq.-1) then 
               LUMI =  d1q*d2b + s1q*s2b + 0*b1q*b2b
            end if 
         else if (idub.eq.1) then 
            if (iq.eq.1) then 
               LUMI = u1q*d2b + c1q*s2b
            else if (iq.eq.-1) then 
               LUMI = d1q*u2b + s1q*c2b
            end if 
         end if 

      else if (ix.eq.22) then 
         
         if (idub.eq.0) then 
            if (iq.eq.1) then 
               LUMI = u1b*u2q + c1b*c2q
            else if (iq.eq.-1) then 
               LUMI=  d1b*d2q + s1b*s2q + 0*b1b*b2q
            end if 
         else if (idub.eq.1) then 
            if (iq.eq.1) then 
               LUMI = d1b*u2q + s1b*c2q
            else if (iq.eq.-1) then 
               LUMI = u1b*d2q + c1b*s2q
            end if 
         end if 

c ---------------------------------------------------------------
c               (q-g)+(q-g) incoming state 
      else if (ix.eq.30) then 

         if (iq.eq.1) then 
            LUMI = u1q*g2x + c1q*g2x
     &           + g1x*u2q + g1x*c2q
         else if (iq.eq.-1) then 
            LUMI = d1q*g2x + s1q*g2x + 0*b1q*g2x
     &           + g1x*d2q + g1x*s2q + 0*g1x*b2q
         end if 
            
      else if (ix.eq.31) then 
         
         if (iq.eq.1) then 
            LUMI = u1q*g2x + c1q*g2x
         else if (iq.eq.-1) then 
            LUMI = d1q*g2x + s1q*g2x + 0*b1q*g2x
         end if 
            
      else if (ix.eq.32) then 

         if (iq.eq.1) then 
            LUMI = g1x*u2q + g1x*c2q
         else if (iq.eq.-1) then 
            LUMI = g1x*d2q + g1x*s2q + 0*g1x*b2q
         end if 
            
c ---------------------------------------------------------------
c               (g-qbar)+(qbar-g);(g-qbar')+(qbar'-g) incoming state 
      else if (ix.eq.40) then 

         if (idub.eq.0) then 
            if (iq.eq.1) then 
               LUMI = u1b*g2x + c1b*g2x
     &              + g1x*u2b + g1x*c2b
            else if (iq.eq.-1) then 
               LUMI = d1b*g2x + s1b*g2x + 0*b1b*g2x
     &              + g1x*d2b + g1x*s2b + 0*g1x*b2b
            end if 

         else if (idub.eq.1) then 
            if (iq.eq.1) then 
               LUMI = d1b*g2x + s1b*g2x + 0*b1b*g2x
     &              + g1x*d2b + g1x*s2b + 0*g1x*b2b
            else if (iq.eq.-1) then 
               LUMI = u1b*g2x + c1b*g2x
     &              + g1x*u2b + g1x*c2b
            end if 
         end if

      else if (ix.eq.41) then 

         if (idub.eq.0) then 
            if (iq.eq.1) then 
               LUMI = g1x*u2b + g1x*c2b
            else if (iq.eq.-1) then 
               LUMI = g1x*d2b + g1x*s2b + 0*g1x*b2b
            end if 
         else if (idub.eq.1) then 
            if (iq.eq.1) then 
               LUMI = g1x*d2b + g1x*s2b + 0*g1x*b2b
            else if (iq.eq.-1) then 
               LUMI=  g1x*u2b + g1x*c2b
            end if 
         end if 

      else if (ix.eq.42) then 

         if (idub.eq.0) then 
            if (iq.eq.1) then 
               LUMI = u1b*g2x + c1b*g2x
            else if (iq.eq.-1) then 
               LUMI = d1b*g2x + s1b*g2x + 0*b1b*g2x
            end if 
         else if (idub.eq.1) then 
            if (iq.eq.1) then 
               LUMI = d1b*g2x + s1b*g2x + 0*b1b*g2x
            else if (iq.eq.-1) then 
               LUMI = u1b*g2x + c1b*g2x
            end if 
         end if 

c -------------------------------------
c               (g-g) incoming state 
      else if (ix.eq.50) then 
         LUMI = g1x * g2x
      else if (ix.eq.51) then 
         LUMI = g1x * g2x /2.D0
      else if (ix.eq.52) then 
         LUMI = g1x * g2x /2.D0

c ---------------------------------------------------------------
c               (b-bbar)+(bbar-b), bottoms only
      else if (ix.eq.60) then 
         LUMI = b1q*b2b + b1b*b2q
      else if (ix.eq.61) then 
         LUMI = b1q*b2b
      else if (ix.eq.62) then 
         LUMI = b1b*b2q

c ---------------------------------------------------------------
c               (b-g)+(b-g) bottoms only
      else if (ix.eq.70) then 
         LUMI = b1q*g2x + g1x*b2q
      else if (ix.eq.71) then 
         LUMI = b1q*g2x
      else if (ix.eq.72) then 
         LUMI = g1x*b2q
            
c ---------------------------------------------------------------
c               (g-bbar)+(bbar-g) bottoms only
      else if (ix.eq.80) then 
         LUMI = b1b*g2x + g1x*b2b
      else if (ix.eq.81) then 
         LUMI = g1x*b2b
      else if (ix.eq.82) then 
         LUMI = b1b*g2x

c ---------------------------------------------------------------
c               (b-b)+(bbar-bbar), bottoms only
      else if (ix.eq.90) then 
         LUMI = b1q*b2q + b1b*b2b
      else if (ix.eq.91) then 
         LUMI = b1q*b2q
      else if (ix.eq.92) then 
         LUMI = b1b*b2b

c ---------------------------------------------------------------
c               (b-qbar)+(bbar-q), one bottom
      else if (ix.eq.100) then 

         if (iq.eq.1) then 
            LUMI = b1q*u2b + b1q*c2b
     &           + b2q*u1b + b2q*c1b
     &           + u1q*b2b + c1q*b2b
     &           + u2q*b1b + c2q*b1b
         else if (iq.eq.-1) then 
            LUMI = b1q*d2b + b1q*s2b
     &           + b2q*d1b + b2q*s1b
     &           + d1q*b2b + s1q*b2b
     &           + d2q*b1b + s2q*b1b
         end if 
            
      else if (ix.eq.101) then 

         if (iq.eq.1) then 
            LUMI = b1q*u2b + b1q*c2b
     &           + b2q*u1b + b2q*c1b
         else if (iq.eq.-1) then 
            LUMI = b1q*d2b + b1q*s2b
     &           + b2q*d1b + b2q*s1b
         end if 

      else if (ix.eq.102) then 
         
         if (iq.eq.1) then 
            LUMI = u1q*b2b + c1q*b2b
     &           + u2q*b1b + c2q*b1b
         else if (iq.eq.-1) then 
            LUMI = d1q*b2b + s1q*b2b
     &           + d2q*b1b + s2q*b1b
         end if 
            
c ---------------------------------------------------------------
c               (b-q), one bottom
      else if (ix.eq.110) then 

         if (iq.eq.1) then 
            LUMI = b1q*u2q + b1q*c2q
     &           + b2q*u1q + b2q*c1q
         else if (iq.eq.-1) then 
            LUMI = b1q*d2q + b1q*s2q
     &           + b2q*d1q + b2q*s1q
         end if 
            
      else if (ix.eq.111) then 

         if (iq.eq.1) then 
            LUMI = b1q*u2q + b1q*c2q
         else if (iq.eq.-1) then 
            LUMI = b1q*d2q + b1q*s2q
         end if 

      else if (ix.eq.112) then 
         
         if (iq.eq.1) then 
            LUMI = b2q*u1q + b2q*c1q
         else if (iq.eq.-1) then 
            LUMI = b2q*d1q + b2q*s1q
         end if 

c ---------------------------------------------------------------
c               (bbar-qbar), one bottom
      else if (ix.eq.120) then 

         if (iq.eq.1) then 
            LUMI = u1b*b2b + c1b*b2b
     &           + u2b*b1b + c2b*b1b
         else if (iq.eq.-1) then 
            LUMI = d1b*b2b + s1b*b2b
     &           + d2b*b1b + s2b*b1b
         end if 
            
      else if (ix.eq.121) then 

         if (iq.eq.1) then 
            LUMI = u2b*b1b + c2b*b1b
         else if (iq.eq.-1) then 
            LUMI = d2b*b1b + s2b*b1b
         end if 
            
      else if (ix.eq.122) then 
         
         if (iq.eq.1) then 
            LUMI = u1b*b2b + c1b*b2b
         else if (iq.eq.-1) then 
            LUMI = d1b*b2b + s1b*b2b
         end if 
            
c -------------------------------------
      else
         print*, " LUMI: ix not set correctly, set to zero ",ix
         LUMI = 0.D0
      end if 

      return 
      end








