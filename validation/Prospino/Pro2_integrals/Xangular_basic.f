cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c all the integrals needed, named like bible, but factor Pi missing    c
c      [but this factor included in the wrapping ANG*** routines]      c
c                                                                      c
c note that the powers in the name for the functions correspond to     c
c  the indices in the bible, negative powers are counted positive !    c
c                                                                      c
c n.b. the cutoff in dum2 makes a big difference,                      c
c      it should not be too small                                      c
c                                                                      c
c added [tp]: A4P2P0, A4P2M1, A4P1M1                                   c
c                                                                      c
c additional imaginary parts of integrals A4** -> I4**                 c
c additional integrals for one dimension X4**                          c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function ANG2_EXT(n1,n2,arg1,arg2)

      implicit none 
      
      integer n1,n2 
      real*8  arg1(1:3),arg2(1:3),pi
      real*8  A4P0P0_NEW
      real*8  A4M2P1_NEW
      real*8  A4M1P0_NEW, A4M1P1_NEW, A4M1P2_NEW
      real*8  A4P1M2_NEW, A4P1M1_NEW, A4P1P0_NEW
      real*8  A4P1P1_NEW, A4P1P2_NEW
      real*8  A4P2M2_NEW, A4P2M1_NEW, A4P2P0_NEW
      real*8  A4P2P1_NEW, A4P2P2_NEW

      call INITIALIZE_ANGULAR_ZERO
      pi = 4.D0 * atan(1.D0)

      if      ((n1.eq. 0).and.(n2.eq. 0)) then 
         ANG2_EXT = pi*A4P0P0_NEW(arg1,arg2)

      else if ((n1.eq. 2).and.(n2.eq.-1)) then 
         ANG2_EXT = pi*A4M2P1_NEW(arg1,arg2)

      else if ((n1.eq. 1).and.(n2.eq. 0)) then 
         ANG2_EXT = pi*A4M1P0_NEW(arg1,arg2)
      else if ((n1.eq. 1).and.(n2.eq.-1)) then 
         ANG2_EXT = pi*A4M1P1_NEW(arg1,arg2)
      else if ((n1.eq. 1).and.(n2.eq.-2)) then 
         ANG2_EXT = pi*A4M1P2_NEW(arg1,arg2)

      else if ((n1.eq.-1).and.(n2.eq. 2)) then 
         ANG2_EXT = pi*A4P1M2_NEW(arg1,arg2)
      else if ((n1.eq.-1).and.(n2.eq. 1)) then 
         ANG2_EXT = pi*A4P1M1_NEW(arg1,arg2)
      else if ((n1.eq.-1).and.(n2.eq. 0)) then 
         ANG2_EXT = pi*A4P1P0_NEW(arg1,arg2)
      else if ((n1.eq.-1).and.(n2.eq.-1)) then 
         ANG2_EXT = pi*A4P1P1_NEW(arg1,arg2)
      else if ((n1.eq.-1).and.(n2.eq.-2)) then 
         ANG2_EXT = pi*A4P1P2_NEW(arg1,arg2)

      else if ((n1.eq.-2).and.(n2.eq. 2)) then 
         ANG2_EXT = pi*A4P2M2_NEW(arg1,arg2)
      else if ((n1.eq.-2).and.(n2.eq. 1)) then 
         ANG2_EXT = pi*A4P2M1_NEW(arg1,arg2)
      else if ((n1.eq.-2).and.(n2.eq. 0)) then 
         ANG2_EXT = pi*A4P2P0_NEW(arg1,arg2)
      else if ((n1.eq.-2).and.(n2.eq.-1)) then 
         ANG2_EXT = pi*A4P2P1_NEW(arg1,arg2)
      else if ((n1.eq.-2).and.(n2.eq.-2)) then 
         ANG2_EXT = pi*A4P2P2_NEW(arg1,arg2)

      else
         print*," ANG2_EXT: function not implemented ",n1,n2
         call HARD_STOP
      end if 

      return
      end

c ----------------------------------------------------------------------
      real*8 function AIM2_EXT(n1,n2,arg1,arg2)

      implicit none 
      
      integer n1,n2 
      real*8  arg1(1:3),arg2(1:3),pi
      real*8  K4M1P1_NEW, K4P1M1_NEW, K4P1P0_NEW, K4P1P1_NEW, K4M2P1_NEW

      call INITIALIZE_ANGULAR_ZERO
      pi = 4.D0 * atan(1.D0)

      if      ((n1.eq. 2).and.(n2.eq.-1)) then 
         AIM2_EXT = pi*K4M2P1_NEW(arg1,arg2)
      else if ((n1.eq. 1).and.(n2.eq.-1)) then 
         AIM2_EXT = pi*K4M1P1_NEW(arg1,arg2)
      else if ((n1.eq.-1).and.(n2.eq. 1)) then 
         AIM2_EXT = pi*K4P1M1_NEW(arg1,arg2)
      else if ((n1.eq.-1).and.(n2.eq. 0)) then 
         AIM2_EXT = pi*K4P1P0_NEW(arg1,arg2)
      else if ((n1.eq.-1).and.(n2.eq.-1)) then 
         AIM2_EXT = pi*K4P1P1_NEW(arg1,arg2)

      else
         print*," AIM2_EXT: function not implemented ",n1,n2
         call HARD_STOP
      end if 

      return
      end

c ----------------------------------------------------------------------
      real*8 function ANG1_EXT(n1,arg_x,arg_y)

      implicit none 

      integer n1
      real*8  arg_x,arg_y,pi
      real*8  X4M2_NEW,X4M1_NEW,X4P0_NEW,X4P1_NEW,X4P2_NEW

      call INITIALIZE_ANGULAR_ZERO
      pi = 4.D0 * atan(1.D0)

      if      (n1.eq. 2) then
         ANG1_EXT = pi*X4M2_NEW(arg_x,arg_y)
      else if (n1.eq. 1) then
         ANG1_EXT = pi*X4M1_NEW(arg_x,arg_y)
      else if (n1.eq. 0) then
         ANG1_EXT = pi*X4P0_NEW(arg_x,arg_y)
      else if (n1.eq.-1) then  
         ANG1_EXT = pi*X4P1_NEW(arg_x,arg_y)
      else if (n1.eq.-2) then  
         ANG1_EXT = pi*X4P2_NEW(arg_x,arg_y)

      else 
         print*," ANG1_EXT: function not implemented ",n1
         call HARD_STOP
      end if 

      return 
      end 

c ----------------------------------------------------------------------
      subroutine INITIALIZE_ANGULAR_ZERO

      implicit none 

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      zero_2 = 1.D-16
      zero_3 = 1.D-10

      return 
      end 

c ----------------------------------------------------------------------
      real*8 function X4P2_NEW(X,Y)
      real*8 X,Y

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if ( (x.gt.0.D0).or.(x**2.lt.y**2) ) then 
         print*," X4P2: function not called correctly ",x,y
         call HARD_STOP 
      end if 

      X4P2_NEW = - X/(X**2-Y**2)**(3.D0/2.D0)

      return 
      end

c ----------------------------------------
      real*8 function X4P1_NEW(X,Y)
      real*8 X,Y

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if ( (x.gt.0.D0).or.(x**2.lt.y**2) ) then 
         print*," X4P1: function not called correctly ",x,y
         call HARD_STOP 
      end if 

      X4P1_NEW = - 1.D0/(X**2-Y**2)**(1.D0/2.D0)

      return 
      end

c ----------------------------------------
      real*8 function X4P0_NEW(X,Y)
      real*8 X,Y

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      X4P0_NEW = 1.D0

      return 
      end

c ----------------------------------------
      real*8 function X4M1_NEW(X,Y)
      real*8 X,Y

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if ( (x.gt.0.D0).or.(x**2.lt.y**2) ) then 
         print*," X4M1: function not called correctly ",x,y
         call HARD_STOP 
      end if 

      X4M1_NEW = X

      return 
      end

c ----------------------------------------
      real*8 function X4M2_NEW(X,Y)
      real*8 X,Y

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if ( (x.gt.0.D0).or.(x**2.lt.y**2) ) then 
         print*," X4M2: function not called correctly ",x,y
         call HARD_STOP 
      end if 

      X4M2_NEW = X**2 + Y**2/2.D0

      return 
      end


c----------------------------------------------------------------------
      real*8 function A4P0P0_NEW(a1,a2)
c      real*8  a,b,c,au,bu,cu
      real*8  a1(1:3),a2(1:3)
      integer n1

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      do n1=1,3
         if (a1(n1).ne.0.D0) print*, " A4P0P0: error zero=",a1(n1),n1
         if (a2(n1).ne.0.D0) print*, " A4P0P0: error zero=",a2(n1),n1
      end do

      if (a1(3).ne.0.D0) print*, " A4P0P0: error a1(3)=",a1(3)
c      a  = a1(1)
c      b  = a1(2)
c      c  = a1(3)
c      au = a2(1)
c      bu = a2(2)
c      cu = a2(3)

      A4P0P0_NEW = 2.D0

      return
      end

c ----------------------------------------
      real*8 function A4P1M1_NEW(a1,a2)
      real*8 a,b,au,bu,a1(1:3),a2(1:3)
      real*8 dum2

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P1M1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)

      dum2 = abs(a+b)
      if (dum2.lt.zero_2) then 
         A4P1M1_NEW = -2.D0*BU/A
      else 
ctp         print*, " A4P1M1: that would be new and not yet checked "
         A4P1M1_NEW = 
     +        2.D0*BU/B 
     +       + (AU*B - A*BU)/B**2 * LOG(abs((A+B)/(A-B)))
      end if 

      return
      end

c ----------------------------------------
      real*8 function A4P2M1_NEW(a1,a2)
      real*8 a,b,au,bu,a1(1:3),a2(1:3)
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P2M1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)

      A4P2M1_NEW = 
     +     2.D0*(AU*B - A*BU)/b/(A**2-B**2)
     +    + BU/B**2 * LOG(abs((A+B)/(A-B)))
      return
      end

c ----------------------------------------
      real*8 function A4P2P0_NEW(a1,a2)
      real*8 a,b,a1(1:3),a2(1:3)
      real*8 dum2
      integer n1

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      do n1=1,3
         if (a2(n1).ne.0.D0) print*, " A4P2P0: error zero=",a2(n1),n1
      end do

      if (a1(3).ne.0.D0) print*, " A4P2P0: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)

      dum2 = abs(a+b)
      if (dum2.lt.zero_2) then 
         A4P2P0_NEW = -1.D0/A**2 
      else 
         A4P2P0_NEW = 2.D0/( A**2 - B**2 )
      end if 

      return
      end

c ----------------------------------------
      real*8 function A4P1P2_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),X,Y
      real*8 dum2,dum3 
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P1P2: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      dum2 = abs(a+b)
      dum3 = abs(au**2-bu**2-cu**2)/au**2
      if ((dum2.lt.zero_2).and.(dum3.lt.zero_3)) then 
         print*, " A4P1P2: error dum3=",dum3 
      else if ((dum2.lt.zero_2).and.(dum3.gt.zero_3)) then 
         X = AU +BU
         Y = BU**2 +CU**2
         A4P1P2_NEW =  1.D0/A/X**2
     +     * ( LOG(abs(X**2/(AU**2 -Y))) +2.D0*(Y +AU*BU)/(AU**2 -Y) )
      else 
         X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
         A4P1P2_NEW = 
     +     (2.D0*A*(BU**2+CU**2) -2.D0*B*AU*BU)/(AU**2-BU**2-CU**2)/X
     +     + B*(B*AU -A*BU)/X**(3.D0/2.D0) 
     +     *LOG(abs((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X))))
      end if
         
      return
      end

c ----------------------------------------
      real*8 function A4P2P1_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),X,Y
      real*8 dum2
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P2P1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      dum2 = abs(a+b)
      if (dum2.lt.zero_2) then 
         X = AU +BU
         Y = BU**2 +CU**2
         A4P2P1_NEW = 1.D0/A**2/X
     +    * ((Y +AU*BU)/X**2*LOG(abs(X**2/(AU**2-Y))) 
     +  - 2.D0*CU**2/X**2 -1.D0)
      else 
         X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
         A4P2P1_NEW = 
     +     2*B*(B*AU -A*BU)/(A**2 -B**2)/X
     +    + (A*(BU**2 +CU**2)-B*AU*BU)/X**(3.D0/2.D0)
     +    *LOG(abs((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X))))
      end if 

      return
      end

c ----------------------------------------
      real*8 function A4P1P1_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),X
      real*8 dum2,dum3 
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P1P1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      dum2 = abs(a+b)
      dum3 = abs(au**2-bu**2-cu**2)/au**2
      if ((dum2.lt.zero_2).and.(dum3.lt.zero_3)) then 
         A4P1P1_NEW = -2.D0/A/(AU +BU)*LOG(abs(2.D0*AU/(AU +BU)))
      else if ((dum2.lt.zero_2).and.(dum3.ge.zero_3)) then 
         A4P1P1_NEW = LOG(abs((AU +BU)**2/(AU**2-BU**2-CU**2)))
     +                /A/(AU +BU)
      else if ((dum2.gt.zero_2).and.(dum3.lt.zero_3)) then 
         print*, " A4P1P1: error dum2, dum3 ",dum2,dum3
         call HARD_STOP
      else 
         X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
         A4P1P1_NEW = 
     +     LOG(abs((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X))))
     +     /SQRT(X) 
      end if 

      return
      end

c ----------------------------------------
      real*8 function A4P1P0_NEW(a1,a2)
      real*8 a,b,a1(1:3),a2(1:3)
      integer n1
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      do n1=1,3
         if (a2(n1).ne.0.D0) print*, " A4P1P0: error zero=",a2(n1),n1
      end do

      if (a1(3).ne.0.D0) print*, " A4P1P0: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)

      A4P1P0_NEW = 1.D0/B*LOG(abs((A +B)/(A-B)))

      return
      end

c ----------------------------------------
      real*8 function A4M1P0_NEW(a1,a2)
      real*8 a,a1(1:3),a2(1:3)
      integer n1
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      do n1=1,3
         if (a2(n1).ne.0.D0) print*, " A4M1P0: error zero=",a2(n1),n1
      end do

      if (a1(3).ne.0.D0) print*, " A4M1P0: error a1(3)=",a1(3)
      a  = a1(1)

      A4M1P0_NEW = 2.D0*A

      return
      end

c ----------------------------------------
      real*8 function A4M1P1_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),Y
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4M1P1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      Y = BU**2 + CU**2
      A4M1P1_NEW = 
     +     2.D0*B*BU/Y
     +    + (A*Y -B*AU*BU)/Y**(3.D0/2.D0)
     +    * LOG(abs((AU +SQRT(Y))/(AU -SQRT(Y))))

      return
      end

c ----------------------------------------
      real*8 function A4M2P1_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),Y
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4M2P1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      Y = BU**2 + CU**2
      A4M2P1_NEW = 
     +     4.D0*A*B*BU/Y
     +    + B**2*AU*(CU**2 -2.D0*BU**2)/Y**2
     +    + ( (A*Y -B*AU*BU)**2/Y**(5.D0/2.D0)
     +    -B**2*CU**2*(AU**2 -BU**2 -CU**2)/2.D0/Y**(5.D0/2.D0) )
     +    * LOG(abs((AU +SQRT(Y))/(AU -SQRT(Y))))

      return
      end

c ----------------------------------------
      real*8 function A4P2M2_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3)
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P2M2: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      A4P2M2_NEW = 
     +     2.D0*(BU**2 -CU**2)/B**2
     +    + 2.D0*(B*AU -A*BU)**2/B**2/(A**2 -B**2)
     +    + (A*CU**2 +2.D0*BU*(B*AU -A*BU))/B**3
     +    * LOG(abs((A + B)/(A-B)))

      return
      end

c ----------------------------------------
      real*8 function A4M1P2_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),Y

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4M1P2: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      Y = BU**2 +CU**2
      
      A4M1P2_NEW = 
     +     2.D0*(A*Y - B*AU*BU)/Y/(AU**2 -Y)
     +    + B*BU/Y**(3.D0/2.D0)
     +     * LOG(abs((AU +SQRT(Y))/(AU -SQRT(Y))))

      return
      end

c ----------------------------------------
      real*8 function A4P0P2_NEW(a1,a2)
      real*8 au,bu,cu,a1(1:3),a2(1:3)
      integer n1
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      do n1=1,3
         if (a1(n1).ne.0.D0) print*, " A4P0P2: error zero=",a1(n1),n1
      end do

      if (a1(3).ne.0.D0) print*, " A4P0P2: error a1(3)=",a1(3)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      A4P0P2_NEW = 2.D0/(AU**2 -BU**2 -CU**2 )

      return
      end

c ----------------------------------------
      real*8 function A4P2P2_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),X,Y
      real*8 dum2
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P2P2: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      dum2 = abs(a+b)
      if (dum2.lt.zero_2) then 
         X = AU +BU
         Y = BU**2 +CU**2
         A4P2P2_NEW = 1.D0/A**2/X**2 * (
     +          (3.D0*CU**2/X**2 + 2.D0*BU/X)
     +           * LOG(abs(X**2/(AU**2 -Y)))
     +          - 8.D0*CU**2/X**2 + 2.D0*Y/(AU**2 -Y) -1.D0 )
      else 
         X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
         Y = BU**2 +CU**2
         A4P2P2_NEW = 
     +     2.D0*B**2/(A**2 -B**2)/X + 2.D0*Y/(AU**2 -Y)/X
     +     - 6.D0*B**2*CU**2/X**2
     +    +( B*BU/X**(3.D0/2.D0)
     +    + 3.D0*B*(B*AU -A*BU)*(A*Y -B*AU*BU)/X**(5.D0/2.D0))
     +    *LOG(abs((A*AU -B*BU + SQRT(X))/(A*AU -B*BU - SQRT(X))))
      end if 

      return
      end

c ----------------------------------------
      real*8 function A4P1M2_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3)
      real*8 dum2

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P1M2: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      dum2 = abs(a+b)
      if (dum2.lt.zero_2) then 
         A4P1M2_NEW = (CU**2 -4.D0*AU*BU -2.D0*BU**2)/A
      else 
         print*, " A4P1M2: integral not implemented "
      endif 

      return
      end

c----------------------------------------------------------------------
C***  HERE WE USE A WIDTH OF THE PARTICLE del=gam^2*m^2 
      real*8 function CBP1P1_NEW(del,a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),del
      real*8 dum2,dum3

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " A4P0P0: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)
 
      dum2 = abs(a+b)
      dum3 = abs(au**2-bu**2-cu**2)/au**2
      if ((dum2.lt.zero_2).and.(dum3.lt.zero_3)) then 
         print*, " CBP1P1: integral not implemented (1) "
      else if ((dum2.lt.zero_2).and.(dum3.ge.zero_3)) then 
         CBP1P1_NEW = LOG(ABS((AU +BU)**2/(AU**2 -BU**2 -CU**2)))/A *
     +           (AU +BU)/((AU +BU)**2 +DEL)
      else 
         print*, " CBP1P1: integral not implemented (2) "
      end if 



      return
      end

c----------------------------------------------------------------------
      real*8 function K4M2P1_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),Y
      COMPLEX*16 KA,KDEL,KONE
      
      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " K4M2P1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      KDEL = (0.D0,1.D-16)
      KONE = (0.D0,-1.D0)
      KA = AU - KDEL

      Y = BU**2 + CU**2
      K4M2P1_NEW = 
     +      ( (A*Y -B*AU*BU)**2/Y**(5.D0/2.D0)
     +    -B**2*CU**2*(AU**2 -BU**2 -CU**2)/2.D0/Y**(5.D0/2.D0) )
     +    * REAL( KONE * LOG( (KA +SQRT(Y))/(KA -SQRT(Y))) )

      return
      end

c ----------------------------------------
      real*8 function K4P1M1_NEW(a1,a2)
      real*8 a,b,au,bu,a1(1:3),a2(1:3)
      COMPLEX*16 KA,KDEL,KONE

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " K4P1M1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)

      KDEL = (0.D0,1.D-16)
      KONE = (0.D0,-1.D0)
      KA = AU - KDEL
      K4P1M1_NEW = 
     +    + (AU*B -A*BU)/B**2
     +    * REAL( KONE * LOG((KA +B)/(KA -B)) )

      return
      end


c ----------------------------------------
      real*8 function K4P1P1_NEW(a1,a2)
C***  THE IMAGINARY PART OF A COMPLEX INTEGRAL
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),X
      COMPLEX*16 KA,KDEL,KONE

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " K4P1P1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      KDEL = (0.D0,1.D-16)
      KONE = (0.D0,-1.D0)
      KA = A - KDEL
      X = (AU*B)**2 +(BU*A)**2 +(CU*A)**2 -(CU*B)**2-2.D0*A*B*AU*BU
      K4P1P1_NEW = REAL( KONE * 
     +     LOG((KA*AU -B*BU + SQRT(X))/(KA*AU -B*BU - SQRT(X)) )
     +     /SQRT(X) )

      return
      end


c ----------------------------------------
      real*8 function K4M1P1_NEW(a1,a2)
      real*8 a,b,au,bu,cu,a1(1:3),a2(1:3),Y
      COMPLEX*16 KAU,KDEL,KONE

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      if (a1(3).ne.0.D0) print*, " K4M1P1: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)
      au = a2(1)
      bu = a2(2)
      cu = a2(3)

      KDEL = (0.D0,1.D-16)
      KONE = (0.D0,-1.D0)
      KAU = AU - KDEL
      Y = BU**2 + CU**2
      K4M1P1_NEW = 
     +    + (A*Y -B*AU*BU)/Y**(3.D0/2.D0)
     +    * REAL( KONE * LOG((KAU +SQRT(Y))/(KAU -SQRT(Y))) )

      return
      end


c ----------------------------------------
      real*8 function K4P1P0_NEW(a1,a2)
      real*8 a,b,a1(1:3),a2(1:3)
      COMPLEX*16 KA,KDEL,KONE
      integer n1

      real*8                      zero_2,zero_3
      common/ANGULAR_BASIC_CUTOFF/zero_2,zero_3

      do n1=1,3
         if (a2(n1).ne.0.D0) print*, " K4P1P0: error zero=",a2(n1),n1
      end do

      if (a1(3).ne.0.D0) print*, " K4P1P0: error a1(3)=",a1(3)
      a  = a1(1)
      b  = a1(2)

      KDEL = (0.D0,1.D-16)
      KONE = (0.D0,-1.D0)
      KA = A - KDEL
      K4P1P0_NEW = 1.D0/B* REAL(KONE * LOG((KA +B)/(KA-B)))

      return
      end



































