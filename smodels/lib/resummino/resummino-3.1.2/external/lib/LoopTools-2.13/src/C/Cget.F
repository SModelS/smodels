* Cget.F
* the three-point tensor coefficients
* this file is part of LoopTools
* improvements by M. Rauch
* last modified 20 Jan 16 th

#include "externals.h"
#include "types.h"

#define npoint 3
#include "defs.h"


	subroutine XCpara(para, p1, p2, p1p2, m1, m2, m3)
	implicit none
	DVAR para(1,*)
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	P(1) = p1
	P(2) = p2
	P(3) = p1p2

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
	end

************************************************************************

	memindex function XCget(p1, p2, p1p2, m1, m2, m3)
	implicit none
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	memindex cacheindex
	external cacheindex, XCcoefx

#ifdef COMPLEXPARA
	memindex Cget
	external Cget
#endif

	DVAR para(1,Pcc)

#ifdef COMPLEXPARA
	if( abs(Im(p1)) + abs(Im(p2)) + abs(Im(p1p2)) .gt. 0 )
     &    print *, "Complex momenta not implemented"
	if( abs(Im(m1)) + abs(Im(m2)) + abs(Im(m3)) .eq. 0 ) then
	  XCget = Cget(p1, p2, p1p2, m1, m2, m3) - offsetC
	  return
	endif
#endif

	call XCpara(para, p1, p2, p1p2, m1, m2, m3)
	XCget = cacheindex(para, Cval(1,0), XCcoefx, RC*Pcc, Ncc, Cno)
	end

************************************************************************

	subroutine XCput(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType res(*)
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	external XCcoefx

	DVAR para(1,Pcc)

#ifdef COMPLEXPARA
	if( abs(Im(p1)) + abs(Im(p2)) + abs(Im(p1p2)) .gt. 0 )
     &    print *, "Complex momenta not implemented"
	if( abs(Im(m1)) + abs(Im(m2)) + abs(Im(m3)) .eq. 0 ) then
	  call Cput(res, p1, p2, p1p2, m1, m2, m3)
	  return
	endif
#endif

	call XCpara(para, p1, p2, p1p2, m1, m2, m3)
	call cachecopy(res, para, Cval(1,0), XCcoefx, RC*Pcc, Ncc, Cno)
	end

************************************************************************

	subroutine XC0nocache(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType res(*)
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	DVAR para(1,Pcc)

#ifdef COMPLEXPARA
	if( abs(Im(p1)) + abs(Im(p2)) + abs(Im(p1p2)) .gt. 0 )
     &    print *, "Complex momenta not implemented"
	if( abs(Im(m1)) + abs(Im(m2)) + abs(Im(m3)) .eq. 0 ) then
	  call C0nocache(res, p1, p2, p1p2, m1, m2, m3)
	  return
	endif
#endif

	call XCpara(para, p1, p2, p1p2, m1, m2, m3)
	call C0func(res, para)
	end

************************************************************************

	ComplexType function XC0i(i, p1, p2, p1p2, m1, m2, m3)
	implicit none
	integer i
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	memindex XCget
	external XCget

	memindex b

	b = XCget(p1, p2, p1p2, m1, m2, m3)
	XC0i = Cval(i+epsi,b)
	end

************************************************************************

	ComplexType function XC0(p1, p2, p1p2, m1, m2, m3)
	implicit none
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	ComplexType XC0i
	external XC0i

	XC0 = XC0i(cc0, p1, p2, p1p2, m1, m2, m3)
	end

************************************************************************

	subroutine XCcoefx(C, para)
	implicit none
	ComplexType C(*)
	DVAR para(1,*)

#include "lt.h"

	memindex b12, b23, b13

	memindex XBget
	external XBget

	b12 = XBget(P(1), M(1), M(2))
	b23 = XBget(P(2), M(2), M(3))
	b13 = XBget(P(3), M(1), M(3))

	call XCcoeff(C, para, Bval(1,b12), Bval(1,b23), Bval(1,b13))
	end

************************************************************************

	subroutine XCcoeff(C, para, B12, B23, B13)
	implicit none
	ComplexType C(*), B12(*), B23(*), B13(*)
	DVAR para(1,*)

#include "lt.h"

	memindex XBget
	external XBget

	DVAR p1, p2, p1p2, m1, m2, m3
	DVAR f1, f2
	QVAR G(2,2)
	ComplexType bsum(0:2), b1sum(0:2), b00sum(0:2), b11sum(0:2)
	ComplexType in(0:2,2)
	logical dump

#ifdef SOLVE_EIGEN
	QVAR Ginv(2,2)
#define IN(i) in(:,i)
#define SOLVE_SETUP XInverse(2, G,2, Ginv,2)
#define SOLVER(b) XSolve(2, G,2, Ginv,2, b)
#else
	integer perm(2)
#define IN(i) in(:,perm(i))
#define SOLVE_SETUP XDecomp(2, G,2, perm)
#define SOLVER(b) XSolve(2, G,2, b)
#endif
#define OUT(i) in(:,i)

	integer finite(0:2)
	data finite /1, 0, 0/

	serial = serial + 1
	dump = ibits(debugkey, DebugC, 1) .ne. 0 .and.
     &    serial .ge. debugfrom .and. serial .le. debugto
	if( dump ) call XDumpPara(3, para, "Ccoeff")

	m1 = M(1)
	m2 = M(2)
	m3 = M(3)
	p1 = P(1)
	p2 = P(2)
	p1p2 = P(3)

	f1 = m2
	  f1 = f1 - m1
	  f1 = f1 - p1
	f2 = m3
	  f2 = f2 - m1
	  f2 = f2 - p1p2

	G(1,1) = 2*p1
	G(2,2) = 2*p1p2
	G(1,2) = p1
	  G(1,2) = G(1,2) + p1p2
	  G(1,2) = G(1,2) - p2
	G(2,1) = G(1,2)
	call SOLVE_SETUP

	bsum = B23(BB0) + B23(BB1)
	b1sum = B23(BB1) + B23(BB11)
	b00sum = B23(BB00) + B23(BB001)
	b11sum = B23(BB11) + B23(BB111)

	call XC0func(C(CC0), para)

	IN(1) = f1*C(CC0) - B23(BB0) + B13(BB0)
	IN(2) = f2*C(CC0) - B23(BB0) + B12(BB0)
	call SOLVER(in)
	C(CC1) = OUT(1)
	C(CC2) = OUT(2)

	C(CC00) = .5D0*(m1*C(CC0) -
     &    .5D0*(f1*C(CC1) + f2*C(CC2) - B23(BB0) - finite))

	IN(1) = f1*C(CC1) + bsum - 2*C(CC00)
	IN(2) = f2*C(CC1) + bsum + B12(BB1)
	call SOLVER(in)
	C(CC11) = OUT(1)
	C(CC12) = OUT(2)

	IN(1) = f1*C(CC2) - B23(BB1) + B13(BB1)
	IN(2) = f2*C(CC2) - B23(BB1) - 2*C(CC00)
	call SOLVER(in)
	C(CC12) = .5D0*(C(CC12) + OUT(1))
	C(CC22) = OUT(2)

c$$$  
c$$$	  MOD MARCEL: Do not need all Coeffs.
c$$$	  

c$$$  C(CC001) = 1/3D0*(m1*C(CC1) -
c$$$     &    .5D0*(f1*C(CC11) + f2*C(CC12) + bsum + finite/3D0))
c$$$	C(CC002) = 1/3D0*(m1*C(CC2) -
c$$$     &    .5D0*(f1*C(CC12) + f2*C(CC22) - B23(BB1) + finite/3D0))
c$$$
c$$$	bsum = bsum + b1sum
c$$$	IN(1) = f1*C(CC11) - bsum - 4*C(CC001)
c$$$	IN(2) = f2*C(CC11) - bsum + B12(BB11)
c$$$	call SOLVER(in)
c$$$	C(CC111) = OUT(1)
c$$$	C(CC112) = OUT(2)
c$$$
c$$$	IN(1) = f1*C(CC22) - B23(BB11) + B13(BB11)
c$$$	IN(2) = f2*C(CC22) - B23(BB11) - 4*C(CC002)
c$$$	call SOLVER(in)
c$$$	C(CC122) = OUT(1)
c$$$	C(CC222) = OUT(2)
c$$$
c$$$	C(CC0000) = 1/4D0*(m1*C(CC00) -
c$$$     &    .5D0*(f1*C(CC001) + f2*C(CC002) - B23(BB00) -
c$$$     &      finite*(m1 + m2 + m3 - .25D0*(p1 + p2 + p1p2))/6D0))
c$$$
c$$$	IN(1) = f1*C(CC001) + b00sum - 2*C(CC0000)
c$$$	IN(2) = f2*C(CC001) + b00sum + B12(BB001)
c$$$	call SOLVER(in)
c$$$	C(CC0011) = OUT(1)
c$$$	C(CC0012) = OUT(2)
c$$$
c$$$	IN(1) = f1*C(CC002) - B23(BB001) + B13(BB001)
c$$$	IN(2) = f2*C(CC002) - B23(BB001) - 2*C(CC0000)
c$$$	call SOLVER(in)
c$$$	C(CC0012) = .5D0*(C(CC0012) + OUT(1))
c$$$	C(CC0022) = OUT(2)
c$$$
c$$$	bsum = bsum + b1sum + b11sum
c$$$	IN(1) = f1*C(CC111) + bsum - 6*C(CC0011)
c$$$	IN(2) = f2*C(CC111) + bsum + B12(BB111)
c$$$	call SOLVER(in)
c$$$	C(CC1111) = OUT(1)
c$$$	C(CC1112) = OUT(2)
c$$$
c$$$	IN(1) = f1*C(CC222) - B23(BB111) + B13(BB111)
c$$$	IN(2) = f2*C(CC222) - B23(BB111) - 6*C(CC0022)
c$$$	call SOLVER(in)
c$$$	C(CC1222) = OUT(1)
c$$$	C(CC2222) = OUT(2)
c$$$
c$$$	IN(1) = f1*C(CC122) + b11sum - 2*C(CC0022)
c$$$	IN(2) = f2*C(CC122) + b11sum - 4*C(CC0012)
c$$$	call SOLVER(in)
c$$$	C(CC1122) = OUT(1)
c$$$	C(CC1222) = .5D0*(C(CC1222) + OUT(2))

	if( dump ) call XDumpCoeff(3, C)
	end

