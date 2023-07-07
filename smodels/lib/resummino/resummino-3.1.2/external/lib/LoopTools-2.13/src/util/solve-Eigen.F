* solve-Eigen.F
* computation of the inverse and solution of a linear system
* by diagonalizing the matrix with the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of LoopTools
* last modified 9 Dec 10 th

#include "externals.h"
#include "types.h"
#include "defs.h"


* A matrix is considered diagonal if the sum of the squares
* of the off-diagonal elements is less than EPS.

#define EPS 2D0**(-102)


************************************************************************
** XEigen diagonalizes a complex symmetric n-by-n matrix.
** Input: n, A = n-by-n matrix
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^T = U A U^-1 with U U^T = 1.

	subroutine XEigen(n, A,ldA, d, U,ldU)
	implicit none
	integer n, ldA, ldU
	QVAR A(ldA,*), U(ldU,*), d(*)

	integer p, q, j, sweep
	QREAL red, off, thresh
	QVAR delta, t, invc, s, x, y
	QVAR ev(2,MAXDIM)

	QREAL sq
	QVAR c
	sq(c) = QRE(c*QCC(c))

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = A(p,p)
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .04D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + sq(A(p,q))
	    enddo
	  enddo
	  if( .not. off .gt. EPS ) return

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      delta = A(p,q)
	      off = sq(delta)
	      if( sweep .gt. 4 .and. off .lt.
     &              EPS*max(sq(ev(2,p)), sq(ev(2,q))) ) then
	        A(p,q) = 0
	      else if( off .gt. thresh ) then
	        x = .5D0*(ev(2,p) - ev(2,q))
	        y = sqrt(x**2 + delta**2)
	        t = x - y
	        s = x + y
	        if( sq(t) .lt. sq(s) ) t = s

	        t = delta/t
	        delta = delta*t
	        ev(1,p) = ev(1,p) + delta
	        ev(2,p) = d(p) + ev(1,p)
	        ev(1,q) = ev(1,q) - delta
	        ev(2,q) = d(q) + ev(1,q)

	        invc = sqrt(t**2 + 1)
	        s = t/invc
	        t = t/(invc + 1)

	        do j = 1, p - 1
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + s*(y - t*x)
	          A(j,q) = y - s*(x + t*y)
	        enddo

	        do j = p + 1, q - 1
	          x = A(p,j)
	          y = A(j,q)
	          A(p,j) = x + s*(y - t*x)
	          A(j,q) = y - s*(x + t*y)
	        enddo

	        do j = q + 1, n
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + s*(y - t*x)
	          A(q,j) = y - s*(x + t*y)
	        enddo

	        A(p,q) = 0

	        do j = 1, n
	          x = U(p,j)
	          y = U(q,j)
	          U(p,j) = x + s*(y - t*x)
	          U(q,j) = y - s*(x + t*y)
	        enddo
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in XEigen"
	end

************************************************************************

	subroutine XDet(n, A,ldA, det)
	implicit none
	integer n, ldA
	QVAR A(ldA,*), det

	QVAR d(MAXDIM), U(MAXDIM,MAXDIM)
	integer p

	call XEigen(n, A,ldA, d, U,MAXDIM)

	det = 1
	do p = 1, n
	  det = det*d(p)
	enddo
	end

************************************************************************
** XInverse forms the (pseudo)inverse of a symmetric n-by-n matrix.
** Input: n, A = n-by-n matrix, symmetric
** (only the upper triangle of A needs to be filled).
** Output: Ainv = (pseudo)inverse of A

	subroutine XInverse(n, A,ldA, Ainv,ldAinv)
	implicit none
	integer n, ldA, ldAinv
	QVAR A(ldA,*), Ainv(ldAinv,*)

	integer p, q, j
	QVAR U(MAXDIM,MAXDIM), d(MAXDIM), t

#if 0
	PRINT *, "SEigen"
	PRINT *, "A11=", A(1,1)
	PRINT *, "A12=", A(1,2)
	PRINT *, "A21=", A(2,1)
	PRINT *, "A22=", A(2,2)
	if( ldA .gt. 2 ) then
	PRINT *, "A13=", A(1,3), A(3,1)
	PRINT *, "A23=", A(2,3), A(3,2)
	PRINT *, "A33=", A(3,3)
	endif
	PRINT *, "-----------"
#endif

	call XEigen(n, A,ldA, d, U,MAXDIM)

* form (pseudo)inverse U^T d^-1 U

	do p = 1, n
	  do q = 1, n
	    t = 0
	    do j = 1, n
	      if( abs(d(j)) .gt. EPS )
     &          t = t + U(j,q)*U(j,p)/d(j)
	    enddo
	    Ainv(q,p) = t
	  enddo
	enddo

#if 0
	PRINT *, "d1=", d(1)
	PRINT *, "d2=", d(2)
	if( ldA .gt. 2 ) PRINT *, "d3=", d(3)
	PRINT *, "-----------"

	PRINT *, "U11=", U(1,1)
	PRINT *, "U12=", U(1,2)
	PRINT *, "U21=", U(2,1)
	PRINT *, "U22=", U(2,2)
	if( ldA .gt. 2 ) then
	PRINT *, "U13=", U(1,3)
	PRINT *, "U23=", U(2,3)
	PRINT *, "U31=", U(3,1)
	PRINT *, "U32=", U(3,2)
	PRINT *, "U33=", U(3,3)
	endif
	PRINT *, "-----------"

	PRINT *, "Ai11=", Ainv(1,1)
	PRINT *, "Ai12=", Ainv(1,2)
	PRINT *, "Ai21=", Ainv(2,1)
	PRINT *, "Ai22=", Ainv(2,2)
	if( ldA .gt. 2 ) then
	PRINT *, "Ai13=", Ainv(1,3)
	PRINT *, "Ai23=", Ainv(2,3)
	PRINT *, "Ai31=", Ainv(3,1)
	PRINT *, "Ai32=", Ainv(3,2)
	PRINT *, "Ai33=", Ainv(3,3)
	endif
	PRINT *, "-----------"
#endif
	end

************************************************************************

	subroutine XSolve(n, A,ldA, Ainv,ldAinv, b)
	implicit none
	integer n, ldA, ldAinv
	QVAR A(ldA,*), Ainv(ldAinv,*)
	ComplexType b(*)

	integer i, j
	ComplexType x(MAXDIM)

#if defined(QUAD) && !defined(COMPLEXPARA)

	QVAR tr, ti

    	do i = 1, n
	  tr = 0
	  ti = 0
	  do j = 1, n
	    tr = tr + Ainv(i,j)*Re(b(j))
	    ti = ti + Ainv(i,j)*Im(b(j))
	  enddo
	  x(i) = ToComplex(tr, ti)
	enddo

#else

	ComplexType delta(MAXDIM), t

	do i = 1, n
	  t = 0
	  do j = 1, n
	    t = t + Ainv(i,j)*b(j)
	  enddo
	  x(i) = t
	enddo

* improve on x
	do i = 1, n
	  t = 0
	  do j = 1, n
	    t = t + A(i,j)*x(j)
	  enddo
	  delta(i) = b(i) - t
	enddo
	do i = 1, n
	  t = 0
	  do j = 1, n
	    t = t + Ainv(i,j)*delta(j)
	  enddo
	  x(i) = x(i) + t
	enddo

#endif

	do i = 1, n
	  b(i) = x(i)
	enddo
	end

