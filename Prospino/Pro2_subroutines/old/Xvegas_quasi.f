c$$$	program test
c$$$
c$$$	implicit none
c$$$
c$$$	integer          ncomp
c$$$	parameter       (ncomp=4)
c$$$
c$$$	integer          ndim, ncall, maxiter, neval, i1
c$$$	double precision result(ncomp)
c$$$	double precision accuracy
c$$$
c$$$	external         func
c$$$
c$$$	accuracy = 1.d-3
c$$$	
c$$$	ndim    = 3
c$$$
c$$$	ncall   = 10000
c$$$	maxiter = 12
c$$$	
c$$$	call vegas(func, ncomp, result, accuracy, ndim,
c$$$     +    ncall, maxiter, neval)
c$$$
c$$$c	print*, " neval  = ",neval
c$$$c	do i1=1,ncomp
c$$$c	   print*, " result(i1) = ",i1,result(i1)
c$$$c	end do
c$$$
c$$$	end 
c$$$
c$$$c -------------------------------------------------------	
c$$$	subroutine func(ndim, x, ncomp, result)
c$$$
c$$$	implicit none
c$$$
c$$$	integer          ndim, ncomp, i1
c$$$	double precision x(ndim), result(ncomp)
c$$$
c$$$c	print*, " input func: ndim = ",ndim
c$$$c	print*, " input func: x    = ",x
c$$$c	print*, " input func: ncomp= ",ncomp
c$$$
c$$$	do i1=1,ncomp
c$$$	   result(1) = 1.d0
c$$$	   result(2) = x(1)
c$$$	   result(3) = x(1)*x(2)
c$$$	   result(4) = x(1)*x(2)*x(3)
c$$$c	   print*, " output func ",i1,result(i1)
c$$$	end do
c$$$
c$$$	return
c$$$	end
c$$$c -------------------------------------------------------	

* vegas.F
* VEGAS Monte Carlo integration of a vector either with ordinary random
* numbers or with quasi-random numbers
* based on code by M. Martinez, J. Illana, J. Bossert, and A. Vicini
* this file is part of FormCalc
* last modified 14 Sep 01 th

* ndim_max is the maximum number of dimensions the integrand may have.
* Note: if your Fortran compiler allows, replace ndim_max by ndim
* (i.e. #define ndim_max ndim) to dimension the arrays dynamically.
*#define ndim_max 3 

* ncomp is the number of components of the integrand vector.
*#define ncomp 2 

* If DEBUG is defined, the error is printed out after each iteration
* so that one can tune the parameters for a particular integral.
ctp #define DEBUG

* RNG determines which random-number generator is used:
* - 1 uses Fortran's ran(...) function,
* - 2 uses a Faure quasi-random sequence,
* - 3 uses a Sobol quasi-random sequence (default).
ctp #define RNG 3


************************************************************************
** vegas integrates a vector in the ndim-dimensional unit hypercube
** using the Monte-Carlo method. The ncomp-dimensional integrand is
** invoked via the subroutine func(ndim, x, ncomp, result). After
** sampling ncall points, the grid is refined. This is called one
** iteration. The iterations loop terminates if the relative error is
** below accuracy, or after maxiter iterations.

** Caution: the refinement of the grid is done only with respect to
** the highest component of the integrand vector.

	subroutine vegas_hahn(func, ncomp, result, accuracy, ndim,
     +    ncall, maxiter, neval)
	implicit none
	integer ncomp, i1
	double precision result(ncomp), accuracy
	integer ndim, ncall, maxiter, neval
	external func

	double precision fun(ncomp), sfun(ncomp), sfun2(ncomp)
	double precision sint(ncomp), sweight(ncomp)
	double precision fun2, sint2, weight
	double precision r, dr, xo, xn, err
	integer iter, call, dim, grid, g, c

	integer rng
	parameter (rng = 1)

	integer ngrid
	parameter (ngrid = 100)

	double precision xi(ngrid, ndim), d(ngrid, ndim)
	double precision x(ndim), imp(ngrid), tmp(ngrid - 1)
	integer pos(ndim)

	if (rng.eq.1) then
	   call inirandom_1(maxiter*ncall, ndim)
	else if (rng.eq.2) then
	   call inirandom_2(maxiter*ncall, ndim)
	else if (rng.eq.3) then
	   call inirandom_3(maxiter*ncall, ndim)
	else
	   print*, " rng wrong, set to 1 "
	   call inirandom_1(maxiter*ncall, ndim)
	end if
	iter = 0
	neval = 0
	do c = 1, ncomp
	  sint(c) = 0
	  sweight(c) = 0
	enddo
	sint2 = 0

* define the initial distribution of intervals
	do grid = 1, ngrid
	  r = dble(grid)/ngrid
	  do dim = 1, ndim
	    xi(grid, dim) = r
	  enddo
	enddo

* iterations loop
1	continue
	iter = iter + 1

* initialize iteration variables
	do c = 1, ncomp
	  sfun(c) = 0
	  sfun2(c) = 0
	enddo
	do dim = 1, ndim
	  do grid = 1, ngrid
	    d(grid, dim) = 0
	  enddo
	enddo

	do call = 1, ncall
	  weight = 1D0/ncall

* compute the point position
	  if (rng.eq.1) then
	     call getrandom_1(x)
	  else if (rng.eq.2) then
	     call getrandom_2(x)
	  else if (rng.eq.3) then
	     call getrandom_3(x)
	  else 
	     print*, " rng wrong, set to 1 "
	     call getrandom_1(x)
	  end if
	  do dim = 1, ndim
	    r = x(dim)*ngrid + 1
	    grid = int(r)
	    xo = 0
	    if(grid .gt. 1) xo = xi(grid - 1, dim)
	    xn = xi(grid, dim) - xo
	    x(dim) = xo + (r - grid)*xn
	    pos(dim) = grid
	    weight = weight*xn*ngrid
	  enddo

* compute the function value
	  call func(ndim, x, ncomp, fun)
	  do c = 1, ncomp
	    fun2 = fun(c)*weight
	    sfun(c) = sfun(c) + fun2
	    fun2 = fun2**2
	    sfun2(c) = sfun2(c) + fun2
	  enddo
	  do dim = 1, ndim
	    d(pos(dim), dim) = d(pos(dim), dim) + fun2
	  enddo
	enddo
	neval = neval + ncall

* compute the integral and error values
	do c = 1, ncomp
	  fun2 = sfun(c)**2
	  r = sfun2(c)*ncall - fun2
	  if(r .ne. 0) then
	    weight = fun2/abs(r)*(ncall - 1)
	    sweight(c) = sweight(c) + weight
	    sint(c) = sint(c) + sfun(c)*weight
	  endif
	  if(sweight(c) .eq. 0) then
	    result(c) = 0
	  else
	    result(c) = sint(c)/sweight(c)
	  endif
	enddo
	sint2 = sint2 + fun2
	err = sqrt(sint2/(sweight(ncomp)*iter))/abs(result(ncomp))
ctp #ifdef DEBUG
ctp	print *, "iteration ", iter, "  error ", err
	if (ncomp.eq.1) then 
	   print ('(A7,I3,A12,G15.6,A8,G13.4)'), 
     &            "  iter=",iter,
     &            "  result=",result,
     &            "  error=",err
	else if (ncomp.eq.2) then 
	   print ('(A7,I3,A12,2(G15.6),A8,G13.4)'), 
     &            "  iter=",iter,
     &            "  result(i)=",result,
     &            "  error=",err
	else if (ncomp.eq.3) then 
	   print ('(A7,I3,A12,3(G15.6),A8,G13.4)'), 
     &            "  iter=",iter,
     &            "  result(i)=",result,
     &            "  error=",err
	else if (ncomp.eq.4) then 
	   print ('(A7,I3,A12,4(G15.6),A8,G13.4)'), 
     &            "  iter=",iter,
     &            "  result(i)=",result,
     &            "  error=",err
	else 
	   do i1=1,ncomp
	      print ('(A7,I3,A12,I3,2x,G15.6,A8,G13.4)'), 
     &            "  iter=",iter,
     &            "  result(i)=",i1,result(i1),
     &            "  error=",err
	   end do
	end if 
ctp #endif
	if(result(ncomp) .eq. 0 .or. err .lt. accuracy) return
	if(iter .gt. maxiter) then
	  print *,
     +      "Warning: VEGAS failed to reach the desired accuracy."
	  print *, "Remaining relative error: ", err
	  return
	endif

* redefine the grid (importance sampling)
* - smooth the f^2 value stored for each interval
	do dim = 1, ndim
	  xo = d(1, dim)
	  xn = d(2, dim)
	  d(1, dim) = .5D0*(xo + xn)
	  x(dim) = d(1, dim)
	  do grid = 2, ngrid - 1
	    r = xo + xn
	    xo = xn
	    xn = d(grid + 1, dim)
	    d(grid, dim) = (r + xn)/3D0
	    x(dim) = x(dim) + d(grid, dim)
	  enddo
	  d(ngrid, dim) = .5D0*(xo + xn)
	  x(dim) = x(dim) + d(ngrid, dim)
	enddo

* - compute the importance function of each interval
	do dim = 1, ndim
	  r = 0
	  do grid = 1, ngrid
	    imp(grid) = 0
	    if(d(grid, dim) .gt. 0) then
	      xo = x(dim)/d(grid, dim)
	      imp(grid) = ((xo - 1)/xo/log(xo))**1.5D0
	    endif
	    r = r + imp(grid)
	  enddo
	  r = r/ngrid

* - redefine the size of each interval
	  dr = 0
	  xn = 0
	  g = 0
	  do grid = 1, ngrid - 1
	    do while(dr .lt. r)
	      g = g + 1
	      dr = dr + imp(g)
	      xo = xn
	      xn = xi(g, dim)
	    enddo
	    dr = dr - r
	    tmp(grid) = xn - (xn - xo)*dr/imp(g)
	  enddo
	  do grid = 1, ngrid - 1
	    xi(grid, dim) = tmp(grid)
	  enddo
	  xi(ngrid, dim) = 1
	enddo

	goto 1
	end


ctp #if RNG == 1

************************************************************************
** inirandom sets up the random-number generator to produce at most
** max dims-dimensional quasi-random vectors

	subroutine inirandom_1(max, dims)
	implicit none
	integer max, dims

	integer ndim
	common /rngdata_1/ ndim

	ndim = dims
	end


************************************************************************
** getrandom is a subtractive Mitchell-Moore random-number generator.
** The algorithm is n(i) = (n(i - 24) - n(i - 55)) mod m, implemented
** as a circular array with n(i + 55) = n(i) and m = 2^30 in this
** version. The array n has been initialized by setting n(i) = i and
** running the algorithm 100,000 times. Code by Ronald Kleiss.

	subroutine getrandom_1(array)
	implicit none
	double precision array(*)

	integer ndim
	common /rngdata_1/ ndim

	integer dim, j, k, l, m, n(55)
	parameter (m = 2**30)
	data k /55/, l /31/
	data n /
     +    980629335, 889272121, 422278310,1042669295, 531256381,
     +    335028099,  47160432, 788808135, 660624592, 793263632,
     +    998900570, 470796980, 327436767, 287473989, 119515078,
     +    575143087, 922274831,  21914605, 923291707, 753782759,
     +    254480986, 816423843, 931542684, 993691006, 343157264,
     +    272972469, 733687879, 468941742, 444207473, 896089285,
     +    629371118, 892845902, 163581912, 861580190,  85601059,
     +    899226806, 438711780, 921057966, 794646776, 417139730,
     +    343610085, 737162282,1024718389,  65196680, 954338580,
     +    642649958, 240238978, 722544540, 281483031,1024570269,
     +    602730138, 915220349, 651571385, 405259519, 145115737 /

	do dim = 1, ndim
	  k = mod(k, 55) + 1
	  l = mod(l, 55) + 1
	  j = n(l) - n(k)
	  if(j .lt. 0) j = j + m
	  n(k) = j
	  array(dim) = dble(j)/m
	enddo
	end

ctp #elif RNG == 2

************************************************************************
** inirandom sets up the random-number generator to produce a Faure
** sequence of at most max dims-dimensional quasi-random vectors.
** Adapted from ACM TOMS algorithm 659, see
** http://www.acm.org/pubs/citations/journals/toms/1988-14-1/p88-bratley

	subroutine inirandom_2(max, dims)
	implicit none
	integer max, dims

	integer coeff(0:19, 0:19)
	integer ndim, prime, nextn, testn, digits
	common /rngdata_2/ coeff, ndim, prime, nextn, testn, digits

	integer i, j, h

	integer primes(40)
	save primes

	data primes /
     +    1, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11,
     +    13, 13, 17, 17, 17, 17, 19, 19,
     +    23, 23, 23, 23, 29, 29, 29, 29,
     +    29, 29, 31, 31, 37, 37, 37, 37,
     +    37, 37, 41, 41, 41 /

	ndim = dims
	prime = primes(dims)
	testn = prime**4
	nextn = testn - 1
	digits = 3

	h = nint(log(dble(max + testn))/log(dble(prime)))
	coeff(0, 0) = 1
	do j = 1, h
	  coeff(j, 0) = 1
	  coeff(j, j) = 1
	enddo
	do j = 1, h
	  do i = j + 1, h
	    coeff(i, j) =
     +        mod(coeff(i - 1, j) + coeff(i - 1, j - 1), prime)
	  enddo
	enddo
	end


************************************************************************
** getrandom generates a vector of random numbers

	subroutine getrandom_2(array)
	implicit none
	double precision array(*)

	integer coeff(0:19, 0:19)
	integer ndim, prime, nextn, testn, digits
	common /rngdata_2/ coeff, ndim, prime, nextn, testn, digits

	integer y(0:19), digit, dim, d, k, p
	double precision r

	p = testn
	k = nextn
	do digit = digits, 0, -1
	  p = p/prime
	  d = mod(k, p)
	  y(digit) = (k - d)/p
	  k = d
	enddo

	r = 0
	do digit = digits, 0, -1
	  r = (r + y(digit))/prime
	enddo
	array(1) = r

	do dim = 2, ndim
	  r = 0
	  p = 1
	  do digit = 0, digits
	    k = 0
	    do d = digit, digits
	      k = k + coeff(d, digit)*y(d)
	    enddo
	    y(digit) = mod(k, prime)
	    p = p*prime
	    r = r + dble(y(digit))/p
	  enddo
	  array(dim) = r
	enddo

	nextn = nextn + 1
	if(nextn .eq. testn) then
	   testn = testn*prime
	   digits = digits + 1
	endif
	end

ctp #else

************************************************************************
** inirandom sets up the random-number generator to produce a Sobol
** sequence of at most max dims-dimensional quasi-random vectors.
** Adapted from ACM TOMS algorithm 659, see
** http://www.acm.org/pubs/citations/journals/toms/1988-14-1/p88-bratley

	subroutine inirandom_3(max, dims)
	implicit none
	integer max, dims

	integer v(40, 30), lastq(40)
	integer ndim, count, norm
	common /rngdata_3/ v, lastq, ndim, count, norm

	integer bits, powers, degree, newv, dim, bit, deg, k

	integer poly(2:40), vinit(2:40, 1:8)
	save poly, vinit

	data poly / 3, 7, 11, 13, 19, 25, 37, 59, 47,
     +    61, 55, 41, 67, 97, 91, 109, 103, 115, 131,
     +    193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
     +    213, 191, 253, 203, 211, 239, 247, 285, 369, 299 /

	data (vinit(dim, 1), dim = 2, 40) / 39*1 /
	data (vinit(dim, 2), dim = 3, 40) /
     +          1, 3, 1, 3, 1, 3, 3, 1,
     +    3, 1, 3, 1, 3, 1, 1, 3, 1, 3,
     +    1, 3, 1, 3, 3, 1, 3, 1, 3, 1,
     +    3, 1, 1, 3, 1, 3, 1, 3, 1, 3 /
	data (vinit(dim, 3), dim = 4, 40) /
     +             7, 5, 1, 3, 3, 7, 5,
     +    5, 7, 7, 1, 3, 3, 7, 5, 1, 1,
     +    5, 3, 3, 1, 7, 5, 1, 3, 3, 7,
     +    5, 1, 1, 5, 7, 7, 5, 1, 3, 3 /
	data (vinit(dim, 4), dim = 6, 40) /
     +                  1, 7, 9, 13, 11,
     +    1, 3, 7, 9, 5, 13, 13, 11, 3, 15,
     +    5, 3, 15, 7, 9, 13, 9, 1, 11, 7,
     +    5, 15, 1, 15, 11, 5, 3, 1, 7, 9 /
	data (vinit(dim, 5), dim = 8, 40) /
     +                            9, 3, 27,
     +    15, 29, 21, 23, 19, 11, 25, 7, 13, 17,
     +    1, 25, 29, 3, 31, 11, 5, 23, 27, 19,
     +    21, 5, 1, 17, 13, 7, 15, 9, 31, 9 /
	data (vinit(dim, 6), dim = 14, 40) /
     +            37, 33, 7, 5, 11, 39, 63,
     +    27, 17, 15, 23, 29, 3, 21, 13, 31, 25,
     +    9, 49, 33, 19, 29, 11, 19, 27, 15, 25 /
	data (vinit(dim, 7), dim = 20, 40) /
     +                                  13,
     +    33, 115, 41, 79, 17, 29, 119, 75, 73, 105,
     +    7, 59, 65, 21, 3, 113, 61, 89, 45, 107 /
	data (vinit(dim, 8), dim = 38, 40) / 7, 23, 39 /

	k = max
	bits = 0
	do while(k .ne. 0)
	  bits = bits + 1
	  k = ishft(k, -1)
	enddo

	do bit = 1, bits
	  v(1, bit) = 1
	enddo

	do dim = 2, dims
	  powers = poly(dim)

	  k = powers
	  degree = -1
	  do while(k .ne. 0)
	    degree = degree + 1
	    k = ishft(k, -1)
	  enddo

	  do bit = 1, degree
	    v(dim, bit) = vinit(dim, bit)
	  enddo

	  do bit = degree + 1, bits
	    newv = v(dim, bit - degree)
	    k = powers
	    do deg = degree, 1, -1
	      if(btest(k, 0))
     +          newv = ieor(newv, ishft(v(dim, bit - deg), deg))
	      k = ishft(k, -1)
	    enddo
	    v(dim, bit) = newv
	  enddo
	enddo

	do bit = 1, bits - 1
	  do dim = 1, dims
	    v(dim, bit) = ishft(v(dim, bit), bits - bit)
	  enddo
	enddo
	norm = ishft(1, bits)

	count = 0
	ndim = dims
	do dim = 1, dims
	  lastq(dim) = 0
	enddo
	end


************************************************************************
** getrandom generates a vector of random numbers

	subroutine getrandom_3(array)
	implicit none
	double precision array(*)

	integer v(40, 30), lastq(40)
	integer ndim, count, norm
	common /rngdata_3/ v, lastq, ndim, count, norm

	integer c, zerobit, dim

	c = count
	zerobit = 1
	do while(btest(c, 0))
	  zerobit = zerobit + 1
	  c = ishft(c, -1)
	enddo

	do dim = 1, ndim
	  lastq(dim) = ieor(lastq(dim), v(dim, zerobit))
	  array(dim) = dble(lastq(dim))/norm
	enddo

	count = count + 1
	end

ctp #endif

