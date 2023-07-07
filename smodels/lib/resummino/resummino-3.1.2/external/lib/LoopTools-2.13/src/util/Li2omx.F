* Li2omx.F
* the dilogarithm function of 1 - x
* this file is part of LoopTools
* last modified 1 Aug 12 th

#include "externals.h"
#include "types.h"
#include "defs.h"


	ComplexType function XLi2omx(x)
	implicit none
	DVAR x

	ComplexType spence
	external spence

	XLi2omx = spence(1, ToComplex(x), 0D0)
	end

************************************************************************
* adapter code for C++

	subroutine XLi2omxsub(res, x)
	implicit none
	ComplexType res
	DVAR x

	ComplexType spence
	external spence

	res = spence(1, ToComplex(x), 0D0)
        end

