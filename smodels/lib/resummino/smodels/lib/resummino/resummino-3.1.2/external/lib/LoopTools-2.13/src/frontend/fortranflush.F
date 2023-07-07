* fortranflush.F
* C-callable function to flush the Fortran stdout (unit 6)
* this file is part of LoopTools
* last modified 6 Oct 09 th


	subroutine fortranflush()
	implicit none

	call flush(6)
	end

