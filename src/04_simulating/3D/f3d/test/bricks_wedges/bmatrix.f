c-------------------------------------------------------------------------------
* RG 30/1/'94: not checked
*
	subroutine BMATRIX (B,SH)
*
C Defines the strain-displacement matrix B
*
c SNGL	implicit real (a-h,o-z)
	implicit double precision (a-h,o-z)
	parameter (NEN=8, NSTR=6, NDOF=3)
	parameter (NEE=NDOF*NEN)
	parameter (ZERO=0.)
c-pass
	dimension B(NSTR,NEE),SH(4,NEN)
*
	k = 0
	do 100 i=1,NEN
	    k = k + 1
	    B(1,k) = SH(1,i)
	    B(2,k) = ZERO
	    B(3,k) = ZERO
	    B(4,k) = SH(2,i)
	    B(5,k) = SH(3,i)
	    B(6,k) = ZERO
	    k = k + 1
	    B(1,k) = ZERO
	    B(2,k) = SH(2,i)
	    B(3,k) = ZERO
	    B(4,k) = SH(1,i)
	    B(5,k) = ZERO
	    B(6,k) = SH(3,i)
	    k = k + 1
	    B(1,k) = ZERO
	    B(2,k) = ZERO
	    B(3,k) = SH(3,i)
	    B(4,k) = ZERO
	    B(5,k) = SH(1,i)
	    B(6,k) = SH(2,i)
100	continue
*
	return
	end
