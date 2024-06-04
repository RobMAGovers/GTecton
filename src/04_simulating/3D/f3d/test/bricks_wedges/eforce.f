C-------------------------------------------------------------------------------
* RG 30/1/'94: not checked
*
	subroutine EFORCE (EVP,P,XL,DUM1,QUAD,ierr)

	USE MODELDEFINITION
	USE CONSTANTS, only: G
*
c Computes the effective forces P at each node
*
c SNGL	implicit real (a-h,o-z)
	implicit double precision (A-H,O-Z)
	parameter (NEN=8, NSD=3, NDOF=3, NSTR=6)
	parameter (NEE=NEN*NDOF)
c-pass
	integer ierr
	logical QUAD
	dimension XL(NSD,NEN),P(NEE),EVP(NSTR)
c-locl
	dimension RG(NEN),SG(NEN),TG(NEN),XS(NSD,NSD),SH(4,NEN),
  >	 B(NSTR,NEE)
	save RG,SG,TG
c-init
	data RG/-1., 1., 1.,-1.,-1., 1., 1.,-1./
data SG/-1.,-1., 1., 1.,-1.,-1., 1., 1./
data TG/-1.,-1.,-1.,-1., 1., 1., 1., 1./
*
	if (IPOINT.eq.1) then
	    NINT = 1
	    W    = 8.
	    G_local    = 0.
	else
	    NINT = 8
	    W    = 1.
!	    G_local = G
	endif
*
	do 200 l=1,NINT
c	    Compute spatial derivatives of shape functions
	    call SHAP30 (RG(l)*G_local,
  >               SG(l)*G_local,
  >               TG(l)*G_local,
  >               XL,DET,SH,XS,QUAD,
  >	     ierr,3)
	    if (ierr.ne.0) return
	    call BMATRIX (B,SH)
	    do 100 i=1,NEE
		EFFE=B(1,i)*EVP(1)+B(2,i)*EVP(2)+B(3,i)*EVP(3)
  >		 +B(4,i)*EVP(4)+B(5,i)*EVP(5)+B(6,i)*EVP(6)
		P(i)=P(i)+W*DET*EFFE
100	    continue
200	continue

	return
	end
