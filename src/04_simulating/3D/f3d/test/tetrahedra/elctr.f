	subroutine ELCTR (XELM,IEN,X,D,LMF,TFAULT,
  >   SKEW,DXE,FACTOR,IERR)

	USE MODELDEFINITION
*
c Calculates (deformed) element center coordinates XELM
*
	implicit none
	integer NSD,NDOF,NEN
	parameter (NSD=3, NDOF=3, NEN=4)
	double precision QUART
	parameter (QUART=0.25d0)
c-pass
	integer IEN,LMF,IERR
	double precision XELM,X,D,TFAULT,SKEW,DXE,FACTOR
	dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*),
  >	 LMF(NDOF,NEN),TFAULT(NDOF,*),
  >	 SKEW(2,*),DXE(NDOF,NEN)
c-comm
	integer NTYPE,NUMEL,NUMAT,IOPT,IPOINT,NUMPR,NUMSLP,
  >	 NUMFN,NPRE,LGDEF,IRESDU,NUMSTR,IGRAV,IVLIM,
  >	 NUMWNK,NSURF,NSED,INCOMP
  >	 NOCOMPR,NSLSKEW

c-locl
	double precision xl,dl
	dimension xl(NSD,NEN),dl(NDOF,NEN)
*
c	localize coordinates
	call LCOORD (X,xl,IEN)
	if (FACTOR.gt.0.) then
	    call LDISP  (dl,D,IEN,NDOF,NEN)
	    call AddFaultDisplacement (dl,LMF,TFAULT,NDOF,NEN)
	    call ADDSNE (dl,DXE,NDOF,NEN)
	    call REZONE (xl,dl,FACTOR)
	endif
	XELM(1) = (xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))*QUART
	XELM(2) = (xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))*QUART
	XELM(3) = (xl(3,1)+xl(3,2)+xl(3,3)+xl(3,4))*QUART
*
	return
	end
