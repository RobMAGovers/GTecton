subroutine LOADLK (B,NSLIP,STN,X,D,SKEW,IEN,LM, &
                   LMX,LMF,TFAULT,MAT,PRPMAT,DXE,OSKEW,ID,IDX,NSTR, &
                   NSD,NDOF,NEN)

USE MODELDEFINITION
use constants,       only: dfmin,dlogmin,dfmax,dlogmax, eps, zero, one, two
use debugmodule,     only: itest, idebug, lustat
use timestepmodule,  only: STEP,ITWORK,IMWORK,IMTWRK

!
! Computes load due to unlocking an interface
!
implicit none

!-pass
integer          :: NSTR,NSD,NDOF,NEN
integer          :: IEN(NEN,*),LM(NDOF,NEN,*),LMX(NDOF,NEN,*), *
                 :: LMF(NDOF,NEN,*), MAT(*), IDX(NDOF,*),NSLIP(5,*)
double precision :: B(*),STN(NSTR,*),X(NSD,*),D(NDOF,*), &
                    SKEW(2,*), TFAULT(NDOF,*),PRPMAT(6,*), &
                    DXE(NDOF,NEN,*),OSKEW(2,*)

!-locl
logical          :: quad
double precision :: PL(24),XL(24),DL(24),BL(3),ROT(3,3)
integer          :: i, nn, n, m
double precision :: thick

do nn=1,NUMSN

	if (NSLIP(3,nn).eq.1 .and. NSLIP(4,nn).eq.1) then
		cycle
	endif

	n = NSLIP(1,nn)
	m = MAT(n)
	THICK = PRPMAT(6,m)
	quad = (IEN(3,n).ne.IEN(4,n))

	call CLEAR (PL,NEE,"PL")
	call LCOORD (X,XL,IEN(1,n),0)

	if (LGDEF.ne.0 .and. (.not.CAUCHY)) then
		call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
		call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN)
		call ADDSNE (DL,DXE(1,1,n),NDOF,NEN)
		call REZONE (XL,DL,ONE)
	endif

	call EFORCE (STN(1,n),PL,XL,THICK,quad,ierr)

	if (ierr.ne.0) then
		write(stderr,10) n
10        format(1x,'unlock2: error in element ',I5)
		call stoper()
	endif

	if (NUMROT.ne.0) then
		if (NSLSKEW.gt.0) then
		    CALL RPFORC (PL,OSKEW,IEN(1,n),NDOF,NEN)
		else
		    CALL RPFORC (PL,SKEW,IEN(1,n),NDOF,NEN)
		endif
	endif
	call ADDFOR (B,PL,LMX(1,1,n),NEE,nn,40)
enddo

do nn=1,NUMSN
	n = IDSLP(nn)
	do i=1,NDOF
	    if (IDX(i,n).gt.0) then
			DIFORC(i,n) = TWO * B(IDX(i,n))
		endif
	enddo
enddo

return
end subroutine
