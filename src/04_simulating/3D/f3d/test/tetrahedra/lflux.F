 subroutine LFLUX (XL,TL,IEN,PRPTEM,FLUX,IERR)
!
! Routine to calculate element heat flux
!
 implicit none
 integer NSD,NEN,NDOF
 parameter (NSD=3,NEN=4,NDOF=3)
!-pass
 integer IEN,IERR
 double precision XL,TL,PRPTEM,FLUX
 dimension XL(NSD,NEN),TL(NEN),PRPTEM(6),FLUX(NSD),IEN(NEN)

!-locl
 double precision DET,SH,GRAD
 dimension GRAD(NDOF),SH(4,NEN)
!
!    Compute temperature gradient
 call SHPTET (0d0,0d0,0d0,XL,DET,SH,IERR,1)
 if (IERR.ne.0) return
!
 GRAD(1)=TL(1)*SH(1,1)+TL(2)*SH(1,2)+TL(3)*SH(1,3)+TL(4)*SH(1,4)
 GRAD(2)=TL(1)*SH(2,1)+TL(2)*SH(2,2)+TL(3)*SH(2,3)+TL(4)*SH(2,4)
 GRAD(3)=TL(1)*SH(3,1)+TL(2)*SH(3,2)+TL(3)*SH(3,3)+TL(4)*SH(3,4)
 FLUX(1) = -PRPTEM(4)*GRAD(1)
 FLUX(2) = -PRPTEM(5)*GRAD(2)
 FLUX(3) = -PRPTEM(6)*GRAD(3)
!
 return
 end
