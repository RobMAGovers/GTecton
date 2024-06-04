 subroutine FCGEOM (XL,LIEN,AREA,XN)

 use constants, only: dfmin,dlogmin,dfmax,dlogmax
!
! Calculates area of an element face, perpendicular to face normal,
! and the unit normal vector.
!
 implicit none
 integer NSD,NEN
 parameter (NSD=3,NEN=4)
 double precision TWO
 parameter (TWO=2d0)
!-pass
 integer LIEN
 double precision XL,AREA,XN
 dimension XL(NSD,NEN),LIEN(3),XN(NSD)

!-locl
 integer i
 double precision XLL
 dimension XLL(NSD,4)
!
!    Localize face coordinates
 do i=1,3
XLL(1,i) = XL(1,LIEN(i))
XLL(2,i) = XL(2,LIEN(i))
XLL(3,i) = XL(3,LIEN(i))
 enddo
!
 XN(1)=(XLL(2,2)-XLL(2,1))*(XLL(3,3)-XLL(3,1)) &
    -(XLL(3,2)-XLL(3,1))*(XLL(2,3)-XLL(2,1))
 XN(2)=(XLL(3,2)-XLL(3,1))*(XLL(1,3)-XLL(1,1)) &
    -(XLL(1,2)-XLL(1,1))*(XLL(3,3)-XLL(3,1))
 XN(3)=(XLL(1,2)-XLL(1,1))*(XLL(2,3)-XLL(2,1)) &
    -(XLL(2,2)-XLL(2,1))*(XLL(1,3)-XLL(1,1))
 AREA=SQRT(XN(1)*XN(1)+XN(2)*XN(2)+XN(3)*XN(3))
 if (ABS(AREA).lt.dfmin) then
write(stderr,1)
 1        format(1x,'FCGEOM: zero element face area')
  return
 endif
 XN(1)=XN(1)/AREA
 XN(2)=XN(2)/AREA
 XN(3)=XN(3)/AREA
 AREA=AREA/TWO
!
 return
 end
