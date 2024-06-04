 subroutine ADDPR (B,X,D,IEN,LM,DUM1,DUM2,IELNO,ISIDE,PRES, &
   LMF,TFAULT,SKEW,DXE)

 USE MODELDEFINITION
 use constants, only: eps,third,half,one

!
! Subroutine to add pressures to load vector
!
 implicit none
 integer NEN,NSD,NDOF
 parameter (NEN=4, NSD=3, NDOF=3)
!-pass
 integer IEN,LM,IELNO,ISIDE,LMF
 double precision B,X,D,DUM1,DUM2,PRES,TFAULT,SKEW,DXE
 dimension IELNO(*),DUM1(*),IEN(NEN,*),X(NSD,*),D(NDOF,*), &
   ISIDE(*),PRES(*),DUM2(*),B(*),LM(NDOF,NEN,*), &
 LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
 DXE(NDOF,NEN,*)

!-locl
 integer k,n,lien,l,ii
 double precision xl,dl,area,xn,f,angle,tmp,rot
 dimension f(NDOF),xl(NSD,NEN),dl(NDOF,NEN),lien(3),xn(NSD), &
   angle(2),tmp(NSD),rot(3,3)
!
 do k=1,IABS(NUMPR)
!        Set element number N
n = IELNO(k)
!        Localize coordinates
call LCOORD (X,xl,IEN(1,n))
    if (LGDEF.ne.0) then
       call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
       call AddFaultDisplacement (dl,LMF(1,1,n),TFAULT,NDOF,NEN)
       call ADDSNE (dl,DXE(1,1,n),NDOF,NEN)
       call REZONE(xl,dl,ONE)
    endif
!        Get local nodal point numbers of the side on which pressure works
call SIDENP(ISIDE(k),.true.,lien)
!        Calculate surface area and unit normal vector
call FCGEOM (xl,lien,area,xn)
!        Define surface mean of SKEW angles
angle(1) = ( -SKEW(1,IEN(lien(1),n)) - &
       SKEW(1,IEN(lien(2),n))-SKEW(1,IEN(lien(3),n)) ) * HALF
angle(2) = ( -SKEW(2,IEN(lien(1),n)) - &
       SKEW(2,IEN(lien(2),n))-SKEW(2,IEN(lien(3),n)) ) * HALF
if (ABS(angle(1)).gt.EPS .or. ABS(angle(2)).gt.EPS) then
       call FORMRT (angle,rot,NDOF)
       call EQUATE (tmp,xn,NDOF)
       call VMPRD (rot,tmp,xn,NDOF,3)
    endif
f(1) = xn(1) * PRES(k) * area * THIRD
    f(2) = xn(2) * PRES(k) * area * THIRD
    f(3) = xn(3) * PRES(k) * area * THIRD
do l=1,NDOF
ii = LM(l,lien(1),n)
if (ii.ne.0) B(ii) = B(ii) + f(l)
ii = LM(l,lien(2),n)
if (ii.ne.0) B(ii) = B(ii) + f(l)
ii = LM(l,lien(3),n)
if (ii.ne.0) B(ii) = B(ii) + f(l)
enddo
 enddo
!
if (iecho.eq.1) write(stdout,1)
 1   format(1x,'Pressure loads are being added')
!
 return
 end
