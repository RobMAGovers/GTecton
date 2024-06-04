subroutine WINKLP (BTOT,DELD,IWELM,IWSIDE,IWTIME,WPRES,IEN,MAT,X,D,PRPMAT,LM, &
  LMF,TFAULT,SKEW,DXE,LMX,NSLIP,DX)
!
! Subroutine to compute Winkler restoring pressures from displacements (in local
! coordinate frame), and add them to the global force vector BTOT
!
USE MATERIALSMODULE, only : LMAT, MaterialIDOfElement
USE MODELDEFINITION
USE ALGEBRA
use constants,       only : zero, eps, half, one, useDXE
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only : elementSlipperyID
#endif

implicit none
integer, parameter :: NEN=4, NSTR=5, NDOF=3, NSD=2

!-pass
integer          :: iwelm, iwside, iwtime, ien, mat, lm, lmf, lmx, nslip
double precision :: wpres, x, prpmat, btot, deld, d, tfault, skew, dxe, dx
dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*),MAT(*),X(NSD,*),PRPMAT(6,*), &
   LM(NDOF,NEN,*),BTOT(*),DELD(NDOF,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*),DX(NDOF,*)

!-locl
double precision :: f(NDOF),fl(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3),tmp(NDOF), &
  anglevec(2),wt,area
integer :: lien(4),WinklerEntry,globalelementnumber,timecode,materialnumber,nodecounter, &
  globalnodenumber,idof,eqno
logical :: quad
!
if (NUMWNK.eq.0) then
    return
endif

do 400 WinklerEntry=1,IABS(NUMWNK)

    globalelementnumber = IABS(IWELM(WinklerEntry))
    timecode = IWTIME(WinklerEntry)
    if (timecode.lt.0) then
        if (timecode.ne.-1.and.NSTEP.ge.-timecode) then
            goto 400
        endif
    elseif (timecode.eq.0) then
        goto 400
    else
        if (NSTEP.lt.(timecode-1)) then
            goto 400
        endif
    endif

    call LCOORD (X,XL,IEN(1,globalelementnumber),0)
 
    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,globalelementnumber),NDOF,NEN)
        call AddFaultDisplacement (DL,globalelementnumber,TFAULT,NDOF,NEN,15)
        call ADDSNE (DL,globalelementnumber,NDOF,NEN,useDXE)
        call REZONE (XL,DL,ONE)
    endif

!   fill LIEN with local side numbers of element side IWSIDE(WinklerEntry)
    quad = IEN(3,globalelementnumber).ne.IEN(4,globalelementnumber)
    call SIDENP(IABS(IWSIDE(WinklerEntry)),quad,lien)
    
!   compute side length or surface area into variable AREA
    call FCGEOM (XL,lien,area,xn,"opn, addwnk")

    if (IWSIDE(WinklerEntry).eq.5 .or. IWSIDE(WinklerEntry).eq.6) then
!       Winkler pressure on top (5) or underside (6) of elements
        if (quad) then
            wt = area/4d0
        else
            wt = area/3d0
        endif
        f(1) = zero
        f(2) = zero
        f(3) = wt*WPRES(WinklerEntry)
    else
!       Winkler pressure perpendicular to one of element sides 1-4
        materialnumber = MaterialIDOfElement(globalelementnumber)
        wt = PRPMAT(6,materialnumber)*area/2D0    ! element thickness * side length / 2
        f(1) = ABS(xn(1)) * WPRES(WinklerEntry) * wt
        f(2) = ABS(xn(2)) * WPRES(WinklerEntry) * wt
        f(3) = zero
    endif

! compute local nodal forces FL from global nodal forces F for first two nodes
    do nodecounter=1,2
        globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
        anglevec(1) = -SKEW(1,globalnodenumber)
        anglevec(2) = 0d0

        if (ABS(anglevec(1)).gt.EPS) then
            call FORMRT (anglevec,rot,NDOF)
            call EQUATE (tmp,f,NDOF)
            call VMPRD (rot,tmp,fl,NDOF,3)
        else
            call EQUATE (fl,f,NDOF)
        endif

!       add local forces to load vector 
        do idof=1,NDOF
            eqno = LM(idof,lien(nodecounter),globalelementnumber)
            if (eqno.gt.0) then
                BTOT(eqno) = BTOT(eqno) - fl(idof)*DELD(idof,globalnodenumber)
            endif
        enddo
    enddo

    nodecounter = 3
    if (lien(nodecounter).ne.0) then
! compute local nodal forces fl from global nodal forces for third node
        globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
        anglevec(1) = -SKEW(1,globalnodenumber)
        anglevec(2) = 0d0

        if (ABS(anglevec(1)).gt.EPS) then
            call FORMRT (anglevec,rot,NDOF)
            call EQUATE (tmp,f,NDOF)
            call VMPRD (rot,tmp,fl,NDOF,3)
        else
            call EQUATE (fl,f,NDOF)
        endif

!       add local forces to load vector 
        do idof=1,NDOF
            eqno = LM(idof,lien(nodecounter),globalelementnumber)
            if (eqno.gt.0) then
                BTOT(eqno) = BTOT(eqno) - fl(idof)*DELD(idof,globalnodenumber)
            endif
        enddo
    endif

    nodecounter = 4
    if (lien(nodecounter).ne.0) then
! compute local nodal forces fl from global nodal forces for fourth node
        globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
        anglevec(1) = -SKEW(1,globalnodenumber)
        anglevec(2) = 0d0

        if (ABS(anglevec(1)).gt.EPS) then
            call FORMRT (anglevec,rot,NDOF)
            call EQUATE (tmp,f,NDOF)
            call VMPRD (rot,tmp,fl,NDOF,3)
        else
            call EQUATE (fl,f,NDOF)
        endif

!       add local forces to load vector 
        do idof=1,NDOF
            eqno = LM(idof,lien(nodecounter),globalelementnumber)
            if (eqno.gt.0) then
                BTOT(eqno) = BTOT(eqno) - fl(idof)*DELD(idof,globalnodenumber)
            endif
        enddo
    endif

400 continue
!
return
end
