subroutine strsar (IEN,X,D,STN,LMF,TFAULT,DXE, &
                   SKEW,MAT,PRPMAT,&
                   headangle,headsize,NEN,NSD,NDOF,NSTR, &
                   DEFORM,deffac,lu,ierr,ibase,istep, &
                   iarcont,verbose,luMon,tecton)
!
! routine for converting stresses to arrows in principal stress directions

USE MODELDEFINITION
USE CONSTANTS,       only: zero, deg2rad, dfmax
use materialsmodule, only: lmat, strinc
use timestepmodule,  only: nintg, maxit
use debugmodule,     only: iecho
use iomodule,        only: ascout, stderr
use postProcessData, only: GMT,connout
use monitorVars,     only: elmon, nnel, nsid, sidemon, monel, monsid
use plotControl,     only: hatangle, hatsize, &
                           Sscale, Ssb, difstn
use errors,          only: ok

implicit none

!-pass
integer            :: NEN,NSD,NDOF,NSTR,lu,ierr,ibase, &
                      istep,iarcont, luMon
double precision   :: deffac, eScale, esb
double precision   :: headangle,headsize

integer            :: ien, lmf, lmx, mat
double precision   :: stn, x, d, tfault, skew, dxe, prpmat
dimension          :: IEN(NEN,*),X(NSD,*),D(NDOF,*),STN(NSTR,*), &
                      LMF(NDOF,NEN,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*), &
                      SKEW(2,*),MAT(*),PRPMAT(6,*)

logical            :: DEFORM,VERBOSE
integer            :: tecton
!-locl
integer, parameter :: MAXRRW=25
real               :: D2REAL
double precision   :: normal
double precision   :: st(6),CENTR(3),arrow(3,MAXRRW),EVAL(3),EVEC(9), &
                      STmin(6),STmax(6),XL(24)

double precision   :: def
integer            :: nn, n, m, i, j
double precision   :: pois, deg
double precision   :: v1, v2
double precision   :: smallevecx, smallevecy
double precision   :: shear

integer            :: narrw

integer, external  :: ipltrw


!-init
data arrow(1,1),narrw/0d0,0/

ierr = OK
def  = ZERO

if (DEFORM) then
    def = deffac
endif
!
if (elmon.and.(tecton.eq.0)) then
    do nn=1,nnel
        n = monel(nn)
        if (difstn) then
!           deviatoric stress
            call LSTRES (st,STN(1,n),NSTR)
            if (strinc) then
                pois = 5d-1
            else
                if (LMAT) then
                    pois = PRPMAT(2,n)
                else
                    pois = PRPMAT(2,MAT(n))
                endif
            endif
            call stndev (st,pois)
        else
            call LSTRES (st,STN(1,n),NSTR)
        endif

        if (verbose) then
            if (difstn) then
                write(luMon,7,advance='no')
 7            format(' Deviatoric')
            endif
            write(luMon,1) n,(ST(i),i=1,NSTR)
 1            format(1x,'stress in element ',I5,': ',6(1PG12.4))
        else
            write(luMon,*) (ST(i),i=1,NSTR)
        endif

    enddo

else if (sidemon.and.(tecton.eq.0)) then

    do nn=1,nsid
        n = monsid(1,nn)
        m = monsid(2,nn)

        if (difstn) then

            call LSTRES (st,STN(1,n),NSTR)

            if (strinc) then
                pois = 5d-1
            else
                if (LMAT) then
                    pois = PRPMAT(2,n)
                else
                    pois = PRPMAT(2,MAT(n))
                endif
            endif

            call stndev (st,pois)
        else
            call LSTRES (ST,STN(1,n),NSTR)
        endif

        call traction (X,IEN(1,n),D,LMF(1,1,n),TFAULT, &
                       n,ST,SKEW,DEFORM,m,normal,shear)

        if (VERBOSE) then
            write(luMon,2) m,n,normal,shear
 2               format(1x,'normal and shear stress on side ',I1, &
             ' of element ',I5,': ',2(1PG12.4))
        else
            write(luMon,*) normal,shear
        endif
    enddo

else if (tecton.ne.0) then

    do nn=1,NUMEL
        call LSTRES (ST,STN(1,nn),NSTR)
        if (tecton.eq.2) then
            ST(1) = 5d-1*(ST(1)+ST(2))
            ST(2) = ST(1)
            do i=3,NSTR
                ST(i) = ZERO
            enddo
        endif
        write(luMon,3) nn,0,1,(ST(i),i=1,NSTR)
 3        format(2I5,I6,6(1PE14.6))
    enddo

    write(luMon,4)
 4        format('end pre-stress data')

else

    hatangle = headangle
    hatsize  = Ssb  *  Sscale  * headsize
!          [km]   = [MPa] x [km/MPa] *  factor

    do nn=1,NSTR
        STmin(nn) =  dfmax
        STmax(nn) = -dfmax
    enddo

    do nn=ibase,NUMEL,istep
        n = nn
!            calculate element center coordinates

        call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
                    SKEW, def,ierr)

        if (ierr.ne.0) then
            call exitp(1)
        endif

        if (difstn) then
            call LSTRES (st,STN(1,n),NSTR)

            if (strinc) then
                pois = 5D-1
            else
                if (LMAT) then
                    pois = PRPMAT(2,n)
                else
                    pois = PRPMAT(2,MAT(n))
                endif
            endif

            call stndev (st,pois)
        else
            call LSTRES (ST,STN(1,n),NSTR)
        endif

        do i=1,NSTR
            if (STmin(i) .gt. ST(i)) then
                STmin(i) =    ST(i)
            endif
            if (STmax(i) .lt. ST(i)) then
                STmax(i) =    ST(i)
            endif
        enddo

!       calculate Eigenvectors
        call Eigen (ST,EVAL,EVEC)

        if (GMT) then
            if (EVAL(1).ge.EVAL(2)) then
                smallevecx = EVEC(4)
                smallevecy = EVEC(5)
                v1 = EVAL(1)
                v2 = EVAL(2)
            else
                smallevecx = EVEC(1)
                smallevecy = EVEC(2)
                v1 = EVAL(2)
                v2 = EVAL(1)
            endif

            deg = atan2(smallevecy,smallevecx)
            deg = 90.0 - deg/deg2rad

            v1 = v1 * Sscale
            v2 = v2 * Sscale

            if (ascout) then
                write(lu,*,err=100) CENTR(1),CENTR(2),v1,v2,deg
            else
                write(lu,err=100) D2REAL(CENTR(1)), &
                   D2REAL(CENTR(2)),D2REAL(v1),D2REAL(v2), &
                   D2REAL(deg)
            endif

        else
            stop "non-GMT plotting not implemented"
        endif

    enddo

    if (iecho.eq.1 .or. &
        iecho.ge.4) then
        call wmnmx(STmin,STmax,NSTR,0)
    endif

endif

return
100  call prerr('stress')
     call exitp(1)

end subroutine

!-----------------------------------------------------------------------

subroutine traction (X,IEN,D,LMF,TFAULT,elemID, &
                     STN,SKEW,DEFORM,ISIDE,SN,ST)
!
use constants, only: useDXE
use modeltopology, only: nsd, ndof, nen, nstr
use debugmodule, only: iecho
use iomodule,        only: stderr

implicit none
!-pass
integer          :: ien, lmf, iside
double precision :: x,d,stn, skew, tfault
dimension        :: IEN(NEN),X(NSD,*),D(NDOF,*),STN(NSTR),SKEW(2,*), &
                    LMF(NDOF,NEN),TFAULT(NDOF,*)

logical          :: DEFORM
integer          :: elemID
double precision :: sn, st, s
!-locl
double precision :: XL(NSD,NEN),DL(NDOF,NEN),XT(2),XN(2),T(2)
integer          :: i,j
!
!    determine unit normal vector and tangential vector
call LCOORD (X,XL,IEN)

if (DEFORM) then
    call LDISP  (DL,D,IEN,NDOF,NEN)
    call AddFaultDisplacement (DL,LMF,TFAULT,NDOF,NEN,130)
    call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
    call REZONE (XL,DL,1d0)
endif

i = ISIDE
j = ISIDE + 1

if (IEN(3).eq.IEN(4)) then
    if (j.gt.3) then
        j=1
    endif
else
    if (j.gt.4) then
        j=1
    endif
endif

!    tangential vector
XT(1) = XL(1,j)-XL(1,i)
XT(2) = XL(2,j)-XL(2,i)

s = SQRT(XT(1)*XT(1)+XT(2)*XT(2))

XT(1) = XT(1)/s
XT(2) = XT(2)/s

!    positively oriented element: outward normal is (y1,-x1)
XN(1) = XT(2)
XN(2) = -XT(1)

!    traction
T(1) = STN(1)*XN(1)+STN(3)*XN(2)
T(2) = STN(3)*XN(1)+STN(2)*XN(2)

!    normal stress
SN = T(1)*XN(1)+T(2)*XN(2)

!    tangential stress
ST = T(1)*XT(1)+T(2)*XT(2)

if (iecho.eq.5) then
    write(stderr,1) XT(1),XT(2),XN(1),XN(2),T(1),T(2)
 1       format(/' stress: unit tangential vector = ',2(1PG14.6)/ &
           ' stress: unit normal vector     = ',2(1PG14.6)/ &
           ' stress: traction vector        = ',2(1PG14.6))
endif

return
end subroutine
