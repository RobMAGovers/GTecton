subroutine gravty (IEN,X,D,MAT,PRPMAT,PRPTEM,T,SURF,GRAV,IXPND, &
                   IYPND,TX,LMTX,LMF,TFAULT,DXE,NXPND,NYPND,xmn,xmx,xstp,lu, &
                   verbose,luMon)

use materialsmodule, only: MaterialIDOfElement
use modeldefinition, only: mode, numel
use constants,       only: re
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: Ygref, g0, iCompens

!    USE MODELDEFINITION

! Routine for calculating gravity signals. Not included as yet: sediments
! and water.

! ICOMPENS=0: Calculate raw gravity data along (erosion) surface
!         =1: Calculate free-air gravity data at Y=0
!         =2: Calculate Bouguer gravity data at Y=0

 implicit none

integer, parameter :: NSD=2,NDOF=2,NEN=4
!-pass
integer            :: lu,luMon,NXPND,NYPND,NSURF
logical            :: verbose

double precision   :: xmn,xmx,xstp

integer            :: ien, mat, ixpnd, iypnd, lmtx, lmf
double precision   :: x, d, prpmat, prptem, t, surf, grav, tx, tfault, dxe

dimension          :: IEN(NEN,*),X(NSD,*),D(NDOF,*),MAT(*),PRPMAT(6,*), &
                      PRPTEM(6,*),T(*),IXPND(2,*),IYPND(2,*),SURF(NSD,*), &
                      GRAV(NDOF),LMTX(NEN,*),TX(*),TFAULT(NDOF,*),LMF(NDOF,NEN,*), &
                      DXE(NDOF,NEN,*)



!-locl
real               :: D2REAL
integer            :: i, numnn, n, nn, m, numbr
double precision   :: xloc, yloc
double precision   :: g, dg, dens, gmag, texp


double precision, external :: gravel, srfhgt

gmag = SQRT(GRAV(1)*GRAV(1)+GRAV(2)*GRAV(2))

!       Calculate number of data points
numbr = INT(ABS((xmx-xmn)/xstp)) + 1

if (numbr.le.0) then
    return
endif

if (icompens.eq.0) then
    write(lu,*) -999,' raw surface gravity anomalies'
else if (icompens.eq.1) then
    write(lu,*) -999,' free-air gravity anomalies at Y=',Ygref
else if (icompens.eq.2) then
    write(lu,*) -999,' Bouguer gravity anomalies at Y=',Ygref
endif

write(lu,*) 2,numbr

do i=1,numbr
    xloc = DBLE(i-1)*xstp + xmn
    if (icompens.ne.2 .and. NSURF.gt.0) then
        yloc = srfhgt(SURF,NSURF,xloc)
    else
!        Bouguer: calculate contribution from below reference level
        yloc = Ygref
    endif

    if (elmon) then
        numnn = nnel
        g = 0d0
    else
        numnn = NUMEL
        g = -g0
    endif

    do nn=1,numnn

        if (elmon) then
            n = monel(nn)
        else
            n = nn
        endif

        m = MaterialIDOfElement(n)

        if (MODE.gt.2) then
            TEXP = PRPTEM(3,m)
        endif

        DENS = PRPMAT(5,m)

        dg = gravel(IEN,X,D,T,IXPND,IYPND,TX,LMTX,LMF(1,1,n), &
        TFAULT,DXE(1,1,n),n,NXPND,NYPND,DENS,TEXP,xloc, &
        yloc,icompens)

        g = g + dg
    enddo

    if (icompens.eq.1) then
        g = g - gmag*(yloc-Ygref)*2d0/RE
    endif

    if (elmon) then
        if (verbose) then
            if (icompens.eq.0) then
                write(luMon,1) n,xloc,yloc,g
 1                   format(1x,'raw gravity contribution of element', &
                 1x,I5,' at (',1PG12.4,',',1PG12.4,') = ', &
                 1PG11.3)
            else if (icompens.eq.1) then
                write(luMon,2) n,xloc,Ygref,g
 2                   format(1x,'free-air gravity contribution of ', &
                 'element ',I5,' at (',1PG12.4,',',1PG12.4, &
                 ') = ',1PG11.3)
            else if (icompens.eq.2) then
                write(luMon,3) n,xloc,Ygref,g
 3                   format(1x,'Bouguer gravity contribution of ', &
                 'element ',I5,' at (',1PG12.4,',',1PG12.4, &
                 ') = ',1PG11.3)
            endif
        else
            if (icompens.eq.0) then
                write(luMon,*) n,xloc,yloc,g
            else
                write(luMon,*) n,xloc,Ygref,g
            endif
        endif
    else
        write(lu,*) D2REAL(xloc),D2REAL(g)
    endif
enddo

return
end subroutine

!-----------------------------------------------------------------------

subroutine gref (IEN,X,D,MAT,PRPMAT,PRPTEM,T,LMF,TFAULT,DXE, &
                 IXPND,IYPND,TX,LMTX,NXPND,NYPND,X0)

USE MODELDEFINITION
use materialsmodule, only: MaterialIDOfElement
use constants,       only: zero
use plotControl,     only: Ygref, g0

! Routine for calculating reference gravity signal at (X0,Ygref).

implicit none

integer, parameter :: NSD=2,NDOF=2,NEN=4
!-pass
integer            :: NXPND,NYPND
double precision   :: X0

integer            :: ien, mat, ixpnd, iypnd, lmtx, lmf
double precision   :: x, d, prpmat, prptem, t, surf, grav, tx, tfault, dxe

dimension          :: IEN(NEN,*),X(NSD,*),D(NDOF,*),MAT(*),PRPMAT(6,*), &
                      PRPTEM(6,*),T(*),IXPND(2,*),IYPND(2,*),LMTX(NEN,*),TX(*), &
                      TFAULT(NDOF,*),LMF(NDOF,NEN,*),DXE(NDOF,NEN,*)

!-locl
integer            :: n, nn, m
double precision   :: dg, dens, texp
double precision, external :: gravel

g0 = zero
do nn=1,NUMEL

    n = nn
    m = MaterialIDOfElement(n)

    if (MODE.gt.2) then
        TEXP = PRPTEM(3,m)
    endif

    DENS = PRPMAT(5,m)

    dg = gravel(IEN,X,D,T,IXPND,IYPND,TX,LMTX,LMF(1,1,n), &
                TFAULT,DXE(1,1,n),n,NXPND,NYPND,DENS,TEXP,X0,Ygref,0)

    g0 = g0 + dg
enddo

return
end subroutine
!-----------------------------------------------------------------------
double precision function gravel (IEN,X,D,T,IXPND,IYPND,TX,LMTX,LMF,TFAULT,DXE, &
                                  NEL,NXPND,NYPND,RHO0,TEXP,X0,Y0,icompens)

USE MODELDEFINITION
use constants,       only: useDXE, zero, one, third
use iomodule,        only: stderr


! Routine to calculate gravity signal in (x0,y0) of a single element.
! When SYMM=.true., it is assumed that the model is symmetric around x=0, and
! the contribution from both the element at x=X and x=-X are added.

implicit none

integer, parameter :: NEN=4,NSD=2,NDOF=2
integer, parameter :: NENP2=NEN+2
!-pass
integer            :: N,NXPND,NYPND

integer            :: ien, ixpnd, iypnd, lmtx, lmf
double precision   :: x, d, t, tx, tfault, dxe

dimension          :: IEN(NEN,*),X(NSD,*),D(NDOF,*),T(*),IXPND(2,*), &
                      IYPND(2,*),TX(*),LMTX(NEN,*),TFAULT(NDOF,*),LMF(NDOF,NEN), &
                      DXE(NDOF,NEN)

integer            :: NEL, icompens
double precision   :: RHO0,TEXP,X0,Y0


!-locl
logical            :: QUAD,SYMM
double precision   :: XL(NSD,NEN),DL(NDOF,NEN),SH(3,NEN), &
                      XS(NSD,NSD),TL(NEN),XE(NSD,NENP2)
double precision   :: sn, tn, T0, grav, dgrav, area, rho
integer            :: i, m, ierr, ip

double precision, external :: Talwani, gxpand

!-init
SYMM = .false.
gravel = zero
sn = zero
tn = zero
T0 = zero

call LCOORD (X,XL,IEN(1,NEL))
call LDISP  (DL,D,IEN(1,NEL),NDOF,NEN)
call AddFaultDisplacement (DL,LMF,TFAULT,NDOF,NEN,157)
call ADDSNE (DL,NEL,NDOF,NEN,useDXE)
call REZONE (XL,DL,1d0)

if (MODE.gt.3) then
    call LTEMP (IEN(1,NEL),T,TL,NEN)
    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
       call LTEMPX (LMTX(1,NEL),IEN(1,NEL),TX,TL,NEN)
    endif
endif

QUAD = (IEN(3,NEL).ne.IEN(4,NEL))
!    Check whether sign of x-components is uniform
m = NEN - 1

if (.not.QUAD) then
    m=m-1
endif

do i=1,m
    if (XL(1,i)*XL(1,i+1).lt.zero) then
        write(stderr,10)
10           format(1x,'gravel: non-uniform x-component signs'/ &
         9x,'x=0 symmetry violated')
        return
    endif
enddo

if (.not.QUAD) then
    tn = -third
endif

call SHAP20 (sn,tn,XL,AREA,SH,XS,QUAD,ierr)

if (ierr.ne.0) then
    write(stderr,20) NEL
20        format(1x,'gravel: error in element ',I5)
    return
endif

!    Element density
if (MODE.gt.3) then
!        Calculate element center temperature
    do i=1,m+1
        T0 = T0 + SH(3,i)*TL(i)
    enddo

    RHO = RHO0 * (one - TEXP * T0)
else
    RHO = RHO0
endif

!    Select part of the element which is below Y0 (icompens=2)
call ELBELW (Y0,XL,QUAD,XE,n,icompens)
if (n.le.2) then
    return
endif

!    Gravity contribution from element
grav  = Talwani(XE,n,X0,Y0,RHO)

!    Add expansion to 400 km depth and to the left and right
dgrav = gxpand(XE,IXPND,IYPND,NEL,NXPND,NYPND,QUAD,SYMM, &
               X0,Y0,RHO)
grav  = grav + dgrav

if (SYMM) then
!        Gravity contribution from other side of symmetry center
    do i=1,n
        XE(1,i) = -XE(1,i)
    enddo

    dgrav = Talwani(XE,n,X0,Y0,RHO)
    grav = grav - dgrav

endif
!
gravel = grav

return
end function
!-----------------------------------------------------------------------
subroutine elbelw (Y0,XL,QUAD,XE,N,icompens)

! ELementBELoW: routine to determine if the element or parts of the element
! have to be taken into account for the gravity calculation.
! The criterion is that only (parts of) elements below
! the horizontal line y=Y0 are to be included. For this reason, nodal point
! coordinates of the part of the element below this line are returned in
! array XE.
! N is the number of nodal points. If N<=2, no gravity calculation is required.

use constants, only: dfmin

implicit none

integer, parameter :: NEN=4, NSD=2
integer, parameter :: NENP2=NEN+2
!-pass
integer            :: N, icompens
double precision   :: Y0
double precision   :: XL(NSD,NEN), XE(NSD,NENP2)
logical            :: QUAD
!-locl
logical            :: PrevBelow, CurBelow
integer            :: inp, nnp, ip

double precision, external :: fjune

!-init
N = 0

if (QUAD) then
    nnp = NEN
else
    nnp = 3
endif

CurBelow = (XL(2,1)-Y0 .lt. -dfmin .or. icompens.ne.2)

if (CurBelow) then
    N = N + 1
    XE(1,N) = XL(1,1)
    XE(2,N) = XL(2,1)
endif

do inp=2,nnp+1
    ip=inp

    if (ip.gt.nnp) then
        ip=1
    endif

    PrevBelow = CurBelow
    CurBelow = (XL(2,ip)-Y0 .lt. -dfmin .or. icompens.ne.2)

    if (CurBelow.and.PrevBelow) then
        if (inp.le.nnp) then
            N = N + 1
            XE(1,N) = XL(1,ip)
            XE(2,N) = XL(2,ip)
        endif
    else if (PrevBelow) then
        N = N + 1
        XE(1,N) = &
           fjune(XL(1,inp-1),XL(2,inp-1),XL(1,ip),XL(2,ip),Y0)
        XE(2,N) = Y0
    else if (CurBelow) then
        N = N + 1
        XE(1,N) = &
       fjune(XL(1,inp-1),XL(2,inp-1),XL(1,ip),XL(2,ip),Y0)
        XE(2,N) = Y0
        if (inp.le.nnp) then
            N = N + 1
            XE(1,N) = XL(1,ip)
            XE(2,N) = XL(2,ip)
        endif
    endif

enddo
return
end subroutine
!-----------------------------------------------------------------------
double precision function fjune(x1,y1,x2,y2,Y0)

! Function written in june 1993, to calculate intercept with line
! through nodal points X with y=Y0

use constants, only: dfmin

implicit none
!-locl
double precision :: x1,y1,x2,y2,Y0
double precision :: a, b, dy

dy = y2 - y1

if (ABS(dy).le.dfmin) then
    fjune = (x1+x2)/2d0
else
    a = (x2-x1)/dy
    b = (x1*y2-x2*y1)/dy
    fjune = a*Y0+b
endif

return
end
!-----------------------------------------------------------------------
double precision function Talwani (X,N,XP,YP,RHO)

USE CONSTANTS, only: grav_const, zero, two

! Function to calculate vertical gravitational attraction at (XP,YP)
! due to a polygon with corner coordinates X and mass density RHO.
! The coordinate system for input data is X=horizontal to the right,
! and Y=vertical upwards. When the element is positively oriented
! (i.e. anti-clockwise), the gravity vector points downward.
!
! Ref: Talwani, M., J.L. Worzel and M. Landisman, Rapid gravity computations
! for two-dimensional bodies with application to the Mendocino submarine
! fracture zone, J. Geoph. Res. 64, 49-59, 1959.

implicit none

integer, parameter :: NSD=2
!-pass
integer            :: N
double precision   :: X(NSD,N)
double precision   :: XP,YP,RHO
!-local
integer            :: i, j
double precision   :: grav
double precision, external :: tallin
double precision   :: Xi, Yi, Xj, Yj

!-init
Talwani = ZERO
grav    = ZERO

if (N.le.1) then
    return
endif

Xj = X(1,1) - XP
Yj = YP - X(2,1)

do j=2,N+1
    Xi = Xj
    Yi = Yj
    if (j.le.N) then
        Xj = X(1,j) - XP
        Yj = YP - X(2,j)
    else
        Xj = X(1,1) - XP
        Yj = YP - X(2,1)
    endif
    grav = grav + tallin (Xi,Yi,Xj,Yj)
enddo

Talwani = TWO * grav_const * RHO * grav

return
end function
!-----------------------------------------------------------------------
double precision function tallin (Xii,Yii,Xjj,Yjj)

USE CONSTANTS, only: pi, zero, half, two,eps
! Contribution of a line piece with uniform y-sign

implicit none

!-pass
double precision :: Xii,Yii,Xjj,Yjj
!-local
double precision :: Xi, Yi, Xj, Yj
double precision :: xrg, dx, dy
double precision :: R2, Ri2, Rj2, pre
double precision :: thetai, thetaj

double precision, external :: elog

Xi = Xii
Yi = Yii
Xj = Xjj
Yj = Yjj

tallin = ZERO

call polar(Xi,Yi,Ri2,thetai)
if (Ri2.le.EPS) then
    return
endif

call polar(Xj,Yj,Rj2,thetaj)

if (Rj2.le.EPS) then
    return
endif

if (ABS(thetai-thetaj).le.EPS) then
    return
endif

if ((ABS(thetai-PI/TWO).le.EPS .or. &
     ABS(thetai+PI/TWO).le.EPS     ).and. &
    (ABS(thetaj-PI/TWO).le.EPS .or. &
     ABS(thetaj+PI/TWO).le.EPS     )     ) then
    return
endif

if (ABS(Xi-Xj).le.EPS) then
    xrg = (Xi*Xi)/(Xj*Xj)*Rj2/Ri2
    xrg = ELOG(xrg)
    tallin = Xi * HALF * xrg
else if (ABS(Yi-Yj).le.EPS) then
    tallin = Yi*(thetaj-thetai)
else if (ABS(Xj).le.EPS) then
    dx   = Xj-Xi
    dy   = Yj-Yi
    R2   = dx*dx + dy*dy
    pre  = (Xj*dx*dy - Yj*dx*dx)/R2
    xrg  = (Yi - Xi*dy/dx)/Ri2**HALF
    xrg  = ELOG(xrg)
    tallin = pre * (thetai-thetaj+dy/dx*xrg)
else if (ABS(Xi).le.EPS) then
    dx   = Xj-Xi
    dy   = Yj-Yi
    R2   = dx*dx + dy*dy
    pre  = (Yj*dx*dx - Xj*dx*dy)/R2
    xrg  = (Yj - Xj*dy/dx)/Rj2**HALF
    xrg  = ELOG(xrg)
    tallin = pre * (thetaj-thetai+dy/dx*xrg)
else
    dx   = Xj-Xi
    dy   = Yj-Yi
    R2   = dx*dx + dy*dy
    pre  = (Xj*dx*dy - Yj*dx*dx)/R2
    xrg  = (Yi-Xi*dy/dx)/(Yj-Xj*dy/dx)*(Rj2/Ri2)**HALF
    xrg  = ELOG(xrg)
    tallin = pre * (thetai-thetaj+dy/dx*xrg)
endif
return
end function

!-----------------------------------------------------------------------

subroutine polar (Xi,Yi,R2,theta)

USE CONSTANTS, only: pi, zero, two, dfmin, half
! Calculates R**2 and THETA [-PI/2,PI/2]

implicit none
!-locl
double precision :: Xi,Yi,R2,theta


if (ABS(Xi).gt.dfmin) then
    theta = ATAN(Yi/ABS(Xi))
else
    if (Yi.lt.-dfmin) then
        theta = -PI * half
    else if (Yi.gt.dfmin) then
        theta = PI * half
    else
        theta = ZERO
    endif
endif
    R2 = Xi*Xi + Yi*Yi
return
end subroutine

!-----------------------------------------------------------------------

double precision function gxpand(XL,IXPND,IYPND,N,NXPND,NYPND,&
                                 QUAD,SYMM,XP,YP, RHO)

use constants, only: zero, two
use plotControl,     only: rhoa, g0

! Routine to calculate gravity contributions from model continuations
! outside actual model.

implicit none

integer, parameter          :: NSD=2,NEN=4
double precision, parameter :: BOTTOM=-400D3
!-pass
integer                     :: N,NXPND,NYPND
logical                     :: QUAD,SYMM

integer                     :: IXPND,IYPND
double precision            :: XL
dimension                      XL(NSD,*),IXPND(2,*),IYPND(2,*)

double precision            :: XP,YP, rho
!-local
integer                     :: nenn, ix, iy, i
integer                     :: inode, jnode
double precision            :: x,h,t,dg, xn1

double precision, external  :: horsht, versht

!-init
NENN = NEN

if (.not.QUAD) then
    NENN = NEN - 1
endif

gxpand = 0d0
ix = 0
iy = 0

!    Check if element number N occurs in IXPND array
i = 1

100 if (i.gt.NXPND) then
    goto 200
endif

if (N.eq.IXPND(2,i)) then
    ix = i
    goto 200
endif

i = i + 1
goto 100

!    Check if element number N occurs in IYPND array
200 i = 1
300 if (i.gt.NYPND) then
    goto 400
endif

if (N.eq.IYPND(2,i)) then
    iy = i
    goto 400
endif

i = i + 1
goto 300

400    if (ix.gt.0 .and. iy.gt.0) then
    inode = IXPND(1,ix)
    jnode = MOD(inode,NENN)+1

    x = (XL(1,inode)+XL(1,jnode))/TWO
    h = (XL(2,inode)+XL(2,jnode))/TWO
    t = ABS(XL(2,inode)-XL(2,jnode))

    dg = horsht(RHO,t,XP-x,h-YP)
    dg = dg + horsht(RHO,t,-XP-x,h-YP)

    inode = IYPND(1,iy)
    jnode = MOD(inode,NENN)+1

    x = MIN(XL(1,inode),XL(1,jnode))
    h = h - t/TWO
    t = ABS(h-BOTTOM)
    h = (h+BOTTOM)/TWO

    dg = dg + horsht(rhoa,t,XP-x,h-YP)

    if (SYMM) then
        dg = dg + horsht(rhoa,t,-XP-x,h-YP)
    endif

    gxpand = dg

else if (ix.gt.0) then

    inode = IXPND(1,ix)
    jnode = MOD(inode,NENN)+1

    x = (XL(1,inode)+XL(1,jnode))/TWO
    h = (XL(2,inode)+XL(2,jnode))/TWO
    t = ABS(XL(2,inode)-XL(2,jnode))
    xn1 = XL(2,jnode) - XL(2,inode)

    write(111,11) n,inode,jnode,t,xn1,h-yp,xp-x
11        format(3i5,4f12.2)

    if (xn1.lt.ZERO) then
       dg = -horsht(RHO,t,-(XP-x),h-YP)
    else
       dg = horsht(RHO,t,XP-x,h-YP)
    endif

    if (SYMM) then
        dg = dg + horsht(RHO,t,-XP-x,h-YP)
    endif

    gxpand = dg
else if (iy.gt.0) then

    inode = IYPND(1,iy)
    jnode = MOD(inode,NENN)+1

    x = (XL(1,inode)+XL(1,jnode))/TWO
    h = (XL(2,inode)+XL(2,jnode))/TWO
    t = ABS(XL(1,inode)-XL(1,jnode))

    dg = versht(rhoa,t,XP-x,h-YP)
    if (SYMM) then
        dg = dg + versht(rhoa,t,-XP-x,h-YP)
    endif

    gxpand = dg
endif

return
end function
!-----------------------------------------------------------------------
double precision function horsht (RHO,T,X,H)

USE CONSTANTS, only: grav_const, pi, zero, two

! Function for calculating gravity at (X,0) due to semi-infinite horizontal
! sheet with thickness T, density RHO and sheet center height H. The sheet
! starts at x=0.
! Telford, W.M., L.P. Geldart, R.E. Sheriff, D.A. Keys, 1976,
! Applied geophysics, 860pp.

implicit none

integer, parameter :: NSD=2, NEN=4
!-pass
double precision   :: RHO,T,X,H
!-locl
double precision   :: XL(NSD,NEN)
double precision, external   :: Talwani

if (ABS(H).gt.10.0d0*T) then
    horsht = -TWO*grav_const*RHO*T*(PI/TWO+ATAN(X/ABS(H)))
    if (H.gt.ZERO) then
        horsht = -horsht
    endif
else if (ABS(H).gt.1d-6) then
    XL(1,1)=5d7
    XL(2,1)=H-T/TWO
    XL(1,2)=5d7
    XL(2,2)=H+T/TWO
    XL(1,3)=ZERO
    XL(2,3)=H+T/TWO
    XL(1,4)=ZERO
    XL(2,4)=H-T/TWO
    horsht = Talwani (XL,NEN,X,ZERO,RHO)
else
    horsht = ZERO
endif

return
end
!-----------------------------------------------------------------------
double precision function versht (RHO,T,X,H)

USE CONSTANTS, only: grav_const, two

! Function for calculating gravity at (X,0) due to semi-infinite vertical
! sheet with thickness T, density RHO and sheet center height H. The sheet
! starts at x=0.
! Telford, W.M., L.P. Geldart, R.E. Sheriff, D.A. Keys, 1976,
! Applied geophysics, 860pp.

implicit none

double precision, parameter :: BOTTOM=-400D3

!-pass
double precision :: RHO,T,X,H
!-local
double precision :: xrg
double precision, external :: eLog

xrg = (BOTTOM*BOTTOM+X*X)/(X*X+H*H)
xrg = ELOG(xrg)
versht = -grav_const*RHO*T*xrg

return
end
!-----------------------------------------------------------------------
double precision function ELOG(xx)

use constants, only: dfmin, dfmax
use iomodule, only: stderr

! calculates natural logarithm of x.

implicit none

double precision :: x, xx


x = xx

if (x.le.dfmin) then
    write(stderr,1) dfmin
 1        format(1x,'ELOG domain error: x<= ',1PE11.5)
    ELOG=-dfmax
else if (x.ge.dfmax) then
    write(stderr,2) dfmax
 2        format(1x,'ELOG domain error: x>= ',1PE11.5)
    ELOG=DLOG(dfmax)
else
    ELOG=DLOG(x)
endif

return
end
