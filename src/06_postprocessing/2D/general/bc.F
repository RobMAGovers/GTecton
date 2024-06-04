subroutine BCWP (X,D,SKEW,DEFORM, DEFFAC,LU,SYSC,SY,mattsym)


USE ALGEBRA,         only: clear, equate, vmprd, formrt
USE MODELDEFINITION, only: numwnk
use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
USE CONSTANTS,       only: pi, deg2rad,eps
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use debugmodule,     only: iecho
use iomodule,        only: ascout, stderr
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact
use postProcessData, only: GMT,connout
use symbols,         only: maxSym
!
! generate symbols element Winkler boundary conditions
!
implicit none
!-pass
logical DEFORM,SYSC
integer LU
double precision DEFFAC,SY

double precision :: x, d, skew
dimension X(NSD,*),D(NDOF,*),SKEW(2,*)
integer mattsym


!-locl
double precision :: xl(2,4),dl(12),xn(2),tmp(2),rot(3,3)
double precision :: angle(2), winkmx, size
integer          :: i, n, LL, j, k
double precision, external :: arctan
integer          :: nsymbl
double precision, allocatable :: symbl(:,:)

allocate(symbl(2, maxSym))


if (NUMWNK.eq.0) then
    if (.not.GMT) then
        if (ascout) then
            write(lu,*) 2,0
        else
            write(lu) 2.0,0.0
        endif
    endif
    write(*,*) "NUMWNK = 0, returning"
    return
endif
!
if (SYSC) then
    winkmx = -1d30
    do k=1,ABS(NUMWNK)
        tmp(1) = ABS(modeldatactx%WPRES(k))
        if (tmp(1).gt.winkmx) then
            winkmx=tmp(1)
        endif
    enddo
    if (winkmx.lt.EPS) then
        if (iecho.ne.0) then
             write(stderr,1)
 1        format(1x,'bc: no non-zero Winkler pressures')
        endif
        write(*,*) "all Winkler pressures very tiny."
        write(*,*) "plnplt does not have small enough font."
        return
    endif
endif
!


do k=1,IABS(NUMWNK)
    n = IABS(modeldatactx%IWELM(k))
    LL = modeldatactx%IWTIME(k)

    if (LL.lt.0) then
        if (LL.ne.-1 .and. NSTEP.ge.-LL) then
            cycle
        endif
    elseif (LL.eq.0) then
        cycle
    else
        if (NSTEP.lt.(LL-1)) then
            cycle
        endif
    endif

    call CLEAR (symbl,2*MAXSYM,"symbl")
    call LCOORD (X,xl,meshdatactx%IEN(1,n))

    if (DEFORM) then
        call LDISP  (dl,D,meshdatactx%IEN(1,n),NDOF,NEN)
        call REZONE (xl,dl,DEFFAC)
    endif

    i = IABS(modeldatactx%IWSIDE(k))
    j = i + 1

    if (meshdatactx%IEN(3,n).ne.meshdatactx%IEN(4,n)) then
        if (j.gt.NEN) then
            j = 1
        endif
    else
        if (j.ge.NEN) then
            j = 1
        endif
    endif

    angle(1) = (-SKEW(1,meshdatactx%IEN(i,n))-SKEW(1,meshdatactx%IEN(j,n)))*5D-1
    angle(2) = 0.
    xn(1) = xl(2,j) - xl(2,i)
    xn(2) = xl(1,i) - xl(1,j)

    if (ABS(angle(1)).gt.EPS) then
        call FORMRT (angle,rot,ndof)
        call EQUATE (tmp,xn,NDOF)
        call VMPRD (rot,tmp,xn,NDOF,3)
    endif

    angle(1) = ARCTAN(xn(1),xn(2))
    xl(1,1) = (xl(1,i)+xl(1,j))*5d-1
    xl(2,1) = (xl(2,i)+xl(2,j))*5d-1


    if (modeldatactx%WPRES(k).gt.EPS) then

        if (GMT) then
            angle(1) = angle(1)/deg2rad - 180.0
            if (SYSC) then
                size = ABS(modeldatactx%WPRES(k))/winkmx
            else
                size = 1.0
            endif

            if (mattsym.eq.1) then
!               cartesian angle anticlockwise from x-axis
                write(LU,*,err=2000) xl(1,1),xl(2,1),size,13,angle(1)+90
            else if (mattsym.eq.2) then
!               geographic: angle clockwise from N
                write(LU,*,err=2000) xl(1,1),xl(2,1),size,13,-angle(1)
            else
                write(LU,*,err=2000) xl(1,1),xl(2,1),angle(1),size
            endif
!             write(*,*,err=2000) xl(1,1),xl(2,1),angle(1),size

 6            format(4(1X,E15.6))
        else
            stop "non-GMT plotting not implemented"
        endif

    else if (modeldatactx%WPRES(k).lt.-EPS) then

        if (GMT) then
            angle(1) = angle(1)/deg2rad + 180.0
            if (SYSC) then
                size = ABS(modeldatactx%WPRES(k))
            else
                size = 1.0
            endif

            if (mattsym.eq.1) then
                write(LU,*,err=2000) xl(1,1),xl(2,1),size,13,angle(1)+90
            else if (mattsym.eq.2) then
                write(LU,*,err=2000) xl(1,1),xl(2,1),size,13,-angle(1)
            else
                write(LU,*,err=2000) xl(1,1),xl(2,1),angle(1),size
            endif
            !write(LU,*,err=2000) xl(1,1),xl(2,1),angle,size
!            write(*,*,err=2000) xl(1,1),xl(2,1),angle,size

        else
            stop "non-GMT plotting not implemented"
        endif

    endif

enddo
!
 return
2000    write(stderr,2001)
 2001    format(1x,'BCWP write error')
 return
 end
!-----------------------------------------------------------------------
subroutine BCSTRS (X,D,SKEW,IEN,IELSTR,ISSIDE,STRS,ISTR,DEFORM, &
                   DEFFAC,LU,SYSC,SY,mattsym)

USE ALGEBRA,         only: formrt, clear, equate, vmprd
USE MODELDEFINITION, only: numstr
USE CONSTANTS,       only: pi, deg2rad,eps
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: ascout, stderr
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact
use postProcessData, only: GMT,connout
use symbols,         only: maxSym
!
! generate symbols stress boundary conditions
!
implicit none
!-pass
 logical DEFORM,SYSC
 integer LU
 double precision DEFFAC,SY

integer :: ielstr, isside, istr, ien
double precision :: x, d, skew, strs

 dimension X(NSD,*),D(NDOF,*),SKEW(2,*),STRS(6,*),IELSTR(*), &
   ISSIDE(*),ISTR(2,*),IEN(NEN,*)

integer mattsym


!-locl
double precision :: xl(2,4),dl(12),xn(2),tmp(3),rot(3,3)
double precision :: EVAL(3),EVEC(9)
integer          :: i, j, k, l, it0, it1, n
double precision :: angle(2), deg, scl, size
integer          :: nsymbl
double precision :: smallevecx, smallevecy, strsmx
double precision :: sxx, sxy, syy
double precision :: v1, v2
double precision, allocatable :: symbl(:,:)

allocate(symbl(2, maxSym))


!
 if (NUMSTR.le.0) then
    if (.not.GMT) then
        if (ascout) then
            write(lu,*) 2,0
        else
            write(lu) 2.0,0.0
        endif
    endif
    return
 endif
!
if (SYSC) then
    strsmx = -1d30
    do l=1,NUMSTR
        do k=1,6
            tmp(1) = ABS(STRS(k,l))
            if (tmp(1).gt.strsmx) strsmx=tmp(1)
        enddo
    enddo
    if (ABS(strsmx).lt.EPS) return
endif
!
do k=1,NUMSTR
    ! get element number, application times
    n = IELSTR(k)
    it0 = ISTR(1,k)
    it1 = ISTR(2,k)

    if (DBLE(NSTEP-it0)*DBLE(NSTEP-it1).le.0d0) then
        call CLEAR (symbl,2*MAXSYM,"symbl")
        call LCOORD (X,xl,IEN(1,n))
        if (DEFORM) then
            call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
            call REZONE (xl,dl,DEFFAC)
        endif
        i = ISSIDE(k)  ! side of element where stress applied
        j = i + 1
        if (IEN(3,n).ne.IEN(4,n)) then
            if (j.gt.NEN) j = 1
        else
            if (j.ge.NEN) j = 1
        endif

        angle(1) = (-SKEW(1,IEN(i,n))-SKEW(1,IEN(j,n)))*5D-1
        xn(1) = xl(2,j) - xl(2,i)
        xn(2) = xl(1,i) - xl(1,j)
        if (ABS(angle(1)).gt.EPS) then
            call FORMRT (angle,rot,ndof)
            call EQUATE (tmp,xn,NDOF)
            call VMPRD (rot,tmp,xn,NDOF,3)
        endif

        sxx = STRS(1,k) / DBLE(it1-it0+1)
        syy = STRS(2,k) / DBLE(it1-it0+1)
        sxy = STRS(4,k) / DBLE(it1-it0+1)
        if (gmt) then
            xl(1,1) = (xl(1,i)+xl(1,j))*5d-1
            xl(2,1) = (xl(2,i)+xl(2,j))*5d-1
            tmp(1) = sxx
            tmp(2) = syy
            tmp(3) = sxy
            call Eigen (tmp,EVAL,EVEC)
            if (EVAl(1).ge.EVAL(2)) then
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

            deg = atan2(smallevecx,smallevecy)
            deg = deg/deg2rad

            if (sysc) then
                scl = 1.0/abs(strsmx)
            else
                if (abs(v1).gt.EPS) then
                    scl = 1.0/abs(v1)
                elseif (abs(v2).gt.EPS) then
                    scl = 1.0/abs(v2)
                else
                    write(0,*) 'bc: no non-zero eigenvalue'
                    return
                endif
            endif
            if (mattsym.eq.1 .or. mattsym.eq.2) then
                if (v1.gt.0.0) then
                    if (mattsym.eq.1) then
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v1*scl,12,deg
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v1*scl,12,deg+180.
                    else if (mattsym.eq.2) then
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v1*scl,12,90.-deg
                    endif
                else
                    if (mattsym.eq.1) then
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v1*scl,11,deg
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v1*scl,11,deg+180.
                    else if (mattsym.eq.2) then
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v1*scl,11,90.-deg
                    endif
                endif
                if (v2.gt.0.0) then
                    if (mattsym.eq.1) then
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v2*scl,12,deg+90.
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v2*scl,12,deg-90.
                    else if (mattsym.eq.2) then
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v2*scl,12,-deg
                    endif
                else
                    if (mattsym.eq.1) then
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v2*scl,11,deg+90.
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v2*scl,11,deg-90.
                    else if (mattsym.eq.2) then
                        write(lu,*,err=2000) xl(1,1),xl(2,1),v2*scl,11,-deg
                    endif
                endif
            else
                write(lu,11,err=2000) xl(1,1),xl(2,1),v1,v2,deg
            endif
11          format(E15.6,1X,E15.6,1X,G12.6,1X,G12.6,1X,E15.6)
        else
            stop "non-GMT plotting not implemented"
        endif
    endif
enddo
!
return
2000    write(stderr,2001)
 2001    format(1x,'BCSTRS write error')
return
end subroutine
!-----------------------------------------------------------------------
subroutine SLPPLT (X,D,SKEW,IDSLP,NSLIP,DEFORM,DEFFAC,LU,SY,mattsym)

USE ALGEBRA,         only: clear
USE MODELDEFINITION, only: numsn, numslp
USE CONSTANTS,       only: pi, deg2rad, zero
use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use postProcessData, only: GMT,connout
use symbols,         only: maxSym
!
! generate symbols for slippery extra dofs
!
implicit none
!-pass
logical DEFORM
integer LU
double precision DEFFAC,SY

integer          :: idslp, nslip
double precision :: x, d, skew
dimension X(NSD,*),D(NDOF,*),IDSLP(*),NSLIP(5,*),SKEW(2,*)
integer :: mattsym



!-locl
double precision :: deg,length
double precision :: xl(2)
integer i, j, n, nsymbl
double precision, allocatable :: symbl(:,:)

allocate(symbl(2, maxSym))

!
if (NUMSN.le.0) then
    return
endif
!
do i=1,NUMSN

    call CLEAR (symbl,2*MAXSYM,"symbl")
    n = modeldatactx%IDSLP(i)
    j = 1

    do while (         j.le.NUMSLP .and. &
              NSLIP(2,j).ne.n)
        j = j + 1
    enddo

    if (GMT) then

        xl(1) = meshdatactx%X(1,n)+DEFFAC*modeldatactx%D(1,n)
        xl(2) = meshdatactx%X(2,n)+DEFFAC*modeldatactx%D(2,n)

        length = (NSLIP(3,j))**2 + (NSLIP(4,j))**2
        length = sqrt(length)

        deg = atan2(DBLE(NSLIP(4,j)),DBLE(NSLIP(3,j)))
        deg = (deg - SKEW(1,n))/deg2rad


        if (mattsym.eq.1) then
            write(LU,*,err=2000) xl(1),xl(2),1,21,deg
            !if (NSLIP(3,j).ne.0.and.NSLIP(4,j).ne.0) then
            !    write(LU,*,err=2000) xl(1),xl(2),1,23,-SKEW(1,n)/deg2rad
            !elseif (NSLIP(3,j).ne.0) then
            !    write(LU,*,err=2000) xl(1),xl(2),1,21,-SKEW(1,n)/deg2rad
            !else
            !    write(LU,*,err=2000) xl(1),xl(2),1,22,-SKEW(1,n)/deg2rad
            !endif
        else if (mattsym.eq.2) then
            write(LU,*,err=2000) xl(1),xl(2),1,21,90.-deg
        else
            ! write an arrow in each direction
            write(LU,6,err=2000) xl(1),xl(2),deg,length
 6          format(E15.6,1X,E15.6,1X,E15.6,1X,E15.6)

            deg = deg + 180.0

            write(LU,6,err=2000) xl(1),xl(2),deg,length
        endif

    else
        stop "non-GMT plotting not implemented"
    endif
enddo
!
return

2000    write(stderr,2001)
 2001    format(1x,'SLPPLT write error')

end subroutine
!-----------------------------------------------------------------------
subroutine NBCPLT (IEN,X,D,IBOND,BOND,SKEW,DEFORM,DEFFAC,LU,SY, &
                   SYSC,iflg,mattsym)

USE CONSTANTS, only: pi, deg2rad, zero,eps
USE MODELDEFINITION
USE MESHDATAMODULE, only: meshdatactx
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact
use postProcessData, only: GMT,connout

!
! generate symbols for nodal boundary conditions
!
!  iflg      what to print
!
!    1        displacement
!    2        velocity
!    3        force
!    4        initial displacement
!    5        strain rate
!    6        all of the above
!
!


implicit none
!-pass
logical DEFORM,SYSC
integer LU,iflg
double precision DEFFAC,SY

integer          :: ien, ibond, mattsym
double precision :: x, d, bond, skew
dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),IBOND(NDOF,*), &
   BOND(NDOF,*),SKEW(2,*)

!-locl
logical          :: nbc
double precision :: length, size
double precision :: xl(2)
double precision :: deg
double precision :: dispmx, forcmx, velomx
integer          :: i,j,n
double precision :: tmp
double precision :: bcmag1,bcmag2
double precision :: angleA, angleB

nbc = .false.
!
!    write(*,*) 'entered NBCPLT with nvertex: ', &
!                meshdatactx%Nvglobal

if (SYSC) then
!        determine maxima to boundary condition values
    dispmx = -1d30
    forcmx = -1d30
    velomx = -1d30

    do j=1,meshdatactx%Nvglobal
        do i=1,2
            tmp = ABS(BOND(i,j))
            if ((IBOND(i,j).eq.1).or.(IBOND(i,j).eq.4.and. &
                NSTEP.eq.0)) then
                if (tmp.gt.dispmx) then 
                    dispmx=tmp
                endif
            elseif (IBOND(i,j).eq.2) then
                if (tmp.gt.velomx) then
                    velomx=tmp
                endif
            elseif (IBOND(i,j).eq.3) then
                if (tmp.gt.forcmx) then
                    forcmx=tmp
                endif
            endif
        enddo
    enddo
endif
!
rewind(lu)

do n=1,meshdatactx%Nvglobal

    if (mattsym.eq.1 .or. mattsym.eq.2) then

        ! Use relevant bc magnitude
        bcmag1 = 0
        bcmag2 = 0
        if (IBOND(1,n).eq.iflg) bcmag1 = BOND(1,n)
        if (IBOND(2,n).eq.iflg) bcmag2 = BOND(2,n)

        ! position of node
        xl(1) = X(1,n)+DEFFAC*D(1,n)
        xl(2) = X(2,n)+DEFFAC*D(2,n)

        ! length of boundary condition
        length = sqrt(bcmag1*bcmag1+bcmag2*bcmag2)

        ! orientation of boundary condition
        deg = atan2(BOND(2,n),BOND(1,n))
        deg = (deg - SKEW(1,n))/deg2rad

        if (mattsym.eq.1) then
            angleA = deg
            angleB = deg+90.
        else 
            angleA = 90.-deg
            angleB = -deg
        endif

        ! Displacement boundary conditions (iflg=1)
        if (iflg.eq.1) then
            if (IBOND(1,n).eq.iflg.and.IBOND(2,n).eq.iflg) then
                if (length.lt.1e-20) then
                    write(LU,*,err=2000) xl(1),xl(2),1,1,0
                elseif (sysc) then
                    write(LU,*,err=2000) xl(1),xl(2),length/dispmx,3,angleA
                else
                    write(LU,*,err=2000) xl(1),xl(2),1,3,angleA
                endif
            elseif (IBOND(1,n).eq.iflg) then
                if (length.lt.1e-20) then
                    write(LU,*,err=2000) xl(1),xl(2),1,2,angleA
                elseif (sysc) then
                    write(LU,*,err=2000) xl(1),xl(2),length/dispmx,3,angleA
                else
                    write(LU,*,err=2000) xl(1),xl(2),1,3,angleA
                endif
            elseif (IBOND(2,n).eq.iflg) then
                if (length.lt.1e-20) then
                    write(LU,*,err=2000) xl(1),xl(2),1,2,angleB
                elseif (sysc) then
                    write(LU,*,err=2000) xl(1),xl(2),length/dispmx,3,angleB
                else
                    write(LU,*,err=2000) xl(1),xl(2),1,3,angleB
                endif
            endif
        elseif (iflg.eq.2) then
            if (IBOND(1,n).eq.iflg.or.IBOND(2,n).eq.iflg) then
                if (sysc) then
                    write(LU,*,err=2000) xl(1),xl(2),length/velomx,4,angleA
                else
                    write(LU,*,err=2000) xl(1),xl(2),1,4,angleA
                endif
            endif
        elseif (iflg.eq.3) then
            if (IBOND(1,n).eq.iflg.or.IBOND(2,n).eq.iflg) then
                if (sysc) then
                    write(LU,*,err=2000) xl(1),xl(2),length/forcmx,5,angleA
                else
                    write(LU,*,err=2000) xl(1),xl(2),1,5,angleA
                endif
            endif
        endif
    else

        if (GMT) then

            xl(1) = X(1,n)+DEFFAC*D(1,n)
            xl(2) = X(2,n)+DEFFAC*D(2,n)

            if (iflg.eq.6) then
                ! all
                length = ZERO
                if (IBOND(1,n).ne.0) then
                    length = length + BOND(1,n)**2
                endif
                if (IBOND(2,n).ne.0) then 
                    length = length + BOND(2,n)**2
                endif
                length = SQRT(length)

                if (IBOND(1,n).ne.0.or.IBOND(2,n).ne.0) then
                    if (length.lt.1e-36) then
                        deg = atan2(DBLE(IBOND(2,n)),DBLE(IBOND(1,n)))
                    else
                        deg = atan2(BOND(2,n),BOND(1,n))
                    endif
                    deg = (deg - SKEW(1,n))/deg2rad

                    if (SYSC) then
                        write(LU,6,err=2000) xl(1),xl(2),deg,length
 6                      format(E15.6,1X,E15.6,1X,E15.6,1X,E15.6)
                    else
                        write(LU,6,err=2000) xl(1),xl(2),deg,1.0
                    endif
                endif
            else
                if (IBOND(1,n).eq.iflg .and. &
                    IBOND(2,n).eq.iflg) then
                    ! both x and y have this boundary condtion
! here we have boundry conditions fixed in multiple directions.
! draw an arrow in each direction, in stead of one arrow in 45/45 direction.
!                length = BOND(1,n)**2 + BOND(2,n)**2
!                length = SQRT(length)

                    if (length.lt.1e-36) then
                        deg = atan2(DBLE(IBOND(2,n)),DBLE(IBOND(1,n)))
                    else
                        deg = atan2(BOND(2,n),BOND(1,n))
                    endif

                    deg = (deg - SKEW(1,n))/deg2rad

                    if (SYSC) then
                        write(LU,6,err=2000) xl(1),xl(2),deg+45.0,length
                        write(LU,6,err=2000) xl(1),xl(2),deg-45.0,length
                    else
                        write(LU,6,err=2000) xl(1),xl(2),deg+45.0,1.0
                        write(LU,6,err=2000) xl(1),xl(2),deg-45.0,1.0
                    endif
                else if (IBOND(1,n).eq.iflg) then
                    ! only x has this boundary condition
                    length = ABS(BOND(1,n))
                    if (length.lt.1e-36) then
                        deg = atan2(ZERO,DBLE(IBOND(1,n)))
                    else
                        deg = atan2(ZERO,BOND(1,n))
                    endif
                    deg = (deg - SKEW(1,n))/deg2rad
                    if (SYSC) then
                        write(LU,6,err=2000) xl(1),xl(2),deg,length
                    else
                        write(LU,6,err=2000) xl(1),xl(2),deg,1.0
                    endif
                else if (IBOND(2,n).eq.iflg) then
                    ! only y has this boundary condition
                    length = ABS(BOND(2,n))
                    if (length.lt.1e-36) then
                        deg = atan2(DBLE(IBOND(2,n)),ZERO)
                    else
                        deg = atan2(BOND(2,n),ZERO)
                    endif
                    deg = (deg - SKEW(1,n))/deg2rad

                    if (SYSC) then
                        write(LU,6,err=2000) xl(1),xl(2),deg,length
                    else
                        write(LU,6,err=2000) xl(1),xl(2),deg,1.0
                    endif
                endif
            endif
        else !GMT
            stop "non-GMT plotting not implemented"
        endif

    endif ! MATTMODE CONDITION

enddo
!
return
2000    write(stderr,2001)
 2001    format(1x,'NBCPLT write error')
end subroutine
!-------------------------------------------------------------------------------
subroutine bcsym1 (IBOND,BOND,IEN,SKEW,X,D,NODE,ICOMP,LU,SIZE, &
                   DEFORM,DEFFAC)

USE ALGEBRA,         only: formrt, clear, EQUATE, vmprd
USE MODELDEFINITION
USE CONSTANTS,       only: pi, deg2rad, zero, one,eps
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact
use symbols,         only: maxSym, mfix, mforc, mslid, mvelo, &
                           syfix, syforc, syslid, syvelo
!
implicit none
!-pass
logical DEFORM
integer NODE,ICOMP,LU
double precision SIZE,DEFFAC

integer   :: ibond, ien
double precision :: bond, d, x, skew
dimension IBOND(NDOF,*),BOND(NDOF,*),IEN(NEN,*),SKEW(2,*), D(NDOF,*),X(NSD,*)

!-locl

logical found
integer ielm,j,k,i,nsymbl
double precision dx1,dy1,dx,dy,s,xinpr1,xinpr2
double precision :: bc(2),rot(3,3),temp(2)
double precision :: xl(2,4),dl(12)
double precision :: angle
double precision, allocatable :: symbl(:,:)

allocate(symbl(2, maxSym))

!
if (ICOMP.le.NDOF .and. &
    IBOND(ICOMP,NODE).eq.0) then
     return
endif
!
ielm = 1
found = .false.
do while (ielm.le.NUMEL .and. (.not.found))
    j = 1
    do while (j.le.NEN .and. (.not.found))
        found = (IEN(j,ielm).eq.NODE)
        if (.not.found) then
            j = j + 1
        endif
    enddo
    if (.not.found) then
        ielm = ielm + 1
    endif
enddo

if (.not.found) then
    write(stderr,10) NODE
10       format(1x,'bc fatal error: node ',I5,' not in mesh')
    call exitp(1)
endif

 k = j + 1
 i = j - 1
 if (IEN(3,ielm).eq.IEN(4,ielm)) then
if (k.gt.3) k = 1
if (i.lt.1) i = 3
 else
if (k.gt.4) k = 1
if (i.lt.1) i = 4
 endif
!
 call LCOORD (X,xl,IEN(1,ielm))
 if (ICOMP.gt.NDOF) goto 100
!
!    decide to which of the two possible face normals the bc vector
!    is closest by setting up normal and bc vectors and taking inner product
 call CLEAR (bc,2,"bc")
 bc(ICOMP) = ONE
 if (SKEW(1,NODE).ne.ZERO) then
call FORMRT (SKEW(1,NODE),rot,ndof)
call EQUATE (temp,bc,2)
call VMPRD (rot,temp,bc,2,3)
 endif
 dx1 = xl(1,j)-xl(1,i)
 dy1 = xl(2,j)-xl(2,i)
 s = SQRT(dx1*dx1+dy1*dy1)
!    normal vector (dy1,-dx1)
 if (s.gt.EPS) then
xinpr1 = (dy1*bc(1)-dx1*bc(2))/s
 else
xinpr1 = ZERO
 endif
 dx = xl(1,k)-xl(1,j)
 dy = xl(2,k)-xl(2,j)
 s = SQRT(dx*dx+dy*dy)
 if (s.gt.EPS) then
xinpr2 = (dy*bc(1)-dx*bc(2))/s
 else
xinpr2 = ZERO
 endif
 if (ABS(xinpr1).ge.ABS(xinpr2)) then
dx = dx1
dy = dy1
 endif
!
 call CLEAR (symbl,2*MAXSYM,"symbl")
 if (IBOND(ICOMP,NODE).eq.1 .or.(IBOND(ICOMP,NODE).eq.4.and. &
   NSTEP.eq.0)) then
if (BOND(ICOMP,NODE).eq.ZERO) then
    call CPSYM (SYSLID,symbl,MSLID)
    nsymbl = MSLID
else
    call CPSYM (SYVELO,symbl,MVELO)
nsymbl = MVELO
endif
 elseif (IBOND(ICOMP,NODE).eq.2) then
    call CPSYM (SYVELO,symbl,MVELO)
nsymbl = MVELO
elseif (IBOND(ICOMP,node).eq.3) then
    call CPSYM (SYFORC,symbl,MFORC)
nsymbl = MFORC
 endif
!
!    Orient bc symbol (initially oriented in +x -direction)
 if (ICOMP.eq.1) then
if (dy.lt.ZERO) then
if (BOND(1,NODE).gt.ZERO) then
    call NEGABC (symbl,nsymbl)
else if (BOND(1,NODE).le.ZERO) then
    call ROTSYM (symbl,nsymbl,1.8d2)
endif
else
if (BOND(1,NODE).lt.ZERO) then
    call NEGABC (symbl,nsymbl)
    call ROTSYM (symbl,nsymbl,1.8d2)
endif
endif
 else
if (dx.gt.ZERO) then
if (BOND(2,NODE).gt.ZERO) then
    call NEGABC (symbl,nsymbl)
else if (BOND(2,NODE).le.ZERO) then
    call ROTSYM (symbl,nsymbl,1.8d2)
endif
else
if (BOND(2,NODE).lt.ZERO) then
    call NEGABC (symbl,nsymbl)
    call ROTSYM (symbl,nsymbl,1.8d2)
endif
endif
call ROTSYM (symbl,nsymbl,90d0)
 endif
!
100    if (ICOMP.gt.NDOF) then
call CLEAR (symbl,2*MAXSYM,"symbl")
if (IBOND(1,NODE).eq.1) then
call CPSYM (SYFIX,symbl,MFIX)
    nsymbl = MFIX
else if (IBOND(1,NODE).eq.2) then
call CPSYM (SYVELO,symbl,MVELO)
    nsymbl = MVELO
else if (IBOND(1,NODE).eq.3) then
call CPSYM (SYFORC,symbl,MFORC)
    nsymbl = MFORC
endif
dx = BOND(1,NODE)
dy = BOND(2,NODE)
angle = ATAN2(dy,dx)/deg2rad
call ROTSYM(symbl,nsymbl,angle)
 endif
!
 if (SKEW(1,NODE).ne.ZERO) call ROTSYM (symbl,nsymbl, &
   -SKEW(1,NODE)/deg2rad)
!
!    Scale symbol
if (DEFORM) then
    call LDISP (dl,D,IEN(1,ielm),NDOF,NEN)
    call REZONE (xl,dl,DEFFAC)
endif

call SYSCAL (xl(1,j),symbl,SIZE)
!
!    Write outcome
 call WRITBC (symbl,nsymbl,LU)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine bcsym2 (IBOND,BOND,IEN,SKEW,X,D,NODE,LU,SIZE, &
   DEFORM,DEFFAC)

USE ALGEBRA
USE MODELDEFINITION
USE CONSTANTS,       only: pi, deg2rad, zero,eps
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact
use symbols,         only: maxSym, mfix, mforc, syfix, syvelo, syforc
!
implicit none
!-pass
 logical DEFORM
 integer NODE,LU
 double precision SIZE,DEFFAC
integer :: ibond, ien
double precision :: bond, skew, d, x
 dimension IBOND(NDOF,*),BOND(NDOF,*),IEN(NEN,*),SKEW(2,*), &
   D(NDOF,*),X(NSD,*)

!-locl

logical found
double precision :: xl(2,4),dl(12),c(2)
integer          :: i, j, ii, jj, ielm, jelm, ni, nj, j2
double precision :: dx, dy, s, angle
integer          :: idummy
double precision :: dummy
integer          :: ierr, nsymbl
double precision, allocatable :: symbl(:,:)

allocate(symbl(2, maxSym))

!
 ielm = 1
 found = .false.
 do while (ielm.le.NUMEL .and. (.not.found))
j = 1
do while (j.le.NEN .and. (.not.found))
found = (IEN(j,ielm).eq.NODE)
if (.not.found) j = j + 1
enddo
if (.not.found) ielm = ielm + 1
 enddo
 if (.not.found) then
write(stderr,10) NODE
10       format(1x,'bc fatal error: node ',I5,' not in mesh')
    call exitp(1)
 endif
!
 call LCOORD (X,xl,IEN(1,ielm))
!
 jelm = ielm
100    jelm = jelm + 1
 found = .false.
 do while (jelm.le.NUMEL .and. (.not.found))
j2 = 1
do while (j2.le.NEN .and. (.not.found))
found = (IEN(j2,jelm).eq.NODE)
if (.not.found) j2 = j2 + 1
enddo
if (.not.found) jelm = jelm + 1
 enddo
 if (found) then
i = 0
ni = 4
nj = 4
if (IEN(3,ielm).eq.IEN(4,ielm)) ni = 3
if (IEN(3,jelm).eq.IEN(4,jelm)) nj = 3
do ii=1,ni
if (ii.ne.j) then
    do jj=1,nj
        if (IEN(ii,ielm).eq.IEN(jj,jelm)) i = ii
    enddo
endif
enddo
if (i.eq.0) goto 100
!
!        set up bc direction vector
dx = xl(1,j) - xl(1,i)
dy = xl(2,j) - xl(2,i)
s = SQRT(dx*dx+dy*dy)
dx = dx / s
dy = dy / s
 else
!        boundary condition on a corner node
call ELCTR (c,IEN(1,ielm),X,dummy,idummy,dummy, &
     idummy,dummy,idummy,dummy,SKEW,ZERO,ierr)
dx = xl(1,j) - c(1)
dy = xl(2,j) - c(2)
s = SQRT(dx*dx+dy*dy)
dx = dx / s
dy = dy / s
 endif
!
 call CLEAR (symbl,2*MAXSYM,"symbl")
 if (IBOND(1,NODE).eq.1 .or.(IBOND(1,NODE).eq.4.and. &
   NSTEP.eq.0)) then
call CPSYM (SYFIX,symbl,MFIX)
nsymbl = MFIX
 elseif (IBOND(1,NODE).eq.2) then
call CPSYM (SYFIX,symbl,MFIX)
nsymbl = MFIX
elseif (IBOND(1,node).eq.3) then
    call CPSYM (SYFORC,symbl,MFORC)
nsymbl = MFORC
 endif
!
!    Orient bc symbol (initially oriented in +x -direction)
 angle = ATAN2(dy,dx)/deg2rad
 call ROTSYM (symbl,nsymbl,angle)
 if (SKEW(1,NODE).ne.ZERO) call ROTSYM (symbl,nsymbl, &
   -SKEW(1,NODE)/deg2rad)
!
!    Scale symbol
 if (DEFORM) then
call LDISP (dl,D,IEN(1,ielm),NDOF,NEN)
call REZONE (xl,dl,DEFFAC)
 endif
call SYSCAL (xl(1,j),symbl,SIZE)
!
!    Write outcome
 call WRITBC (symbl,nsymbl,LU)
!
 return
 end
!-------------------------------------------------------------------------------
!subroutine slipperybc(IEN,X,SKEW,NSLIP)
!implicit none
! USE CONSTANTS, only: pi, deg2rad, zero,eps
! USE !!!MODELDEFINITION
! USE MESHDATAMODULE, only: meshdatactx
!return
!end
!-------------------------------------------------------------------------------
subroutine CPSYM (SYMBL1,SYMBL2,N)

USE ALGEBRA, only: clear
use symbols, only: maxSym
!
implicit none
!-pass
integer :: n
double precision :: SYMBL1,SYMBL2(2,maxSym)
dimension SYMBL1(2,*)

!-local
integer :: j
!
call CLEAR (SYMBL2,2*MAXSYM,"symbl2")

do j=1,N
    SYMBL2(1,j) = SYMBL1(1,j)
    SYMBL2(2,j) = SYMBL1(2,j)
enddo
!
 return
 end
!-------------------------------------------------------------------------------
subroutine NEGABC (SYMBL,MAXSYM)

use constants, only: one

implicit none
!-pass
integer          :: maxsym
double precision :: SYMBL(2,MAXSYM)
!-local
integer          :: m, n, i, j

m = 1
n = INT(SYMBL(2,m))

do while (n.gt.0)
    do i=1,n
        j = m + i
        SYMBL(1,j) = SYMBL(1,j) - ONE
    enddo

    m = m + n + 1
    n = 0

    if (m.lt.MAXSYM) then
        n = INT(SYMBL(2,m))
    endif
enddo
!
return
end
!-------------------------------------------------------------------------------
 subroutine FLIPBC (SYMBL,MAXSYM)

implicit none
!-pass
integer          :: maxsym
double precision :: SYMBL(2,MAXSYM)
!-local
integer          :: m, n, i, j

m = 1
n = INT(SYMBL(2,m))
do while (n.gt.0)
    do i=1,n
        j = m + i
        SYMBL(1,j) = -SYMBL(1,j)
    enddo
    m = m + n + 1
    n = 0
    if (m.lt.MAXSYM) then
        n = INT(SYMBL(2,m))
    endif
enddo
!
return
end
!-------------------------------------------------------------------------------
 subroutine ROTSYM (SYMBL,MAXSYM,angle)
 USE CONSTANTS, only: deg2rad
!
! rotates symbol over angle, around (0.,0.)
!
implicit none
!-pass
integer          :: maxsym
double precision :: SYMBL(2,MAXSYM)
!-local
integer          :: m, n, i, j
double precision :: rangl, ca, sa, x, y, angle


rangl = angle * deg2rad
ca = COS(rangl)
sa = SIN(rangl)
m = 1
n = INT(SYMBL(2,m))

do while (n.gt.0)
    do i=1,n
        j = i + m
        if (j.gt.MAXSYM) then
            stop 'ROTSYM: overflow'
        endif
        X = SYMBL(1,j)
        Y = SYMBL(2,j)
        SYMBL(1,j) = X*ca-Y*sa
        SYMBL(2,j) = X*sa+Y*ca
    enddo
    m = m + n + 1
    n = 0
    if (m.lt.MAXSYM) then
        n = INT(SYMBL(2,m))
    endif
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine SYSCAL (XL,SYMBL,SYSIZE)

! moves and rescales the symbol

use symbols, only: maxSym

implicit none
!-pass
double precision SYSIZE

double precision :: symbl, xl
dimension SYMBL(2,*),XL(*)


!-local
integer          :: m, n, i, j,k

!
k = 1
m = INT(SYMBL(2,k))
do while (m.gt.0)
    do i=1,m
        j = i + k
        SYMBL(1,j) = XL(1) + SYMBL(1,j)*SYSIZE
        SYMBL(2,j) = XL(2) + SYMBL(2,j)*SYSIZE
    enddo
    k = k + m + 1
    m = 0
    if (k.lt.MAXSYM) then
        m = INT(SYMBL(2,k))
    endif
enddo
!
 return
 end
!-------------------------------------------------------------------------------
! MATT ADDITION
subroutine BCDEF(mattsym)
implicit none
integer mattsym
logical firstentry
data firstentry/.true./
save firstentry
! if (firstentry) write(0,'(100a)') 'GMT custom symbols by Matt Herman'
firstentry = .false.
    write(0,*) 'Writing boundary condition GMT macro file "bcs.def"'
    open(unit=511,file='bcs.def',status='unknown')
    if (mattsym.eq.1) write(511,'(A)') 'N: 2 or'
    if (mattsym.eq.2) write(511,'(A)') 'N: 2 oa'
    write(511,'(A)') '$2 R'
! Symbol 1: fixed node in 2D boundary condition (fixed triangle)
    write(511,'(A)') 'if $1 == 1 then {'
    write(511,'(A)') '    0 0 M'
    write(511,'(A)') '    -0.787273 0.454545 D'
    write(511,'(A)') '    -0.787273 0.272727 D'
    write(511,'(A)') '    -1 0.363636 D'
    write(511,'(A)') '    -0.787273 0.272727 D'
    write(511,'(A)') '    -0.787273 0.0909091 D'
    write(511,'(A)') '    -1 0.181818 D'
    write(511,'(A)') '    -0.787273 0.0909091 D'
    write(511,'(A)') '    -0.787273 -0.0909091 D'
    write(511,'(A)') '    -1 0 D'
    write(511,'(A)') '    -0.787273 -0.0909091 D'
    write(511,'(A)') '    -0.787273 -0.272727 D'
    write(511,'(A)') '    -1 -0.181818 D'
    write(511,'(A)') '    -0.787273 -0.272727 D'
    write(511,'(A)') '    -0.787273 -0.454545 D'
    write(511,'(A)') '    0 0 D'
! Symbol 2: fixed node in 1D boundary condition (roller triangle)
    write(511,'(A)') '} elseif $1 == 2 then {'
    write(511,'(A)') '    0 0 M'
    write(511,'(A)') '    -0.787273 0.454545 D'
    write(511,'(A)') '    -0.787273 -0.454545 D'
    write(511,'(A)') '    0 0 M'
    write(511,'(A)') '    -0.890000 0.227273 0.20 c -G-'
    write(511,'(A)') '    -0.890000 -0.227273 0.20 c -G-'
! Symbol 3: non-zero nodal displacement boundary condition (skinny arrow)
    write(511,'(A)') '} elseif $1 == 3 then {'
    write(511,'(A)') '    -1.00 0.01 M -Gblack' ! for GMT 5.2.1; in GMT 5.4, use "-G+p" to set fill color to pen color
    write(511,'(A)') '    -0.10 0.01 D'
    write(511,'(A)') '    -0.10 -0.01 D'
    write(511,'(A)') '    -1.00 -0.01 D'
    write(511,'(A)') '    -1.00 0.01 D'
    write(511,'(A)') '    0.00 0.00 M -Gblack' ! for GMT 5.2.1; in GMT 5.4, use "-G+p" to set fill color to pen color
    write(511,'(A)') '    -0.20 0.10 D'
    write(511,'(A)') '    -0.10 0.00 D'
    write(511,'(A)') '    -0.20 -0.10 D'
    write(511,'(A)') '    0.00 0.00 D'
! Symbol 4: non-zero nodal velocity boundary condition (hollow skinny arrow)
    write(511,'(A)') '} elseif $1 == 4 then {'
    write(511,'(A)') '    -1.00 0.03 M'
    write(511,'(A)') '    -0.12 0.03 D'
    write(511,'(A)') '    -0.31 0.13 D'
    !write(511,'(A)') '    S'  ! in GMT 5.4, can use S to stop drawing line segment
    write(511,'(A)') '    -0.22 0.13 D'
    write(511,'(A)') '    0.00 0.00 D'
    write(511,'(A)') '    -0.22 -0.13 D'
    !write(511,'(A)') '    S'  ! in GMT 5.4, can use S to stop drawing line segment
    write(511,'(A)') '    -0.31 -0.13 D'
    write(511,'(A)') '    -0.12 -0.03 D'
    write(511,'(A)') '    -1.00 -0.03 D'
    !write(511,'(A)') '    S'  ! in GMT 5.4, can use S to stop drawing line segment
! Symbol 5: nodal force (arrow with letter "F" on it)
    write(511,'(A)') '} elseif $1 == 5 then {'
    write(511,'(A)') '    -1.00 0.01 M -Gblack' ! for GMT 5.2.1; in GMT 5.4, use "-G+p" to set fill color to pen color
    write(511,'(A)') '    -0.10 0.01 D'
    write(511,'(A)') '    -0.10 -0.01 D'
    write(511,'(A)') '    -1.00 -0.01 D'
    write(511,'(A)') '    -1.00 0.01 D'
    write(511,'(A)') '    0.00 0.00 M -Gblack' ! for GMT 5.2.1; in GMT 5.4, use "-G+p" to set fill color to pen color
    write(511,'(A)') '    -0.20 0.10 D'
    write(511,'(A)') '    -0.10 0.00 D'
    write(511,'(A)') '    -0.20 -0.10 D'
    write(511,'(A)') '    0.00 0.00 D'
    write(511,'(A)') '    -0.30 -0.16 M -G-'
    write(511,'(A)') '    -0.30 -0.06 D'
    write(511,'(A)') '    -0.40 -0.06 D'
    write(511,'(A)') '    -0.40 -0.16 D'
    write(511,'(A)') '    -0.40 -0.06 D'
    write(511,'(A)') '    -0.50 -0.06 D'
    write(511,'(A)') '    -0.30 -0.06 D'
! Symbol 11: element stress boundary condition, compressive/negative (fat arrow)
    write(511,'(A)') '} elseif $1 == 11 then {'
    write(511,'(A)') '    -0.50  0.08 M -Gblack' ! for GMT 5.2.1; in GMT 5.4, use "-G+p" to set fill color to pen color
    write(511,'(A)') '    -0.18  0.08 D'
    write(511,'(A)') '    -0.18  0.16 D'
    write(511,'(A)') '     0.00  0.00 D'
    write(511,'(A)') '    -0.18 -0.16 D'
    write(511,'(A)') '    -0.18 -0.08 D'
    write(511,'(A)') '    -0.50 -0.08 D'
    write(511,'(A)') '    -0.50  0.08 D'
! Symbol 12: element stress boundary condition, dilational/positive (fat arrow)
    write(511,'(A)') '} elseif $1 == 12 then {'
    write(511,'(A)') '    0.00  0.08 M -Gblack' ! for GMT 5.2.1; in GMT 5.4, use "-G+p" to set fill color to pen color
    write(511,'(A)') '    0.32  0.08 D'
    write(511,'(A)') '    0.32  0.16 D'
    write(511,'(A)') '    0.50  0.00 D'
    write(511,'(A)') '    0.32 -0.16 D'
    write(511,'(A)') '    0.32 -0.08 D'
    write(511,'(A)') '    0.00 -0.08 D'
    write(511,'(A)') '    0.00  0.08 D'
! Symbol 13: element Winkler pressure (spring)
    write(511,'(A)') '} elseif $1 == 13  then {'
    write(511,'(A)') '    0 0 M'
    write(511,'(A)') '    0 -0.1 D'
    write(511,'(A)') '    -0.20 -0.16 D'
    write(511,'(A)') '    0.20 -0.28 D'
    write(511,'(A)') '    -0.20 -0.40 D'
    write(511,'(A)') '    0.20 -0.52 D'
    write(511,'(A)') '    -0.20 -0.64 D'
    write(511,'(A)') '    0 -0.70 D'
    write(511,'(A)') '    0 -0.85 D'
    write(511,'(A)') '    S'
    write(511,'(A)') '    -0.20 -0.85 M'
    write(511,'(A)') '    0.20 -0.85 D'
    write(511,'(A)') '    S'
    write(511,'(A)') '    -0.05 -0.85 M'
    write(511,'(A)') '    -0.20 -1.00 D'
    write(511,'(A)') '    S'
    write(511,'(A)') '    0.05 -0.85 M'
    write(511,'(A)') '    -0.10 -1.00 D'
    write(511,'(A)') '    S'
    write(511,'(A)') '    0.15 -0.85 M'
    write(511,'(A)') '    0.00 -1.00 D'
    write(511,'(A)') '    S'
! Symbol 21: element slippery node in x direction (circle offset parallel to diameter)
    write(511,'(A)') '} elseif $1 == 21 then {'
    write(511,'(A)') '    0 0 M'
    write(511,'(A)') '    -0.20 0.0 0.5 180 0 A'
    write(511,'(A)') '    0.20 0.0 0.5 -180 0 A'
! Symbol 22: element slippery node in y direction (circle offset perpendicular to diameter)
    write(511,'(A)') '} elseif $1 == 22 then {'
    write(511,'(A)') '    0 0.08 M'
    write(511,'(A)') '    0.0 0.08 0.5 180 0 A'
    write(511,'(A)') '    0 -0.08 M'
    write(511,'(A)') '    0.0 -0.08 0.5 -180 0 A'
! Symbol 23: element slippery node in both directions (circle offset obliquely)
    write(511,'(A)') '} elseif $1 == 23 then {'
    write(511,'(A)') '    -0.20 0.08 M'
    write(511,'(A)') '    -0.20 0.08 0.5 180 0 A'
    write(511,'(A)') '    0.20 -0.08 M'
    write(511,'(A)') '    0.20 -0.08 0.5 -180 0 A'
! All others: (diamond)
    write(511,'(A)') '} else {'
    write(511,'(A)') '    0 0 0.5 d'
    write(511,'(A)') '}'
close(511)
return
end
