#ifndef EXPERIMENTAL_NORATESTORAGE
subroutine trian (IEN,X,D,STN,STR,RATE,LMF,LMX, &
                  TFAULT,MAT,PRPMAT,PRPLAS,DXE,SKEW, &
                  iflag,DEFORM,deffac,triaout)
#else
subroutine trian (IEN,X,D,STN,STR,LMF,LMX, &
                  TFAULT,MAT,PRPMAT,PRPLAS,DXE,SKEW, &
                  iflag,DEFORM,deffac,triaout)
#endif

USE ALGEBRA
USE MODELDEFINITION
use materialsmodule, only: MaterialIDOfElement, effsign
use debugmodule,     only: iecho
use modeltopology,   only: NEN,NSD,NDOF,NSTR
use iomodule,        only: stderr
use plotControl,     only: ratlog, viscli, vislog
! invokes "triangulate" to generate smoothly interpolated triangles
#ifdef EXPERIMENTAL_NORATESTORAGE
use modeldatamodule, only: modeldatactx
use  meshdatamodule, only:  meshdatactx
#endif

implicit none
!-pass



integer iflag

integer :: ien, lmf, lmx, mat
double precision :: x, d, stn, str, tfault, prpmat, prplas, skew, dxe


#ifndef EXPERIMENTAL_NORATESTORAGE
dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),STN(NSTR,*),STR(NSTR,*), &
   LMF(NDOF,NEN,*),LMX(NDOF,NEN,*),TFAULT(NDOF,*), &
   MAT(*),PRPMAT(6,*),PRPLAS(9,*),RATE(NSTR,*),SKEW(2,*), &
 DXE(NDOF,NEN,*)
double precision :: rate
#else
dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),STN(NSTR,*),STR(NSTR,*), &
   LMF(NDOF,NEN,*),LMX(NDOF,NEN,*),TFAULT(NDOF,*), &
   MAT(*),PRPMAT(6,*),PRPLAS(9,*),SKEW(2,*), &
 DXE(NDOF,NEN,*)
#endif

double precision :: deffac
logical DEFORM
character triaout*(*)
!-locl
#include "lenchar.i"
 integer csystem
 character(len=LENCHAR) cmd
double precision :: CENTR(3),EE(6),BETA(6),BETB(6),XL(24)
 integer, external :: lnblk, nextlu
double precision, external :: effee, effstn, etalin

double precision :: anpwr, eff, emhu, def
integer          :: i, n, nn, m, lu, ierr

#ifdef EXPERIMENTAL_NORATESTORAGE
double precision :: strainRate(NSTR)
#endif

!-init
def  = 0d0

if (DEFORM) then
    def = deffac
endif

if (IDIG(ICVIS,4,1).ne.0) then
    call CLEAR (BETB,NSTR,"BETB")
endif

lu = nextlu(0)
call openf(lu,'/tmp/tria.dat','unknown')

do nn=1,NUMEL

#ifdef EXPERIMENTAL_NORATESTORAGE
    call StrainrateOfElem (meshdatactx%X,&
                           modeldatactx%D,&
                           modeldatactx%DELD, &
                            strainRate,&
                            meshdatactx%IEN,&
                           modeldatactx%LMF, &
                           modeldatactx%TFAULT,&
                           modeldatactx%DFAULT,&
                           modeldatactx%LMX, &
                           modeldatactx%DX, &
                           modeldatactx%SKEW,&
                           modeldatactx%DXE, &
                            nn)
#endif


    n = nn
    m = MaterialIDOfElement(n)
    if (iflag.eq.0) then
!        contour effective stress
        eff = effstn (STN(1,n),PRPMAT(2,m),effsign)
    else if (iflag.eq.1) then
!        contour effective strain
        if (ratlog) then
            eff = effee (STR(1,n),PRPMAT(2,m),.false.)
            eff = TENLOG(ABS(eff))
        else
            eff = effee (STR(1,n),PRPMAT(2,m),effsign)
        endif
    else if (iflag.eq.2) then
!        contour effective total rate
#ifdef EXPERIMENTAL_NORATESTORAGE
        if (ratlog) then
            eff = effee (strainRate(1),PRPMAT(2,m),.false.)
            eff = TENLOG(ABS(eff))
        else
            eff = effee (strainRate(1),PRPMAT(2,m),effsign)
        endif
#else
        if (ratlog) then
            eff = effee (RATE(1,n),PRPMAT(2,m),.false.)
            eff = TENLOG(ABS(eff))
        else
            eff = effee (RATE(1,n),PRPMAT(2,m),effsign)
        endif
#endif
    else if (iflag.eq.3) then
!        contour effective viscous rate
        call FORMBT (STN(1,n),BETA,BETB,PRPMAT(1,m),PRPLAS(1,m))
        BETA(3) = BETA(3)*5D-1
        BETB(3) = BETB(3)*5D-1
        if (IDIG(ICVIS,4,2).ne.0) then
            do i=1,NSTR
                BETA(i) = BETA(i) - BETB(i)
            enddo
        endif
        if (ratlog) then
            eff = effee (BETA,5D-1,.false.)
            eff = TENLOG(eff)
        else
            eff = effee (BETA,5D-1,effsign)
        endif
    else if (iflag.eq.4) then
        if (IDIG(ICVIS,4,2).ne.0) then
!            contour effective plastic rate
            call FORMBT (STN(1,n),BETA,BETB,PRPMAT(1,m), &
                         PRPLAS(1,m))
            BETA(3) = BETA(3)*5D-1
            BETB(3) = BETB(3)*5D-1
            if (ratlog) then
                eff = effee (BETB,PRPMAT(2,m),.false.)
                eff = TENLOG(eff)
            else
                eff = effee (BETB,PRPMAT(2,m),effsign)
            endif
        else
            eff = 0d0
        endif
    else if (iflag.eq.5) then
!        contour effective viscosity
        emhu  = PRPMAT(3,m)
        anpwr = PRPMAT(4,m)
        if (viscli) then
            eff = effstn (STN(1,n),5D-1,.false.)
            eff = etaLin (eff,emhu,anpwr)
        else
            eff = EMHU
        endif
        if (vislog)  then
            eff = TENLOG(eff)
        endif
    else
        write(stderr,1)
 1        format(1x,'trian: fatal iflag error')
        call exitp(1)
    endif

    call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
                SKEW,def, ierr)

    if (ierr.ne.0) then
        call exitp(1)
    endif

    write(lu,*) CENTR(1),CENTR(2),eff

enddo

call closef(lu)

cmd = 'triangulate -p /tmp/tria.dat > '// &
        triaout(1:lnblk(triaout))

if (iecho.eq.1) then
    write(stderr,2,advance='no')
 2        format(1x,'performing triangulation ... ')
endif

ierr = csystem(cmd(1:lnblk(cmd)))
if (ierr.ne.0) then
    call exitp(1)
endif

cmd  = 'rm -f /tmp/tria.dat'
ierr = csystem(cmd(1:lnblk(cmd)))
if (ierr.ne.0) then
    call exitp(1)
endif

if (iecho.eq.1) then
    write(stderr,3)
 3        format(1x,'done')
endif

return
end
