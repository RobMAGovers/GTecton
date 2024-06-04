subroutine fluxar (D,IEN,X,FLUX,LMF,TFAULT,DXE, &
                   SKEW,deffac,headsize,headangle,Fsb,Fscale,ierr, &
                   luMon,NDOF,NEN,NSD,ibase,istep,iarcont,lu,DEFORM,VERBOSE)

USE CONSTANTS,       only: zero, deg2rad
use modeldefinition, only: numel
use iomodule,        only: ascout, stderr
use postProcessData, only: GMT,connout
use monitorVars,     only: elmon, monel, nnel
use errors,          only: exec, ok
!
! routine for converting heat fluxes to arrows
!

implicit none
!-pass
logical          :: DEFORM,VERBOSE
integer          :: luMon,NDOF,NEN,NSD,ibase,istep,iarcont,lu
double precision :: deffac, headsize,headangle, Fsb,Fscale


integer          :: ien, lmf
double precision :: d, x, flux, skew, tfault, dxe
dimension D(NDOF,*),IEN(NEN,*),X(NSD,*),FLUX(NSD,*),SKEW(2,*), &
   LMF(NDOF,NEN,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*)





!-locl
integer, parameter :: MAXRRW=25
integer            :: nn,n,narrw,i,j, ierr
real               :: D2REAL
double precision   :: def, deg
double precision   :: smallevecx
double precision   :: smallevecy
double precision   :: hatangle, hatsize
double precision   :: v1, v2
double precision   :: FLUXL(3),EVAL(3),EVEC(3,3),CENTR(3),arrow(3,MAXRRW)

integer, external  :: ipltrw

!-init

data arrow(1,1),narrw/0d0,0/

ierr  = OK
def   = ZERO

if (DEFORM) then
    def = deffac
endif
!

if (elmon) then
    do nn=1,nnel
        n = monel(nn)
        call LOCFLX (FLUXL,FLUX(1,n),NSD)
        if (verbose) then
            write(luMon,1) n,(FLUXL(i),i=1,NSD)
 1               format(1x,'heat flux in element ', &
             I5,': ',3(1PG12.4))
        else
            write(luMon,*) (FLUXL(i),i=1,NSD)
        endif
    enddo
    goto 1000
endif

hatangle = headangle
hatsize  = Fsb  *  Fscale  * headsize

do nn=ibase,NUMEL,istep
    n = nn
!        calculate element center coordinates
    call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
                SKEW,def, ierr)
    if (ierr.ne.0) then
        write(stderr,2) n
 2        format(1x,'flux: error in element ',I5)
        ierr = EXEC
        goto 1000
    endif
    call LOCFLX (FLUXL,FLUX(1,n),NSD)
    call EigVec (FLUXL,EVAL,EVEC,NSD)
!        calculate arrow coordinates
    ierr = ipltrw (EVAL,EVEC,CENTR,Fscale,hatsize,hatangle, &
                   arrow,narrw,MAXRRW,iarcont,.false.)

    if (ierr.ne.OK) then
        write(stderr,4) n
 4        format(1x,'flux: error in element ',I5)
        ierr = EXEC
        goto 1000
    endif

    if (GMT) then
        if (EVAL(1).ge.EVAL(2)) then
            smallevecx = EVEC(2,1)
            smallevecy = EVEC(2,2)
            v1 = EVAL(1)
            v2 = EVAL(2)
        else
            smallevecx = EVEC(1,1)
            smallevecy = EVEC(1,2)
            v1 = EVAL(2)
            v2 = EVAL(1)
        endif

        deg = atan2(smallevecy,smallevecx)
! todo is this right?
        deg = deg/deg2rad
        v1 = v1*Fscale
        v2 = v2*Fscale

        if (ascout) then
            write(lu,*,err=2000) CENTR(1),CENTR(2),v1,v2,deg
        else
            write(lu,err=2000) D2REAL(CENTR(1)), &
                   D2REAL(CENTR(2)),D2REAL(v1),D2REAL(v2),D2REAL(deg)
        endif

    else
        stop "non-GMT plotting not implemented"
    endif

enddo

1000    return

2000    call prerr('flux')

ierr = EXEC
goto 1000
end
!-------------------------------------------------------------------------------
subroutine LOCFLX (FLUXL,FLUX,NSD)

! routine to localize fluxes

implicit none
!-pass
integer NSD, i
double precision :: FLUXL(NSD),FLUX(NSD)

do i=1,NSD
    FLUXL(i) = FLUX(i)
enddo

return
end
