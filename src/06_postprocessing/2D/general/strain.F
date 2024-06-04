subroutine strnar (STR,IEN,X,D,LMF,LMX,TFAULT, &
                   DXE,SKEW,MAT,PRPMAT,Escale,Esb,&
                   headangle,headsize, &
                   NEN,NSD,NDOF,NSTR,DEFORM,deffac,&
                   lu,ierr,ibase,istep,iarcont, &
                   verbose,luMon,sdev)

USE ALGEBRA
USE MODELDEFINITION
USE CONSTANTS,       only: zero, deg2rad, dfmax
use materialsmodule, only: matdatactx, lmat
use meshdatamodule,  only: meshdatactx
use debugmodule,     only: iecho
use iomodule,        only: ascout, stderr
use postProcessData, only: GMT,connout
use monitorVars,     only: elmon, monel, nnel
use errors,          only: ok
!
! routine for converting strains to arrows in principal stress directions
!
implicit none
!-pass
integer            :: NEN,NSD,NDOF,NSTR,lu,ierr,ibase,istep, &
                      iarcont,luMon
double precision   :: deffac, eScale, esb
double precision   :: headangle,headsize

integer            :: ien, lmf, lmx, mat
double precision   :: str, x, d, tfault, skew, dxe, prpmat

dimension          :: STR(NSTR,*),IEN(NEN,*),X(NSD,*),D(NDOF,*), &
                      LMF(NDOF,NEN,*),LMX(NDOF,NEN,*),TFAULT(NDOF,*), &
                      SKEW(2,*),DXE(NDOF,NEN,*),MAT(*),PRPMAT(6,*)

logical            :: DEFORM,verbose,sdev


!-locl
integer, parameter :: MAXRRW=25
real               :: D2REAL
double precision   :: EVAL(3),EVEC(9),CENTR(3),arrow(3,MAXRRW), &
                      EEmax(6),EEmin(6),sl(6)
double precision   :: def
double precision   :: hatangle, hatsize
integer            :: nn, n, i, j
double precision   :: pois, deg
double precision   :: v1, v2
double precision   :: smallevecx, smallevecy

integer            :: narrw

integer, external  :: ipltrw


!-init
! todo, narrw is fucky
data arrow(1,1),narrw/0d0,0/

ierr = OK
def  = ZERO


if (DEFORM) then
    def = deffac
endif

if (elmon) then
    do nn=1,nnel

        n = monel(nn)
        call EQUATE (sl,STR(1:NSTR,n),NSTR)
    
        if (sdev) then
            if (LMAT) then
                pois = matdatactx%PRPMAT(2,n)
            else
                pois = matdatactx%PRPMAT(2,meshdatactx%MAT(n))
            endif
            call STRDEV (sl,pois)
        endif

        if (verbose) then
            write(luMon,1) n,(sl(i),i=1,NSTR)
 1            format(1x,'strain in element ',I5,': ',6(1PG12.4))
        else
            write(luMon,*) (sl(i),i=1,NSTR)
        endif

    enddo
else ! no elmon
    hatangle = headangle
    hatsize  = Esb  *  Escale  * headsize
!        [km]    = [%]  x  [km/%]  x  factor

    do i=1,NSTR
        EEmin(i) =  dfmax
        EEmax(i) = -dfmax
    enddo

    do nn=ibase,NUMEL,istep
        n = nn
!             calculate element center coordinates

        call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
                        SKEW, def,ierr)

        if (ierr.ne.OK) then
            return
        endif

        call EQUATE (sl,STR(1:NSTR,n),NSTR)

        if (sdev) then
            if (LMAT) then
                pois = matdatactx%PRPMAT(2,n)
            else
                pois = matdatactx%PRPMAT(2,meshdatactx%MAT(n))
            endif
            call STRDEV (sl,pois)
        endif


        do i=1,NSTR
            if (EEmin(i) .gt. sl(i)) then
                EEmin(i)    = sl(i)
            endif
            if (EEmax(i) .lt. sl(i)) then
                EEmax(i)     =sl(i)
            endif
        enddo

!       calculate Eigenvectors
        call Eigen (sl,EVAL,EVEC)

!       calculate arrow coordinates
        ierr = ipltrw (EVAL,EVEC,CENTR,Escale, &
                       hatsize,hatangle,arrow,narrw, &
                       MAXRRW,iarcont,.true.)

        if (ierr.ne.OK) then
            write(stderr,2) n
 2            format(1x,'strain: error in element ',I5)
             return
        endif

        if (GMT) then
            if (EVAL(1).ge.EVAL(2)) then
                smallevecx = EVEC(4)
                smallevecy = EVEC(5)
                v1         = EVAL(1)
                v2         = EVAL(2)
            else
                smallevecx = EVEC(1)
                smallevecy = EVEC(2)
                v1         = EVAL(2)
                v2         = EVAL(1)
            endif


            deg = atan2(smallevecy,smallevecx)
            deg = 90.0 - deg/deg2rad
            v1  = v1 * Escale
            v2  = v2 * Escale

            if (ascout) then
                write(lu,*,err=600) CENTR(1),CENTR(2),v1,v2,deg
            else
                write(lu,err=600) D2REAL(CENTR(1)), &
                       D2REAL(CENTR(2)),D2REAL(v1),D2REAL(v2), &
                       D2REAL(deg)
            endif

        else
            stop "non-GMT plotting not implemented"
        endif

    enddo

    if (iecho.eq.1.or.iecho.ge.4) then
        call wmnmx (EEmin,EEmax,NSTR,1)
    endif

endif

500    return

600 call prerr('strain')
    call exitp(1)

end subroutine
