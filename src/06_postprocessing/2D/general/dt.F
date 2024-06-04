subroutine DtChk (STN,MAT,PRPMAT,IJSUR,XPT,IEN, &
NEN,NSTR,NSUR,NPLDOF,NPLT,verbose,IFLAG,luMon,CST)

use timestepmodule,  only: DELTP, alfap
use materialsmodule, only: matdatactx, &
                           MaterialIDOfElement
use modeldatamodule, only: modeldatactx
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel

!    USE MODELDEFINITION

! Calculate time step characteristics:
!     IFLAG = 0, stabilitity ratio
!     IFLAG = 1, oscillation ratio
!     IFLAG = 2, relative modal accuracy

!    XPT (x-y-f, number of nodes on a triangle, number of triangles)
!    NTRIA number of triangles that can be generated from one element

implicit none

integer, parameter :: NTRIA=2
!-pass
integer            :: NSUR,NPLDOF,NPLT,NSTR,NEN,luMon, iFlag

integer            :: mat, ien, ijsur
double precision   :: stn, prpmat, xpt
dimension          :: STN(NSTR,*),MAT(*),PRPMAT(6,*), &
                      IEN(NEN,*),IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*)

logical verbose,CST

!-locl
integer            :: nn,n,m,NUMEL,NUMTR,ntr,mm, i
double precision   :: lambda
double precision   :: ITRIA(NTRIA,2)
logical            :: QUAD

double precision   :: E, pois, emhu, anpwr, dT

double precision, external :: EigMax
real,             external :: D2REAL
double precision, external :: Stabil, Oscill, Accur

if (elmon) then
    if (verbose) then
        write(luMon,1) DELTP
 1       format(1x,'Actual time step in FE-calculation Dt = ', &
     1PG15.6)
    endif

    do nn=1,nnel
        n      = monel(nn)
        m      = MaterialIDOfElement(n)
        E      = matdatactx%PRPMAT(1,m)
        POIS   = matdatactx%PRPMAT(2,m)
        EMHU   = matdatactx%PRPMAT(3,m)
        ANPWR  = matdatactx%PRPMAT(4,m)
        lambda = EigMax (modeldatactx%STN(1,n),E,POIS,EMHU,ANPWR,ALFAP)

        if (IFLAG.eq.0) then
            dT = Stabil (ALFAP,lambda,DELTP)
            if (verbose) then
                write(luMon,2) n,D2REAL(dT)
 2                   format(1x,'element ',I5, &
                 ': stability requires 1/(',1PG16.6,') Dt ')
            endif
        else if (IFLAG.eq.1) then
            dT = Oscill (ALFAP,lambda,DELTP)
            if (verbose) then
                write(luMon,4) n,D2REAL(dT)
 4                   format(1x,'element ',I5, &
                 ': oscillation limit requires 1/(',1PG16.6, &
                 ') Dt')
            endif
        else if (IFLAG.eq.2) then
            dT = Accur (ALFAP,lambda,DELTP)
            if (verbose) then
                write(luMon,6) n,D2REAL(dT)
 6                   format(1x,'element ',I5, &
                 ': relative modal error = ',1PG16.6)
            endif
        else
            write(stderr,8)
 8            format(/1x,'DtChk: fatal IFLAG error')
            call exitp(1)
        endif

        if (.not.verbose) then
            write(luMon,*) D2REAL(dT)
        endif
    enddo
else
    NUMTR = 0
    do nn=1,NUMEL
        n      = nn
        QUAD   = (IEN(3,n).ne.IEN(4,n))
        m      = MaterialIDOfElement(n)

        E      = matdatactx%PRPMAT(1,m)
        POIS   = matdatactx%PRPMAT(2,m)
        EMHU   = matdatactx%PRPMAT(3,m)
        ANPWR  = matdatactx%PRPMAT(4,m)
        lambda = EigMax(modeldatactx%STN(1,n),E,POIS,EMHU,ANPWR,ALFAP)

        if (IFLAG.eq.0) then
            dT = Stabil (ALFAP,lambda,DELTP)
        endif
        if (IFLAG.eq.1) then
            dT = Oscill (ALFAP,lambda,DELTP)
        endif
        if (IFLAG.eq.2) then
            dT = Accur  (ALFAP,lambda,DELTP)
        endif

        if (CST) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = dT
            XPT(3,2,NUMTR) = dT
            XPT(3,3,NUMTR) = dT
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = dT
                XPT(3,2,NUMTR) = dT
                XPT(3,3,NUMTR) = dT
            endif
        else
!                determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                number of local triangles NTR
            if (ntr.gt.0) then
                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = dT
                enddo
!                loop triangles
                do mm=1,ntr
!                loop surrounding elements
                    do i=1,2
                        n      = ITRIA(i,mm)
                        m      = MaterialIDOfElement(n)

                        E      = matdatactx%PRPMAT(1,m)
                        POIS   = matdatactx%PRPMAT(2,m)
                        EMHU   = matdatactx%PRPMAT(3,m)
                        ANPWR  = matdatactx%PRPMAT(4,m)
                        lambda = EigMax (modeldatactx%STN(1,n),E,POIS, &
                                         EMHU,ANPWR,ALFAP)

                            if (IFLAG.eq.0) then
                           dT = Stabil (ALFAP,lambda,DELTP)
                        endif
                        if (IFLAG.eq.1) then
                           dT = Oscill (ALFAP,lambda,DELTP)
                        endif
                        if (IFLAG.eq.2) then
                           dT = Accur  (ALFAP,lambda,DELTP)
                        endif

                        XPT(3,i+1,NUMTR+mm) = dT
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif
return
end subroutine
