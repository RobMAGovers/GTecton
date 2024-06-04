subroutine efrate (IEN,STN,PRPMAT,PRPLAS,MAT,IJSUR, &
                   XPT,RATE,ratlog,verbose,iflag,luMon,CST)

USE ALGEBRA,         only: tenlog, idig, clear
USE MODELDEFINITION, only: NUMEL, icvis, numtr
use materialsmodule, only: MaterialIDOfElement, effsign
use constants,       only: dfmin, dfmax
use monitorVars,     only: elmon, nnel, monel

!
! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element
!
! iflag=10-19 total strain rate
! iflag=20-29 viscous strain rate
! iflag=30-39 plastic strain rate
! iflag last digit = 0: effective value
!       last digit = 1: Exx
!       last digit = 2: Eyy
!       last digit = 3: Exy
!       last digit = 4: Exz
!       last digit = 5: Eyz
!
implicit none

integer, parameter :: NPLT=3, NPLDOF=3, NEN=4, NSUR=4, NSTR=5
!-pass
integer            :: iflag,luMon
logical            :: ratlog,verbose,CST

integer            :: ien, mat, ijsur
double precision   :: stn, prpmat, prplas, rate, xpt
dimension          :: IEN(NEN,*),STN(NSTR,*),PRPMAT(6,*),PRPLAS(9,*), &
                      MAT(*),IJSUR(NSUR,*),RATE(NSTR,*),XPT(NPLDOF,NPLT,*)

!-locl
integer, parameter :: NTRIA=2
logical            :: QUAD
integer            :: ITRIA(NTRIA,2)
double precision   :: EE(NSTR),BETA(NSTR),BETB(NSTR)
character(len=25)  :: comment
integer            :: n, nn, m, mm, ntr, i, j, l
double precision   :: er

double precision, external ::  effee



if (IDIG(ICVIS,4,1).ne.0) then
    call CLEAR (BETB,NSTR,"BETB")
endif
!
if (elmon) then
    do nn=1,nnel
        n = monel(nn)
        m = MaterialIDOfElement(n)
        if (iflag.ge.20) then
            call FORMBT (STN(1,n),BETA,BETB,PRPMAT(1,m), &
                         PRPLAS(1,m))
!           convert engineering strain rate to regular one
            BETA(3) = BETA(3)*5D-1
            BETB(3) = BETB(3)*5D-1
            if (iflag.ge.20.and.iflag.le.29) then
                do i=1,NSTR
                    EE(i) = BETA(i) - BETB(i)
                 enddo
            else
                call LSTRES (EE,BETB,NSTR)
            endif
        else
!           total
            call LSTRES(EE,RATE(1,n),NSTR)
        endif

        if (MOD(iflag,10).eq.0) then
            ER = effee (EE,PRPMAT(2,m),effsign)
        else
            ER = EE(MOD(iflag,10))
        endif

        if (ratlog) then
            ER = ABS(ER)
            if (ER.lt.dfmin) then
                ER=dfmin
            endif
            ER = TENLOG(ER)
            comment(1:4) = 'LOG '
            l = 5
        else
            l = 1
        endif

        if (iflag.ge.10.and.iflag.le.19) then
            comment(l:l+6) = 'total '
            l = l + 6
        else if (iflag.ge.20.and.iflag.le.29) then
            comment(l:l+8) = 'viscous '
            l = l + 8
        else if (iflag.ge.30.and.iflag.le.39) then
            comment(l:l+8) = 'plastic '
            l = l + 8
        endif

        if (MOD(iflag,10).eq.0) then
            comment(l:l+10) = 'effective '
            l = l + 9
        else if (MOD(iflag,10).eq.1) then
            comment(l:l+3) = 'xx-'
            l = l + 2
        else if (MOD(iflag,10).eq.2) then
            comment(l:l+3) = 'yy-'
            l = l + 2
        else if (MOD(iflag,10).eq.3) then
            comment(l:l+3) = 'xy-'
            l = l + 2
        else if (MOD(iflag,10).eq.4) then
            comment(l:l+3) = 'xz-'
            l = l + 2
        else if (MOD(iflag,10).eq.5) then
            comment(l:l+3) = 'yz-'
            l = l + 2
        endif

        if (verbose) then
            write(luMon,2,advance='no') comment(1:l)
 2            format(1x,a)
            write(luMon,3) n,ER
 3            format('strain rate in element ',I5,' = ',1PG15.6)
        else
            write(luMon,*) ER
        endif

    enddo
else
    NUMTR = 0
    do nn=1,NUMEL
        n = nn

        m = MaterialIDOfElement(n)

        QUAD=(IEN(3,n).ne.IEN(4,n))

        if (CST) then
            if (iflag.ge.20) then
                call FORMBT (STN(1,n),BETA,BETB,PRPMAT(1,m), &
                                PRPLAS(1,m))
                BETA(3) = BETA(3)*5D-1
                BETB(3) = BETB(3)*5D-1
                if (iflag.ge.20.and.iflag.le.29) then
                    do i=1,NSTR
                        EE(i) = BETA(i) - BETB(i)
                    enddo
                else
                    call LSTRES (EE,BETB,NSTR)
                endif
            else
                call LSTRES (EE,RATE(1,n),NSTR)
            endif

            if (MOD(iflag,10).eq.0) then
                ER = effee (EE,PRPMAT(2,m),effsign)
            else
                ER = EE(MOD(iflag,10))
            endif

            if (ratlog) then
                ER = ABS(ER)
                if (ER.lt.dfmin) then
                    ER=dfmin
                endif
                ER = TENLOG(ER)
            endif

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = ER
            XPT(3,2,NUMTR) = ER
            XPT(3,3,NUMTR) = ER
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = ER
                XPT(3,2,NUMTR) = ER
                XPT(3,3,NUMTR) = ER
            endif
        else
!            determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!            number of local triangles NTR
            if (ntr.gt.0) then

                if (iflag.ge.20) then
                    call FORMBT (STN(1,n),BETA,BETB,PRPMAT(1,m), &
                                 PRPLAS(1,m))
                    BETA(3) = BETA(3)*5D-1
                    BETB(3) = BETB(3)*5D-1
                    if (iflag.ge.20.and.iflag.le.29) then
                        do i=1,NSTR
                            EE(i) = BETA(i) - BETB(i)
                        enddo
                    else
                        call LSTRES (EE,BETB,NSTR)
                    endif
                else
                    call LSTRES (EE,RATE(1,n),NSTR)
                endif

                if (MOD(iflag,10).eq.0) then
                    ER = effee (EE,PRPMAT(2,m),effsign)
                else
                    ER = EE(MOD(iflag,10))
                endif

                if (ratlog) then
                    ER = ABS(ER)
                if (ER.lt.dfmin) ER=dfmin
                    ER = TENLOG(ER)
                endif

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = ER
                enddo
!                loop triangles
                do mm=1,ntr
!                loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)

                        m = MaterialIDOfElement(n)

                        if (iflag.ge.20) then
                            call FORMBT (STN(1,n),BETA,BETB, &
                                         PRPMAT(1,m),PRPLAS(1,m))
                            BETA(3) = BETA(3)*5D-1
                            BETB(3) = BETB(3)*5D-1
                            if (iflag.ge.20.and.iflag.le.29) then
                                do j=1,NSTR
                                    EE(j) = BETA(j) - BETB(j)
                                enddo
                            else
                                call LSTRES (EE,BETB,NSTR)
                            endif
                        else
                            call LSTRES (EE,RATE(1,n),NSTR)
                        endif

                        if (MOD(iflag,10).eq.0) then
                            ER = effee (EE,PRPMAT(2,m),effsign)
                        else
                            ER = EE(MOD(iflag,10))
                        endif

                        if (ratlog) then
                            ER = ABS(ER)
                            if (ER.lt.dfmin) then
                                ER=dfmin
                            endif
                            ER = TENLOG(ER)
                        endif

                        XPT(3,i+1,NUMTR+mm) = ER
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif
!
return
end subroutine

!-------------------------------------------------------------------------------

function ShHeat (STN,RATE,POIS)

use constants, only: two

! calculates maximum shear heat production rate per unit volume

implicit none

integer, parameter :: NSTR=5
!-pass
double precision   :: STN(NSTR),RATE(NSTR)
!-locl
double precision   :: dstn(NSTR),ee(NSTR)
double precision   :: shHeat, pois

call LSTRES(ee,RATE,NSTR)
call strdev (ee,pois)
call LSTRES(dstn,STN,NSTR)
call stndev (dstn,pois)

ShHeat = dstn(1)*ee(1) + &
         dstn(2)*ee(2) + &
    TWO*(dstn(3)*ee(3) + &
         dstn(4)*ee(4) + &
         dstn(5)*ee(5))

return
end function

!-------------------------------------------------------------------------------

subroutine efstrn (STR,MAT,PRPMAT,IEN,XPT,IJSUR,verbose,luMon, &
                   CST,IFLAG)

USE MODELDEFINITION, only: NUMEL, NUMTR
use materialsmodule, only: MaterialIDOfElement, &
                           matdatactx, &
                           effsign
use monitorVars,     only: elmon, nnel, monel

!
!    IFLAG=0: effective
!    IFLAG=1: xx
!    IFLAG=2: yy
!    IFLAG=3: xy
!    IFLAG=4: xz
!    IFLAG=5: yz
!    IFLAG=6: zz
!
implicit none

integer, parameter :: NSTR=5,NEN=4,NSUR=4,NPLT=3,NPLDOF=3
!-pass
integer            :: luMon,IFLAG
logical            :: verbose,CST

integer            :: ien, ijsur, mat
double precision   :: str, prpmat, xpt
dimension          :: STR(NSTR,*),MAT(*),PRPMAT(6,*),IEN(NEN,*), &
                      IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*)
!-locl
integer, parameter :: NTRIA=2
logical            :: quad
integer            :: itria(ntria,2)
integer            :: n, nn, m, mm, ntr, i
double precision   :: pois, se

character(len=2)   :: type
dimension          :: type(6)
save type
data type/'xx','yy','xy','xz','yz','zz'/

double precision, external ::  effee


if (elmon) then
    do nn=1,nnel

        n    = monel(nn)
        m    = MaterialIDOfElement(n)
        pois = matdatactx%PRPMAT(2,m)

        if (IFLAG.eq.0) then
            SE = EFFEE (STR(1,n),pois,effsign)
        else if (IFLAG.ge.1.and.IFLAG.le.5) then
            SE = STR(IFLAG,n)
        else
            SE = 0d0
        endif

        if (verbose) then
            if (IFLAG.eq.0) then
                write(luMon,2) n,SE
 2                   format(1x,'Effective strain in',1x, &
                 'element ',I5,' = ',1PG15.6)
            else
                write(luMon,3) type(IFLAG),n,SE
 3                   format(1x,A2,'-strain in element ',I5,' = ', &
                 1PG15.6)
            endif
        else
            write(luMon,*) SE
        endif

    enddo

else
    NUMTR = 0
    do nn=1,NUMEL

        n    = nn
        quad = (IEN(3,n).ne.IEN(4,n))
        m    = MaterialIDOfElement(n)
        pois = matdatactx%PRPMAT(2,m)

        if (IFLAG.eq.0) then
            SE = EFFEE (STR(1,n),pois,effsign)
        else if (IFLAG.ge.1.and.IFLAG.le.5) then
            SE = STR(IFLAG,n)
        else
            SE = 0d0
        endif

        if (CST) then

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = SE
            XPT(3,2,NUMTR) = SE
            XPT(3,3,NUMTR) = SE
            if (quad) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = SE
                XPT(3,2,NUMTR) = SE
                XPT(3,3,NUMTR) = SE
            endif

        else
!           determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,quad)
!           number of local triangles NTR
            if (ntr.gt.0) then
                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = SE
                enddo
!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                     do i=1,2
                        n    = ITRIA(i,mm)
                        m    = MaterialIDOfElement(n)
                        pois = matdatactx%PRPMAT(2,m)

                        if (IFLAG.eq.0) then
                            SE = EFFEE (STR(1,n),pois,effsign)
                        else if (IFLAG.ge.1.and.IFLAG.le.5) then
                            SE = STR(1,n)
                        else
                            SE = 0d0
                        endif

                        XPT(3,i+1,NUMTR+mm) = SE
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine

!-------------------------------------------------------------------------------

subroutine VolumetricStrain (IEN,IJSUR,XPT,MAT,PRPMAT,STR,verbose,luMon, &
                             CST,vrate)

USE MODELDEFINITION, only: NUMEL, NUMTR
use modeldatamodule, only: modeldatactx
use materialsmodule, only: MaterialIDOfElement, &
                           matdatactx
use monitorVars,     only: elmon, nnel, monel

implicit none
!

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
!-pass

!-pass
integer            :: luMon
logical            :: verbose,CST,vrate

integer            :: ien, ijsur, mat
double precision   :: str, prpmat, xpt
dimension          :: STR(NSTR,*),MAT(*),PRPMAT(6,*),IEN(NEN,*), &
                      IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*)





!-local
integer, parameter :: NTRIA=2
logical            :: quad
integer            :: itria(ntria,2)
integer            :: n, nn, m, mm, ntr, i
double precision   :: pois, se

double precision, external :: dilate


if (elmon) then
    do nn=1,nnel
        n    = monel(nn)
        m    = MaterialIDOfElement(n)
        pois = matdatactx%PRPMAT(2,m)
        SE   = dilate (modeldatactx%STR(1,n),pois)
        if (verbose) then
            if (vrate) then
                write(luMon,1) n,SE
 1                   format(1x,'Volumetric strain rate in',1x, &
                 'element ',I5,' = ',1PG15.6)
            else
                write(luMon,2) n,SE
 2                   format(1x,'Volumetric strain in ',1x, &
                 'element ',I5,' = ',1PG15.6)
            endif
        else
            write(luMon,*) SE
        endif
    enddo
else
    NUMTR = 0
    do nn=1,NUMEL

        n    = nn
        m    = MaterialIDOfElement(n)
        pois = matdatactx%PRPMAT(2,m)
        SE   = dilate (modeldatactx%STR(1,n),pois)
        QUAD = (IEN(3,n).ne.IEN(4,n))

        if (CST) then

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = SE
            XPT(3,2,NUMTR) = SE
            XPT(3,3,NUMTR) = SE
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = SE
                XPT(3,2,NUMTR) = SE
                XPT(3,3,NUMTR) = SE
            endif

        else
!           determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!           number of local triangles NTR
            if (ntr.gt.0) then

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = SE
                enddo

                do mm=1,ntr
                    do i=1,2

                        n    = ITRIA(i,mm)
                        m    = MaterialIDOfElement(n)
                        pois = matdatactx%PRPMAT(2,m)
                        SE   = dilate (modeldatactx%STR(1,n),pois)

                        XPT(3,i+1,NUMTR+mm) = SE
                    enddo
                enddo

                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine efstrs (IEN,STN,IJSUR,XPT,MAT,PRPMAT,verbose,luMon, &
                   CST,iflag)

USE MODELDEFINITION, only: NUMEL, NUMTR
use materialsmodule, only: MaterialIDOfElement, &
                           matdatactx, &
                           effsign, strinc
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel

! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element


implicit none
!

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
integer, parameter :: NTRIA=2

!-pass

!-pass
integer            :: luMon,IFLAG
logical            :: verbose,CST

integer            :: ien, ijsur, mat, itria
double precision   :: stn, prpmat, xpt
dimension          :: STN(NSTR,*),MAT(*),PRPMAT(6,*),IEN(NEN,*), &
                      IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), ITRIA(NTRIA,2)

!-local
logical            :: quad
integer            :: n, nn, m, mm, ntr, i
double precision   :: pois, se

double precision, external :: dilate, effstn



if (iflag.eq.4.or.iflag.eq.5) then
    write(stderr,*) 'S1 and S3 contouring not implemented'
    return
endif
!
if (elmon) then
    do nn=1,nnel
        n = monel(nn)
        if (iflag.eq.1) then
            SE = STN(1,n)
        else if (iflag.eq.2) then
            SE = STN(2,n)
        else if (iflag.eq.3) then
            SE = STN(3,n)
        else
            m = MaterialIDOfElement(n)
            if (strinc) then
                pois = 5d-1
            else
                pois = matdatactx%PRPMAT(2,m)
            endif
            SE = effstn (STN(1,n),pois,effsign)
        endif

        if (verbose) then
                write(luMon,1) n,SE
 1               format(1x,'Effective stress in element ',I5,' = ', &
             1PG15.6)
        else
            write(luMon,*) SE
        endif
    enddo
else
    NUMTR = 0
    do nn=1,NUMEL

        n = nn
        QUAD=(IEN(3,n).ne.IEN(4,n))

        if (CST) then

            if (iflag.eq.1) then
                SE = STN(1,n)
            else if (iflag.eq.2) then
                SE = STN(2,n)
            else if (iflag.eq.3) then
                SE = STN(3,n)
            else
                m = MaterialIDOfElement(n)
                if (strinc) then
                    pois = 5d-1
                else
                    pois = matdatactx%PRPMAT(2,m)
                endif
                SE = effstn (STN(1,n),pois,effsign)
            endif

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = SE
            XPT(3,2,NUMTR) = SE
            XPT(3,3,NUMTR) = SE
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = SE
                XPT(3,2,NUMTR) = SE
                XPT(3,3,NUMTR) = SE
            endif
        else
!                determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                number of local triangles NTR
            if (ntr.gt.0) then

                if (iflag.eq.1) then
                    SE = STN(1,n)
                else if (iflag.eq.2) then
                    SE = STN(2,n)
                else if (iflag.eq.3) then
                    SE = STN(3,n)
                else
                    m = MaterialIDOfElement(n)
                    if (strinc) then
                        pois = 5d-1
                    else
                        pois = matdatactx%PRPMAT(2,m)
                    endif
                    SE = effstn (STN(1,n),pois,effsign)
                endif

                do mm=1,ntr
                   XPT(3,1,NUMTR+mm) = SE
                enddo

!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)

                        if (iflag.eq.1) then
                            SE = STN(1,n)
                        else if (iflag.eq.2) then
                            SE = STN(2,n)
                        else if (iflag.eq.3) then
                            SE = STN(3,n)
                        else
                            m = MaterialIDOfElement(n)
                            if (strinc) then
                                pois = 5d-1
                            else
                                pois = matdatactx%PRPMAT(2,m)
                            endif
                            SE = effstn (STN(1,n),pois, &
                                         effsign)
                        endif

                        XPT(3,i+1,NUMTR+mm) = SE
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine frheat (IEN,STN,RATE,IJSUR,XPT,MAT,PRPMAT,verbose, &
                   luMon,CST)

USE MODELDEFINITION, only: NUMEL, NUMTR
use monitorVars,     only: elmon, nnel, monel
use materialsmodule, only: lmat

implicit none
!

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
integer, parameter :: NTRIA=2

!-pass

!-pass
integer            :: luMon
logical            :: verbose,CST

integer            :: ien, ijsur, mat, itria
double precision   :: stn, prpmat, xpt, rate
dimension          :: STN(NSTR,*),MAT(*),PRPMAT(6,*),IEN(NEN,*), &
                      IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), ITRIA(NTRIA,2), &
                      RATE(NSTR,*)
!-local
logical            :: quad
integer            :: n, nn, mm, ntr, i
double precision   :: pois, heat

double precision, external :: shHeat


if (elmon) then
    do nn=1,nnel

        n = monel(nn)

        if (LMAT) then
            pois = PRPMAT(2,n)
        else
            pois = PRPMAT(2,MAT(n))
        endif

        HEAT = ShHeat (STN(1,n),RATE(1,n),pois)

        if (verbose) then
            write(luMon,1) n,HEAT
 1               format(1x,'Heat production rate per unit volume in', &
             1x,'element ',I5,' = ',1PG15.6)
        else
            write(luMon,*) HEAT
        endif

    enddo
else
    NUMTR = 0
    do nn=1,NUMEL

        n    = nn
        QUAD = (IEN(3,n).ne.IEN(4,n))

        if (LMAT) then
            pois = PRPMAT(2,n)
        else
            pois = PRPMAT(2,MAT(n))
        endif

        HEAT = ShHeat (STN(1,n),RATE(1,n),pois)

        if (CST) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = HEAT
            XPT(3,2,NUMTR) = HEAT
            XPT(3,3,NUMTR) = HEAT
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = HEAT
                XPT(3,2,NUMTR) = HEAT
                XPT(3,3,NUMTR) = HEAT
            endif
        else
!           determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!           number of local triangles NTR
            if (ntr.gt.0) then

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = HEAT
                enddo

!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)

                        if (LMAT) then
                            pois = PRPMAT(2,n)
                        else
                            pois = PRPMAT(2,MAT(n))
                        endif

                        HEAT = ShHeat (STN(1,n),RATE(1,n),pois)
                        XPT(3,i+1,NUMTR+mm) = HEAT

                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine presur (MAT,PRPMAT,IEN,STN,IJSUR,XPT,verbose,luMon, &
                   CST)

USE MODELDEFINITION, only: NUMEL, NUMTR
use materialsmodule, only: MaterialIDOfElement, &
                           matdatactx, strinc
use constants,       only: three
use monitorVars,     only: elmon, nnel, monel

implicit none
!

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
integer, parameter :: NTRIA=2
!-pass

!-pass
integer            :: luMon
logical            :: verbose,CST

integer            :: ien, ijsur, mat
double precision   :: stn, prpmat, xpt
dimension          :: STN(NSTR,*),MAT(*),PRPMAT(6,*),IEN(NEN,*), &
                      IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*)
!-local
logical            :: quad
integer            :: n, nn, m, mm, ntr, i
double precision   :: pois, se
integer            :: ITRIA(NTRIA,2)

double precision, external :: shHeat, aj1


if (elmon) then
    do nn=1,nnel
        n = monel(nn)
        m = MaterialIDOfElement(n)

        if (strinc) then
            pois = 5d-1
        else
            pois = matdatactx%PRPMAT(2,m)
        endif

        SE = AJ1(STN(1,n),pois) / three

        if (verbose) then
            write(luMon,1) n,SE
 1               format(1x,'Pressure in element ',I5,' = ', &
             1PG15.6)
        else
            write(luMon,*) SE
        endif

    enddo

else
    NUMTR = 0
    do nn=1,NUMEL

        n    = nn
        QUAD =(IEN(3,n).ne.IEN(4,n))
        m    = MaterialIDOfElement(n)

        if (strinc) then
            pois = 5d-1
        else
            pois = matdatactx%PRPMAT(2,m)
        endif

        SE = AJ1(STN(1,n),pois) / three

        if (CST) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = SE
            XPT(3,2,NUMTR) = SE
            XPT(3,3,NUMTR) = SE
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = SE
                XPT(3,2,NUMTR) = SE
                XPT(3,3,NUMTR) = SE
            endif
        else
!                determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                number of local triangles NTR
            if (ntr.gt.0) then
                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = SE
                enddo

!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2

                        n = ITRIA(i,mm)
                        m = MaterialIDOfElement(n)

                        if (strinc) then
                            pois = 5d-1
                        else
                            pois = matdatactx%PRPMAT(2,m)
                        endif

                        SE  = AJ1(STN(1,n),pois) / three
                        XPT(3,i+1,NUMTR+mm) = SE
                    enddo
                enddo

                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine

!-------------------------------------------------------------------------------

subroutine efflux (IEN,FLUX,IJSUR,XPT,verbose,luMon,NUMEL,CST)

USE ALGEBRA,         only: dot
USE MODELDEFINITION, only: NUMTR
use monitorVars,     only: elmon, nnel, monel

! Subroutine to calculate heat flux magnitudes.
! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

implicit none
!

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
!-pass

!-pass
integer            :: luMon, numel
logical            :: verbose,CST

integer            :: ien, ijsur
double precision   :: flux, xpt
dimension          :: IJSUR(NSUR,*),FLUX(NSD,*),XPT(NPLDOF,NPLT,*), &
                      IEN(NEN,*)



!-local
integer, parameter :: NTRIA=2
integer            :: ITRIA(NTRIA,2)
logical            :: quad
integer            :: n, nn, mm, ntr, i
double precision   :: ef



if (elmon) then
    do nn=1,nnel

        n  = monel(nn)
        ef = DOT (FLUX(1,n),FLUX(1,n),NSD)
        ef = SQRT(ef)

        if (verbose) then
            write(luMon,1) n,ef
 1               format(1x,'Heat flux magnitude in element ',I5, &
             ' = ',1PG15.6)
        else
            write(luMon,*) ef
        endif

    enddo
else
    NUMTR = 0
    do nn=1,NUMEL

        n    = nn
        QUAD = (IEN(3,n).ne.IEN(4,n))

        if (CST) then

            ef = DOT (FLUX(1,n),FLUX(1,n),NSD)
            ef = SQRT(ef)

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = ef
            XPT(3,2,NUMTR) = ef
            XPT(3,3,NUMTR) = ef
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = ef
                XPT(3,2,NUMTR) = ef
                XPT(3,3,NUMTR) = ef
            endif

        else
!                determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                number of local triangles NTR
            if (ntr.gt.0) then

                ef = DOT (FLUX(1,n),FLUX(1,n),NSD)
                ef = SQRT(ef)

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = ef
                enddo

!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2

                        n  = ITRIA(i,mm)
                        ef = DOT (FLUX(1,n),FLUX(1,n),NSD)
                        ef = SQRT(ef)

                        XPT(3,i+1,NUMTR+mm) = ef
                    enddo
                enddo

                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine effyld (STN,IJSUR,XPT,MAT,verbose, luMon,CST)

USE ALGEBRA,         only: dot, clear
USE MODELDEFINITION, only: NUMEL, NUMTR
use materialsmodule, only: MaterialIDOfElement, &
                           matdatactx
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel

! for implementing strain hardening, store full betb array
! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

implicit none

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
integer, parameter :: NTRIA=2

!-pass
integer            :: luMon
logical            :: verbose,CST

integer            :: ijsur, mat
double precision   :: stn, xpt
dimension          :: STN(NSTR,*),MAT(*), &
                      IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*)
!-local
logical            :: quad
integer            :: n, nn, m, mm, ntr, i
double precision   :: STNHRD
double precision   :: BETB(NSTR)
integer            :: ITRIA(NTRIA,2)
double precision   :: fNorm

double precision, external :: shHeat, aj1

!-init
call clear(BETB,NSTR,"BETB")

if (elmon) then
    do nn=1,nnel

        n      = monel(nn)
        m      = MaterialIDOfElement(n)
        STNHRD = matdatactx%PRPLAS(7,m)

        if (STNHRD.gt.1e-6) then
            write(stderr,1)
 1           format(1x,'effyld: strain hardening calculation not',1x, &
         'implemented')
        endif

        call calyld(STN(1,n),matdatactx%PRPLAS(1,m),matdatactx%PRPMAT(1,m),BETB,FNORM)

        if (verbose) then
            write(luMon,2) n,FNORM
 2               format(1x,'normalized yield stress in element ',I5, &
             ' = ',1PG15.6)
        else
            write(luMon,*) FNORM
        endif

    enddo
else
    NUMTR = 0
    do nn=1,NUMEL

        n      = nn
        m      = MaterialIDOfElement(n)
        STNHRD = matdatactx%PRPLAS(7,m)

        if (STNHRD.gt.1e-6) then
            write(stderr,1)
        endif

        QUAD=(meshdatactx%IEN(3,n).ne.meshdatactx%IEN(4,n))

        if (CST) then
            call calyld (STN(1,n),matdatactx%PRPLAS(1,m),matdatactx%PRPMAT(1,m), &
                         BETB,FNORM)

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = FNORM
            XPT(3,2,NUMTR) = FNORM
            XPT(3,3,NUMTR) = FNORM
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = FNORM
                XPT(3,2,NUMTR) = FNORM
                XPT(3,3,NUMTR) = FNORM
            endif
         else
!           determine triangle element indices
               call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!           number of local triangles NTR
            if (ntr.gt.0) then
                call calyld (STN(1,n),matdatactx%PRPLAS(1,m),matdatactx%PRPMAT(1,m), &
                             BETB,FNORM)
                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = FNORM
                enddo

!               loop triangles
                 do mm=1,ntr
!                   loop surrounding elements
                     do i=1,2
                        n = ITRIA(i,mm)
                        m = MaterialIDOfElement(n)
                        if (STNHRD.gt.1e-6) then
                            write(stderr,1)
                        endif

                        call calyld (modeldatactx%STN(1,n),matdatactx%PRPLAS(1,m), &
                                     matdatactx%PRPMAT(1,m),BETB,FNORM)

                        XPT(3,i+1,NUMTR+mm) = FNORM
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine
!-------------------------------------------------------------------------------

subroutine calyld (STN,PRPLAS,PRPMAT,BETB,FNORM)
!
! Computes FNORM
!
implicit none
!-pass

double precision   :: STN,PRPLAS,PRPMAT,BETB,FNORM
dimension          :: STN(*),PRPLAS(*),PRPMAT(*),BETB(*)

double precision   :: uniyld, pstrs

! compute scalar plastic stress PSTRS
call YIELDS (STN,PRPLAS,PRPMAT,pstrs)

! compute FNORM
call INYLD (STN,PRPLAS,PRPMAT,pstrs,BETB,FNORM,uniyld)

if (FNORM.le.0.0) then
    FNORM = 0.0
endif

return
end subroutine

!-------------------------------------------------------------------------------
subroutine ctrtemp (IEN,X,D,LMF,TFAULT, &
                    T,TX,LMTX,IJSUR,XPT,SKEW,verbose,luMon,CST)

USE MODELDEFINITION, only: NUMEL, NUMTR
use constants,       only: useDXE, one
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
!-pass
logical            :: verbose,CST
integer            :: luMon

integer            :: ijsur, ien, lmf, lmtx
double precision   :: xpt, x, d, tfault, t, tx, skew
dimension          :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*),&
                      TFAULT(NDOF,*), T(*),TX(*),LMTX(NEN,*),SKEW(2,*)


!-locl
logical            :: QUAD
double precision   :: XL(24),DL(24)
integer            :: ITRIA(NTRIA,2)
integer            :: n, nn, mm, ntr, i, ierr
double precision   :: Tc

double precision, external :: Tcentr

ierr = 0

if (elmon) then
    do nn=1,nnel

        n = monel(nn)

        call LCOORD (X,XL,IEN(1,n))
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,170)
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
        call REZONE (XL,DL,ONE)

        Tc = Tcentr (T,TX,XL,IEN(1,n),LMTX(1,n),ierr)

        if (ierr.ne.0) then
            write(stderr,1) n
 1            format(1x,'ctrtem: error in element ',I5)
            call exitp(1)
        endif

        if (verbose) then
            write(luMon,2) n,Tc
 2               format(1x,'temperature in element ',I5, &
             ' = ',1PG15.6)
        else
            write(luMon,*) Tc
        endif

    enddo
else
    NUMTR = 0
    do nn=1,NUMEL

        n = nn

        call LCOORD (X,XL,IEN(1,n))
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,171)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)
        call REZONE (XL,DL,ONE)

        Tc = Tcentr (T,TX,XL,IEN(1,n),LMTX(1,n),ierr)

        if (ierr.ne.0) then
            write(stderr,1) n
            call exitp(1)
        endif

        QUAD=(IEN(3,n).ne.IEN(4,n))

        if (CST) then

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = Tc
            XPT(3,2,NUMTR) = Tc
            XPT(3,3,NUMTR) = Tc
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = Tc
                XPT(3,2,NUMTR) = Tc
                XPT(3,3,NUMTR) = Tc
            endif

        else
!                determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                number of local triangles NTR
            if (ntr.gt.0) then

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = Tc
                enddo

!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)

                        call LCOORD (X,XL,IEN(1,n))
                        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
                        call AddFaultDisplacement (DL,n,TFAULT,NDOF, NEN, 172)
                        call ADDSNE (DL,n,NDOF,NEN,useDXE)
                        call REZONE (XL,DL,ONE)

                        Tc = Tcentr(T,TX,XL,IEN(1,n),LMTX(1,n), &
                                    ierr)

                        if (ierr.ne.0) then
                            write(stderr,1) n
                            call exitp(1)
                        endif

                        XPT(3,i+1,NUMTR+mm) = Tc
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif
return
end subroutine

!-------------------------------------------------------------------------------

subroutine plbeta (IEN,IJSUR,XPT,STR,verbose,luMon,ix,CST)

USE MODELDEFINITION, only: NUMEL, NUMTR
use constants,       only: one
use monitorVars,     only: elmon, nnel, monel

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
!-pass
logical            :: verbose,CST
integer            :: ix, luMon

integer            :: ijsur, ien
double precision   :: xpt, str
dimension          :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      IEN(NEN,*), STR(NSTR,*)





!-locl
logical            :: QUAD
integer            :: ITRIA(NTRIA,2)
double precision   :: se
integer            :: n, nn, mm, ntr, i


if (elmon) then
    do nn=1,nnel
        n = monel(nn)
        SE = STR(ix,n) + ONE
        if (verbose) then
            write(luMon,1) n,SE
 1        format(1x,'beta in',1x,'element ',I5,' = ',1PG15.6)
        endif
    enddo
else
    NUMTR = 0
    do nn=1,NUMEL

        n    = nn
        quad = (IEN(3,n).ne.IEN(4,n))
        SE   = STR(ix,n) + ONE

        if (CST) then

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = SE
            XPT(3,2,NUMTR) = SE
            XPT(3,3,NUMTR) = SE
            if (quad) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = SE
                XPT(3,2,NUMTR) = SE
                XPT(3,3,NUMTR) = SE
            endif
         else
!                determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!            number of local triangles NTR
            if (ntr.gt.0) then

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = SE
                enddo

                do mm=1,ntr
                    do i=1,2
                        n   = ITRIA(i,mm)
                        SE  = STR(ix,n) + ONE

                        XPT(3,i+1,NUMTR+mm) = SE
                    enddo
                enddo

                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine MassConserv (X,D,IEN,MAT,PRPMAT,PRPTEM,XPT,IJSUR, &
                        T,T0,TX,STN,LMTX,verbose,luMon,CST)

USE MODELDEFINITION, only: NUMTR, NUMEL
use monitorVars,     only: elmon, nnel, monel

! Routine to check accuracy of solutions, based upon element volume changes
! that may be expected from thermal expansion and isotropic stress (pressure).

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5


!-pass
double precision   :: X,D,PRPMAT,PRPTEM,XPT,T,T0, stn, tx
integer            :: IEN,MAT,IJSUR,luMon
logical            :: verbose,CST
integer            :: LMTX

dimension          :: X(NSD,*),D(NDOF,*),IEN(NEN,*),MAT(*),PRPMAT(6,*), &
                      PRPTEM(6,*),XPT(NPLDOF,NPLT,*),IJSUR(NSUR,*),T(*),T0(*), &
                      STN(NSTR,*),TX(*),LMTX(NEN,*)
!-locl
double precision   :: error
integer            :: nn,n,mm,i,j,ntr,iopt
logical            :: QUAD
integer            :: ITRIA(NTRIA,2)

double precision, external :: dMass

if (elmon) then
    do nn=1,nnel
        n     = monel(nn)
        error = dMass (n,X,D,IEN,MAT,PRPMAT,PRPTEM,T,T0,TX, &
                       STN,LMTX,IOPT)
        if (verbose) then
            write(luMon,10) n,error
10               format(1x,'mass increase in element ',I5, &
             ': ',1PG12.4,' %')
        else
            write(luMon,*) error
        endif
    enddo
else
    NUMTR = 0
    do nn=1,NUMEL
        n    = nn
        QUAD = (IEN(3,n).ne.IEN(4,n))
        if (CST) then
            error = dMass (n,X,D,IEN,MAT,PRPMAT,PRPTEM,T,T0,TX, &
                           STN,LMTX,IOPT)
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = error
            XPT(3,2,NUMTR) = error
            XPT(3,3,NUMTR) = error
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = error
                XPT(3,2,NUMTR) = error
                XPT(3,3,NUMTR) = error
            endif
        else
!            determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                   number of local triangles NTR
            if (ntr.gt.0) then
                error = dMass (n,X,D,IEN,MAT,PRPMAT,PRPTEM,T,T0, &
                               TX,STN,LMTX,IOPT)
                j = NUMTR
                do mm=1,ntr
                    j = j + 1
                    XPT(3,1,j) = error
                enddo
!               loop triangles
                j = NUMTR
                do mm=1,ntr
!                   loop surrounding elements
                    j = j + 1
                    do i=1,2
                        n     = ITRIA(i,mm)
                        error = dMass (n,X,D,IEN,MAT,PRPMAT, &
                                       PRPTEM,T,T0,TX,STN,LMTX,IOPT)
                        XPT(3,i+1,j) = error
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif
return
end subroutine
!-----------------------------------------------------------------------
double precision function ELMVOL (XL,QUAD)

! Calculates element volume
use modeltopology, only: nsd, nen

implicit none

!-pass
logical          :: QUAD
double precision :: XL(NSD,NEN)
!-locl
double precision :: a(nsd),b(nsd),area
integer          :: i

do i=1,NSD
    a(i) = XL(i,2) - XL(i,1)
    b(i) = XL(i,3) - XL(i,2)
enddo

area = 5d-1*(a(1)*b(2)-a(2)*b(1))

if (QUAD) then
    do i=1,NSD
        a(i) = XL(i,4) - XL(i,3)
        b(i) = XL(i,1) - XL(i,4)
    enddo
    area = area+5d-1*(a(1)*b(2)-a(2)*b(1))
endif

ELMVOL = area
return
end function

!-----------------------------------------------------------------------

double precision function dMass (n,X,D,IEN,MAT,PRPMAT,PRPTEM, &
                                 T,T0,TX,STN,LMTX, IOPT)

use materialsmodule, only: MaterialIDOfElement, &
                           matdatactx
use constants,       only: dfmin
use iomodule,        only: stderr

! Calculates relative mass change percentage


implicit none

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=3,NSD=2,NEN=4,NSUR=4,NSTR=5
!-pass
double precision   :: X,D,PRPMAT,PRPTEM,T,TX,T0,STN
integer            :: n,IEN,MAT,IOPT, lmtx
dimension          :: X(NSD,*),D(NDOF,*),IEN(NEN,*),MAT(*),PRPMAT(6,*), &
                      PRPTEM(6,*),T(*),T0(*),TX(*),LMTX(NEN,*),STN(NSTR,*)

!-locl
double precision   :: DENSIT,TEXP,THICK,V0,V1,TC0,TC1,M0,M1, &
                      RHO1,STNDEN,E,POIS
logical            :: QUAD
integer            :: m, ierr
double precision   :: XL(24),DL(24)

double precision, external :: ELMVOL, TCENTR

!-init
dMass  = 0d0
ierr   = 0

m      = MaterialIDOfElement(n)
E      = matdatactx%PRPMAT(1,m)
POIS   = matdatactx%PRPMAT(2,m)
DENSIT = matdatactx%PRPMAT(5,m)
THICK  = matdatactx%PRPMAT(6,m)
TEXP   = matdatactx%PRPTEM(3,m)
QUAD   = (IEN(3,n).ne.IEN(4,n))

call LCOORD (X,XL,IEN(1,n))


!   Calculate initial volume V0 and temperature TC0
V0  = ELMVOL (XL,QUAD) * THICK
TC0 = TCENTR (T0,TX,XL,IEN(1,n),LMTX(1,n),ierr)

if (ierr.ne.0) then
    return
endif

! Assume zero pressure at t=0
M0 = DENSIT * V0

if (M0.lt.dfmin) then
    write(stderr,10) n
10        format(1x,'dMass: element ',I5,' has initial mass <= 0')
    return
endif

!    Calculate current volume V1 and temperature TC1
call LDISP (DL,D,IEN(1,n),NDOF,NEN)
call REZONE (XL,DL,1d0)

V1 = ELMVOL (XL,QUAD) * THICK
TC1 = TCENTR (T,TX,XL,IEN(1,n),LMTX(1,n),ierr)

if (ierr.ne.0) then
    return
endif

RHO1 = DENSIT * (1d0 - TEXP*(TC1-TC0)) * &
       STNDEN (E,POIS,STN(1,n),IOPT)
M1   = RHO1 * V1

if (M1.lt.dfmin) then
    write(stderr,20) n
20        format(1x,'dMass: element ',I5,' has mass <= 0')
    return
endif
!    Relative mass change
dMass = ABS(M1-M0)/M0 * 1d2

return
end function

!-----------------------------------------------------------------------

double precision function STNDEN (E,POIS,STN,IOPT)

use constants, only: one, two

! Calculates density

implicit none
!-pass
double precision :: E,POIS
integer          :: IOPT
double precision :: STN
dimension        :: STN(*)
!-locl
double precision :: DILET

if (IOPT.eq.0) then
    ! Plane strain
    DILET = (one+POIS)*(one-two*POIS)/E
else
    ! Plane stress
    DILET = (one-two*POIS)/E
endif

DILET = DILET * (STN(1)+STN(2))
STNDEN = one / (one+DILET)

return
end function

!-------------------------------------------------------------------------------

subroutine TDNSTY (MAT,IEN,X,D,LMF,TFAULT,T,TX,LMTX,PRPMAT, &
                   PRPTEM,IJSUR,XPT,TOUT,verbose,luMon,IERR)

USE MODELDEFINITION, only: LGDEF, NUMEL, NUMTR
use materialsmodule, only: MaterialIDOfElement
use iomodule,        only: stderr
use constants,       only: one
use monitorVars,     only: elmon, nnel, monel
use errors,          only: exec

! computes density
!
implicit none


integer, parameter :: NEN=4,NSD=2,NDOF=3
integer, parameter :: NPLT=3,NPLDOF=3,NSUR=4,NTRIA=2
!-pass
logical            :: TOUT,verbose
integer            :: MAT,IEN,LMF,LMTX,IJSUR,luMon,IERR
double precision   :: X,D,TFAULT,T,TX,PRPMAT,PRPTEM,XPT
dimension          :: MAT(*),IEN(NEN,*),X(NSD,*),D(NDOF,*), &
                      LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
                      T(*),TX(*),LMTX(NEN,*),PRPMAT(6,*),PRPTEM(6,*), &
                      IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*)


!-locl
integer            :: nn,n,m
double precision   :: densit,texp,t0,rho
double precision   :: xl(NSD,NEN),dl(NSD,NEN)

double precision, external :: Tcentr
!
if (elmon) then
    do nn=1,nnel

        n = monel(nn)
        m = MaterialIDOfElement(n)

        densit = PRPMAT(5,m)

        if (TOUT) then
            texp   = PRPTEM(3,m)
            call LCOORD (X,xl,IEN(1,n))

            if (LGDEF.ne.0) then
                call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
                call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,173)
                call ADDSNE (dl,n,NDOF,NEN,0)
                call REZONE (xl,dl,ONE)
            endif

            t0 = Tcentr (T,TX,xl,IEN(1,n),LMTX(1,n),ierr)

            if (ierr.ne.0) then
                write(stderr,1) n
 1                format(1x,'density: error in element ',I8)
                IERR = EXEC
                return
            endif

            rho = densit*(ONE-texp*t0)
        else
            rho = densit
        endif

        if (verbose) then
            write(luMon,2) n,rho
 2            format(1x,'Density in element ',I8,' = ',1PG15.6)
        else
            write(luMon,*) rho
        endif

    enddo
else
    NUMTR = 0

    do nn=1,NUMEL
        n = nn
        m = MaterialIDOfElement(n)
        densit = PRPMAT(5,m)

        if (TOUT) then
            texp   = PRPTEM(3,m)
            call LCOORD (X,xl,IEN(1,n))

            if (LGDEF.ne.0) then
                call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
                call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,174)
                call ADDSNE (dl,n,NDOF,NEN,0)
                call REZONE (xl,dl,ONE)
            endif

            t0 = Tcentr (T,TX,xl,IEN(1,n),LMTX(1,n),ierr)

            if (ierr.ne.0) then
                write(stderr,1) n
                IERR = EXEC
                return
            endif

            rho = densit*(ONE-texp*t0)
        else
            rho = densit
        endif

        NUMTR = NUMTR + 1
        XPT(3,1,NUMTR) = rho
        XPT(3,2,NUMTR) = rho
        XPT(3,3,NUMTR) = rho
        if (IEN(3,n).ne.IEN(4,n)) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = rho
            XPT(3,2,NUMTR) = rho
            XPT(3,3,NUMTR) = rho
        endif
    enddo
endif

return
end subroutine
