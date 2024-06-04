subroutine efrate (IEN,STN,PRPMAT,PRPLAS,MAT,IJSUR, &
                   XPT,RATE,ratlog,verbose,iflag,luMon,CST)

USE ALGEBRA
USE MODELDEFINITION, only: NUMEL, NUMTR
use materialsmodule, only: MaterialIDOfElement, effsign
use constants,       only: dfmin
use modeldefinition, only: icvis
use timestepmodule,  only: nintg, maxit
use materialsmodule, only: lmat
use monitorVars,     only: elmon, nnel, monel


! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element
!
! irate=10-19 total strain rate
! irate=20-29 viscous strain rate
! irate=30-39 plastic strain rate
! irate last digit =0: effective value
!       last digit =1: Exx
!       last digit =2: Eyy
!       last digit =3: Exy
!       last digit =4: S1
!       last digit =5: S3

implicit none

integer, parameter :: NPLT=3, NPLDOF=3, NEN=4, NSUR=4, NSTR=4
!-pass
logical            :: ratlog,verbose,CST
integer            :: ien, mat, ijsur
double precision   :: stn, prpmat, prplas, rate, xpt
dimension IEN(NEN,*),STN(NSTR,*),PRPMAT(6,*),PRPLAS(9,*), &
   MAT(*),IJSUR(NSUR,*),RATE(NSTR,*),XPT(NPLDOF,NPLT,*)

integer            :: iflag, lumon
!-locl
integer, parameter :: NTRIA=2
logical            :: QUAD
integer            :: ITRIA(NTRIA,2)
double precision   :: EE(NSTR),BETA(NSTR),BETB(NSTR)
character(len=25)  :: comment
integer            :: nn, n, m, l, i, j, mm, ntr
double precision   :: ER

double precision, external :: effee

!-init
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
!            convert engineering strain rate to regular one
            BETA(4) = BETA(4)*5D-1
            BETB(4) = BETB(4)*5D-1
            if (iflag.ge.20.and.iflag.le.29) then
                do i=1,NSTR
                    EE(i) = BETA(i) - BETB(i)
                enddo
            else
                call LSTRES (EE,BETB,NSTR)
            endif
        else
!            total
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
            comment(l:l+3) = 'rr-'
            l = l + 2
        else if (MOD(iflag,10).eq.2) then
            comment(l:l+3) = 'zz-'
            l = l + 2
        else if (MOD(iflag,10).eq.3) then
            comment(l:l+3) = 'tt-'
            l = l + 2
        else if (MOD(iflag,10).eq.4) then
            comment(l:l+3) = 'rz-'
            l = l + 2
        endif

        if (verbose) then
            write(luMon,2,advance='no') comment(1:l)
 2               format(1x,a)
            write(luMon,3) n,ER
 3               format('strain rate in element ',I5,' = ',1PG15.6)
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
                BETA(4) = BETA(4)*5D-1
                BETB(4) = BETB(4)*5D-1
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
                    BETA(4) = BETA(4)*5D-1
                    BETB(4) = BETB(4)*5D-1
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
                            BETA(4) = BETA(4)*5D-1
                            BETB(4) = BETB(4)*5D-1
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
double precision function ShHeat (STN,RATE,POIS)

use constants, only : two

! calculates maximum shear heat production rate per unit volume

implicit none

integer, parameter :: NSTR=4
!-pass
double precision   :: STN(NSTR),RATE(NSTR)
!-locl
double precision   :: dstn(NSTR),ee(NSTR)
double precision   :: pois
!

call LSTRES(ee,RATE,NSTR)
call strdev (ee,pois)
call LSTRES(dstn,STN,NSTR)
call stndev (dstn,pois)

ShHeat = dstn(1)*ee(1) + dstn(2)*ee(2) + TWO*dstn(3)*ee(3)
!
return
end function
!-------------------------------------------------------------------------------
subroutine efstrn (STR,MAT,PRPMAT,IEN,XPT,IJSUR,verbose,luMon, &
                   CST,IFLAG)

USE MODELDEFINITION, only: NUMEL, NUMTR
use monitorVars,     only: elmon, nnel, monel
use materialsModule, only: effsign

!
!       IFLAG=0: effective
!       IFLAG=1: rr
!       IFLAG=2: zz
!       IFLAG=3: rz
!       IFLAG=4: rt
!       IFLAG=5: zt
!       IFLAG=6: tt
!

implicit none

integer, parameter :: NSTR=4
integer, parameter :: NEN=4
integer, parameter :: NSUR=4
integer, parameter :: NPLT=3
integer, parameter :: NPLDOF=3


!-pass
integer          :: luMon
logical          :: verbose,CST

integer          :: mat, prpmat, ien, ijsur
double precision :: str, xpt
dimension STR(NSTR,*),MAT(*),PRPMAT(6,*),IEN(NEN,*), &
   IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*)

integer          :: iflag

!-locl
integer, parameter :: NTRIA=2
logical          :: quad

integer  itria(ntria,2)
character(len=2)type
dimension type(6)
save type
data type/'rr','zz','rz','rt','zt','tt'/

integer          :: n, nn, i, mm, ntr
double precision :: SE

double precision, external :: effee

if (elmon) then
    do nn=1,nnel
        n = monel(nn)

        if (IFLAG.eq.0) then
            SE = EFFEE (STR(1,n),0d0,effsign)
        else if (IFLAG.eq.1) then
            SE = STR(1,n)
        else if (IFLAG.eq.2) then
            SE = STR(2,n)
        else if (IFLAG.eq.3) then
            SE = STR(4,n)
        else if (IFLAG.eq.4) then
            SE = 0d0
        else if (IFLAG.eq.5) then
            SE = 0d0
        else if (IFLAG.eq.6) then
            SE = STR(3,n)
        endif


        if (verbose) then
            if (IFLAG.eq.0) then
                write(luMon,2) n,SE
 2                   format(1x,'Effective strain in element ',I5, &
                 ' = ',1PG15.6)
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
        n = nn

        if (IFLAG.eq.0) then
            SE = EFFEE (STR(1,n),0d0,effsign)
        else if (IFLAG.eq.1) then
            SE = STR(1,n)
        else if (IFLAG.eq.2) then
            SE = STR(2,n)
        else if (IFLAG.eq.3) then
            SE = STR(4,n)
        else if (IFLAG.eq.4) then
            SE = 0d0
        else if (IFLAG.eq.5) then
            SE = 0d0
        else if (IFLAG.eq.6) then
            SE = STR(3,n)
        endif

        quad = (IEN(3,n).ne.IEN(4,n))

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
!            determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,quad)
!            number of local triangles NTR
            if (ntr.gt.0) then

                if (IFLAG.eq.0) then
                    SE = EFFEE (STR(1,n),0d0,effsign)
                else if (IFLAG.eq.1) then
                    SE = STR(1,n)
                else if (IFLAG.eq.2) then
                    SE = STR(2,n)
                else if (IFLAG.eq.3) then
                    SE = STR(4,n)
                else if (IFLAG.eq.4) then
                    SE = 0d0
                else if (IFLAG.eq.5) then
                    SE = 0d0
                else if (IFLAG.eq.6) then
                    SE = STR(3,n)
                endif

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = SE
                enddo

!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)
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
use materialsmodule, only: matdatactx, &
                           MaterialIDOfElement
use monitorVars,     only: elmon, nnel, monel

!
implicit none

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=2,NSD=2,NEN=4,NSUR=4,NSTR=4
!-pass
integer            :: luMon
logical            :: verbose,vrate,CST
integer            :: ien, ijsur, mat
double precision   :: xpt, prpmat, str
dimension IEN(NEN,*),IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*),MAT(*), &
   PRPMAT(6,*),STR(NSTR,*)
!-locl
integer, parameter :: NTRIA=2
integer            :: ITRIA(NTRIA,2)
double precision   :: EE(NSTR),XL(24)
logical            :: QUAD




integer            :: n, nn, i, m, mm, ntr
double precision   :: SE, pois
double precision, external :: dilate


if (elmon) then
    do nn=1,nnel

        n = monel(nn)
        SE = dilate (modeldatactx%STR(1,n),0d0)

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
        n = nn
        SE = dilate (modeldatactx%STR(1,n),0d0)
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
!            determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!            number of local triangles NTR

            if (ntr.gt.0) then

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = SE
                enddo

                do mm=1,ntr
                    do i=1,2
                        n = ITRIA(i,mm)
                        m = MaterialIDOfElement(n)
                        pois = matdatactx%PRPMAT(2,m)
                        SE = dilate (modeldatactx%STR(1,n),pois)
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
subroutine efstrs (IEN,STN,IJSUR,XPT,MAT,PRPMAT,verbose,luMon, &
                   CST,iflag)

use constants,       only: half
USE MODELDEFINITION, only: NUMEL, NUMTR
use modeltopology,   only: nsd, ndof, nstr, nen
use monitorVars,     only: elmon, nnel, monel
use materialsModule, only: effsign

! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element
! iflag=1 Srr
! iflag=2 Szz
! iflag=3 Srz
! iflag=4 S1
! iflag=5 S3
! iflag=6 Stt

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3
integer, parameter :: NPLDOF=3
integer, parameter :: NSUR=4

!-pass
logical            :: CST
integer            :: ien, ijsur, itria, mat
double precision   :: stn, xpt, prpmat
dimension IEN(NEN,*),STN(NSTR,*),IJSUR(NSUR,*), &
 XPT(NPLDOF,NPLT,*),ITRIA(NTRIA,2),MAT(*),PRPMAT(6,*)
integer            :: iflag
integer            :: luMon

!-locl
logical          :: verbose,QUAD
integer          :: n, nn, i, mm
double precision :: SE
integer          :: ntr

double precision, external :: EFFSTN

if (elmon) then
    do nn=1,nnel
        n = monel(nn)


        if (iflag.eq.1) then
            SE = STN(1,n)
        else if (iflag.eq.2) then
            SE = STN(2,n)
        else if (iflag.eq.3) then
            SE = STN(4,n)
        else if (iflag.eq.6) then
            SE = STN(3,n)
        else if (iflag.eq.4) then
            SE = 0d0
        else if (iflag.eq.5) then
            SE = 0d0
        else
!                   effective deviatoric stress
            SE = EFFSTN (STN(1,n),0d0,effsign)
        endif




        if (verbose) then
            if (iflag.eq.1) then
                write(*,1) 'rr',n,SE
 1                   format(1x,A2,'-component of stress in element ',I12,' = ',1PG15.6)
            else if (iflag.eq.2) then
                write(*,1) 'zz',n,SE
            else if (iflag.eq.3) then
                write(*,1) 'rz',n,SE
            else if (iflag.eq.6) then
                write(*,1) 'tt',n,SE
            else if (iflag.eq.4) then
                write(*,1) 'S1',n,SE
            else if (iflag.eq.5) then
                write(*,1) 'S3',n,SE
            else
                write(*,2) n,SE
 2                   format(1x,'Effective stress in element ',I12,' = ',1PG15.6)
            endif



        else
            write(*,*) SE
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
            else if (iflag.eq.6) then
                SE = STN(3,n)
            else if (iflag.eq.3) then
                SE = STN(4,n)
            else if (iflag.eq.4) then
                SE = 0d0
            else if (iflag.eq.5) then
                SE = 0d0
            else
                SE = EFFSTN (STN(1,n),0d0,effsign)
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
!                   determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                   number of local triangles NTR
            if (ntr.gt.0) then
                if (iflag.eq.1) then
                    SE = STN(1,n)
                else if (iflag.eq.2) then
                    SE = STN(2,n)
                else if (iflag.eq.6) then
                    SE = STN(3,n)
                else if (iflag.eq.3) then
                    SE = STN(4,n)
                else if (iflag.eq.4) then
                    SE = 0d0
                else if (iflag.eq.5) then
                    SE = 0d0
                else
                    SE = EFFSTN (STN(1,n),0d0,effsign)
                endif


                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = SE
                enddo
!                       loop triangles
                do mm=1,ntr
!                           loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)

                        if (iflag.eq.1) then
                            SE = STN(1,n)
                        else if (iflag.eq.2) then
                            SE = STN(2,n)
                        else if (iflag.eq.6) then
                            SE = STN(3,n)
                        else if (iflag.eq.3) then
                            SE = STN(4,n)
                        else if (iflag.eq.4) then
                            SE = 0d0
                        else if (iflag.eq.5) then
                            SE = 0d0
                        else
                            SE = EFFSTN (STN(1,n),0d0,effsign)
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
use materialsmodule, only: matdatactx, lmat
use monitorVars,     only: elmon, nnel, monel

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NDOF=2,NSD=2,NEN=4,NSUR=4,NSTR=4
!-pass
logical            :: CST, verbose
integer            :: luMon

integer            :: ien, ijsur, mat, itria
double precision   :: stn, prpmat, xpt, rate
 dimension IEN(NEN,*),STN(NSTR,*),IJSUR(NSUR,*),MAT(*), &
   PRPMAT(6,*),XPT(NPLDOF,NPLT,*),ITRIA(NTRIA,2),RATE(NSTR,*)

!-locl
logical            :: QUAD
integer            :: n, nn, mm, i, ntr
double precision   :: heat, pois

double precision, external   :: ShHeat

if (elmon) then
    do nn=1,nnel
        n = monel(nn)

        if (LMAT) then
            pois = matdatactx%PRPMAT(2,n)
        else
            pois = matdatactx%PRPMAT(2,MAT(n))
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
        n = nn
        QUAD=(IEN(3,n).ne.IEN(4,n))

        if (LMAT) then
            pois = matdatactx%PRPMAT(2,n)
        else
            pois = matdatactx%PRPMAT(2,MAT(n))
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
!                determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                number of local triangles NTR
            if (ntr.gt.0) then
                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = HEAT
                enddo
!                loop triangles
                do mm=1,ntr
!                    loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)
                        if (LMAT) then
                            pois = matdatactx%PRPMAT(2,n)
                        else
                            pois = matdatactx%PRPMAT(2,MAT(n))
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
use constants,       only: three
use monitorVars,     only: elmon, nnel, monel

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NDOF=2,NSD=2,NEN=4,NSUR=4,NSTR=4

!-pass
logical            :: CST
integer            :: luMon

integer            :: ien, ijsur, itria, mat
double precision   :: stn, xpt, prpmat
 dimension IEN(NEN,*),STN(NSTR,*),IJSUR(NSUR,*), &
   XPT(NPLDOF,NPLT,*),ITRIA(NTRIA,2),MAT(*),PRPMAT(6,*)

!-locl
logical            :: verbose,QUAD
integer            :: n, nn, mm, i, ntr

double precision   :: SE

double precision, external :: AJ1

if (elmon) then
    do nn=1,nnel
        n = monel(nn)

        SE = AJ1(STN(1,n),0d0) / three

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
        n = nn

        QUAD=(IEN(3,n).ne.IEN(4,n))
        SE = AJ1(STN(1,n),0d0) / three

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
!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)

                        SE  = AJ1(STN(1,n),0d0) / three
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
subroutine efflux (IEN,FLUX,IJSUR,XPT,verbose,luMon,NUMEL,CST)

USE ALGEBRA,         only: DOT
USE MODELDEFINITION, only: NUMTR
use meshdatamodule,  only: meshdatactx
use monitorVars,     only: elmon, nnel, monel

! Subroutine to calculate heat flux magnitudes.
! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NDOF=2,NSD=2,NEN=4,NSUR=4,NSTR=4

!-pass
integer            :: luMon, numel
logical            :: verbose,CST

integer            :: ijsur, ien
double precision   :: flux, xpt
dimension IJSUR(NSUR,*),FLUX(NSD,*),XPT(NPLDOF,NPLT,*), &
   IEN(NEN,*)


!-locl
integer            :: ITRIA(NTRIA,2)
integer            :: n, nn, mm, i, ntr
logical            :: QUAD
double precision   :: ef


if (elmon) then

    do nn=1,meshdatactx%neglobal

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

    do nn=1,meshdatactx%neglobal
        n = nn
        QUAD=(IEN(3,n).ne.IEN(4,n))

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
        else ! CST == .true.
!           determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)

!           number of local triangles NTR
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
                        n = ITRIA(i,mm)
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
subroutine effyld (IJSUR,XPT,MAT,verbose, luMon,CST)

USE ALGEBRA
USE MODELDEFINITION, only: NUMEL, NUMTR
use materialsmodule, only: MaterialIDOfElement, &
                             matdatactx
use meshdatamodule,  only:  meshdatactx
use modeldatamodule, only: modeldatactx
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel

! for implementing strain hardening, store full betb array
! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NDOF=2,NSD=2,NEN=4,NSUR=4,NSTR=4

!-pass
logical            :: CST
integer            :: luMon

integer            :: ijsur, itria, mat
double precision   :: xpt
dimension IJSUR(NSUR,*), &
   XPT(NPLDOF,NPLT,*),ITRIA(NTRIA,2),MAT(*)

!-locl
double precision   :: BETB(NSTR)
logical            :: verbose,QUAD
integer            :: n, nn, m, mm, i, ntr
double precision   :: STNHRD, fNorm

!-init
call clear(BETB,NSTR,"BETB")

if (elmon) then
    do nn=1,nnel

        n = monel(nn)
        m = MaterialIDOfElement(n)

        STNHRD = matdatactx%PRPLAS(7,m)
        if (STNHRD.gt.1e-6) then
            write(stderr,1)
 1           format(1x,'effyld: strain hardening calculation not',1x, &
         'implemented')
        endif

        call calyld(modeldatactx%STN(1,n), &
                      matdatactx%PRPLAS(1,m), &
                      matdatactx%PRPMAT(1,m), &
                               BETB,FNORM)
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
        n = nn
        m = MaterialIDOfElement(n)

        STNHRD = matdatactx%PRPLAS(7,m)
        if (STNHRD.gt.1e-6) then
            write(stderr,1)
        endif

        QUAD=(meshdatactx%IEN(3,n).ne.meshdatactx%IEN(4,n))
        if (CST) then
            call calyld(modeldatactx%STN(1,n), &
                          matdatactx%PRPLAS(1,m), &
                          matdatactx%PRPMAT(1,m), &
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
                call calyld(modeldatactx%STN(1,n), &
                              matdatactx%PRPLAS(1,m), &
                              matdatactx%PRPMAT(1,m), &
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

                        call calyld(modeldatactx%STN(1,n), &
                                      matdatactx%PRPLAS(1,m), &
                                      matdatactx%PRPMAT(1,m), &
                                           BETB,FNORM)

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
double precision :: STN,PRPLAS,PRPMAT,BETB,FNORM
dimension STN(*),PRPLAS(*),PRPMAT(*),BETB(*)
!-local
double precision :: pstrs, uniyld

! compute scalar plastic stress PSTRS
call YIELDS (STN,PRPLAS,PRPMAT,pstrs)
! compute FNORM
call INYLD (STN,PRPLAS,PRPMAT,pstrs,BETB,FNORM,uniyld)

if (FNORM.le.0.0) then
    FNORM = 0.0
endif

return
end

!-------------------------------------------------------------------------------
subroutine ctrtemp (IEN,X,D,LMF,TFAULT, &
                    T,TX,LMTX,IJSUR,XPT,SKEW,verbose,luMon,CST)

USE MODELDEFINITION, only: NUMEL, NUMTR
use iomodule,        only: stderr
use constants,       only: one
use monitorVars,     only: elmon, nnel, monel

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NDOF=2,NSD=2,NEN=4,NSUR=4,NSTR=4

!-pass
logical            :: verbose,CST
integer            :: luMon

integer            :: ijsur, ien, lmf, lmtx
double precision   :: xpt, x, d, tfault, t, tx, skew
dimension IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
   IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   T(*),TX(*),LMTX(NEN,*),SKEW(2,*)


!-locl
logical QUAD
integer            :: ITRIA(NTRIA,2)
double precision   :: XL(24), DL(24)
integer            :: n, nn, m, mm, ntr, i, ierr
double precision   :: Tc

double precision, external :: Tcentr

ierr = 0

if (elmon) then
    do nn=1,nnel
        n = monel(nn)

        call LCOORD (X,XL,IEN(1,n))
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,118)
        call ADDSNE (DL,n,NDOF,NEN,0)
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
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,119)
        call ADDSNE (DL,n,NDOF,NEN,0)
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
!                loop triangles
                do mm=1,ntr
!                    loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)
                        call LCOORD (X,XL,IEN(1,n))
                        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
                        call AddFaultDisplacement (DL,n,TFAULT,NDOF, &
                                                   NEN,120)
                        call ADDSNE (DL,n,NDOF,NEN,0)
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
integer, parameter :: NPLT=3,NPLDOF=3,NEN=4,NSUR=4,NSTR=4

!-pass
integer            :: luMon,ix
logical            :: verbose,CST

integer            :: ijsur, ien
double precision   :: xpt, str
dimension IEN(NEN,*),IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
   STR(NSTR,*)



!-locl
integer            :: ITRIA(NTRIA,2)
logical            :: quad
integer            :: n, nn, m, mm, ntr, i, ierr
double precision   :: SE

!
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

        n = nn
        quad = (IEN(3,n).ne.IEN(4,n))
        SE = STR(ix,n) + ONE

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
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!           number of local triangles NTR
            if (ntr.gt.0) then

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = SE
                enddo

                do mm=1,ntr
                    do i=1,2
                        n = ITRIA(i,mm)
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
use modeltopology,   only: ndof, nsd, nen, nsur, nstr
use monitorVars,     only: elmon, nnel, monel

!    USE MATERIALSMODULE, only: LMAT

! Routine to check accuracy of solutions, based upon element volume changes
! that may be expected from thermal expansion and isotropic stress (pressure).

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3
integer, parameter :: NPLDOF=3

!-pass
double precision   :: X,D,PRPMAT,PRPTEM,XPT,T,T0, stn, tx
integer            :: IEN,MAT,IJSUR,luMon
logical            :: verbose,CST
integer :: LMTX

dimension X(NSD,*),D(NDOF,*),IEN(NEN,*),MAT(*),PRPMAT(6,*), &
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
        n = monel(nn)
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
        n = nn
        QUAD=(IEN(3,n).ne.IEN(4,n))
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
!           determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!           number of local triangles NTR
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
                        n = ITRIA(i,mm)
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
end
!-----------------------------------------------------------------------
function ELMVOL (XL,QUAD)

! Calculates element volume

implicit none

integer, parameter :: NSD=2,NEN=4
!-pass
double precision   :: ELMVOL
logical            :: QUAD
double precision   :: XL(NSD,NEN)
!-locl
double precision   :: area
double precision   :: a(NSD),b(NSD)
integer            :: i

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
end

!-----------------------------------------------------------------------

double precision function dMass (n,X,D,IEN,MAT,PRPMAT,PRPTEM,T,T0,TX, &
                                 STN,LMTX, IOPT)

USE MODELDEFINITION, only: NUMEL, NUMTR
use materialsmodule, only: matdatactx, &
                           MaterialIDOfElement
use constants,       only: dfmin
use materialsmodule, only: lmat
use timestepmodule,  only: nintg, maxit
use iomodule,        only: stderr

! Calculates relative mass change percentage

implicit none

integer, parameter :: NPLT=3,NPLDOF=3,NDOF=2,NSD=2,NEN=4,NSUR=4,NSTR=4
!-pass
double precision   :: X,D,PRPMAT,PRPTEM,T,TX,T0,STN
integer            :: n,IEN,MAT,IOPT, lmtx
dimension          :: X(NSD,*),D(NDOF,*),IEN(NEN,*),MAT(*),PRPMAT(6,*), &
                      PRPTEM(6,*),T(*),T0(*),TX(*),LMTX(NEN,*),STN(NSTR,*)

!-locl
double precision   :: DENSIT,TEXP,THICK,V0,V1,TC0,TC1,M0,M1, &
                      TCENTR,ELMVOL,RHO1,STNDEN,E,POIS
logical            :: QUAD
integer            :: m, ierr
double precision   :: XL(24),DL(24)

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

!    Calculate initial volume V0 and temperature TC0
V0 = ELMVOL (XL,QUAD) * THICK
TC0 = TCENTR (T0,TX,XL,IEN(1,n),LMTX(1,n),ierr)

if (ierr.ne.0) then
    return
endif

!    Assume zero pressure at t=0
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
M1 = RHO1 * V1

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

integer, parameter :: NEN=4,NSD=2,NDOF=2
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
integer            :: nn,n,m,ntr
double precision   :: densit,texp,t0,rho
integer            :: ITRIA(NTRIA,2)
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
                call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,121)
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
                call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,122)
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
