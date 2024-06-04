 subroutine criteria (IEN,IJSUR,XPT,verbose, &
 luMon,CST,xn,X,D,ISIDE,LMF,DXE,SKEW, &
 TFAULT,side)
! routine that use a criteria to find the best triangles

USE MODELDEFINITION, only: numel
USE MODELTOPOLOGY,   only: NSD, NDOF
use constants,       only: useDXE, one
use modeldefinition, only: lgdef
use modeltopology,   only: nen, nstr
use monitorVars,     only: elmon, nnel, monel

implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3,NPLDOF=3,NSUR=4
!-pass
logical CST

integer :: ien, ijsur, itria, iside, lmf
double precision :: d, xpt, x, skew
dimension IEN(NEN,*),IJSUR(NSUR,*),D(NDOF,*), &
  XPT(NPLDOF,NPLT,*), &
  ITRIA(NTRIA,2),X(NSD,*),ISIDE(*),LMF(NDOF,NEN,*),SKEW(2,*)


double precision ::  XL(NSD,NEN),xn(NSD)

double precision :: dxe, tfault
dimension :: DXE(NDOF,NEN,*),TFAULT(NDOF,*)

integer   :: lien(2),side(4)
integer   :: lumon


!-locl
    logical verbose, QUAD
double precision :: DL(NDOF,NEN)

integer :: i, mm, n, nn, numtr, ntr
double precision :: ar

double precision, external :: Shcriteria

if (elmon) then
    do nn=1,nnel
        n = monel(nn)
          QUAD= IEN(3,n).ne.IEN(4,n)
           call LCOORD (X,XL,IEN(1,n))
        if (LGDEF.ne.0) then
            call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,132)
               call ADDSNE (DL,n,NDOF,NEN,useDXE)
               call REZONE (XL,DL,1d0)
        endif
           AR = Shcriteria (QUAD,XL,xn,side,ISIDE)
        if (verbose) then
            write(luMon,1) n,AR
 1               format(1x,'AR', &
             1x,'element ',I5,' = ',1PG15.6)
        else
            write(luMon,*) AR
        endif
   enddo
else
    NUMTR = 0
    do nn=1,NUMEL
        n = nn
        QUAD=(IEN(3,n).ne.IEN(4,n))
           call LCOORD (X,XL,IEN(1,n))
           if (LGDEF.ne.0) then
               call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
               call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,133)
               call ADDSNE (DL,n,NDOF,NEN,useDXE)
            call REZONE (XL,DL,ONE)
          endif
         AR = Shcriteria (QUAD,XL,xn,side,ISIDE)
        if (CST) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = AR
            XPT(3,2,NUMTR) = AR
            XPT(3,3,NUMTR) = AR
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = AR
                XPT(3,2,NUMTR) = AR
                XPT(3,3,NUMTR) = AR
            endif
        else
!                   determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                   number of local triangles NTR
            if (ntr.gt.0) then
                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = AR
                enddo

!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2
                        n = ITRIA(i,mm)
                        QUAD=IEN(3,n).ne.IEN(4,n)
                        call LCOORD (X,XL,IEN(1,n))
                        if (LGDEF.ne.0) then
                            call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
                            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,134)
                            call ADDSNE (DL,n,NDOF,NEN,useDXE)
                            call REZONE (XL,DL,ONE)
                        endif
                        AR = Shcriteria (QUAD,XL,xn,side,ISIDE)
                        XPT(3,i+1,NUMTR+mm) = AR
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif
!
return
end
!-------------------------------------------------------------
double precision function Shcriteria (QUAD,XL,xn,side,ISIDE)

USE MODELDEFINITION
use constants,       only: two
use modeltopology,   only: nstr, nen, nsd

implicit none
!-pass
double precision :: STN(NSTR),RATE(NSTR)
logical          :: QUAD
integer          :: side, iside

!-locl
integer          :: lien
double precision :: xn, xl
dimension lien(2),xn(NSD)
dimension  XL(NSD,NEN),side(4),ISIDE(*)

double precision :: AMINAREA,AMAXAREA
!-locl
integer :: k, kmax
double precision :: area

if (QUAD) then
    kmax=4
else
    kmax=3
endif

do k=1,kmax
    call SIDENP(k,QUAD,lien)
    call FCGEOM(XL,lien,area,xn)
!      write(stderr,37)
    side(k)= area
enddo

if (QUAD) then
    AMAXAREA=max(side(1),side(2),side(3),side(4))
    AMINAREA=min(side(1),side(2),side(3),side(4))
else
    AMAXAREA=max(side(1),side(2),side(3))
    AMINAREA=min(side(1),side(2),side(3))
endif

Shcriteria = AMINAREA/AMAXAREA
!
return
end
!------------------------------------------------------------------------

#ifdef EXPERIMENTAL_NORATESTORAGE


 subroutine ratear (IEN,X,D,STN,PRPMAT,PRPLAS,MAT,LMF, &
   TFAULT,DXE,SKEW,Rscale,Rsb,hatangle, &
   hatsize,NEN,NSD,NDOF,NSTR,DEFORM,deffac,lu,ierr,ratlog,ibase, &
   istep,iarcont,verbose,irate,luMon,ratdev)

USE ALGEBRA
USE MODELDEFINITION
USE CONSTANTS,       only: zero, deg2rad, dfmin, dfmax
use materialsmodule, only: MaterialIDOfElement, &
                           lmat
use timestepmodule,  only: nintg, maxit
use debugmodule,     only: iecho
use modeldatamodule, only: modeldatactx
use postProcessData, only: GMT,connout

!
! routine for converting strain rates to arrows in principal stress directions
!
! irate=10-19 total strain rate
! irate=20-29 viscous strain rate
! irate=30-39 plastic strain rate
!

implicit none

!-pass
 integer NEN,NSD,NDOF,NSTR,lu,ierr,ibase,istep,iarcont, &
   irate,luMon

integer          :: ien, mat, lmf
double precision :: x,d, stn, prpmat, prplas, tfault, dxe, skew

dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),STN(NSTR,*), &
 PRPMAT(6,*),MAT(*),PRPLAS(9,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),DXE(NDOF,NEN,*),SKEW(2,*)

logical          :: DEFORM,ratlog,verbose,ratdev

double precision :: rsb, rscale, deffac, pois
double precision :: smallevecx, smallevecy
double precision :: v1, v2





!-locl
integer, parameter :: MAXRRW=25
 real D2REAL

double precision :: strainRate(nstr)

double precision :: ee, centr, arrow, eval, evec, eemax, eemin, beta, betb

 dimension ee(6),CENTR(3),arrow(3,MAXRRW),EVAL(3), &
 EVEC(9),EEmax(6),EEmin(6),BETA(6),BETB(6)


integer           :: n,m,i,j,nn
integer           :: itype

double precision  :: hatangle, hatsize ! of arrows
double precision  :: def, deg

integer           :: narrw ! used as loop index, hence must be integer

integer, external :: ipltrw

!-init


data arrow(1,1), narrw/0d0, 0/

ierr = OK
def  = ZERO

if (DEFORM) then
    def = deffac
endif

if (IDIG(ICVIS,4,1).ne.0) then
    call CLEAR (BETB,NSTR,"BETB")
endif


itype = INT(irate/10)-1


write(*,*) "entered ratear"

if (elmon) then
    do nn=1,nnel
        n = monel(nn)
        m = MaterialIDOfElement(n)

        call StrainrateOfElem(X,&
                              D,&
                              modeldatactx%DELD,&
                              strainRate,IEN,LMF,TFAULT,&
                              modeldatactx%DFAULT,modeldatactx%LMX, &
                              modeldatactx%DX,SKEW,modeldatactx%DXE, n)

        if (itype.ne.0) then
            call FORMBT(STN(1,n),BETA,BETB,PRPMAT(1,m), &
                       PRPLAS(1,m))
            BETA(3) = BETA(3)*5D-1
            BETB(3) = BETB(3)*5D-1
            if (itype.eq.1) then
                do i=1,NSTR
                    ee(i) = BETA(i)-BETB(i)
                enddo
            else
                call LSTRES(ee,BETB,NSTR)
            endif
        else
            if (ratdev) then
                call LSTRES(ee,strainRate(1:nstr),NSTR)
                pois = PRPMAT(2,m)
                call strdev (ee,pois)
            else
                call LSTRES(EE,strainRate(1:nstr),NSTR)
            endif
        endif

        if (ratlog) then
            do i=1,NSTR
                if (ABS(EE(i)).lt.dfmin) then
                    EE(i) = dfmin
                endif
            enddo
            if (verbose) then
                if (itype.eq.0) then
                    write(luMon,2) n, &
             (TENLOG(ABS(EE(i))),i=1,NSTR)
 2                       format(1x,'LOG total rates in element ',I5, &
                     ': ',6(1PG12.4))
                else if (itype.eq.1) then
                    write(luMon,3) n, &
             (TENLOG(ABS(EE(i))),i=1,NSTR)
 3                       format(1x,'LOG visc. rates in element ',I5, &
                     ': ',6(1PG12.4))
                   else if (itype.eq.2) then
                    write(luMon,4) n, &
             (TENLOG(ABS(EE(i))),i=1,NSTR)
 4                       format(1x,'LOG plast. rates in element ',I5, &
                     ': ',6(1PG12.4))
                endif
            else
                write(luMon,*) (TENLOG(ABS(EE(i))),i=1,NSTR)
            endif
        else
            if (verbose) then
                if (itype.eq.0) then
                    write(luMon,6) n,(EE(i),i=1,NSTR)
 6                       format(1x,'total rates in element ',I5, &
                     ': ',6(1PG12.4))
                else if (itype.eq.1) then
                    write(luMon,7) n,(ee(i),i=1,NSTR)
 7                       format(1x,'viscous rates in element ',I5, &
                     ': ',6(1PG12.4))
                else if (itype.eq.2) then
                    write(luMon,8) n,(ee(i),i=1,NSTR)
 8                       format(1x,'plastic rates in element ',I5, &
                     ': ',6(1PG12.4))
                endif
            else
                write(luMon,*) (ee(i),i=1,NSTR)
            endif
        endif
    enddo
else
    ! we do not monitor elements, just output all of them.

!    hatangle = hatangle
    hatsize  = Rsb  *  Rscale  * hatsize
!   [km]   = [/s] x [km.s] *  factor

    do i=1,NSTR
        EEmin(i) =  dfmax
        EEmax(i) = -dfmax
    enddo

    do nn=ibase,NUMEL,istep
        n = nn
        m = MaterialIDOfElement(n)

        call StrainrateOfElem(X,&
                              D,&
                              modeldatactx%DELD,&
                              strainRate,IEN,LMF,TFAULT,&
                              modeldatactx%DFAULT,modeldatactx%LMX, &
                              modeldatactx%DX,SKEW,modeldatactx%DXE, n)

        write(*,*) "new ratear: ", n, strainRate


!             calculate element center coordinates
        call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
                 SKEW,def,ierr)
        if (ierr.ne.0) then
            write(stderr,9) n
 9            format(1x,'ratear: error in element ',I5)
            call stoper()
        endif

        write(*,*) "centr: ", centr
        write(*,*) "itype: ", itype

        if (itype.ne.0) then
            call FORMBT(STN(1,n),BETA,BETB,PRPMAT(1,m), PRPLAS(1,m))
            BETA(3) = BETA(3)*5D-1
            BETB(3) = BETB(3)*5D-1
            if (itype.eq.1) then
                do i=1,NSTR
                    ee(i) = BETA(i)-BETB(i)
                enddo
            else
                call LSTRES(ee,BETB,NSTR)
            endif
    
        else
            write(*,*) "ratdev", ratdev
            if (ratdev) then
                call LSTRES(ee,strainRate(1:nstr),NSTR)
                pois = PRPMAT(2,m)
                call strdev (ee,pois)
            else
                call LSTRES(ee,strainRate(1:nstr),NSTR)
            endif
        endif

        write(*,*) "ee:", ee

        do i=1,NSTR
            if (EEmin(i).gt.ee(i)) then
                EEmin(i)=ee(i)
            endif
            if (EEmax(i).lt.ee(i)) then
                EEmax(i)=ee(i)
            endif
        enddo

        call Eigen (ee,EVAL,EVEC)
!            calculate arrow coordinates
        ierr = ipltrw (EVAL,EVEC,CENTR,Rscale,hatsize,hatangle, &
                 arrow,narrw,MAXRRW,iarcont,.true.)
        if (ierr.ne.OK) then
            write(stderr,10) n
10            format(1x,'rate: error in element ',I5)
            goto 50
        endif

        write(*,*) "EVAL", EVAL
        write(*,*) "EVEC", EVEC

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
! TODO: deg/deg2rad seems odd. Is this correct?
            deg = 90.0 - deg/deg2rad
            v1 = v1*Rscale
            v2 = v2*Rscale

            write(*,*) "new ratear writes: ", CENTR(1),CENTR(2),v1,v2,deg


            if (ascout) then
                write(lu,*,err=100) CENTR(1),CENTR(2),v1,v2,deg
            else
                write(lu,err=100) &
                   D2REAL(CENTR(1)), &
                   D2REAL(CENTR(2)), &
                   D2REAL(v1), &
                   D2REAL(v2), &
                   D2REAL(deg)
            endif
        else
            stop "non-GMT plotting not implemented"
        endif

    enddo ! loop over the elements

    if (iecho.eq.1 .or. iecho.ge.4) then
        call wmnmx(EEmin,EEmax,NSTR,2)
    endif

endif

50    return
!
100    call prerr('rate')
 call exitp(1)
end

#else

 subroutine ratear (IEN,X,D,STN,PRPMAT,PRPLAS,MAT,RATE,LMF, &
   TFAULT,DXE,SKEW,Rscale,Rsb,hatangle, &
   hatsize,NEN,NSD,NDOF,NSTR,DEFORM,deffac,lu,ierr,ratlog,ibase, &
   istep,iarcont,verbose,irate,luMon,ratdev)

USE ALGEBRA
USE MODELDEFINITION
USE CONSTANTS,       only: zero, deg2rad, dfmin, dfmax
use materialsmodule, only: MaterialIDOfElement, &
                           lmat
use timestepmodule,  only: nintg, maxit
use debugmodule,     only: iecho
use iomodule,        only: ascout, stderr
use postProcessData, only: GMT,connout
use monitorVars,     only: elmon, nnel, monel
use errors,          only: ok
!
! routine for converting strain rates to arrows in principal stress directions
!
! irate=10-19 total strain rate
! irate=20-29 viscous strain rate
! irate=30-39 plastic strain rate
!

implicit none

!-pass
 integer NEN,NSD,NDOF,NSTR,lu,ierr,ibase,istep,iarcont, &
   irate,luMon

integer          :: ien, mat, lmf
double precision :: x,d, stn, prpmat, rate, prplas, tfault, dxe, skew

dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),STN(NSTR,*), &
 PRPMAT(6,*),MAT(*),RATE(NSTR,*),PRPLAS(9,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),DXE(NDOF,NEN,*),SKEW(2,*)

logical          :: DEFORM,ratlog,verbose,ratdev

double precision :: rsb, rscale, deffac, pois
double precision :: smallevecx, smallevecy
double precision :: v1, v2





!-locl
integer, parameter :: MAXRRW=25
 real D2REAL

double precision :: ee, centr, arrow, eval, evec, eemax, eemin, beta, betb

 dimension ee(6),CENTR(3),arrow(3,MAXRRW),EVAL(3), &
 EVEC(9),EEmax(6),EEmin(6),BETA(6),BETB(6)


integer           :: n,m,i,j,nn
integer           :: itype

double precision  :: hatangle, hatsize ! of arrows
double precision  :: def, deg

integer           :: narrw ! used as loop index, hence must be integer

integer, external :: ipltrw

!-init


data arrow(1,1), narrw/0d0, 0/

ierr = OK
def  = ZERO

if (DEFORM) then
    def = deffac
endif

if (IDIG(ICVIS,4,1).ne.0) then
    call CLEAR (BETB,NSTR,"BETB")
endif


itype = INT(irate/10)-1
!

if (elmon) then
    do nn=1,nnel
       n = monel(nn)
        m = MaterialIDOfElement(n)
        if (itype.ne.0) then
            call FORMBT(STN(1,n),BETA,BETB,PRPMAT(1,m), &
                       PRPLAS(1,m))
            BETA(3) = BETA(3)*5D-1
            BETB(3) = BETB(3)*5D-1
            if (itype.eq.1) then
                do i=1,NSTR
                    ee(i) = BETA(i)-BETB(i)
                enddo
            else
                call LSTRES(ee,BETB,NSTR)
            endif
        else
            if (ratdev) then
                call LSTRES(ee,RATE(1,n),NSTR)
                pois = PRPMAT(2,m)
                call strdev (ee,pois)
            else
                call LSTRES(EE,RATE(1,n),NSTR)
            endif
        endif

        if (ratlog) then
            do i=1,NSTR
                if (ABS(EE(i)).lt.dfmin) then
                    EE(i) = dfmin
                endif
            enddo
            if (verbose) then
                if (itype.eq.0) then
                    write(luMon,2) n, &
             (TENLOG(ABS(EE(i))),i=1,NSTR)
 2                       format(1x,'LOG total rates in element ',I5, &
                     ': ',6(1PG12.4))
                else if (itype.eq.1) then
                    write(luMon,3) n, &
             (TENLOG(ABS(EE(i))),i=1,NSTR)
 3                       format(1x,'LOG visc. rates in element ',I5, &
                     ': ',6(1PG12.4))
                   else if (itype.eq.2) then
                    write(luMon,4) n, &
             (TENLOG(ABS(EE(i))),i=1,NSTR)
 4                       format(1x,'LOG plast. rates in element ',I5, &
                     ': ',6(1PG12.4))
                endif
            else
                write(luMon,*) (TENLOG(ABS(EE(i))),i=1,NSTR)
            endif
        else
            if (verbose) then
                if (itype.eq.0) then
                    write(luMon,6) n,(EE(i),i=1,NSTR)
 6                       format(1x,'total rates in element ',I5, &
                     ': ',6(1PG12.4))
                else if (itype.eq.1) then
                    write(luMon,7) n,(ee(i),i=1,NSTR)
 7                       format(1x,'viscous rates in element ',I5, &
                     ': ',6(1PG12.4))
                else if (itype.eq.2) then
                    write(luMon,8) n,(ee(i),i=1,NSTR)
 8                       format(1x,'plastic rates in element ',I5, &
                     ': ',6(1PG12.4))
                endif
            else
                write(luMon,*) (ee(i),i=1,NSTR)
            endif
        endif
    enddo
else
    ! we do not monitor elements, just output all of them.

!    hatangle = hatangle
    hatsize  = Rsb  *  Rscale  * hatsize
!   [km]   = [/s] x [km.s] *  factor

    do i=1,NSTR
        EEmin(i) =  dfmax
        EEmax(i) = -dfmax
    enddo

    do nn=ibase,NUMEL,istep
        n = nn
        m = MaterialIDOfElement(n)


!        write(*,*) "old ratear: ", n, RATE(1:nstr,n)


!             calculate element center coordinates
        call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
                 SKEW,def,ierr)
        if (ierr.ne.0) then
            write(stderr,9) n
 9            format(1x,'ratear: error in element ',I5)
            call stoper()
        endif

!        write(*,*) "centr: ", centr
!        write(*,*) "itype: ", itype


        if (itype.ne.0) then
            call FORMBT(STN(1,n),BETA,BETB,PRPMAT(1,m), &
           PRPLAS(1,m))
            BETA(3) = BETA(3)*5D-1
            BETB(3) = BETB(3)*5D-1
            if (itype.eq.1) then
                do i=1,NSTR
                    ee(i) = BETA(i)-BETB(i)
                enddo
            else
                call LSTRES(ee,BETB,NSTR)
            endif

        else
!            write(*,*) "ratdev", ratdev
            if (ratdev) then
                call LSTRES(ee,RATE(1,n),NSTR)
                pois = PRPMAT(2,m)
                call strdev (ee,pois)
            else
                call LSTRES(ee,RATE(1,n),NSTR)
            endif
        endif

!        write(*,*) "ee:", ee

        do i=1,NSTR
            if (EEmin(i).gt.ee(i)) then
                EEmin(i)=ee(i)
            endif
            if (EEmax(i).lt.ee(i)) then
                EEmax(i)=ee(i)
            endif
        enddo

        call Eigen (ee,EVAL,EVEC)
!            calculate arrow coordinates
        ierr = ipltrw (EVAL,EVEC,CENTR,Rscale,hatsize,hatangle, &
                 arrow,narrw,MAXRRW,iarcont,.true.)
        if (ierr.ne.OK) then
            write(stderr,10) n
10            format(1x,'rate: error in element ',I5)
            goto 50
        endif

!        write(*,*) "EVAL", EVAL 
!        write(*,*) "EVEC", EVEC


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
! TODO: deg/deg2rad seems odd. Is this correct?
            deg = 90.0 - deg/deg2rad
            v1 = v1*Rscale
            v2 = v2*Rscale

!            write(*,*) "old ratear writes: ", CENTR(1),CENTR(2),v1,v2,deg


            if (ascout) then
                write(lu,*,err=100) CENTR(1),CENTR(2),v1,v2,deg
            else
                write(lu,err=100) &
                   D2REAL(CENTR(1)), &
                   D2REAL(CENTR(2)), &
                   D2REAL(v1), &
                   D2REAL(v2), &
                   D2REAL(deg)
            endif
        else
            stop "non-GMT plotting not implemented"
        endif

    enddo ! loop over the elements

    if (iecho.eq.1 .or. iecho.ge.4) then
        call wmnmx(EEmin,EEmax,NSTR,2)
    endif

endif

50    return
!
100    call prerr('rate')
 call exitp(1)
end



#endif



!-------------------------------------------------------------------------------
subroutine sclbar (x,y,sb,scale,headangle,headsize,hor,iarcont, &
                   lu,ierr,both)

USE ALGEBRA
use iomodule,        only: ascout, stderr
use errors,          only: ok

! plots a scale bar at (x,y)

implicit none
!-pass
integer lu,ierr,iarcont
logical hor,both
double precision x, y, headangle, headsize, sb, scale


!-locl
integer, parameter :: MAXRRW=25
double precision :: CENTR(3),arrow(3,MAXRRW),EVAL(3),EVEC(3,3)
double precision :: hatangle, hatsize

 real D2REAL

!-local
!integer :: ierr
integer :: i,j
integer :: narrw
integer, external :: ipltrw

!-init
 data arrow(1,1),narrw/0d0,0/

 ierr     = OK
 hatangle = headangle
 hatsize  = sb  *  scale  * headsize
 call clear(EVEC,9,"EVEC")
 call clear(EVAL,3,"EVAL")

!    calculate arrow coordinates
 CENTR(1) = x
 CENTR(2) = y
 CENTR(3) = 0d0
 EVEC(1,1) = 1d0
 EVEC(2,2) = 1d0
 EVEC(3,3) = 1d0
 if (hor) then
EVAL(1) = sb
if (both) then
EVAL(2) = sb*0.1d0
else
EVAL(2) = 0.d0
endif
 else
if (both) then
EVAL(1) = sb*0.1d0
else
EVAL(1) = 0.d0
endif
EVAL(2) = sb
 endif
 ierr = ipltrw (EVAL,EVEC,CENTR,scale,hatsize,hatangle, &
 arrow,narrw,MAXRRW,iarcont,both)
 if (ierr.ne.OK) then
write(stderr,1)
 1        format(1x,'scalebar error')
 else
if (narrw.gt.0) then
if (ascout) then
    write(lu,2,err=100) 4e0,REAL(narrw)
    write(lu,2,err=100) &
         ((arrow(i,j),i=1,3),j=1,narrw)
 2            format(1x,3(1PG12.6))
else
    write(lu,err=100) REAL(4.0),REAL(narrw)
    write(lu,err=100) &
         ((D2REAL(arrow(i,j)),i=1,3),j=1,narrw)
    endif
endif
 endif
 return
100    call prerr('sclbar')
 call exitp(1)
 end
!-------------------------------------------------------------------------------
double precision function etaLin (effstr,emhu,anpwr)

USE ALGEBRA
use constants, only: dfmin,dlogmin
! calculates linear viscosity equivalent

implicit none
!-pass
double precision emhu,anpwr,effstr
!-locl
double precision :: eta, xpow
double precision :: x, y


eta    = 1D50
xpow   = anpwr - 1d0

if (xpow.gt.1E-6) then
    if (effstr.gt.dfmin) then
        x = TENLOG(effstr)
        x = xpow * x
        if (x.gt.dlogmin) then
            y = TENLOG(emhu)
            y = anpwr * y
            eta = y - x
            eta = 1D1 ** eta
        endif
    endif
else
    eta = emhu
endif

etaLin = eta

return
end
!-------------------------------------------------------------------------------
subroutine grid (IEN,X,D,XPT,LMF,TFAULT,DXE, &
                 SKEW,NUMBRS,DEFORM,deffac,lu)

USE ALGEBRA
USE MODELDEFINITION
use constants,       only: useDXE, zero
use modeltopology,   only: NEN,NSD,NDOF
use iomodule,        only: ascout, stderr
use postProcessData, only: GMT,connout
use monitorVars,     only: elmon, nnel, monel, npmon, nnp, monnp


! subroutine to write element coordinates to output file
! NUMBRS=.true.  : P: 'pt'  commands, GMT: psxy with elm numbers
! NUMBRS=.false. : P: datafile, GMT: psxy multisegment
! DEFORM=.true.  : deformed mesh
! DEFORM=.false. : original coordinates

 implicit none
!-pass
 integer lu
double precision :: deffac

integer          :: ien, lmf
double precision :: x, d, centr, xpt, tfault, dxe, skew
 dimension X(NSD,*),D(NDOF,*),IEN(NEN,*),CENTR(3),XPT(3,*), &
   LMF(NDOF,NEN,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*), &
   SKEW(2,*)

 logical NUMBRS,DEFORM



!-locl
double precision :: XL(24),DL(24)
logical QUAD
integer i, n, nn, ntot, nxpt, ierr
double precision :: def

!-init
 if (DEFORM) then
def = deffac
 else
def = ZERO
call CLEAR(DL,24,"DL")
 endif

 ntot = NUMEL

 if (elmon) then
ntot=nnel
 endif

 if (NUMBRS) then

if (.not.GMT) then
     write(lu,*) '** ',ntot,' elements'
endif

do nn=1,ntot
    n = nn
    if (elmon) then
        n=monel(nn)
    endif

    call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
                SKEW,def,ierr)

    if (ierr.ne.0) then
        call exitp(1)
    endif

    if (GMT) then
        write(lu,*) CENTR(1),CENTR(2),n
    else
        stop "non-GMT plotting not implemented"
    endif

enddo
 else
NXPT = 0

do nn=1,ntot

    n = nn

    if (elmon) then
        n=monel(nn)
    endif
    
    call LCOORD (X,XL,IEN(1,n))

    QUAD = (IEN(3,n).ne.IEN(4,n))

    if (DEFORM) then
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,135)
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
    endif


    call wrelms (XL,DL,XPT,QUAD,DEFORM,deffac,NXPT)
    
enddo

! itype 4 -> contour data
call wGMTfile (XPT,NXPT,4,lu,ascout,stderr,ierr)

 endif
!
 return
 end
!-------------------------------------------------------------------------------
subroutine OutSrf (ISURF,SURF,SURF1,X,D,NSD,NDOF,lu,luMon, &
                   verbose,sedsav)

USE ALGEBRA
USE MODELDEFINITION
use monitorVars,     only: nnp, npmon, monnp

implicit none
!-pass
integer NSD,NDOF,lu,luMon

integer :: isurf
double precision :: x, d, surf, surf1
dimension X(NSD,*),D(NDOF,*),ISURF(*),SURF(NSD,*),SURF1(NSD,*)

logical sedsav,verbose

!-locl
double precision :: SURFL(3)
integer :: i, j, n, nn

if (NSURF.le.0) then
    return
endif

if (npmon) then
    do i=1,NSURF
        n = ISURF(i)

!       Occurs in monitor list?
        j = 1
100     if (j.gt.nnp) then
            cycle
        endif

        nn = monnp(j)
        if (n.eq.nn) then
            goto 200
        endif
        j = j + 1
        goto 100

200     if (sedsav) then
            call EQUATE (SURFL,SURF1(1:NSD,i),NSD)
            call UPDATE (SURFL,X(1,n),NSD)
            call UPDATE (SURFL,D(1,n),NSD)
        else
            call EQUATE (SURFL,SURF(1:NSD,i),NSD)
        endif

        if (verbose) then
            write(luMon,10) n,(SURFL(j),j=1,NSD)
10               format(1x,'Erosion surface node ',I5,': ', &
             3(1PG14.6,3X))
        else
            write(luMon,*) (SURFL(j),j=1,NSD)
        endif

    enddo
else
    write(lu,*) 2,NSURF
    do i=1,NSURF
        n = ISURF(i)
        if (sedsav) then
            call EQUATE (SURFL,SURF1(1:NSD,i),NSD)
            call UPDATE (SURFL,X(1,n),NSD)
            call UPDATE (SURFL,D(1,n),NSD)
        else
            call EQUATE (SURFL,SURF(1:NSD,i),NSD)
        endif
        write(lu,*) (SURFL(j),j=1,NSD)
    enddo
endif

return
end subroutine
!-----------------------------------------------------------------------
real function D2REAL (x)

double precision x
double precision, parameter :: FLMIN=1.17549E-38, FLMAX=3.40282E+38

if (x.ge.FLMAX) then
    D2REAL = FLMAX
else if (x.le.-FLMAX) then
    D2REAL = -FLMAX
else if (x.gt.0. .and. x.le.FLMIN) then
    D2REAL = FLMIN
else if (x.lt.0. .and. x.ge.-FLMIN) then
    D2REAL = -FLMIN
else
    D2REAL = REAL(x)
endif

return
end
!-------------------------------------------------------------------------------
subroutine conn (lu)
!
! subroutine to write connection file for plotting with gmt to output file
!
use iomodule,        only: stderr
USE MODELDEFINITION, only: numel
use postProcessData, only: GMT,connout
use monitorVars,     only: elmon

 implicit none
!-pass
 integer lu


!-local
integer n, j, ntot

!-init
 ntot = NUMEL
 if (elmon) then
write (stderr,10)
10        format(1x,'connection file not relevant in monitor mode')
return
 endif
!
 if (.not.GMT) then
write (stderr,20)
20       format(1x,'connection file not relevant for plotting with', &
     1x,'"P"')
return
 endif
!
 j = 0
 do n=1,ntot
    write (lu,*,err=100) j , j+1 , j+2
    j = j+3
 enddo
!
 return
100    call prerr("conn")
 return
 end
