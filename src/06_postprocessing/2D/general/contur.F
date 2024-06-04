subroutine cntdis (IEN,X,D,DELD,NSLIP,DXE,DELDXE,NFAULT,LMF, &
                   TFAULT,DFAULT,SKEW,IMPRINT,ITPRINT, &
                   MAXSTP,DELT,XPT,NDOF,NUMTR,DEFORM,VERBOSE,LUMON,DEFFAC)
! cntdis ?= contour displacement

! this is used do either print the effective velocity,
! of the effective displacement


USE ALGEBRA
USE CONSTANTS,       only: zero, deg2rad, useDELDXE, useDXE, &
                           disp, velo, &
                           total, normal, slip, split, &
                           coordx, coordy, coordz, coordr, effective, &
                           isAxiPlt
USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELDEFINITION, only: NTYPE,NUMEL,NUMAT,IOPT,IPOINT,NUMPR,NUMSLP, &
                           NUMFN,NPRE,LGDEF,IRESDU,NUMSTR,IGRAV,IVLIM, &
                           NUMWNK,NSURF,NSED,INCOMP,NOCOMPR,NSLSKEW,ishell,radius

use debugmodule,     only: debug
use modeltopology,   only: nen, nsd
use monitorVars,     only: nnp, npmon, monnp
use plotControl,     only: jShell, printThisDisp
use documentation,   only: plotVariety
!
! processes scalar quantities for displacements or velocities
!
implicit none
integer, parameter :: NPLT=3, NPLDOF=3, NSUR=4
!-pass
logical DEFORM,VERBOSE
integer NDOF,LUMON
double precision :: deffac

integer          :: ien, nslip, lmf, imprint, itprint, maxstp
double precision :: x, d, deld, dxe, deldxe, nfault, tfault, dfault, skew, delt, xpt
dimension IEN(NEN,*),X(NSD,*), &
 D(NDOF,*),DELD(NDOF,*), &
 NSLIP(5,*),DXE(NDOF,NEN,*),DELDXE(NDOF,NEN,*), &
 NFAULT(3,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),DFAULT(NDOF,*), &
 SKEW(2,*), &
 IMPRINT(*),ITPRINT(*),MAXSTP(*),DELT(*), &
 XPT(NPLDOF,NPLT,*)
!-locl
integer :: n, nn, i, j, k
logical processed,firstentry
character(len=35) string
double precision :: xl(NDOF,NEN), dl(NSD,NEN), val(NEN)
save firstentry
integer :: kk,l,ll,mm,ne
integer :: numtr
double precision :: de, def
integer :: ifn, isn
double precision :: ylat

!-init
data firstentry/.true./
NUMTR  = 0


def = ZERO
if (DEFORM) then
    def = deffac
endif

firstentry = .false.

if (.not.npmon) then
    ! no monitoring of nodes, print the appropriate displacement
    ! of all the nodes so that a contour plot can be made.

    do nn=1,meshdatactx%Neglobal
        n = nn

        call LCOORD(X,xl,IEN(1,n))

        if (def.gt.0d0) then
            call LDISP (dl,D,IEN(1,n),NDOF,NEN)
            call ADDSNE (dl,n,NDOF,NEN,useDXE)
            call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,141)
            if (JSHELL.ne.0) then
                call SPHDISP (xl,dl)
            endif
            call REZONE(xl,dl,def)
        endif

        call CLEAR (dl,NDOF*NEN,"dl")
        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq.    normal) then
            if (printThisDisp%dispOrVelo .eq. disp) then
                call LDISP (dl,D,IEN(1,n),NDOF,NEN)
            else
                call LDISP (dl,DELD,IEN(1,n),NDOF,NEN)
            endif
        endif

        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq. slip) then
            if (printThisDisp%dispOrVelo .eq. disp) then
               call ADDSNE (dl,n,NDOF,NEN,useDXE)
            else
               call ADDSNE (dl,n,NDOF,NEN,useDELDXE)
            endif
        endif

        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq. split) then
            if (printThisDisp%dispOrVelo .eq. disp) then
                call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,142)
            else
                call AddFaultDisplacement (dl,n,DFAULT,NDOF,NEN,143)
            endif
        endif

        if (printThisDisp%coord .eq. effective) then
            do i=1,NEN
!                k = 1+(i-1)*NDOF
!                val(i) = SQRT(DOT(dl(k),dl(k),NDOF))
                val(i) = norm2(dl(:,i))
				

!                if (printThisDisp%dispType .eq. total) then
!                    if (abs(dl(k)).ge.abs(dl(k+1))) then
!                        if (dl(k).lt.0.) val(i) = -val(i)
!                    else
!                        if (dl(k+1).lt.0.) val(i) = -val(i)
!                    endif
!                endif

            enddo
        else if (printThisDisp%coord .eq. coordx .or. &
                 printThisDisp%coord .eq. coordr) then
            do i=1,NEN
!                k = 1+(i-1)*NDOF
!                val(i) = dl(k)
                val(i) = dl(1,i)
            enddo
        else if (printThisDisp%coord .eq. coordy) then
            do i=1,NEN
!                k = 2+(i-1)*NDOF
!                val(i) = dl(k)
                val(i) = dl(2,i)
            enddo
        else if (printThisDisp%coord .eq. coordz) then
            do i=1,NEN
				if (plotVariety.eq.isAxiplt) then
                    ! the only 2D variety with a z-component
					val(i) = dl(2,i)
                else
                    val(i) = ZERO ! no 3D should happen here. (maybe opn?)
                endif
            enddo
        endif

        ! elmscalar put the coordinates in xl and the values in val
        ! in the XPT array, so that the XPT array can be sent to 
        ! subroutine writria that write the data to output file.
        call elmscalar (IEN(1,n),xl,val,XPT,NUMTR)
    enddo 
else 

    ! Nodes are being monitored.
    ! We are no longer printing an actual contour plot,
    ! but just the list of nodes.
    ! However when monitored nodes are slippery nodes, or faulted,
    ! and the slippery of faulted displacement is requested,
    ! they have to be printed multiple times for each element.

    do nn=1,nnp
        n = monnp(nn)
        processed = .false.



        if (NUMSLP.gt.0.and.(printThisDisp%dispType .eq. total .or. &
                             printThisDisp%dispType .eq. slip)) then

            ! Slippery displacement must be added.
            ! Loop through the slippery array to study
            ! it's connectivity print a line for each element.

            do isn=1,NUMSLP
                if (NSLIP(2,isn).ne.n) then
                    cycle
                endif
                processed = .true.
                ne = NSLIP(1,isn)
                if (printThisDisp%dispType .eq. total) then
                          if (printThisDisp%dispOrVelo .eq. disp) then
                       call LDISP (dl,D,IEN(1,ne),NDOF,NEN)
                    else
                       call LDISP (dl,DELD,IEN(1,ne),NDOF,NEN)
                    endif
                else
                    call CLEAR (dl,NDOF*NEN,"dl")
                endif

                if (printThisDisp%dispOrVelo .eq. disp) then
                   call ADDSNE (dl,ne,NDOF,NEN,useDXE)
                else
                   call ADDSNE (dl,ne,NDOF,NEN,useDELDXE)
                endif


                if (printThisDisp%coord .eq. effective) then
                    do i=1,NEN
!                        k = 1+(i-1)*NDOF
!                        val(i) = SQRT(DOT(dl(k),dl(k),NDOF))
                         val(i) = norm2(dl(:,i))
                    enddo
                else if (printThisDisp%coord .eq. coordx .or. &
                         printThisDisp%coord .eq. coordr) then
                    do i=1,NEN
!                        k = 1+(i-1)*NDOF
!                        val(i) = dl(k)
                         val(i) = dl(1,i)
                    enddo
                else if (printThisDisp%coord .eq. coordy) then
                    do i=1,NEN
!                        k = 2+(i-1)*NDOF
!                        val(i) = dl(k)
                         val(i) = dl(2,i)
                    enddo
                else if (printThisDisp%coord .eq. coordz) then
                    do i=1,NEN
                        if (NDOF.eq.3) then
!                            k = 3+(i-1)*NDOF
!                            val(i) = dl(k)
                            val(i) = dl(3,i)
                        else
                            val(i) = ZERO
                        endif
                    enddo
                endif
                k = 1
                do while (k.le.NEN.and.IEN(k,ne).ne.n)
                    k = k + 1
                enddo
                if (VERBOSE) then
                    if (printThisDisp%dispType .eq. total) then
                        string = 'total '
                        l = 7
                    else if (printThisDisp%dispType .eq. slip) then
                        string = 'slippery '
                        l = 10
                    endif

                    if (printThisDisp%coord .eq. effective) then
                        string(l:l+13) = 'magnitude of '
                        l = l + 13
                    else if (printThisDisp%coord .eq. coordx .or. &
                             printThisDisp%coord .eq. coordr) then
                        string(l:l+2) = 'x-'
                        l = l + 2
                    else if (printThisDisp%coord .eq. coordy) then
                        string(l:l+2) = 'y-'
                        l = l + 2
                    else if (printThisDisp%coord .eq. coordz) then
                        string(l:l+2) = 'z-'
                        l = l + 2
                    endif

                    if (printThisDisp%dispOrVelo .eq. disp) then
                        string(l:l+12) = 'displacement'
                        l = l + 12
                    else
                        string(l:l+8) = 'velocity'
                        l = l + 8
                    endif

                    write(luMon,3,advance='no') string(1:len(trim(string)))
 3                       format(1x,a)
                    write(luMon,4) n,ne,val(k)
 4                       format(' of node ',I12,' in element ',I12, &
                     ' = ',1PG15.6)
                else
                    write(luMon,*) val(k)
                endif
            enddo
        endif

        if (processed) then 
            cycle
        endif


        if (NUMFN.gt.0.and. (printThisDisp%dispType .eq. total .or. &
                             printThisDisp%dispType .eq. split)) then


            ! Faulted displacement must be added.
            ! Loop through the faults array to study        
            ! it's connectivity print a line for each element.

            if (debug) then
                write(*,*) "entered faulted block"
            endif

            do ifn=1,NUMFN
                if (NFAULT(2,ifn).ne.n) then
                    cycle
                endif
                processed = .true.
                ne = NFAULT(1,ifn)
                if (printThisDisp%dispType .eq. total) then
                    if (printThisDisp%dispOrVelo .eq. disp) then
                       call LDISP (dl,D,IEN(1,ne),NDOF,NEN)
                    else
                       call LDISP (dl,DELD,IEN(1,ne),NDOF,NEN)
                    endif
                else
                    call CLEAR (dl,NDOF*NEN,"dl")
                endif

                      if (printThisDisp%dispOrVelo .eq. disp) then
                    call AddFaultDisplacement (dl,ne,TFAULT,NDOF,NEN,144)
                else
                    call AddFaultDisplacement (dl,ne,DFAULT,NDOF,NEN,145)
                endif

                if (printThisDisp%coord .eq. effective) then
                    do i=1,NEN
!                        k = 1+(i-1)*NDOF
!                        val(i) = SQRT(DOT(dl(k),dl(k),NDOF))
                        val(i) = norm2(dl(:,i))
 
                   enddo
                else if (printThisDisp%coord .eq. coordx) then
                    do i=1,NEN
!                        k = 1+(i-1)*NDOF
!                        val(i) = dl(k)
                         val(i) = dl(1,i)

                    enddo
                else if (printThisDisp%coord .eq. coordy) then
                    do i=1,NEN
!                        k = 2+(i-1)*NDOF
!                        val(i) = dl(k)
                        val(i) = dl(2,i)
                    enddo
                else if (printThisDisp%coord .eq. coordz) then
                    do i=1,NEN
                        if (NDOF.eq.3) then
!                            k = 3+(i-1)*NDOF
!                            val(i) = dl(k)
                            val(i) = dl(3,i)

                        else
                            val(i) = ZERO
                        endif
                    enddo
                endif

                k = 1

                do while (k.le.NEN.and.IEN(k,ne).ne.n)
                    k = k + 1
                enddo

                if (VERBOSE) then

                    if (printThisDisp%dispType .eq. total) then
                        string = 'total '
                        l = 7
                    else if (printThisDisp%dispType .eq. split) then
                        string = 'faulted '
                        l = 9
                    endif

                    if (printThisDisp%coord .eq. effective) then
                        string(l:l+13) = 'magnitude of '
                        l = l + 13
                    else if (printThisDisp%coord .eq. coordx .or. &
                             printThisDisp%coord .eq. coordr) then
                        string(l:l+2) = 'x-'
                        l = l + 2
                    else if (printThisDisp%coord .eq. coordy) then
                        string(l:l+2) = 'y-'
                        l = l + 2
                    else if (printThisDisp%coord .eq. coordz) then
                        string(l:l+2) = 'z-'
                        l = l + 2
                    endif

                    if (printThisDisp%dispOrVelo .eq. disp) then
                        string(l:l+12) = 'displacement'
                        l = l + 12
                    else
                        string(l:l+8) = 'velocity'
                        l = l + 8
                    endif

                    write(luMon,3) string(1:l)
                    write(luMon,4) n,ne,val(k)
                else
                    write(luMon,*) val(k)
                endif

                if (debug) then
                    write(*,*) "leaving faulted block"
                endif

            enddo

        endif

        if (processed) then
            cycle
        endif

        if (printThisDisp%dispType .ne. normal .and. &
            printThisDisp%dispType .ne. total) then
            cycle
        endif

        if (printThisDisp%dispOrVelo .eq. disp) then
            do i=1,NDOF
                val(i) = D(i,n)
            enddo
        else
            do i=1,NDOF
                val(i) = DELD(i,n)
            enddo
        endif

        if (printThisDisp%coord .eq. effective) then
            de = SQRT(DOT(val,val,NDOF))
        else if (printThisDisp%coord .eq. coordx .or. &
                 printThisDisp%coord .eq. coordr) then
            de = val(1)
        else if (printThisDisp%coord .eq. coordy) then
            de = val(2)
        else if (printThisDisp%coord .eq. coordz) then
            if (plotVariety .eq. isAxiPlt) then
                de = val(2)
            else
                write(*,*) "wanting to print z coordinate, maar there is only x and y"
                stop "Leaving GTECTON..."
 !                if (NDOF.eq.3) then
!                    de = val(3)
!                else
!                    de = ZERO
!                endif
            endif


        endif

        if (VERBOSE) then
    
              if (printThisDisp%dispType .eq. total) then
                string = 'total '
                l = 7
            else if (printThisDisp%dispType .eq. normal) then 
                string = 'average '
                l = 9
            endif

            if (printThisDisp%coord .eq. effective) then
                string(l:l+13) = 'magnitude of '
                l = l + 13
            else if (printThisDisp%coord .eq. coordx) then
                string(l:l+2) = 'x-'
                l = l + 2
            else if (printThisDisp%coord .eq. coordy) then
                string(l:l+2) = 'y-'
                l = l + 2
            else if (printThisDisp%coord .eq. coordz) then
                string(l:l+2) = 'z-'
                l = l + 2
            endif

 
            if (printThisDisp%dispOrVelo .eq. disp) then
                string(l:l+12) = 'displacement'
                l = l + 12
            else
                string(l:l+8) = 'velocity'
                l = l + 8
            endif

            write(luMon,3) string(1:l)
            write(luMon,5) n,de
 5               format(" of node ",I5," = ",1PG15.6)
        else
            write(luMon,*) de
        endif

    enddo
endif
!

return
end
!-------------------------------------------------------------------------------
subroutine cntemp (IEN,X,D,XPT,LMF,TFAULT,DXE, &
   T,TLST,LMTX,TX,TLSTX,SKEW,NEN,NSD,NDOF,NUMTR,DEFORM, &
   NPLDOF,NPLT,verbose,iflag,luMon,deffac)

USE MODELDEFINITION, only: NUMEL, NUMFN, NUMSLP
use monitorVars,     only: nnp, monnp, npmon
use plotControl,     only: jShell, iRelt

!
! produce contour bodies for temperatures
! IFLAG = 0,    absolute temperatures
! IFLAG = 1,    incremental temperatures
! IFLAG = 2,    relative temperatures
!
implicit none
!-pass
integer NEN,NSD,NDOF,NPLDOF,NPLT,iflag, lumon, numtr
double precision :: deffac

integer          :: ien, lmf, lmtx
double precision :: x, d, t, tlst, xpt, tfault, dxe, tx, tlstx, skew
dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),T(*),TLST(*), &
 XPT(NPLDOF,NPLT,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   DXE(NDOF,NEN,*),LMTX(NEN,*),TX(*),TLSTX(*),SKEW(2,*)

logical DEFORM,verbose

!-locl
double precision :: XL(24),DL(24),TL(8)
integer          :: n, nn, i, k
double precision :: effD, def
!-init

NUMTR  = 0
def = 0d0

if (DEFORM) then
    def = deffac
endif

if (.not.npmon) then
    do nn=1,NUMEL

        n = nn
        call LCOORD(X,XL,IEN(1,n))

        if (DEFORM) then
            call LDISP (DL,D,IEN(1,n),NDOF,NEN)
            call ADDSNE (DL,n,NDOF,NEN,0)
            call AddFaultDisplacement(DL,n,TFAULT,NDOF,NEN,146)
            if (JSHELL.ne.0) call SPHDISP(XL,DL)
            call REZONE (XL,DL,def)
        endif

        call LTEMP (IEN(1,n),T,TL,NEN)

        if (iflag.ne.0) then
            do i=1,NEN
                k = IEN(i,n)
                TL(i) = TL(i) - TLST(k)
            enddo
        endif

        call elmscalar (IEN(1,n),XL,TL,XPT,NUMTR)

    enddo
else

    do nn=1,nnp
        n = monnp(nn)

        if (iflag.eq.0) then
            effD = T(n)
        else if (iflag.eq.1) then
            effD = T(n) - TLST(n)
        else if (iflag.eq.2) then
            effD = T(n) - TLST(n)
        endif

        if (verbose) then
            if (iflag.eq.0) write(luMon,5) n,effD
 5              format(1x,'temperature at node ',I5, &
                ' = ',1PG15.6)
            if (iflag.eq.1) write(luMon,6) n,effD
 6              format(1x,'temperature increment at node ',I5, &
                ' = ',1PG15.6)
            if (iflag.eq.2) write(luMon,7) irelt,n,effD
 7              format(1x,'temperature change since t=',i5,1x, &
                'at node ',I5,' = ',1PG15.6)
        else
            write(luMon,*) effD
        endif

    enddo
endif

return
end
!-------------------------------------------------------------------------------
subroutine elmscalar (IEN,XL,VAL,XPT,NUMTR)

USE ALGEBRA
use constants, only: zero

!
! contours nodal scalar quantities per element
!
implicit none

integer, parameter :: NPLT=3
integer, parameter :: NPLDOF=3
integer, parameter :: NSD=2
integer, parameter :: NEN=4
!-pass
integer            :: NUMTR
integer            :: ien
double precision   :: xpt, xl, val
dimension          :: XPT(NPLDOF,NPLT,*),XL(NSD,NEN),VAL(NEN),IEN(NEN)
!-locl
logical            :: QUAD
double precision   :: XC(NSD),SH(3,NEN),XS(NSD,NSD)
integer            :: i, j, ierr
double precision   :: det, valc
!-init
QUAD = (IEN(3).ne.IEN(4))
!
if (QUAD) then
!-   calculate extra point
    call SHAP20 (zero,zero,XL,DET,SH,XS,QUAD,ierr)
    call CLEAR (XC,NSD,"XC")
    VALC = zero
    do j=1,NEN
        do i=1,NSD
            XC(i) = XC(i) + SH(3,j)*XL(i,j)
        enddo
        VALC = VALC + SH(3,j)*VAL(j)
    enddo

    NUMTR = NUMTR + 1

    XPT(1,1,NUMTR) = XL(1,1)
    XPT(2,1,NUMTR) = XL(2,1)
    XPT(3,1,NUMTR) = VAL(1)
    XPT(1,2,NUMTR) = XL(1,2)
    XPT(2,2,NUMTR) = XL(2,2)
    XPT(3,2,NUMTR) = VAL(2)
    XPT(1,3,NUMTR) = XC(1)
    XPT(2,3,NUMTR) = XC(2)
    XPT(3,3,NUMTR) = VALC

    NUMTR = NUMTR + 1

    XPT(1,1,NUMTR) = XL(1,2)
    XPT(2,1,NUMTR) = XL(2,2)
    XPT(3,1,NUMTR) = VAL(2)
    XPT(1,2,NUMTR) = XC(1)
    XPT(2,2,NUMTR) = XC(2)
    XPT(3,2,NUMTR) = VALC
    XPT(1,3,NUMTR) = XL(1,3)
    XPT(2,3,NUMTR) = XL(2,3)
    XPT(3,3,NUMTR) = VAL(3)

    NUMTR = NUMTR + 1

    XPT(1,1,NUMTR) = XL(1,3)
    XPT(2,1,NUMTR) = XL(2,3)
    XPT(3,1,NUMTR) = VAL(3)
    XPT(1,2,NUMTR) = XL(1,4)
    XPT(2,2,NUMTR) = XL(2,4)
    XPT(3,2,NUMTR) = VAL(4)
    XPT(1,3,NUMTR) = XC(1)
    XPT(2,3,NUMTR) = XC(2)
    XPT(3,3,NUMTR) = VALC

    NUMTR = NUMTR + 1

    XPT(1,1,NUMTR) = XL(1,4)
    XPT(2,1,NUMTR) = XL(2,4)
    XPT(3,1,NUMTR) = VAL(4)
    XPT(1,2,NUMTR) = XC(1)
    XPT(2,2,NUMTR) = XC(2)
    XPT(3,2,NUMTR) = VALC
    XPT(1,3,NUMTR) = XL(1,1)
    XPT(2,3,NUMTR) = XL(2,1)
    XPT(3,3,NUMTR) = VAL(1)

 else

    NUMTR = NUMTR + 1

    XPT(1,1,NUMTR) = XL(1,1)
    XPT(2,1,NUMTR) = XL(2,1)
    XPT(3,1,NUMTR) = VAL(1)
    XPT(1,2,NUMTR) = XL(1,2)
    XPT(2,2,NUMTR) = XL(2,2)
    XPT(3,2,NUMTR) = VAL(2)
    XPT(1,3,NUMTR) = XL(1,3)
    XPT(2,3,NUMTR) = XL(2,3)
    XPT(3,3,NUMTR) = VAL(3)

endif

return
end
