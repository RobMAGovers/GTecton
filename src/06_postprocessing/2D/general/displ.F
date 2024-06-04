subroutine dispar (X,D,DELD,NSLIP,DXE,DELDXE, &
                   NFAULT,LMF,TFAULT,DFAULT,IEN,SKEW,IMPRINT,ITPRINT, &
                   MAXSTP,DELT,DEFFAC,HEADSIZE,HEADANGLE,DSB,DSCALE,IERR, &
                   LUMON,IBASE,ISTEP,IARCONT,LU,DEFORM,VERBOSE)

USE MODELDEFINITION
use modeldatamodule, only: modeldatactx
USE MESHDATAMODULE
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use debugmodule,     only: iecho
use monitorVars,     only: nnp, monnp, npmon
use errors,          only: ok
use constants,       only: disp, velo, &
                           total, normal, slip, split, &
                           coordx, coordy, coordz, coordr, effective
use plotControl,     only: printThisDisp
!
implicit none
!-pass
integer          :: IERR,LUMON,IBASE,ISTEP,IARCONT,LU
logical          :: DEFORM,VERBOSE
double precision :: deffac, HEADSIZE,HEADANGLE, dsb, dscale

integer          :: ien, lmf, nslip, nfault, imprint, itprint, maxstp
double precision :: d, dxe, skew, tfault, deld, deldxe, x, dfault, delt
dimension IEN(NEN,*),D(NDOF,*),DXE(NDOF,NEN,*),SKEW(2,*), &
 LMF(NDOF,NEN,*),TFAULT(NDOF,*),DELD(NDOF,*),DELDXE(NDOF,NEN,*), &
 X(NSD,*),DFAULT(NDOF,*),NSLIP(5,*),NFAULT(3,*), &
 IMPRINT(*),ITPRINT(*),MAXSTP(*),DELT(*)




!-locl
logical slippery,faulted,firstentry
save firstentry
integer          :: i, luo
integer          :: n, nn, nbase, npnum, nstep
double precision :: def
double precision :: hatangle, hatsize
!-init


data firstentry/.true./

ierr   = OK
def    = 0d0

if (DEFORM) then
    def = deffac
endif

if (firstentry) then
  firstentry = .false.
endif

if (NPMON) then
    npnum = nnp
    luo   = LUMON
    nstep = 1
    nbase = 1
else
    npnum = meshdatactx%Nvglobal
    luo   = LU
    nbase = IBASE
    nstep = ISTEP
endif

hatangle = HEADANGLE
hatsize  = DSB  *  DSCALE  * HEADSIZE

! go write the file for the arrows, with starting locations, and distances
! in x and y.

! nstep, because it is possible to not print every arrow, but only the one of for every
! nstep-th vertex.


do nn=nbase,npnum,nstep

    n = nn

    if (NPMON) then
        n = monnp(nn)
    endif

!   slippery node displacements

    if (NUMSLPglobal.gt.0.and.(printThisDisp%dispType .eq. total .or. &
                               printThisDisp%dispType .eq. slip)) then
        slippery = .true.

        call DISP2 (modeldatactx%NSLIP,5,NUMSLP,n,slippery,def, &
                    luo,hatsize,Dscale,hatangle,IARCONT,VERBOSE,IEN,X,D, &
                    DXE,SKEW,LMF,TFAULT,DFAULT,DELD, &
                    DELDXE,2,nn,ierr)

        if (IERR.ne.OK) then
            return
        endif
    else
        slippery = .false.
    endif

    if (slippery) then
        cycle
    endif

!   faulted node displacements

    if (NUMFNglobal.gt.0.and. (printThisDisp%dispType .eq. total .or. &
                               printThisDisp%dispType .eq. split)) then
        faulted = .true.

        call DISP2 (NFAULT,3,NUMFN,n,faulted, &
                    def,luo,hatsize,Dscale,hatangle,IARCONT,VERBOSE, &
                    IEN,X,D,DXE,SKEW,LMF,TFAULT, &
                    DFAULT,DELD,DELDXE,1,nn,IERR)

        if (IERR.ne.OK) return
    else
        faulted = .false.
    endif

    if (faulted) then
        cycle
    endif



!   regular displacements
    if (printThisDisp%dispType .eq. total .or. &
        printThisDisp%dispType .eq. normal) then

        call DISP2 (NFAULT,3,NUMFN,n, &
                    faulted,def,luo,hatsize,Dscale, &
                    hatangle,IARCONT,VERBOSE,IEN,X,D,DXE, &
                    SKEW,LMF,TFAULT,DFAULT,DELD,DELDXE,3,nn,IERR)


        if (IERR.ne.OK) then
            return
        endif
    endif

enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine DISP2 (NSlipOrFault,NDIM,NUMSF,NODE,ELMBASED,DEF, &
 LUO,HATSIZE,SCALE,HATANGLE,IARCONT,VERBOSE,IEN,X,D,DXE, &
 SKEW,LMF,TFAULT,DFAULT,DELD,DELDXE,callerID,objectID,IERR)

USE ALGEBRA
USE MODELDEFINITION
USE CONSTANTS,       only: deg2rad, useDXE, useDELDXE, &
                           disp, velo, &
                           total, normal, slip, split, &
                           coordx, coordy, coordz, coordr, effective

use meshdatamodule
use modeldatamodule, only: modeldatactx
use debugmodule,     only: iecho
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: ascout, stderr
use postProcessData, only: GMT,connout
use monitorVars,     only: npmon
use errors,          only: exec, ok
use plotControl,     only: jShell, printThisDisp

!
! processes vector quantities for displacement or velocity
!
implicit none
!-pass
integer :: NSlipOrFault, iarcont, ierr
double precision DEF,hatsize,scale,hatangle
logical ELMBASED,VERBOSE
integer NDIM,NUMSF,NODE,luo
! this sub is both called for slippery and for faulted nodes, hence NSlipOrFault

integer          :: ien,lmf
double precision :: d,dxe,skew,tfault,dfault,deld,deldxe,x
dimension NSlipOrFault(NDIM,*)
dimension IEN(NEN,*),D(NDOF,meshdatactx%nvglobal),DXE(NDOF,NEN,*), &
 SKEW(2,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),DFAULT(NDOF,*), &
 DELD(NDOF,*),DELDXE(NDOF,NEN,*),X(NSD,*)
integer :: objectID ! index of the calling vertex/element
integer :: callerID ! 1 = faulted; 2 = slippery, 3 = normal
!-locl
integer, parameter :: MAXRRW=25

double precision  :: length

double precision  :: datavec(3)
double precision  :: eval(3)
double precision  :: evec(3,3)
double precision  :: arrow(3,MAXRRW)
double precision  :: dl(12)
double precision  :: xl(8)
double precision  :: locvec(2)

real, external :: D2REAL
character(len=22)string

integer :: narrw
data arrow(1,1),narrw /0d0,0/

double precision :: deg, ylat

integer :: i, isn, j, kk, l, ll, local, ne

integer, external :: ipltrw
!


! ******************* IF ELMBASED == true ******************
! *** displacements of the elements (and not that of the vertices
! *** is computed
! *********************************************************



if (.not.ELMBASED) then
    goto 300
endif

ELMBASED = .false.


do isn=1,NUMSF ! NUMSF = nr of faulted nodes

    if (NSlipOrFault(2,isn).ne.NODE) then
        cycle
    endif

    ELMBASED = .true.
    ne = NSlipOrFault(1,isn)
    call CLEAR(dl,12,"dl")


!   filling dl to get values for the displacement.
    if (printThisDisp%dispOrVelo .eq. disp) then

        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq. normal) then
            call LDISP (dl,D,IEN(1,ne),NDOF,NEN)
        endif

        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq. slip) then
            call ADDSNE (dl,ne,NDOF,NEN,useDXE)
        endif

        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq. split) then
            call AddFaultDisplacement (dl,ne,TFAULT,NDOF,NEN, 184)
        endif

    else
        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq. normal) then
            call LDISP (dl,DELD,IEN(1,ne),NDOF,NEN)
        endif
        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq. slip) then
            call ADDSNE (dl,ne,NDOF,NEN,useDELDXE)
        endif
        if (printThisDisp%dispType .eq. total .or. &
            printThisDisp%dispType .eq. split) then
            call AddFaultDisplacement (dl,ne,DFAULT,NDOF,NEN, 185)
        endif
    endif

!   identify local node number of NODE in element NE
    do kk=1,NEN
        local = kk
        if (IEN(local,ne).eq.NODE) then
            goto 100
        endif
    enddo

!   set data (displacent/velocity) vector
100 do l=1,NDOF
        ll = l + (local-1)*NDOF
        datavec(l) = dl(ll)
    enddo



    if (npmon) then
        if (VERBOSE) then
            if (printThisDisp%dispType .eq. total) then
                string = 'total '
                l = 7
            else if (printThisDisp%dispType .eq. normal) then
                string = 'average '
                l = 9
            else if (printThisDisp%dispType .eq. slip) then
                string = 'slippery '
                l = 10
            else if (printThisDisp%dispType .eq. split) then
                string = 'faulted '
                l = 9
            endif

            if (printThisDisp%dispOrVelo .eq. disp) then
                string(l:l+12) = 'displacement'
                l = l + 12
            else
                string(l:l+8) = 'velocity'
                l = l + 8
            endif

            write(LUO,1,advance='no') string(1:len(trim(string)))
 1               format(1x,a)

            write(LUO,2) NODE,ne,(datavec(l),l=1,NDOF)
 2               format(' of node ',I12,' in element ',I12,': ', &
             3(1PG15.6))

        else
            write(LUO,*) (datavec(l),l=1,NDOF)
        endif
        cycle
    endif

!   define location vector

    call LCOORD (X,XL,IEN(1,ne))
    call LDISP (dl,D,IEN(1,ne),NDOF,NEN)
    call ADDSNE (dl,ne,NDOF,NEN,useDXE)


    call AddFaultDisplacement (dl,ne,TFAULT,NDOF,NEN, 186)

    ll = (local-1)*NDOF

    if (JSHELL.eq.0) then
        locvec(1) = xl(ll+1) + DEF*dl(ll+1)
        locvec(2) = xl(ll+2) + DEF*dl(ll+2)
    else
        ylat = xl(ll+2)*DEG2RAD
        locvec(1) = xl(ll+1) + DEF*dl(ll+1) / &
                    (DEG2RAD*RADIUS*COS(ylat))
        locvec(2) = xl(ll+2) + DEF*dl(ll+2) / &
                    (DEG2RAD*RADIUS)
    endif

    call EigVec (datavec,eval,evec,NDOF)

    ierr = ipltrw (eval,evec,locvec,SCALE,HATSIZE,HATANGLE, &
     arrow,narrw,MAXRRW,IARCONT,.false.)

    if (ierr.ne.OK) then
        write(stderr,3) NODE
        write(*,3) NODE

 3           format(1x,'disp2: error on node ',I5)
        ierr = EXEC
        return
    endif

    if (GMT) then

          length = datavec(1)**2.0 + datavec(2)**2.0
          length = sqrt(length)*SCALE

          deg = atan2(datavec(2),datavec(1))
          deg = deg/DEG2RAD

          if (ascout) then
            write(LUO,6,err=2000)     locvec(1), &
                                      locvec(2), &
                                      deg, &
                                      length
 6                   format(E15.6,1X,E15.6,1X,E15.6,1X,E15.6)


          else
            write(LUO,err=2000) D2REAL(locvec(1)), &
                                D2REAL(locvec(2)), &
                                D2REAL(deg), &
                                D2REAL(length)
          endif

    else
        stop "non-GMT plotting not implemented"
    endif

    return

enddo
return
!

! ******************* IF ELMBASED == false ******************
! *** displacements of the vertices (and not that of the elements
! *** is computed
! *********************************************************

300     if (.true.) then
!            write(*,*) 'entering vertex displacements', npmon
endif


if (npmon) then

    do l=1,NDOF                
        if (printThisDisp%dispOrVelo .eq. disp) then
            datavec(l)=modeldatactx%D(l,NODE)
        else
            datavec(l)=modeldatactx%DELD(l,NODE)
        endif
    enddo

    if (VERBOSE) then
        if (printThisDisp%dispType .eq. total) then
            string = 'total '
            l = 7
        else if (printThisDisp%dispType .eq. normal) then
            string = 'average '
            l = 9
        else if (printThisDisp%dispType .eq. slip) then
            string = 'slippery '
            l = 10
        else if (printThisDisp%dispType .eq. split) then
            string = 'faulted '
            l = 9
        endif

        if (printThisDisp%dispOrVelo .eq. disp) then
            string(l:l+12) = 'displacement'
            l = l + 12
        else
            string(l:l+8) = 'velocity'
            l = l + 8
        endif

        write(LUO,1,advance="no") string(1:l)

        write(LUO,5) NODE,(datavec(l),l=1,NDOF)

 5           format(" of node ",I5,": ",3(1PG15.6))
    else
        write(LUO,*) (datavec(l),l=1,NDOF)
    endif

    return
endif ! monitoring on.


do l=1,NDOF
   if (printThisDisp%dispOrVelo .eq. disp) then
       ! displacement
       datavec(l)=D(l,NODE)
   else
       ! velocity (DELta D)
       datavec(l)=DELD(l,NODE)
   endif
enddo



call EigVec (datavec,eval,evec,NDOF)



if (JSHELL.eq.0) then
    locvec(1) = X(1,NODE) + DEF*D(1,NODE)
    locvec(2) = X(2,NODE) + DEF*D(2,NODE)
else
    ylat = X(2,NODE)*DEG2RAD

    locvec(1) = X(1,NODE) + DEF*D(1,NODE) / &
     (DEG2RAD*RADIUS*COS(ylat))

    locvec(2) = X(2,NODE) + DEF*D(2,NODE) / &
     (DEG2RAD*RADIUS)

endif

ierr = ipltrw (eval,evec,locvec,SCALE,HATSIZE,HATANGLE,arrow, &
               narrw,MAXRRW,iarcont,.false.)

if (ierr.ne.OK) then
    write(stderr,3) NODE
    ierr = EXEC
    return
endif


if (GMT) then

    length = datavec(1)**2+datavec(2)**2
    length = sqrt(length)*SCALE
    deg = atan2(datavec(2),datavec(1))
    deg = deg/DEG2RAD

    ! actually write the file.
    if (ascout) then
        write(LUO,6,err=2000)     locvec(1), &
                                  locvec(2), &
                                  deg, &
                                  length


    else
        write(LUO,err=2000) D2REAL(locvec(1)), &
                            D2REAL(locvec(2)), &
                            D2REAL(deg), &
                            D2REAL(length)
    endif
else
    stop "non-GMT plotting not implemented"
endif

return
!
2000    call prerr('displ2')
ierr = EXEC
return
end subroutine
