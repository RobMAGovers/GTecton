subroutine xgood (cmd,ncmd,ierr)

USE MODELDEFINITION, only: NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use iomodule
use postProcessData, only: GMT, connout, CST, deffac, deformed, &
                           setria, triang
use monitorVars,     only: elmon
use plotcontrol,     only: nplt
use errors,          only: ok, syntax

! Routine that X-ecutes the "good" command to find the best elements.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)
!-locl
#include "lenchar.i"
character(len=LENCHAR) :: file
integer :: ISIDE
double precision :: side
double precision :: xn

!external
integer, external :: lnblk
!-init
ierr = OK
file = ' '

if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10       format(1x,'Usage: good[lement] file')
    ierr = SYNTAX
    goto 1000
endif

if (ncmd.ge.2) then
    file = cmd(2)
endif

if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,20)
20       format(1x,'good: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif
!
if (triang) then

    write(stderr,51)
51       format(1x,'goodelement: no triangulation available')
    goto 1000
endif


if (setria .and. .not.elmon) then
!           calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
                 plotdatactx%IJSUR,&
                 plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)

    setria = .false.

    if (NUMTR.le.0) then
        write(stderr,60)
60           format(1x,'goodelement: no contour triangles generated')
        goto 1000
    endif

    if (iecho.ge.4) then 
        call tstria (plotdatactx%XPT,NPLT,NUMTR,'trel.dat')
    endif
endif

call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000

call criteria (  meshdatactx%IEN, &
                 plotdatactx%IJSUR,&
                 plotdatactx%XPT, &
                 verbose,luMon,CST,xn, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
 ISIDE,modeldatactx%LMF, &
 modeldatactx%DXE,modeldatactx%SKEW,modeldatactx%tfault,side)

if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 1, ierr)
endif

if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end
!------------------------------------------------------------------------------

subroutine xfemi (cmd,ncmd,ierr)

USE ALGEBRA
USE MODELDATAMODULE
USE MESHDATAMODULE
use modeldefinition, only: ISHELL,MODE,NUMFN
use debugmodule,     only: iecho, parallellog, itest
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule
use postProcessData, only: setria,CST
use plotcontrol,     only: T0calc, jShell, CRDMIN,CRDMAX
use errors,          only: exec

! Routine that X-ecutes the "femi" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
integer fnblnk,lusv
integer, external :: lnblk,nextlu

if (ncmd.gt.1) then
    FEIN = cmd(2)
    call shvar (FEIN)
endif

if (ncmd.gt.2) then
    TEMPIN = cmd(3)
    call shvar (TEMPIN)
endif

if (iecho.eq.1) then
    write(stderr,1) 'reading "',FEIN(fnblnk(FEIN):lnblk(FEIN)),'"'
1   format(1x,80A)
endif

if (ifeout.ne.0) then
!   save file pointer/logical unit of FEDSK
    lusv = luin
    luin = nextlu(20)
endif

! establish the basic geometric dimensions and element type
call INTCON (NSD,NDOF,NSTR,NEN)

! read in system and mechanical data
if (iecho.eq.1 .and. MODE.gt.2) then 
    write(stderr,1) 'reading "',TEMPIN(fnblnk(TEMPIN):lnblk(TEMPIN)),'"'
endif

! read TECIN.DAT
call READFE(NSD,NDOF,NSTR,NEN,ierr)

if (ierr.ne.0) then
    ierr = EXEC
    return
endif

if (itest.eq.1) then
    call ELMTEST (meshdatactx%X, meshdatactx%IEN)
    if (iecho.eq.1) write(*,*) 'elements tested OK'
endif

! recompute differential displacements and test input
if (NUMFN.gt.0) then
    call faultequalize (modeldatactx%NFAULT,modeldatactx%FAULT,NUMFN,NDOF)
endif

! read in thermal data
call RDTHRM (T0calc,ierr)
if (ierr.ne.0) then
    ierr = EXEC
    return
endif

! transfer material properties to elements
call ELMAT ()

ifein  = 1
setria = .true.
if (ifeout.ne.0) then
    luin = lusv
endif
 
if (ISHELL.ne.0) then
    JSHELL = 1
    ISHELL = 0
endif

call coordsMinMax (NSD,CRDMIN,CRDMAX)

if (iecho.eq.2) then
    call tstinp (0)
endif

if (iecho.eq.8) then
    call ParallelLog("elastc", "Finished allocation")
endif

! generate element ij-arrays
if (.not.CST) then
    call elmsys (ierr)
endif

if (iecho.eq.8) then
    call ParallelLog("xfemi", "finished")
endif

return
end
!-------------------------------------------------------------------------------
subroutine xfemo (cmd,ncmd,ierr)

use filemodule,      only: BaseFedsk, &
                           index_filename
use timestepmodule,  only: maxitime, nstep
use debugmodule,     only: iecho, parallelLog
use iomodule
use monitorVars,     only: elmon, npmon
use errors,          only: exec, ok, syntax

!
! Routine that X-ecutes the "femo" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) filtyp
integer openform,ifopen

integer :: i, ireq
integer :: jfile, jform

integer, external :: lnblk, ifbopen ! from tools
logical there
!-init

if (iecho.eq.8) then
    call parallelLog("xfemo","entered sub")
endif

!
if (ifein.eq.0) then
    write(stderr,1)
1   format(1x,'femout: first read FE-input file')
    ierr = EXEC
    goto 200
 endif

2   format(1x,'femo[ut] [file] [asc[ii]|bin[ary]]')

jform = 0
jfile = 0
i = 2

! loop through the command line arguments
do while (i.le.ncmd) 
    if (cmd(i)(1:3).eq.'asc') then
        if (jform.eq.0) then
            ascin = .true.
            jform = 1
        else
            write(stderr,2)
            ierr = SYNTAX
            goto 200
        endif
    else if (cmd(i)(1:3).eq.'bin') then
        if (jform.eq.0) then
            ascin = .false.
            jform = 1
        else
            write(stderr,2)
            ierr = SYNTAX
            goto 200
        endif
    else
        if (jfile.eq.0) then
            BaseFEDSK = cmd(i)
            jfile = 1
            call shvar (BaseFEDSK)
            fedskrt = lnblk(Basefedsk)
            if (Basefedsk(fedskrt:fedskrt).eq.'/') then
                Basefedsk = basefedsk(1:fedskrt)//'FEDSK.DAT'
            endif
        else
            write(stderr,2)
            ierr = SYNTAX
            goto 200
        endif
    endif
    i = i + 1
enddo

! set fedsk name to iteration 0
FEDSK = index_filename(BaseFedsk, maxitime, 0, .false.)
inquire (file=FEDSK,exist=there)
if (.not. there) then
    write(stderr,15) TRIM(FEDSK)
15  format('femout: FEM output file does not exist: ',256A)
    ierr = SYNTAX
    goto 200
endif

if (ifeout.ne.0) then
    if (iecho.eq.2) write(stderr,*) 're-opening output file'
    call closef(luin)
endif

if (jform.eq.0) then
    ierr = openform(luin,FEDSK,ascin)
    if (ierr.ne.OK) then
        goto 200
    endif
else
    if (ascin) then
        ierr = ifopen  (luin,FEDSK,'old')
        if (ierr.ne.0) then
            ierr = EXEC
            goto 200
        endif
    else
        ierr = ifbopen (luin,FEDSK,'old')
        if (ierr.ne.0) then
            ierr = EXEC
            goto 200
        endif
    endif
endif

if (iecho.eq.1) then
    filtyp = 'binary'
    if (ascin) then
        filtyp = 'ASCII'
    endif

    write(stderr,3) filtyp(1:6), &
     '"',FEDSK(1:lnblk(FEDSK)),'"'
3   format(1x,'opened ',A6,' file ',80A)
endif
!
if (iecho.eq.8) then
    call parallelLog("xfemo","reading fedsk")
endif

ireq = NSTEP
call FOUTRD (ireq,ierr)

if (iecho.eq.8) then
    call parallelLog("xfemo","leaving sub")
endif

200    return
end
!-------------------------------------------------------------------------------
subroutine xwfemi (cmd,ncmd,ierr)

! Routine that X-ecutes the "wf" command.
use iomodule
use errors,          only: exec, syntax


implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)


!-locl
integer :: lu
integer, external :: nextlu
#include "lenchar.i"
!

write(*,*) "warning; xwfemi stil has to be validated"

if (ifein.eq.0) then
    write(stderr,1)
 1       format(1x,'wfemi: first read FE-input file')
    ierr = EXEC
    goto 200
endif

if (ncmd.lt.2) then
    write(stderr,2)
 2       format(1x,'wf[emi] file')
    ierr = SYNTAX
    goto 200
endif

lu = nextlu(0)
call openf (lu,cmd(2),'unknown')
call writfe(ierr,lu)
call closef (lu)

200     return
end
!-------------------------------------------------------------------------------

subroutine xError (cmd,ncmd,ierr)

! Compute the norm of the a posteriori error vector following Verfurth J. Comp. Appl. Math. 1994:
! for each element,
! COmpute Sum over F . n
! With F the displacement gradient tensor,
! and the n the outgoing normal vector of every side,
! 
! check the difference with the incoming, summed over the neighbors.

use meshdatamodule,  only: meshdatactx, &
                           elementNeighborTableComplete, &
                           BuildElementNeighbourTable2D
use modeldatamodule, only: modeldatactx
use algebra,         only: vmprd, vectorLength, areaOfTriangleIn3D
use modeltopology,   only: ndof, nsd, nen
use constants,       only: useDXE
use plotControl,     only: jShell


implicit none



integer ncmd,ierr
character cmd*(*)
dimension cmd(*)

integer :: iElement, iSide, iPoint, neighborID
integer :: fileID
double precision, allocatable :: SumFnDiff(:,:)
logical :: elemIsQuadrilateral

double precision :: BP(ndof, nsd) ! ndof and nsd are not the same in axi, opn(?)
!double precision :: BP(meshdatactx%neglobal,ndof, nsd) ! ndof and nsd are not the same in axi, opn(?)

double precision :: sideVec(nsd), normalVec(nsd)
double precision :: xl(nsd,nen), dl(nsd,nen)

double precision :: contributionVec(nsd), contributionNorm

! for debugging:
logical          :: debugThis
integer          :: thisPoint

! get filename from  command line args
if (ncmd.ne.2) then
    write(*,*) "usage: error [filename]"
else
    fileID = 382
    open(unit=fileID, file=cmd(2))
endif

if (.not. elementNeighborTableComplete) then
    call BuildElementNeighbourTable2D()
endif

allocate(SumFnDiff(2,meshdatactx%neglobal))

do iElement = 1, meshdatactx%neglobal

    ! get local coordinates and displacments
    call LCOORD (meshdatactx%X,XL,meshdatactx%IEN(1,iElement), iElement)
    elemIsQuadrilateral = (meshdatactx%IEN(3,iElement).ne.meshdatactx%IEN(4,iElement))

    if (JSHELL.ne.0) then
        call SHLCRD (XL,elemIsQuadrilateral)
    endif

    call LDISP (DL,modeldatactx%D,meshdatactx%IEN(1,iElement),NDOF,NEN)
    call AddFaultDisplacement (DL,iElement,modeldatactx%TFAULT,NDOF,NEN,267)
    call ADDSNE   (DL,iElement,NDOF,NEN,useDXE)
!    call REZONE   (XL,DL,1d0)
    call BPMATRIX (XL,DL,elemIsQuadrilateral,BP,ierr)

    do iSide = 1,3 
        if (meshdatactx%ElementNeighbours(iElement, iSide) .eq. 0) then
            ! the element is on edge of the domain, it has no
            ! neighbor on this side. Skip it.
        else
            ! there is a neighbor here. 
            ! add the contribution of the normal vector, and subtract is from the neighbor.

            if (iSide.eq.1) then ! side from point 1 to 2
                sideVec = xl(:,2) - xl(:,1)
                call getNormal2D (sideVec, normalVec)

            else if (iSide.eq.2) then ! side from point 2 to 3
                sideVec = xl(:,3) - xl(:,2)
                call getNormal2D (sideVec, normalVec)

            else ! (iSide.eq.3) ! side from point 3 to 1
                sideVec = xl(:,1) - xl(:,3)
                call getNormal2D (sideVec, normalVec)

            endif

            call VMPRD(BP, normalVec, contributionVec, 2, 2)
                          
            SumFnDiff(:,iElement) = SumFnDiff(:,iElement) + contributionVec
            neighborID = meshdatactx%ElementNeighbours(iElement, iSide)
            SumFnDiff(:,neighborID) = SumFnDiff(:,neighborID) + contributionVec

        endif
    enddo

enddo

do iElement = 1, meshdatactx%neglobal
    call LCOORD (meshdatactx%x,xl,meshdatactx%IEN(:,iElement), iElement)
    do iPoint = 1, 3
        write(fileID, *) xl(1,iPoint), xl(2,iPoint), vectorLength(SumFnDiff(:,iElement),nsd)
    enddo
enddo
close(fileID)


! possibly print BP matrix values
debugThis = .false.
if (debugThis) then

    open(unit=411, file="BP11.dat")
    open(unit=412, file="BP12.dat")
    open(unit=421, file="BP21.dat")
    open(unit=422, file="BP22.dat")

    do iElement = 1, meshdatactx%neglobal

        call LCOORD (meshdatactx%X,XL,meshdatactx%IEN(1,iElement), iElement)
        elemIsQuadrilateral = (meshdatactx%IEN(3,iElement).ne.meshdatactx%IEN(4,iElement))
        if (JSHELL.ne.0) then
            call SHLCRD (XL,elemIsQuadrilateral)
        endif
        call LDISP (DL,modeldatactx%D,meshdatactx%IEN(1,iElement),NDOF,NEN)
        call AddFaultDisplacement (DL,iElement,modeldatactx%TFAULT,NDOF,NEN,267)
        call ADDSNE   (DL,iElement,NDOF,NEN,useDXE)
    !    call REZONE   (XL,DL,1d0)

        call BPMATRIX (XL,DL,elemIsQuadrilateral,BP,ierr)

        do iPoint = 1,3
            thisPoint = meshdatactx%IEN(iPoint,iElement)
            write(411,*) meshdatactx%x(:, thisPoint), BP(1,1)
            write(412,*) meshdatactx%x(:, thisPoint), BP(1,2)
            write(421,*) meshdatactx%x(:, thisPoint), BP(2,1)
            write(422,*) meshdatactx%x(:, thisPoint), BP(2,2)
        enddo

    enddo

    close(411)
    close(412) 
    close(421) 
    close(422) 

endif


deallocate(SumFnDiff)

end subroutine

subroutine getNormal2D (vec, normal)
! small helper for xError above.
! Take for granted that the sequence of points
! follows the element counterclockwise.
! Otherwise the elm check in tecin, during femi, 
! would have given an error because the determinant would
! have become negative.
use algebra, only: vectorLength
implicit none
double precision :: vec(2), normal(2)

! normals point outside of the triangle
normal(1) =  vec(2)
normal(2) = -vec(1)

! normalize the normal
normal = normal / vectorLength(normal,2)

end subroutine

!-------------------------------------------------------------------------------

subroutine xpartition (cmd,ncmd,ierr)
! for now, this only works on the initial mesh


use meshdatamodule, only: meshdatactx


implicit none

! ! for append and ifein

!-pass  
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)
!-local
integer iElem, iNode, nodeID
integer partitionFileID
integer partitionID

integer, external :: nextlu
double precision :: x, y


!if (ifein.eq.0) then
!    write(stderr,*) 'partition: first read FE-input file'
!    ierr = EXEC
!    return
!endif

!if (ncmd.lt.2) then
!    write(stderr,2)
! 2       format(1x,'part[ition] file')
!    ierr = SYNTAX
!    return
!endif



partitionFileID = nextlu(0)
open(unit=partitionFileID, file=cmd(2))

do iElem=1,meshdatactx%neglobal
    do iNode=1,3
        nodeID = meshdatactx%ien(iNode, iElem)
        x = meshdatactx%X(1,nodeID)
        y = meshdatactx%X(2,nodeID)
        partitionID = meshdatactx%e2p(iElem)

        write(partitionFileID,*) x, y, partitionID

    enddo
enddo

call closef(partitionFileID)


end subroutine

!-------------------------------------------------------------------------------

subroutine xInvariant(cmd,ncmd,ierr)

use modeldefinition, only: numel
USE MODELDATAMODULE, only: modeldatactx
USE MESHDATAMODULE,  only: meshdatactx
use materialsmodule, only: matdatactx, &
                           MaterialIDOfElement

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)

!local
#include "lenchar.i"
character(len=LENCHAR) fileName
integer :: invariantType
integer :: iElem
integer, external :: lnblk
double precision :: invariant

double precision :: pois
double precision :: stress(3)

double precision, external :: aj1, aj2, aj3

if (cmd(2)(1:1).eq.'1') then
    invariantType = 1
else if (cmd(2)(1:1).eq.'2') then
    invariantType = 2
else if (cmd(2)(1:1).eq.'3') then
    invariantType = 3
endif

fileName = cmd(3)


write(*,*) "ready with", invariantType, trim(fileName), meshdatactx%neglobal, numel

open(unit=48, file=fileName)



do iElem = 1, meshdatactx%neglobal
    pois    = matdatactx%PRPMAT(2,MaterialIDOfElement(iElem)) ! 2 is always POIS.
    stress  = modeldatactx%stn(:,iElem)

    if (invariantType .eq. 1) then
        invariant = AJ1(stress, pois)
    else if (invariantType .eq. 2) then
        invariant = AJ2(stress, pois)
    else if (invariantType .eq. 3) then
        invariant = AJ3(stress, pois)
    endif

    write(48, *) meshdatactx%X(:,meshdatactx%ien(1, iElem)), invariant
    write(48, *) meshdatactx%X(:,meshdatactx%ien(2, iElem)), invariant
    write(48, *) meshdatactx%X(:,meshdatactx%ien(3, iElem)), invariant

enddo

close(48)

end subroutine

!-------------------------------------------------------------------------------





subroutine xgrid (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule
use postProcessData, only: deffac,setria,deformed
use plotcontrol,     only: grnum
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "grid" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk
integer :: i
!-init
i     = 2
file = ' '

10   format(1x,'Usage: gr[id] file [[no]numbers]')

! initialise grnum
grnum =.false.

100     if (i.le.ncmd) then
    if (cmd(i)(1:2).eq.'nu') then
        grnum =.true.
    else if (cmd(i)(1:4).eq.'nonu') then
        grnum =.false.
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
    endif
    i = i + 1
    goto 100
endif

if (lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'grid: first read FE-input file')
    ierr = EXEC
    goto 1000
endif
!
call OfOpen (ierr,file,ascout.or.grnum,.false.,append,.false., &
             lufo,luMon)

if (ierr.ne.OK) goto 1000

call grid (      meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
 plotdatactx%XPT,modeldatactx%LMF,modeldatactx%tfault, &
 modeldatactx%DXE,modeldatactx%SKEW, &
 grnum,deformed,deffac,lufo)

setria = .true.

call OfClose (ierr,file,.false.,lufo,luMon)

1000    return
end
!-------------------------------------------------------------------------------
subroutine xnp (cmd,ncmd,ierr)

USE MODELDATAMODULE, only: modeldatactx
USE MESHDATAMODULE,  only: meshdatactx
use debugmodule,     only: iecho, parallelLog
use iomodule
use postProcessData, only: deffac,deformed
use monitorVars,     only: npmon
use plotcontrol,     only: npnum
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "np" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)


  
! for /PLTVAR/ -> deformed,deffac etc



!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk
logical total,tecton
save total
integer :: i
!-init
data total /.true./
i     = 2
file = ' '
tecton = .false.

if (iecho.eq.8) then
    call parallelLog("xnp", "entered sub")
endif


10   format(1x,'Usage: np file [[no]numbers] [tecton]',1x, &
 '[[av]erage|[to]tal] [raw]')

100     if (i.le.ncmd) then
    if (cmd(i)(1:3).eq.'raw') then
!                write (*,*) 'raw detected'
        call npraw(cmd(2))
        return
    else if (cmd(i)(1:2).eq.'nu') then
        npnum = .true.
    else if (cmd(i)(1:4).eq.'nonu') then
        npnum = .false.
    else if (cmd(i)(1:2).eq.'av') then
        total = .false.
    else if (cmd(i)(1:2).eq.'to') then
        total = .true.
    else if (cmd(i).eq.'tecton') then
        tecton = .true.
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
    endif
    i = i + 1
    goto 100
endif
if (lnblk(file).le.0 .and. (.not.npmon)) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'np: first read FE-input file')
    ierr = EXEC
    goto 1000
endif

call OfOpen(ierr,file,.true.,npmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000

if (iecho.eq.8) then
    call parallelLog("xnp", "calling nodals")
endif

call nodals (    meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
 modeldatactx%NSLIP,modeldatactx%DXE,npnum, &
 deformed,deffac,lufo,luMon,verbose,total,tecton)

if (iecho.eq.8) then
    call parallelLog("xnp", "leaving sub")
endif


call OfClose (ierr,file,npmon,lufo,luMon)

1000    return
end
!-------------------------------------------------------------------------------
 subroutine npraw (dumpfile)

USE MESHDATAMODULE, only: meshdatactx

#include "lenchar.i"
character(len=LENCHAR) dumpfile

integer :: i

!        write(*,*) 'writing coordinates to file: ', dumpfile

open(unit=123, file=dumpfile)
do i=1,meshdatactx%nvglobal
    write(123,*) meshdatactx%X(1,i), meshdatactx%X(2,i)
enddo
close(123)

 end subroutine
!-------------------------------------------------------------------------------


subroutine xstress (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NPRE
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: deffac,deformed
use monitorVars,     only: elmon, sidemon
use plotControl,     only: difstn, &
                           hatangle,hatsize, &
                           iarcont, ielbas, ielstp
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "stress" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer tecton
integer, external :: lnblk
integer :: i
save tecton
data tecton /0/

if ((.not.elmon).and.(.not.sidemon).and. ncmd.lt.2) then
    write(stderr,10)
10       format(1x,'Usage: stre[ss] file [[no[dev]] [tecton]')
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'stress: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30       format(1x,'stress: no stresses calculated')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (NPRE.eq.0) then
        write(stderr,40)
40           format(1x,'stress: first read FE-output file')
        ierr = EXEC
        goto 1000
    else
        write(stderr,50)
50           format(1x,'writing pre-stress INPUT data')
    endif
else
    if (.not.Mout) then
        write(stderr,60) NSTEP
60           format(1x,'stress: no mechanical output at itime = ',I5)
        ierr = EXEC
        goto 1000
    endif
endif

i = 2
file = ' '
100     if (ncmd.ge.i) then
    if (cmd(i)(1:5).eq.'nodev') then
        difstn = .false.
    elseif (cmd(i)(1:3).eq.'dev') then
        difstn = .true.
    elseif (cmd(i).eq.'tecton') then
        tecton = 1
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (tecton.gt.0 .and. (.not.difstn)) tecton = 2

call OfOpen(ierr,file,ascout,elmon.or.sidemon.or.(tecton.gt.0), &
 append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000

call strsar (    meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
                 modeldatactx%STN,&
                 modeldatactx%LMF,modeldatactx%tfault, &
 modeldatactx%DXE,modeldatactx%SKEW, &
                 meshdatactx%mat, &
 matdatactx%prpmat,hatangle,hatsize,NEN,NSD,NDOF,NSTR,deformed, &
 deffac,lufo,ierr,ielbas,ielstp,iarcont,verbose,luMon,tecton)

if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon.or.sidemon.or.(tecton.gt.0),lufo, &
 luMon)

1000    return
end
!-------------------------------------------------------------------------------
subroutine xcf (cmd,ncmd,ierr)

use iomodule
use errors,          only: exec, ok, syntax, prgerr

! Routine that X-ecutes the "cf" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
logical there
integer, external :: lnblk, nextlu, ifopen

!-init
ierr = OK

if (ncmd.lt.2) then
    write(stderr,1)
 1       format(1x,'Usage: cf file')
    ierr = SYNTAX
    goto 1000
endif

file = cmd(2)
call shvar(file)
inquire (file=file,exist=there)

if (.not.there) then
    write(stderr,2) 'cf: file does not exist: "', &
     cmd(2)(1:lnblk(cmd(2))),'"'
 2       format(1x,80A)
    ierr = EXEC
    goto 1000
endif

if (ilu+1.gt.MaxStack) then
    write(stderr,3) MaxStack
 3       format(1x,'cf: nesting depth limit = ',I2)
    ierr = EXEC
    goto 1000
endif

if (lucmd(ilu).ne.stdin) then
    fstack(ilu) = fcmd
else
    fstack(ilu) = 'stdin'
endif

ilu = ilu + 1
lucmd(ilu) = nextlu(0)
ierr = ifopen(lucmd(ilu),file,'old')

if (ierr.eq.0) then
    if (lnblk(file).gt.LEN(fcmd)) then
        write(stderr,4)
 4           format(1x,'cf: LEN(fcmd) insufficient')
        ierr = PRGERR
        return
    endif
    fcmd = file(1:lnblk(file))
else
    ierr = EXEC
    ilu  = ilu - 1
    fcmd = fstack(ilu)
    if (fcmd.eq.'stdin') fcmd = ' '
endif

1000    return
end
!-------------------------------------------------------------------------------
subroutine xstrain (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: deffac,deformed
use monitorVars,     only: elmon
! todo, all this stuff should be directly included in sub strnar,
!       and not passed along
use plotControl,     only: Escale,Esb,&
                           hatangle,hatsize, &
                           ielbas,ielstp,iarcont, sdev
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "strain" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk, nextlu, ifopen
integer :: i
!-init
ierr = OK
file = ' '

10   format(1x,'Usage: stra[in] file')

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'strain: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'strain: no strains calculated')
    else
        write(stderr,40)
40           format(1x,'strain: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'strain: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif

i = 2
100     if (i.le.ncmd) then
    if (cmd(i)(1:3).eq.'dev') then
        sdev = .true.
    else if (cmd(i)(1:5).eq.'nodev') then
        sdev = .false.
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
    endif
    i = i + 1
    goto 100
endif

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) goto 1000



call strnar (    modeldatactx%STR, &
                 meshdatactx%IEN,  &
                 meshdatactx%X,    &
                 modeldatactx%D,   &
 modeldatactx%LMF,modeldatactx%LMX, &
 modeldatactx%tfault,modeldatactx%DXE, &
 modeldatactx%SKEW,      meshdatactx%mat, &
 matdatactx%prpmat,Escale,Esb,hatangle,hatsize,NEN, &
 NSD,NDOF,NSTR,deformed,deffac,lufo,ierr,ielbas,ielstp,iarcont, &
 verbose,luMon,sdev)

if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end
!-------------------------------------------------------------------------------
subroutine xrate (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: deffac,deformed
use monitorVars,     only: elmon
! todo include this only in ratear
use plotControl,     only: Rscale,Rsb,iRate, ratlog, &
                           hatangle,hatsize, &
                           ielbas,ielstp,iarcont, sdev
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "rate" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
logical ratdev
save ratdev
integer, external :: lnblk
integer :: i
logical :: totalRate
!-init
data ratdev/.false./
ierr = OK
totalRate = .false.

if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10       format(1x,'Usage: ra[te] file [[no]log]',1x, &
     '[to[tal]|[vi[scous]|pl[astic]] [no[dev]]')
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'rate: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'rate: no rates calculated')
    else
        write(stderr,40)
40           format(1x,'rate: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (NSTEP.le.0 .or. (.not.Mout)) then
    write(stderr,50) NSTEP
50       format(1x,'rate: no strain rates output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif

i    = 2
file = ' '
100     if (ncmd.ge.i) then
    if (cmd(i)(1:2).eq.'lo') then
        ratlog = .true.
    elseif (cmd(i)(1:4).eq.'nolo') then
        ratlog = .false.
    else if (cmd(i)(1:2).eq.'to') then
        totalRate = .true.
        irate  = MOD(irate,10)+10
    else if (cmd(i)(1:2).eq.'vi') then
        irate  = MOD(irate,10)+20
    else if (cmd(i)(1:2).eq.'pl') then
        irate  = MOD(irate,10)+30
    else if (cmd(i)(1:3).eq.'dev') then
        ratdev = .true.
    else if (cmd(i)(1:5).eq.'nodev') then
        ratdev = .false.
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
!
call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000
!

call ratear (     meshdatactx%IEN, &
                  meshdatactx%X,   &
                 modeldatactx%D,  &
                 modeldatactx%STN,&
                   matdatactx%prpmat, &
                        matdatactx%PRPLAS,&
                 meshdatactx%mat, &
#ifndef EXPERIMENTAL_NORATESTORAGE
                modeldatactx%rate, &
#endif
            modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
                 modeldatactx%SKEW,Rscale,Rsb, &
 hatangle,hatsize,NEN,NSD,NDOF,NSTR,deformed,deffac,lufo,ierr, &
 ratlog,ielbas,ielstp,iarcont,verbose,irate,luMon,ratdev,totalRate)


if (ierr.ne.OK) goto 1000
!
call OfClose (ierr,file,elmon,lufo,luMon)
!
1000    return
end
!-------------------------------------------------------------------------------

subroutine xSlipSides(cmd,ncmd,ierr)

USE MODELDATAMODULE, only: modeldatactx
USE MESHDATAMODULE,  only: meshdatactx
use modeldefinition, only: NUMSLPglobal

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)

! write the elements to file, and give the +1, -1, or o for slippery elements

integer, allocatable :: slipSides(:)
character(len=255)   :: WriteToFileName='<undef>'
integer, parameter   :: writeToFileHandle=83
integer              :: elemID, sideID, iElem, iSlip

if (ncmd.ne.2) then
    write(0,1) 
1   format("usage: slipsides filename")
    goto 100
else
    WriteToFileName = trim(cmd(2))
endif
if (WriteToFileName .eq. '<undef>') then
    write(0,1) "usage: slipsides filename"
    goto 100
endif

allocate(slipSides(meshdatactx%neglobal))
slipSides = 0

if (NUMSLPglobal .gt. 0) then
    ! step through the slippery array
    do iSlip = 1, NUMSLPglobal
        elemID = modeldatactx%NSLIPglobal(1,iSlip)
        sideID = modeldatactx%NSLIPglobal(3,iSlip)

        if (slipSides(elemID).eq.0) then
            ! this element has not yet another slippery node to it.
            slipSides(elemID) = sideID
        elseif (slipSides(elemID).eq.-sideID) then
            ! there is a conflict in this elements.
            slipSides(elemID) = 2
        elseif (slipSides(elemID).eq.2) then
            ! this node has already got a prior conflict recorded. Do nothing.
        elseif (slipSides(elemID).eq.sideID) then
            ! do nothing. This element already has thegood side
!           write(*,*) "Something is fishy with the weight of the slippery element"
!           write(*,*) "Should be 1, 0, or -1, but is", slipSides(elemID)
        endif
    enddo
endif

! write them to file:
open(unit=writeToFileHandle, file=WriteToFileName)

do iElem=1,meshdatactx%neglobal
    write(writeToFileHandle,*) meshdatactx%x(1,meshdatactx%ien(1,iElem)), &
                               meshdatactx%x(2,meshdatactx%ien(1,iElem)), slipSides(iElem)
    write(writeToFileHandle,*) meshdatactx%x(1,meshdatactx%ien(2,iElem)), &
                               meshdatactx%x(2,meshdatactx%ien(2,iElem)), slipSides(iElem)
    write(writeToFileHandle,*) meshdatactx%x(1,meshdatactx%ien(3,iElem)), &
                               meshdatactx%x(2,meshdatactx%ien(3,iElem)), slipSides(iElem)
enddo

close(writeToFileHandle)


100 continue
end subroutine


!-------------------------------------------------------------------------------


subroutine xestress (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule
use documentation,   only: plotVariety
use constants,       only: isPlnplt, isAxiplt, isOpnplt
use timestepmodule,  only: nstep
use postProcessData, only: GMT, connout, cst, deffac, deformed, setria, triang
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "effstre" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk
integer :: i, iflag
!-init
ierr = OK
file = ' '

if ((.not.elmon) .and. ncmd.lt.2) then
    if (plotVariety.eq.isPlnplt) then
        write(stderr,10)
    else if (plotVariety.eq.isAxiplt) then
        write(stderr,11)
    else if (plotVariety.eq.isOpnplt) then
        write(stderr,12)
    else
        write(stderr,*) "Sorry, xestress did not recognize the plot variety"
        write(stderr,*) "This should not happen. Please contact model support"
    endif

10       format(1x,'Usage: estre[ss] [e|xx|yy|xy|s1|s3] file')
11       format(1x,'Usage: estre[ss] [e|rr|zz|rz|rt|zt|tt] file')
12       format(1x,'Usage of xestress in opnplt not yet defined')

    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'effstress: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'effstress: no stresses calculated')
    else
        write(stderr,40)
40           format(1x,'effstress: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'effstress: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,55)
55       format(1x,'effstress: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif
!
i    = 2
file = ' '
iflag = 0
100     if (ncmd.ge.i) then

    if (plotVariety.eq.isPlnplt) then

        if (cmd(i).eq.'e') then
            iflag = 0
        else if (cmd(i).eq.'xx') then
            iflag = 1
        else if (cmd(i).eq.'yy') then
            iflag = 2
        else if (cmd(i).eq.'xy') then
            iflag = 3
        else if (cmd(i).eq.'s1') then
            iflag = 4
        else if (cmd(i).eq.'s3') then
            iflag = 5
        else
            if (lnblk(file).gt.0) then
                write(stderr,10)
                ierr = SYNTAX
                goto 1000
            endif
            file = cmd(i)
        endif
        i = i + 1

    else if (plotVariety.eq.isAxiplt) then

        if (cmd(i).eq.'e') then
            iflag = 0
        else if (cmd(i).eq.'rr') then
            iflag = 1
        else if (cmd(i).eq.'zz') then
            iflag = 2
        else if (cmd(i).eq.'rz') then
            iflag = 3
        else if (cmd(i).eq.'rt') then
            write(*,*) "r theta stress not yet implemented in axisymmetric code"
        else if (cmd(i).eq.'zt') then
            write(*,*) "z theta stress not yet implemented in axisymmetric code"
        else if (cmd(i).eq.'tt') then
            iflag = 6
        else if (cmd(i).eq.'s1') then
            iflag = 4
        else if (cmd(i).eq.'s3') then
            iflag = 5
        else
            if (lnblk(file).gt.0) then
                write(stderr,10)
                ierr = SYNTAX
                goto 1000
            endif
            file = cmd(i)
        endif
        i = i + 1

    endif
    goto 100

endif

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
!

if (triang) then
!           invoke "triangulate" to generate smoothly interpolated triangles
    call trian ( meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,  &
                modeldatactx%STN, &
                modeldatactx%STR, &
#ifdef EXPERIMENTAL_NORATESTORAGE
     modeldatactx%LMF,modeldatactx%LMX, &
#else
     modeldatactx%rate,modeldatactx%LMF,modeldatactx%LMX, &
#endif
     modeldatactx%tfault, &
                 meshdatactx%mat, & !meshdatactx%mat,
     matdatactx%prpmat, &
     matdatactx%PRPLAS,modeldatactx%DXE,modeldatactx%SKEW,0,deformed, &
     deffac,file)
    goto 1000
endif


if (setria .and. .not.elmon) then
!           calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
    setria = .false.


    if (NUMTR.le.0) then
        write(stderr,60)
60           format(1x,'effstress: no contour triangles generated')
        goto 1000
    endif

    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif

call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) goto 1000

call EFSTRS (meshdatactx%IEN, &
             modeldatactx%STN, &
             plotdatactx%IJSUR,plotdatactx%XPT, &
             meshdatactx%mat, &
             matdatactx%prpmat,verbose,luMon,CST,iflag)

if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 2,ierr)
endif

if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end
!-------------------------------------------------------------------------------
subroutine xestrain (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule
use documentation,   only: plotVariety
use constants,       only: isPlnplt, isAxiplt, isOpnplt
use timestepmodule,  only: nstep
use postProcessData, only: GMT,    connout, &
                           deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "effstra/estra" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk
integer :: i, iflag
!-init
ierr = OK
i    = 2
file = ' '
iflag = 0
!
!        write(*,*) 'entered xestra in xeclib, with NUMTR: ', NUMTR

100     if (ncmd.ge.i) then

    if (plotVariety.eq.isPlnplt) then

        if      (cmd(i).eq.'e') then
            iflag = 0
        else if (cmd(i).eq.'xx') then
            iflag = 1
        else if (cmd(i).eq.'yy') then
            iflag = 2
        else if (cmd(i).eq.'xy') then
            iflag = 3
        else if (cmd(i).eq.'xz') then
            iflag = 4
        else if (cmd(i).eq.'yz') then
            iflag = 5
        else if (cmd(i).eq.'zz') then
            iflag = 6
        else
            if (lnblk(file).gt.0) then
                write(stderr,10)
10               format(1x,'Usage: estra[in] [e|xx|yy|zz|xy|xz|yz]')
                ierr = SYNTAX
                goto 1000
            endif
            file = cmd(i)
        endif

    else if (plotVariety.eq.isAxiplt) then

        if      (cmd(i).eq.'e') then
            iflag = 0
        else if (cmd(i).eq.'rr') then
            iflag = 1
        else if (cmd(i).eq.'zz') then
            iflag = 2
        else if (cmd(i).eq.'rz') then
            iflag = 3
        else if (cmd(i).eq.'rt') then
            iflag = 4
        else if (cmd(i).eq.'zt') then
            iflag = 5
        else if (cmd(i).eq.'tt') then
            iflag = 6
        else
            if (lnblk(file).gt.0) then
                write(stderr,11)
11               format(1x,'Usage: estra[in] [e|rr|zz|rz|rt|zt|tt]')
                ierr = SYNTAX
                goto 1000
            endif
            file = cmd(i)
        endif

    else if (plotVariety.eq.isOpnplt) then
        write(*,*) "Sorry, Estra for opn not yet fixed"
    else
        write(*,*) "Sorry, Did not recognize plot variety."
        write(*,*) "This should not happen. Please contact model support"
    endif

    i = i + 1

    goto 100
endif

!        write(*,*) 'iflag set to: ', iflag

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'effstrain: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'effstrain: no strains calculated')
    else
        write(stderr,40)
40           format(1x,'effstrain: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'effstrain: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,55)
55       format(1x,'effstrain: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif
!


if (triang) then
!            write(*,*) 'doing trian'
!           invoke "triangulate" to generate smoothly interpolated triangles
    call trian (meshdatactx%IEN, &
                meshdatactx%X,  &
               modeldatactx%D, &
               modeldatactx%STN, &
               modeldatactx%STR, &
#ifdef EXPERIMENTAL_NORATESTORAGE
     modeldatactx%LMF,modeldatactx%LMX, &
#else
     modeldatactx%rate,modeldatactx%LMF,modeldatactx%LMX, &
#endif
     modeldatactx%tfault,&
                meshdatactx%mat, &
                 matdatactx%PRPMAT, & 
     matdatactx%PRPLAS,modeldatactx%DXE,modeldatactx%SKEW,1,deformed, &
     deffac,file)
    goto 1000
endif



if (setria .and. .not.elmon) then
!           calculate triangle coordinates in XPT
!            write(*,*) 'doing tricrd'
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
    setria = .false.

!            write(*,*) 'xestra says: NUMTR = ', NUMTR
    

    if (NUMTR.le.0) then
        write(stderr,60)
60           format(1x,'effstrain: no contour triangles generated')
        goto 1000
    endif

    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif

call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) goto 1000

call efstrn (modeldatactx%STR, & 
              meshdatactx%mat, & 
               matdatactx%prpmat, &
              meshdatactx%IEN, & 
              plotdatactx%XPT, &
              plotdatactx%IJSUR,verbose,luMon,CST,iflag)

if (.not.elmon) then
   call writria (plotdatactx%XPT,NUMTR,lufo,3,ierr)
endif

if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end
!-----------------------------------------------------------------------
subroutine xrot (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: GMT,    connout, &
                           deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "rot" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
!-init
ierr = OK
file = ' '
!
100     if (ncmd.ge.2) file = cmd(2)
!
if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
10       format(1x,'Usage: rot[ation] file')
    ierr = SYNTAX
    goto 1000
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'rotation: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'rotation: no mechanical solution calculated')
    else
        write(stderr,40)
40           format(1x,'rotation: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'rotation: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,55)
55       format(1x,'rot: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif
!
if (triang) then
    goto 1000
endif
!
if (setria .and. .not.elmon) then
!           calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
    setria = .false.
    if (NUMTR.le.0) then
        write(stderr,60)
60           format(1x,'rotation: no contour triangles generated')
        goto 1000
    endif
    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif
!
call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000
!


call INFROT (plotdatactx%IJSUR,plotdatactx%XPT, &
                 meshdatactx%IEN, &
                 meshdatactx%X,  &
                 modeldatactx%D, &
      modeldatactx%LMF, &
 modeldatactx%tfault, &
 modeldatactx%SKEW,verbose,luMon,CST)


!
if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 4,ierr)
endif
if (ierr.ne.OK) goto 1000
!
call OfClose (ierr,file,elmon,lufo,luMon)
!
1000    return
end
!-----------------------------------------------------------------------
subroutine xrotr (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: GMT,    connout, &
                           deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "rotr" command. Cmoputer the rotation *rate*

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
!-init
ierr = OK
file = ' '
!
100     if (ncmd.ge.2) then
    file = cmd(2)
endif
!
if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
10       format(1x,'Usage: rotr[ate] file')
    ierr = SYNTAX
    goto 1000
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'rotrate: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'rotrate: no mechanical solution calculated')
    else
        write(stderr,40)
40      format(1x,'rotrate: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50  format(1x,'rotrate: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif

if (NSTEP.eq.0) then
    write(stderr,52) 
52  format(1x,'rotrate: no rotation rates at itime = 0')
    ierr = EXEC
    goto 1000
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,55)
55  format(1x,'rotrate: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif
!
if (triang) then
    goto 1000
endif
!
if (setria .and. .not.elmon) then

!           calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
    setria = .false.
    if (NUMTR.le.0) then
        write(stderr,60)
60      format(1x,'rotrate: no contour triangles generated')
        goto 1000
    endif
    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif
!
call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) then
    goto 1000
endif
!

call INFROTR (plotdatactx%IJSUR,&
              plotdatactx%XPT, &    ! right rotation rate in here
              meshdatactx%IEN, &
              meshdatactx%X,  &
              modeldatactx%D, &
              modeldatactx%DELD, &
              modeldatactx%LMF, &
              modeldatactx%tfault, &
              modeldatactx%SKEW,&
              verbose,luMon,CST)
!
if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 4,ierr)
endif

if (ierr.ne.OK) then
    goto 1000
endif
!
call OfClose (ierr,file,elmon,lufo,luMon)
!
1000    return
end
!-------------------------------------------------------------------------------
subroutine xeffr (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR, NTYPE
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use materialsmodule, only: matdatactx
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: GMT,    connout, &
                           deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt, irate, ratlog
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "era" command.
! for effective strain rate

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
integer :: i
!-init
ierr = OK

10   format(1x,'Usage: era[te] file [[no]log]',1x, &
 '[to[tal]|[vi[scous]|pl[astic]] [e|xx|yy|xy|xz|yz]')
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'erate: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'erate: no strain rates calculated')
    else
        write(stderr,40)
40           format(1x,'erate: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout .or. NSTEP.le.0) then
    write(stderr,50) NSTEP
50       format(1x,'erate: no strain rate output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,55)
55       format(1x,'erate: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif
!
i    = 2
file = ' '
100     if (ncmd.ge.i) then
    if (cmd(i)(1:2).eq.'lo') then
        ratlog = .true.
    elseif (cmd(i)(1:4).eq.'nolo') then
        ratlog = .false.
    else if (cmd(i)(1:2).eq.'to') then
        irate  = MOD(irate,10)+10
    else if (cmd(i)(1:2).eq.'vi') then
        irate  = MOD(irate,10)+20
    else if (cmd(i)(1:2).eq.'pl') then
        irate  = MOD(irate,10)+30
    else if (cmd(i).eq.'e') then
        irate = INT(irate/10)*10
    elseif (cmd(i).eq.'xx') then
        irate = INT(irate/10)*10 + 1
    else if (cmd(i).eq.'yy') then
        irate = INT(irate/10)*10 + 2
    else if (cmd(i).eq.'xy') then
        irate = INT(irate/10)*10 + 3
    else if (cmd(i).eq.'xz'.and.NTYPE.eq.2) then
        irate = INT(irate/10)*10 + 4
    else if (cmd(i).eq.'yz'.and.NTYPE.eq.2) then
        irate = INT(irate/10)*10 + 5
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (triang) goto 1000

if (setria .and. .not.elmon) then
!   calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
                 plotdatactx%IJSUR,plotdatactx%XPT, &
                 modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
                 modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
    setria = .false.
    if (NUMTR.le.0) then
        write(stderr,60)
60           format(1x,'erate: no contour triangles generated')
        goto 1000
    endif
    if (iecho.ge.4) then
       call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif

call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000

if (iecho.ge.4) write(*,*) 'calling efrate with irate=',irate

call efrate (meshdatactx%IEN, &!meshdatactx%ien &
 modeldatactx%STN,matdatactx%prpmat,matdatactx%PRPLAS,meshdatactx%mat, &
 plotdatactx%IJSUR,plotdatactx%XPT,&
#ifndef EXPERIMENTAL_NORATESTORAGE
modeldatactx%rate,&
#endif
ratlog,verbose,irate,luMon,CST)

if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 5,ierr)
endif


if (ierr.ne.OK) then
    goto 1000
endif

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end
!-------------------------------------------------------------------------------
subroutine xte (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, ITMODE, NUMTR
use modeldatamodule
use meshdatamodule
use debugmodule,     only: iecho, parallellog
use iomodule
use modeltopology,   only: ndof, nen, nsd
use timestepmodule,  only: nstep
use postProcessData, only: GMT, connout, deffac, deformed, setria
use monitorVars,     only: elmon, npmon
use plotControl,     only: itemp, nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "te" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
logical tecton
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk
integer :: i

!-init
file = ' '
tecton = .false.
i = 2

if (iecho.eq.8) then
    call ParallelLog("xtemp", "Entered")
endif


10   format(1x,'Usage: te[mperature] [tecton] file')

100     if (i.le.ncmd) then
    if (cmd(i).eq.'tecton') then
        tecton = .true.
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
    endif
    i = i + 1
    goto 100
endif

if ((.not.npmon) .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'temp: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (ITMODE.ne.0 .and. MODE.ge.3) then
        write(stderr,30)
30           format(1x,'temperatures: write temperature INPUT data')
    else
        if (MODE.ne.5 .and. MODE.ne.6) then
            write(stderr,40)
40               format(1x,'temp: no temperatures calculated')
            ierr = EXEC
            goto 1000
        else
            write(stderr,50)
50               format(1x,'temp: first read FE output file')
            ierr = EXEC
            goto 1000
        endif
    endif
else
    if (MODE.ne.5 .and. MODE.ne.6) then
        write(stderr,40)
        ierr = EXEC
        goto 1000
    endif
    if (.not.Tout) then
        write(stderr,60) NSTEP
60           format(1x,'temp: no thermal output at itime = ',I5)
        ierr = EXEC
        goto 1000
    endif
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,65)
65       format(1x,'temp: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif
!
if (tecton) then
    call OfOpen(ierr,file,.true.,npmon,append,.false.,lufo, &
     luMon)
    if (ierr.ne.OK) goto 1000
    call tecout (modeldatactx%T,lufo,luMon)
    goto 1000
endif

call OfOpen(ierr,file,ascout,npmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000

setria = .true.



if (itemp.eq.0) then

    if (iecho.eq.8) then
        call ParallelLog("xtemp", "Reading temperature type 0")
    endif

    call cntemp (meshdatactx%IEN, &
                 meshdatactx%X,  &
                 modeldatactx%D, &
     plotdatactx%XPT,modeldatactx%LMF, &
     modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%T,modeldatactx%TLST,modeldatactx%LMTX,modeldatactx%TX, &
     modeldatactx%TLSTX,modeldatactx%SKEW,NEN,NSD,NDOF,NUMTR,deformed,plotDOFperElem, &
     NPLT,verbose,0,luMon,deffac)

else if (itemp.eq.1) then

    if (iecho.eq.8) then
        call ParallelLog("xtemp", "Reading temperature type 1")
    endif

    call cntemp (meshdatactx%IEN, &
                 meshdatactx%X,  &
                 modeldatactx%D, &
     plotdatactx%XPT,modeldatactx%LMF, &
     modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%T,modeldatactx%TLST,modeldatactx%LMTX,modeldatactx%TX, &
     modeldatactx%TLSTX,modeldatactx%SKEW,NEN,NSD,NDOF,NUMTR,deformed,plotDOFperElem, &
     NPLT,verbose,1,luMon,deffac)

else if (itemp.eq.2) then

    if (iecho.eq.8) then
        call ParallelLog("xtemp", "Reading temperature type 2")
    endif

    call cntemp (meshdatactx%IEN, &
                 meshdatactx%X,  &
                 modeldatactx%D, &
     plotdatactx%XPT,modeldatactx%LMF, &
     modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%T,modeldatactx%T0,modeldatactx%LMTX,modeldatactx%TX, &
     modeldatactx%TLSTX,modeldatactx%SKEW,NEN,NSD,NDOF,NUMTR,deformed,plotDOFperElem, &
     NPLT,verbose,2,luMon,deffac)

else
    write(stderr,70)
70       format(1x,'xte: fatal itemp error')
    call exitp(1)
endif


if (iecho.eq.8) then
    call ParallelLog("xtemp", "Finished reading temperature")
endif


if (iecho.eq.8) then
    call ParallelLog("xtemp", "writing and testing triangles")
endif

if (.not.npmon) then

    call writria (plotdatactx%XPT,NUMTR,lufo,6,ierr)

    if (ierr.ne.OK) then
        goto 1000
    endif

    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR,'trel.dat')
    endif

endif

if (iecho.eq.8) then
    call ParallelLog("xtemp", "writing and testing triangles done")
endif



call OfClose (ierr,file,npmon,lufo,luMon)

1000 if (iecho.eq.8) then
    call ParallelLog("xtemp", "Done")
endif
   return
end
!-------------------------------------------------------------------------------
subroutine xmo (cmd,ncmd,ierr)

  ! Routine that X-ecutes the "monitor" command.

  use modeldefinition, only: NUMEL
  use meshdatamodule,  only: meshdatactx
  use modeltopology,   only: NEN
  use iomodule
  use errors,          only: exec, syntax

  implicit none

  !-pass
  integer ncmd,ierr
  character cmd*(*)
  dimension cmd(*)


  if (ncmd.lt.3) then
      write(stderr,10)
10    format(1x,'Usage: mo[nitor] flag number[s]'/ &
       5X,'with flag one of: n[ode], e[lement], s[side]')
      ierr = SYNTAX
      goto 1000
  endif

  if (ifein.eq.0) then
      write(stderr,20)
20    format(1x,'monitor: read FE input file first')
      ierr = EXEC
      goto 1000
  endif

  call setmon (cmd,ncmd,meshdatactx%ien,meshdatactx%neglobal,meshdatactx%nvglobal,ierr)

1000 return

end subroutine xmo
!-------------------------------------------------------------------------------
subroutine xnomo (cmd,ncmd,ierr)

  use iomodule,        only: stderr
  use monitorVars,     only: elmon, nnel, npmon, nnp, &
                             dealloc_node_monitor_data,dealloc_elem_monitor_data

  use errors,          only: syntax

  ! Routine that X-ecutes the "nomo" command.

  implicit none
  !-pass
  integer ncmd,ierr
  character cmd*(*)
  dimension cmd(*)
  
  if (ncmd.lt.2) then
      write(stderr,10)
10    format(1x,'Usage: nomo[nitor] flag'/5X,'with flag one of: n[ode], e[lement]')
      ierr = SYNTAX
      return
  endif
  
  if (cmd(2)(1:1).eq.'e') then
      elmon = .false.
      nnel  = 0
      call dealloc_elem_monitor_data
  elseif (cmd(2)(1:1).eq.'n') then
      npmon = .false.
      nnp   = 0
      call dealloc_node_monitor_data
  else
      write(stderr,10)
      ierr = 10
  endif
  
  return

end subroutine xnomo
!-------------------------------------------------------------------------------
subroutine xeflux (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: GMT, connout, &
                           deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "eflux" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)




!-locl
#include "lenchar.i"

 character(len=LENCHAR) file         
 integer, external :: lnblk
!-init
 ierr = OK
 file = ' '

 if ((.not.elmon) .and. ncmd.lt.2) then
write(stderr,10)
10        format(1x,'Usage: efl[ux] file')
ierr = SYNTAX
goto 1000
 endif

 if (ncmd.ge.2) file = cmd(2)

 if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'eflux: read FE input and output files first')
ierr = EXEC
goto 1000
 endif

 if (MODE.ne.5 .and. MODE.ne.6) then
write(stderr,30)
30        format(1x,'eflux: no thermal calculations performed')
ierr = EXEC
goto 1000
 endif

 if (ifeout.eq.0) then
write(stderr,40)
40        format(1x,'eflux: first read FE output file')
ierr = EXEC
goto 1000
 endif

 if (.not.Tout) then
write(stderr,50) NSTEP
50       format(1x,'eflux: no heat flux output at itime = ', &
     I5)
ierr = EXEC
goto 1000
 endif
!
 if (GMT .and. .not.connout .and. .not.elmon) then
write(stderr,55)
55        format(1x,'eflux: first write connectivity file')
ierr = SYNTAX
goto 1000
 endif

 if (triang) then
write(stderr,'("eflux: triangulation not yet available")')
ierr = EXEC
goto 1000
 endif

 if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
           plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
setria = .false.


if (NUMTR.le.0) then
    write(stderr,60)
60        format(1x,'eflux: no contour triangles generated')
    goto 1000
endif

if (iecho.ge.4) then
    call tstria (plotdatactx%XPT,NPLT,meshdatactx%neglobal, 'trel.dat')
endif

 endif


 call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000

 call efflux (meshdatactx%IEN, &
         modeldatactx%Thermalflux, &
         plotdatactx%IJSUR,plotdatactx%XPT,verbose, &
   luMon,meshdatactx%neglobal,CST)


 if (.not.elmon) then
call writria (plotdatactx%XPT,NUMTR,lufo, 7,ierr)
 endif

 if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
 end

!-------------------------------------------------------------------------------
subroutine xflux (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: deffac,deformed
use monitorVars,     only: elmon
use plotControl    ! only: a lot. todo: put dircetly in fluxar
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "flux" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
 integer, external :: lnblk
!-init
 ierr = OK
 file = ' '

 if ((.not.elmon) .and. ncmd.lt.2) then
write(stderr,10)
10        format(1x,'Usage: fl[ux] file')
ierr = SYNTAX
goto 1000
 endif

 if (ncmd.ge.2) file=cmd(2)

 if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'flux: read FE input and output files first')
ierr = EXEC
goto 1000
 endif

 if (ifeout.eq.0) then
if (MODE.ne.5 .and. MODE.ne.6) then
write(stderr,30)
30        format(1x,'flux: no heat fluxes calculated')
else
write(stderr,40)
40        format(1x,'flux: first read FE output file')
endif
ierr = EXEC
goto 1000
 endif

 if (MODE.ne.5 .and. MODE.ne.6) then
write(stderr,30)
ierr = EXEC
goto 1000
 endif

 if (.not.Tout) then
write(stderr,50) NSTEP
50       format(1x,'flux: no heat flux output at itime = ', &
     I5)
ierr = EXEC
goto 1000
 endif

 call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000
!
!    TODO check, is IPD/IPIEN/IPX sequence correct? should not be IPIEN/IPX/IPD/
!    like everywhere else? 

 call fluxar (modeldatactx%D,  &
         meshdatactx%IEN, &
         meshdatactx%X,   &
         modeldatactx%Thermalflux, &
         modeldatactx%LMF, &
   modeldatactx%tfault,modeldatactx%DXE, &
   modeldatactx%SKEW,deffac,hatsize,hatangle,Fsb,Fscale,ierr,luMon,NDOF, &
   NEN,NSD,ielbas,ielstp,iarcont,lufo,deformed,verbose)
 if (ierr.ne.OK) goto 1000
!
call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end
!-------------------------------------------------------------------------------
subroutine xedisp (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use modeltopology,   only: NDOF
use iomodule
use documentation,   only: plotVariety
use constants,       only: isPlnplt, isAxiplt, isOpnplt, &
                           disp, velo, &
                           normal, total, slip, split, &
                           coordx, coordy, coordz, coordr, effective
use timestepmodule,  only: outputcontroldatactx, timestepdatactx, nstep
use postProcessData, only: GMT, connout, deffac, deformed, setria
use monitorVars,     only: npmon, nnp, monnp, elmon
use plotControl,     only: nplt, printThisDisp
use errors,          only: exec, ok, syntax
                      

! Routine that executes the "edisp" command. (or edi)

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
integer :: i
logical :: plotX, plotY
double precision :: totalDisp
integer :: nodeID
!-init

plotX = .false.
plotY = .false.

ierr = OK

if ((.not.npmon) .and. ncmd.lt.2) then
    call xedispUsage(1)

    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'edispl: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30       format(1x,'edispl: no mechanical calculations performed')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    write(stderr,40)
40       format(1x,'edispl: first read FE output file')
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'edispl: no displacement output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,55)
55       format(1x,'edispl: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif

i    = 2
file = ' '

! set defaults, possible overwriting setting from a   
! previous call to xvelo, xdispl, xevel, xedisp
printThisDisp%dispOrVelo = disp
printThisDisp%dispType   = total
printThisDisp%coord      = effective

100     if (i.le.ncmd) then
    if (cmd(i).eq.'total') then
        printThisDisp%dispType = total
    elseif (cmd(i).eq.'normal') then
           printThisDisp%dispType = normal
    elseif (cmd(i).eq.'slip') then
        printThisDisp%dispType = slip
    elseif (cmd(i).eq.'split') then
        printThisDisp%dispType = split
    elseif (cmd(i).eq.'e') then
        printThisDisp%coord = effective
    elseif (cmd(i).eq.'x' .and. plotVariety.eq.isPlnplt) then
        printThisDisp%coord = coordx
        plotX = .true.
    elseif (cmd(i).eq.'r' .and. plotVariety.eq.isAxiplt) then
        printThisDisp%coord = coordr
        plotX = .true.
    elseif (cmd(i).eq.'y' .and. plotVariety.eq.isPlnplt) then
        printThisDisp%coord = coordy
        plotY = .true.
    elseif (cmd(i).eq.'z' .and. plotVariety.eq.isAxiplt) then
        printThisDisp%coord = coordz
        plotY = .true.
    elseif (cmd(i).eq.'z') then
        printThisDisp%coord = coordz

    else
        if (lnblk(file).gt.0) then
            call xedispUsage(2)

            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (.not.npmon .and. lnblk(file).le.0) then
    call xedispUsage(3)

    ierr = SYNTAX
    goto 1000
endif
!
!        write(*,*) 'opening file: ', trim(file)

call OfOpen(ierr,file,ascout,npmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) then
    goto 1000
endif
!
setria = .true.


! cntdis takes care of ptinting the displacement,
! taking slippery and faulted nodes into account.
call cntdis (meshdatactx%IEN, &
             meshdatactx%X,   &
            modeldatactx%D,   &
            modeldatactx%DELD,&
 modeldatactx%NSLIP, &
 modeldatactx%DXE,modeldatactx%DX,modeldatactx%nfault,modeldatactx%LMF,&
 modeldatactx%tfault,modeldatactx%dfault, &
 modeldatactx%SKEW,outputcontroldatactx%IMPRINT, &
 outputcontroldatactx%ITPRINT, timestepdatactx%MAXSTP, &
 timestepdatactx%DELT,plotdatactx%XPT,NDOF,NUMTR,deformed, &
 verbose,luMon,deffac)





if (.not. npmon) then
    ! plot the values of plotdatactx%XPT to file
    call writria (plotdatactx%XPT,NUMTR,lufo,8,ierr)
   
    if (ierr.ne.OK) then 
        goto 1000
    endif

    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR,'trel.dat')
    endif

endif
call OfClose (ierr,file,npmon,lufo,luMon)

1000    return
end
!-------------------------------------------------------------------------------
subroutine xedispUsage(tag)

use documentation,       only: plotVariety
use constants,           only: isPlnplt, isAxiplt, isOpnplt
use iomodule

implicit none


integer :: tag


10   format(1x,'Usage: edi[sp] file [total|normal|slip|split] [e|x|y|z]')
11   format(1x,'Usage: edi[sp] file [total|normal|slip|split] [e|r|z]')
12   format(1x,'Usage: edi commands of opn not yet specified')


if      (plotVariety.eq.isPlnplt) then
    write(stderr,10)
else if (plotVariety.eq.isAxiplt) then
    write(stderr,11)
else if (plotVariety.eq.isOpnplt) then
    write(stderr,12)
else
    write(stderr,*) "xedisp does not recognize plotVariety: ", plotVariety
    write(stderr,*) "This should not happen. Please contact model support."
endif

write(stderr,*) "debug-tag: ", tag


end subroutine

!-------------------------------------------------------------------------------



subroutine xevel (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use modeltopology,   only: NDOF
use iomodule
use documentation,   only: plotVariety
use constants,       only: isPlnplt, isAxiplt, isOpnplt, &
                           disp, velo, &
                                    normal, total, slip,    split, &
                                    coordx, coordy, coordz, coordr, effective

use timestepmodule,  only: timestepdatactx, outputcontroldatactx, nstep
use postProcessData, only: GMT, connout, &
                           deffac,setria,deformed
use monitorVars,     only: elmon, npmon
use plotControl,     only: nplt, printThisDisp
use errors,          only: exec, ok, syntax


!
! Routine that X-ecutes the "evel" command. (Effective velocity)
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
integer :: i
!-init
ierr = OK
!
10   format(1x,'Usage: eve[locity] file [total|normal|slip|split] [e|x|y]')
11   format(1x,'Usage: eve[locity] file [total|normal|slip|split] [e|r|z]')
12   format(1x,'Usage: edi commands of opn not yet specified')

if ((.not.npmon) .and. ncmd.lt.2) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'evel: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif
!
if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30       format(1x,'evel: no mechanical calculations performed')
    ierr = EXEC
    goto 1000
endif
!
if (ifeout.eq.0) then
    write(stderr,40)
40       format(1x,'evel: first read FE output file')
    ierr = EXEC
    goto 1000
endif
!
if (.not.Mout .or. NSTEP.le.0) then
    write(stderr,50) NSTEP
50       format(1x,'evel: no velocity output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,55)
55       format(1x,'evel: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif

i    = 2
file = ' '

! set defaults, possible overwriting setting from a   
! previous call to xvelo, xdispl, xevel, xedisp
printThisDisp%dispOrVelo = velo
printThisDisp%dispType   = total
printThisDisp%coord      = effective



100     if (i.le.ncmd) then
    if (cmd(i).eq.'total') then
        printThisDisp%dispType = total
    elseif (cmd(i).eq.'normal') then
        printThisDisp%dispType = normal
    elseif (cmd(i).eq.'slip') then
        printThisDisp%dispType = slip
    elseif (cmd(i).eq.'split') then
        printThisDisp%dispType = split
    elseif (cmd(i).eq.'e') then
        printThisDisp%dispType = total
    elseif (cmd(i).eq.'x' .and. plotVariety.eq.isPlnplt) then
        printThisDisp%coord = coordx
    elseif (cmd(i).eq.'r' .and. plotVariety.eq.isAxiplt) then
           printThisDisp%coord = coordr

    elseif (cmd(i).eq.'y' .and. plotVariety.eq.isPlnplt) then
           printThisDisp%coord = coordy
    elseif (cmd(i).eq.'z' .and. plotVariety.eq.isAxiplt) then
           printThisDisp%coord = coordz

    elseif (cmd(i).eq.'z') then
           printThisDisp%coord = coordz


    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (.not.npmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
!
call OfOpen(ierr,file,ascout,npmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000
!
setria = .true.
!

call cntdis (meshdatactx%IEN, &
             meshdatactx%X,   &
            modeldatactx%D,   &
            modeldatactx%DELD,&
 modeldatactx%NSLIP, &
 modeldatactx%DXE,modeldatactx%DX,modeldatactx%nfault,modeldatactx%LMF,&
 modeldatactx%tfault,modeldatactx%dfault, &
 modeldatactx%SKEW,outputcontroldatactx%IMPRINT, &
 outputcontroldatactx%ITPRINT,timestepdatactx%MAXSTP,&
 timestepdatactx%DELT,plotdatactx%XPT,NDOF,NUMTR,deformed, &
 verbose,luMon,deffac)

if (.not.npmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo,9,ierr)
    if (ierr.ne.OK) then
        goto 1000
    endif
    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif
call OfClose (ierr,file,npmon,lufo,luMon)
!
1000    return
end
!-------------------------------------------------------------------------------
subroutine xdispl (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule
use timestepmodule,  only: outputcontroldatactx, &
                                timestepdatactx, &
                           nstep
use postProcessData, only: deffac,deformed
use monitorVars,     only: npmon
use plotControl ! ,      only: a lot, todo do from dispar
use errors,          only: exec, ok, syntax
use constants,       only: total, normal, slip, split, &
                           disp, effective

! Routine that X-ecutes the "displ" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
integer :: dumpfile
integer :: i
!-init
ierr = OK


!        meshdatactx%x contains the initial positions of the vertices
!        these have to be compared to positions in de fedsk file, to
!        compute the displacement.        

if ((.not.npmon) .and. ncmd.lt.2) then
    write(stderr,10)
10       format(1x,'Usage: di[splacements] [normal|slip|split|', &
     'total] file')
    ierr = SYNTAX
    goto 1000
endif

!       if (ncmd.ge.2) file=cmd(2)

if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'displ: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'displ: no displacements calculated')
    else
        write(stderr,40)
40           format(1x,'displ: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'displ: no displacements output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif

dumpfile = 0
i    = 2
file = ' '

! set defaults, possible overwriting setting from a   
! previous call to xvelo, xdispl, xevel, xedisp
printThisDisp%dispOrVelo = disp
printThisDisp%dispType   = total
printThisDisp%coord      = effective


100     if (i.le.ncmd) then

    if (cmd(i).eq.'total') then
        printThisDisp%dispType = total
    elseif (cmd(i).eq.'normal') then
           printThisDisp%dispType = normal
    elseif (cmd(i).eq.'slip') then
           printThisDisp%dispType = slip
    elseif (cmd(i).eq.'split') then
           printThisDisp%dispType = split
    elseif (cmd(i).eq.'e' .or. \
            cmd(i).eq.'x' .or. \
            cmd(i).eq.'y' .or. \
            cmd(i).eq.'z' .or. \
            cmd(i).eq.'r') then
        write(stderr,60) cmd(i)(1:1)
60      format(1x,'displ: with the "',A1,'" argument you ', &
         'probably intended to select a component.'/1x, &
         'You should use the "edisp" command for this')
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (.not.npmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

call OfOpen (ierr,file,ascout,npmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) then
      goto 1000
endif


call dispar ( meshdatactx%X,   &
                 modeldatactx%D,   &
                 modeldatactx%DELD,&
                 modeldatactx%NSLIP,&
     modeldatactx%DXE, &
     modeldatactx%DX, &
     modeldatactx%nfault,modeldatactx%LMF,modeldatactx%tfault,modeldatactx%dfault,&
             meshdatactx%IEN, &
     modeldatactx%SKEW, &
     outputcontroldatactx%IMPRINT, outputcontroldatactx%ITPRINT, &
     timestepdatactx%MAXSTP,timestepdatactx%DELT,deffac,hatsize, &
     hatangle,Dsb,Dscale,ierr,LuMon,inpbas,inpstp,iarcont, &
     lufo,deformed,verbose)



if (ierr.ne.OK) then
    goto 1000
endif

call OfClose (ierr,file,npmon,lufo,luMon)

1000    return
end
!-------------------------------------------------------------------------------
subroutine xvelo (cmd,ncmd,ierr)
!
USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule
use timestepmodule,  only: outputcontroldatactx, &
                                timestepdatactx, &
                           nstep
use postProcessData, only: deffac,deformed
use monitorVars,     only: npmon
use plotControl  ! ,a lot. todo : do from dispar directly
use errors,          only: exec, ok, syntax
use constants,       only: disp, velo, &
                           normal, total, slip, split, &
                           coordx, coordy, coordz, coordr, effective


! Routine that X-ecutes the "velocity" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
integer :: i
!-init
ierr = OK
!
if ((.not.npmon) .and. ncmd.lt.2) then
    write(stderr,10)
10       format(1x,'Usage: ve[locity] [normal|slip|split|', &
     'total] [e|x|y|z] file')
    ierr = SYNTAX
    goto 1000
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'velocity: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'velocity: no displacements calculated')
    else
        write(stderr,40)
40           format(1x,'velocity: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout .or. NSTEP.le.0) then
    write(stderr,50) NSTEP
50       format(1x,'velocity: no velocities output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
i    = 2
file = ' '

! set defaults, possible overwriting setting from a   
! previous call to xvelo, xdispl, xevel, xedisp
printThisDisp%dispOrVelo = velo
printThisDisp%dispType   = total
printThisDisp%coord      = effective



100     if (i.le.ncmd) then
    if (cmd(i).eq.'total') then
        printThisDisp%dispType = total
    elseif (cmd(i).eq.'normal') then
           printThisDisp%dispType = normal
    elseif (cmd(i).eq.'slip') then
        printThisDisp%dispType = slip
    elseif (cmd(i).eq.'split') then
        printThisDisp%dispType = split 
    elseif (cmd(i).eq.'e' .or. \
            cmd(i).eq.'x' .or. \
            cmd(i).eq.'y' .or. \
            cmd(i).eq.'z' .or. \
            cmd(i).eq.'r') then

        write(stderr,60) cmd(i)(1:1)
60      format(1x,'velo: with the "',A1,'" argument you ', &
         'probably intended to select a component.'/1x, &
         'You should use the "evelo" command for this')

    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (.not.npmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
!
call OfOpen (ierr,file,ascout,npmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) goto 1000
!

call dispar ( meshdatactx%X,  &
             modeldatactx%D,  &
             modeldatactx%DELD,&
 modeldatactx%NSLIP,modeldatactx%DXE, &
 modeldatactx%DX, &
 modeldatactx%nfault,modeldatactx%LMF,modeldatactx%tfault,modeldatactx%dfault, &
              meshdatactx%IEN, &
 modeldatactx%SKEW, &
 outputcontroldatactx%IMPRINT, outputcontroldatactx%ITPRINT, &
 timestepdatactx%MAXSTP,timestepdatactx%DELT,deffac,hatsize, &
 hatangle,Vsb,Velsc,ierr,LuMon,inpbas,inpstp,iarcont,lufo, &
 deformed,verbose)
!
! 


call OfClose (ierr,file,npmon,lufo,luMon)
!
1000    return
end
!-------------------------------------------------------------------------------
subroutine xvi (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR, icvis, ivlim, nsur
USE MODELDATAMODULE
USE MESHDATAMODULE
use algebra,         only: idig
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: GMT, connout, &
                           deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl ! a lot todo
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "viscosity" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
integer :: i
!-init
ierr = OK

10   format(1x, &
 'Usage: vi[scosity] file [[no]log] [no[nu[mbers]]]',1x, &
 '[no[li[near]]]')
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'viscosity: read FE input file first')
    ierr = EXEC
    goto 1000
endif

if (ICVIS.ne.0) then
    if (ifeout.eq.0) then
        write(stderr,30)
30            format(1x,'WARNING: writing viscosity input data')
    else
        if (MODE.ne.2 .and. MODE.ne.6) then
            write(stderr,40)
40            format(1x,'viscosity: no viscosities calculated')
            ierr = EXEC
            goto 1000
        else if (.not.Mout) then
            write(stderr,50) NSTEP
50               format(1x,'viscosity: no viscosity output at itime', &
             ' = ',I5)
            ierr = EXEC
            goto 1000
        endif
    endif
endif
!
i    = 2
file = ' '

100    if (ncmd.ge.i) then
    if (cmd(i)(1:2).eq.'lo') then
        vislog = .true.
    elseif (cmd(i)(1:4).eq.'nolo') then
        vislog = .false.
    else if (cmd(i)(1:2).eq.'nu') then
        viscnu = .true.
    else if (cmd(i)(1:4).eq.'nonu') then
        viscnu = .false.
    else if (cmd(i)(1:2).eq.'li') then
        viscli = .true.
    else if (cmd(i)(1:4).eq.'noli') then
        viscli = .false.
    else if (cmd(i)(1:2).eq.'ra') then
        visraw = .true.
    else if (cmd(i)(1:4).eq.'nora') then
        visraw = .false.
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif
 
if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (viscli .and. NSTEP.le.0) then
    write(stderr,51) NSTEP
51       format(1x,'viscosity: no linearized viscosity at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
if (GMT .and. .not.connout .and. .not.elmon) then
    write(stderr,55)
55        format(1x,'viscosity: first write connectivity file')
    ierr = SYNTAX
    goto 1000
endif
!
if (triang) then
!        invoke "triangulate" to generate smoothly interpolated triangles
    call trian (meshdatactx%IEN, &
                meshdatactx%X,   &
           modeldatactx%D,  &
           modeldatactx%STN, &
           modeldatactx%STR, &
#ifdef EXPERIMENTAL_NORATESTORAGE
     modeldatactx%LMF,modeldatactx%LMX, &
#else
     modeldatactx%rate,modeldatactx%LMF,modeldatactx%LMX, &
#endif
       modeldatactx%tfault,meshdatactx%mat,matdatactx%prpmat, &
       matdatactx%PRPLAS,modeldatactx%DXE,modeldatactx%SKEW,5,deformed, &
       deffac,file)
    goto 1000
endif
!
if (setria .and. .not.elmon .and. .not.viscnu) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
    setria = .false.

    if (NUMTR.le.0) then
        write(stderr,60)
60        format(1x,'viscosity: no contour triangles generated')
        goto 1000
    endif

    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif
!
call OfOpen (ierr,file,ascout,viscnu.or.elmon,append,triang, &
             lufo,luMon)

if (ierr.ne.OK) then
    goto 1000
endif
!
if (visraw .and. IVLIM.ne.0 .and. IDIG(ICVIS,4,1).ne.0) then
!        Original viscosities (before viscosity limits were applied)
    call setvisc (meshdatactx%IEN, &
                  meshdatactx%X,   &
                 modeldatactx%D,   &
     modeldatactx%T,modeldatactx%TX, &
     modeldatactx%STN,modeldatactx%STR, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%LMTX,modeldatactx%SKEW,NSD,NSTR,NDOF, &
       NEN)
endif

if (visraw .and. IVLIM.ne.0 .and. IDIG(ICVIS,4,3).ne.0) then
    call vweak(modeldatactx%STR,NSTR)
endif
!
if (viscnu.or.elmon) then
    call nummat (meshdatactx%IEN, &
             meshdatactx%X,   &
            modeldatactx%D,  &
     meshdatactx%MAT, matdatactx%prpmat, &
       modeldatactx%STN,modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
       modeldatactx%SKEW,NEN,NSD,NDOF,NSTR,deformed, &
       deffac,lufo,ierr,luMon,verbose)

    if (ierr.eq.OK) then
        call OfClose (ierr,file,elmon,lufo,luMon)
    endif
    goto 1000
endif
!
call comat (meshdatactx%mat,plotdatactx%IJSUR,plotdatactx%XPT,meshdatactx%IEN,modeldatactx%STN, &
            matdatactx%prpmat,NSUR,NPLDOF,NPLT,NEN,NSTR,CST)

if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 10,ierr)
endif

if (ierr.eq.OK) then
    call OfClose (ierr,file,elmon,lufo,luMon)
endif
!
1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xma (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR, icvis, nSur
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use algebra,         only: idig
use debugmodule,     only: iecho
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule
use postProcessData, only: GMT, connout, &
                           deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: npldof, nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "materials" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
#include "lenchar.i"
 character(len=LENCHAR) file
 integer, external :: lnblk
integer :: i, iflag
!-init
 ierr = OK

 if (ifein.eq.0) then
write(stderr,10)
10        format(1x,'materials: read FE input file first')
ierr = EXEC
goto 1000
 endif

 i     = 2
 file  = ' '
 iflag = 0
100   if (ncmd.ge.i) then
! 0 MAT
! 1-6 PRPMAT
! 7-12 PRPCRP
! 13-21 PRPLAS
! 22-27 PRPTEM
if (cmd(i)(1:2).eq.'nu') then
    iflag = 0
elseif (cmd(i).eq.'E') then
    iflag = 1
else if (cmd(i).eq.'v') then
    iflag = 2
else if (cmd(i)(1:2).eq.'vi') then
    iflag = 3
else if (cmd(i)(1:2).eq.'po') then
    iflag = 4
else if (cmd(i)(1:2).eq.'de') then
    iflag = 5
else if (cmd(i)(1:2).eq.'th') then
    iflag = 6
else if (cmd(i).eq.'Qpl') then
    if (IDIG(ICVIS,4,1).ne.0) then
        iflag = 6+1
    else
        write(stderr,20) 'Qpl',ICVIS
20            format(1x,'mat: ',I3,' not available for ICVIS=',I1)
       ierr = SYNTAX
       return
    endif
else if (cmd(i).eq.'Apl') then
    if (IDIG(ICVIS,4,1).ne.0) then
        iflag = 6+2
    else
        write(stderr,20) 'Apl',ICVIS
        ierr = SYNTAX
        return
    endif
else if (cmd(i).eq.'npl') then
    if (IDIG(ICVIS,4,1).ne.0) then
        iflag = 6+3
    else
        write(stderr,20) 'npl',ICVIS
        ierr = SYNTAX
        return
    endif
else if (cmd(i).eq.'Vpl') then
    if (IDIG(ICVIS,4,1).ne.0) then
        iflag = 6+6
    else
        write(stderr,20) 'Vpl',ICVIS
        ierr = SYNTAX
        return
    endif
else if (cmd(i)(1:2).eq.'he') then
    if (MODE.ge.3) then
        iflag = 21+1
    else
        write(stderr,30) 'heat ',MODE
30            format(1x,'mat: ',I5,' not available for MODE=',I1)
        ierr = SYNTAX
    return
    endif
else if (cmd(i).eq.'Cp') then
    if (MODE.ge.3) then
        iflag = 21+2
    else
        write(stderr,30) 'Cp   ',MODE
        ierr = SYNTAX
        return
    endif
else if (cmd(i)(1:2).eq.'al') then
    if (MODE.ge.3) then
        iflag = 21+3
    else
        write(stderr,30) 'alpha',MODE
        ierr = SYNTAX
        return
    endif
else if (cmd(i).eq.'Kx') then
if (MODE.ge.3) then
    iflag = 21+4
else
    write(stderr,30) 'Kx   ',MODE
       ierr = SYNTAX
    return
endif
else if (cmd(i).eq.'Ky') then
if (MODE.ge.3) then
    iflag = 21+5
else
    write(stderr,30) 'Ky   ',MODE
       ierr = SYNTAX
    return
endif
else if (cmd(i)(1:2).eq.'pl') then
if (IDIG(ICVIS,4,2).ne.0) then
    iflag = 12+1
else
    write(stderr,20) 'pl ',ICVIS
       ierr = SYNTAX
    return
endif
else if (cmd(i)(1:2).eq.'fr') then
if (IDIG(ICVIS,4,2).ne.0) then
    iflag = 12+2
else
    write(stderr,20) 'fr ',ICVIS
       ierr = SYNTAX
    return
endif
else if (cmd(i)(1:2).eq.'un') then
if (IDIG(ICVIS,4,2).ne.0) then
    iflag = 12+3
else
    write(stderr,20) 'un ',ICVIS
       ierr = SYNTAX
    return
endif
else if (cmd(i)(1:2).eq.'fl') then
if (IDIG(ICVIS,4,2).ne.0) then
    iflag = 12+4
else
    write(stderr,20) 'fl ',ICVIS
       ierr = SYNTAX
    return
endif
else if (cmd(i)(1:3).eq.'yld') then
if (IDIG(ICVIS,4,2).ne.0) then
    iflag = 12+5
else
    write(stderr,20) 'yld',ICVIS
       ierr = SYNTAX
    return
endif
else if (cmd(i).eq.'nflow') then
if (IDIG(ICVIS,4,2).ne.0) then
    iflag = 12+6
else
    write(stderr,20) 'nfl',ICVIS
       ierr = SYNTAX
    return
endif
else if (cmd(i)(1:2).eq.'ha') then
if (IDIG(ICVIS,4,2).ne.0) then
    iflag = 12+7
else
    write(stderr,20) 'ha ',ICVIS
       ierr = SYNTAX
    return
endif
else
    if (lnblk(file).gt.0) then
        write(stderr,40)
40            format(1x,'Usage: ma[terial] [flag] file')
        ierr = SYNTAX
        goto 1000
    endif

file = cmd(i)
endif

i = i + 1
goto 100
 endif
!
 if (GMT .and. .not.connout .and. .not.elmon) then
write(stderr,55)
55        format(1x,'material: first write connectivity file')
ierr = SYNTAX
goto 1000
 endif
!
 if (triang) then
!        invoke "triangulate" to generate smoothly interpolated triangles
call trian (meshdatactx%IEN, &
            meshdatactx%X,   &
           modeldatactx%D,   &
           modeldatactx%STN, &
           modeldatactx%STR, &
#ifdef EXPERIMENTAL_NORATESTORAGE
     modeldatactx%LMF,modeldatactx%LMX, &
#else
     modeldatactx%rate,modeldatactx%LMF,modeldatactx%LMX, &
#endif       
       modeldatactx%tfault, &
       meshdatactx%mat, matdatactx%PRPMAT, & !meshdatactx%mat,matdatactx%prpmat, &
       matdatactx%PRPLAS,modeldatactx%DXE,modeldatactx%SKEW,5,deformed, &
       deffac,file)
goto 1000
 endif
!
 if (setria .and. .not.elmon)then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                 modeldatactx%D,  &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
  setria = .false.

if (NUMTR.le.0) then
    write(stderr,60)
60        format(1x,'viscosity: no contour triangles generated')
    goto 1000
endif

if (iecho.ge.4) then
    call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
endif

 endif
!
 call OfOpen (ierr,file,ascout,elmon,append,triang,lufo,luMon)
 if (ierr.ne.OK) then
goto 1000
 endif

  


 if (elmon) then
   call matplt(meshdatactx%mat,&
               matdatactx%prpmat,&
               matdatactx%prpcrp,&
               matdatactx%PRPLAS,&
               matdatactx%prptem, &
               luMon,iflag,verbose)
 else
   if (iflag.ge.0.and.iflag.le.6) then
! general properties
call cmat (plotdatactx%IJSUR,plotdatactx%XPT, & 
           matdatactx%PRPMAT, & !matdatactx%prpmat
           6, &
           iflag,NSUR,NPLDOF,NPLT,NEN,NSTR,CST)

   else if (iflag.ge.7.and.iflag.le.12) then
! creep properties
call cmat (plotdatactx%IJSUR,plotdatactx%XPT, &
           matdatactx%PRPCRP,6, &
           iflag-6,NSUR,NPLDOF,NPLT,NEN,NSTR,CST)

   else if (iflag.ge.13.and.iflag.le.21) then

call cmat (plotdatactx%IJSUR,plotdatactx%XPT, &
           matdatactx%PRPLAS,9, &
           iflag-12,NSUR,NPLDOF,NPLT,NEN,NSTR,CST)

   else if (iflag.ge.22.and.iflag.le.27) then

call cmat (plotdatactx%IJSUR,plotdatactx%XPT, &
           matdatactx%PRPTEM,6, &
           iflag-21,NSUR,NPLDOF,NPLT,NEN,NSTR,CST)
   else
 stop 'program error'
   endif
 endif


 if (.not.elmon) then
call writria (plotdatactx%XPT,NUMTR,lufo, 11,ierr)
 endif

 if (ierr.eq.OK) then
call OfClose (ierr,file,elmon,lufo,luMon)
 endif
!
1000    return
 end
!-------------------------------------------------------------------------------
subroutine xnbc (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule
use postProcessData, only: GMT,deffac, deformed
use plotControl,     only: mattsym, readsy, sysc, sysi
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "nbc" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)




!-locl
#include "lenchar.i"
 integer iflg
character(len=LENCHAR) file         
 integer, external :: lnblk
 save iflg
integer :: i
!-init
 data iflg/6/
 ierr = OK
 file = ' '
 i = 2
!
 if (GMT) then
100        if (i.le.ncmd) then
if (cmd(i).eq.'d') then
    iflg = 1
else if (cmd(i).eq.'v') then
    iflg = 2
else if (cmd(i).eq.'f') then
    iflg = 3
else if (cmd(i).eq.'i') then
    iflg = 4
else if (cmd(i).eq.'r') then
    iflg = 5
else if (cmd(i).eq.'a') then
    iflg = 6
else
    if (lnblk(file).le.0) then
        file = cmd(i)
    else
        write(stderr,10)
10            format(1x,'Usage: nbc [a|d|v|f|i|r|s] file')
    ierr = SYNTAX
    return
    endif
endif

i = i + 1
goto 100
endif

if (lnblk(file).le.0) then
    write(stderr,10)
ierr = SYNTAX
return
endif
if (ifein.eq.0) then
    write(stderr,20)
20            format(1x,'nbc: read FE input file first')
    ierr = EXEC
    return
endif
!
call OfOpen (ierr,file,.true.,.false.,append,.false.,lufo, &
       luMon)
if (ierr.ne.OK) return
 else
if (ncmd.lt.2) then
    write(stderr,30)
30            format(1x,'Usage: nbc file')
    ierr = SYNTAX
    return
endif
file = cmd(2)
!
if (ifein.eq.0) then
    write(stderr,20)
    ierr = EXEC
    return
endif
!
call OfOpen (ierr,file,ascout,.false.,append,.false.,lufo, &
       luMon)
if (ierr.ne.OK) return
!
if (readsy) then
    call rdalsy ()
    readsy = .false.
endif
!
iflg = 6
 endif
!

if (mattsym.eq.1 .or. mattsym.eq.2) call BCDEF(mattsym)

!subroutine NBCPLT (IEN,X,D,IBOND,BOND,SKEW,DEFORM,DEFFAC,LU,SY, &
!   SYSC,iflg)


!    nodal boundary condition symbols
 call NBCPLT (meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
        modeldatactx%IBOND, &
        modeldatactx%BOND, &
   modeldatactx%SKEW,deformed,deffac,lufo,sysi,sysc,iflg,mattsym)
!
call OfClose (ierr,file,.false.,lufo,luMon)
!
 return
 end
!-------------------------------------------------------------------------------
subroutine xslip (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule
use postProcessData, only: GMT
use postProcessData, only: deffac,deformed
use plotControl,     only: mattsym, sysi
use errors,          only: exec, ok, syntax

!
! Routine that X-ecutes the "slip" command.
!
 implicit none
!-pass
 integer ncmd,ierr,iflg
 character cmd*(*)
 dimension cmd(*)




!-locl
#include "lenchar.i"
 character(len=LENCHAR) file
 integer, external :: lnblk
!-init
 ierr = OK
 file = ' '
!
 if (ncmd.ge.2) file = cmd(2)
!
 if (lnblk(file).le.0) then
write(stderr,10)
10        format(1x,'Usage: slip[pery] file')
ierr = SYNTAX
goto 1000
 endif
!
if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'xslip: read FE input file first')
ierr = EXEC
return
 endif
!
 if (GMT) then
call OfOpen (ierr,file,.true.,.false.,append,.false.,lufo, &
       luMon)
 else
            stop "non-GMT plotting not implemented"
 endif

 if (ierr.ne.OK) return

!    slippery dofs symbols
 call SLPPLT (meshdatactx%X,   &
        modeldatactx%D,   &
 modeldatactx%SKEW,modeldatactx%IDSLP,modeldatactx%NSLIP, &
   deformed,deffac,lufo,sysi,mattsym)
!
call OfClose (ierr,file,.false.,lufo,luMon)
!
1000    return
 end

!-------------------------------------------------------------------------------
subroutine xebc (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule
use postProcessData, only: GMT, deffac, deformed
use plotControl,     only: mattsym, sysc, sysi
use errors,          only: exec, ok, syntax

!
! Routine that X-ecutes the "ebc" command.
!
 implicit none
!-pass
 integer ncmd,ierr,iflg
 character cmd*(*)
 dimension cmd(*)




!-locl
#include "lenchar.i"
 character(len=LENCHAR) file,errm
 integer, external :: lnblk
integer :: i
!-init
 file = ' '
 iflg = -1
!
 i = 2
 if (GMT) then
errm = 'Usage: ebc w|s file'
 else
    stop "non-GMT plotting not implemented"
endif

!
 do while (i.le.ncmd)
if (cmd(i).eq.'w') then
iflg = 1
else if (cmd(i).eq.'s') then
iflg = 2
else
if (lnblk(file).le.0) then
    file = cmd(i)
    else
    write(stderr,'(1X,19a)') errm(1:lnblk(errm))
    ierr = SYNTAX
    return
    endif
endif
i = i + 1
 enddo
 if (lnblk(file).le.0) then
write(stderr,'(1X,19a)') errm(1:lnblk(errm))
ierr = SYNTAX
return
 endif
 if (GMT .and. iflg.eq.-1) then
write(stderr,'(1X,19a)') errm(1:lnblk(errm))
ierr = SYNTAX
return
 endif
 if (.not.GMT) iflg = -1


if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'ebc: read FE input file first')
ierr = EXEC
return
 endif


 if (GMT) then
call OfOpen (ierr,file,.true.,.false.,append,.false.,lufo, &
       luMon)
 else
    stop "non-GMT plotting not implemented"
endif

 if (ierr.ne.OK) return

if (mattsym.eq.1 .or. mattsym.eq.2) call BCDEF(mattsym)

!write(*,*) "xebc has lufo: ", lufo

 if (iflg.eq.1 .or. iflg.eq.-1) then
!        Winkler pressure symbols

call BCWP (meshdatactx%X,   &
           modeldatactx%D,   &
     modeldatactx%SKEW,deformed,deffac,lufo,sysc, &
       sysi,mattsym)

 endif
 if (iflg.eq.2 .or. iflg.eq.-1) then
!        stress boundary condition symbols
call BCSTRS (meshdatactx%X,       &
             modeldatactx%D,      &
             modeldatactx%SKEW,   &
             meshdatactx%IEN,     &
             modeldatactx%IELSTR, &
             modeldatactx%ISSIDE, &
             modeldatactx%STRS,   &
             modeldatactx%ISTR,   &
             deformed,            &
             deffac,              &
             lufo,                &
             sysc,                &
             sysi,                &
             mattsym)

endif
!
call OfClose (ierr,file,.false.,lufo,luMon)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine xsc (cmd,ncmd,ierr)

use iomodule
use plotControl   ! only: almost everything
use errors,          only: ok, syntax


! Routine that X-ecutes the "scalebar" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)

#include "lenchar.i"



!-locl
 character(len=LENCHAR) file,string,unit
 logical both,numeric
 integer, external :: lnblk, fnblnk
double precision :: xsb, ysb, sb, scale
double precision, external :: chreal
!-init
 ierr = OK
 file = ' '

10    format(1x,'Usage: sc[alebar] flag x y file [h[or]|v[ert]]')
 if (ncmd.lt.5) then
write(stderr,10)
ierr = SYNTAX
goto 1000
 endif

 both = .true.
 if (cmd(2)(1:4).eq.'stre') then
sb    = Ssb
scale = Sscale
unit  = Sunit
 else if (cmd(2)(1:4).eq.'stra') then
sb    = Esb
scale = Escale
unit  = Eunit
 else if (cmd(2)(1:2).eq.'ra') then
sb    = Rsb
scale = Rscale
unit  = Runit
 else if (cmd(2)(1:2).eq.'di') then
sb    = Dsb
scale = Dscale
both  = .false.
unit  = Dunit
 else if (cmd(2)(1:2).eq.'fl') then
sb    = Fsb
scale = Fscale
both  = .false.
unit  = Funit
 else if (cmd(2)(1:2).eq.'ve') then
sb    = Vsb
scale = Velsc
both  = .false.
unit  = Vunit
 else
write(stderr,20)
20       format(1x,'scalebar: flag should be one of :'/1x, &
     '"stre[ss]", "stra[in]", "ra[te]", "di[spl]", "fl[ux]"',1x, &
     '"ve[locity]"')
ierr = SYNTAX
goto 1000
 endif

 if (numeric(cmd(3)).and.numeric(cmd(4))) then
xsb = chreal(cmd(3))
ysb = chreal(cmd(4))
 else
write(stderr,10)
ierr = SYNTAX
goto 1000
 endif

 if (ncmd.gt.5) then
if (cmd(6)(1:1).eq.'h') then
sbhor = .true.
elseif (cmd(6)(1:1).eq.'v') then
sbhor = .false.
else
write(stderr,10)
ierr = SYNTAX
goto 1000
endif
 endif

!    use 5th argument name as a seed for two file names
 file = cmd(5)(fnblnk(cmd(5)):lnblk(cmd(5))) // '.dat'
 call OfOpen (ierr,file,ascout,.false.,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000

 call sclbar (xsb,ysb,sb,scale,hatangle,hatsize,sbhor,iarcont, &
   lufo,ierr,both)
 if (ierr.ne.OK) goto 1000

 call OfClose (ierr,file,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000

 file = cmd(5)(fnblnk(cmd(5)):lnblk(cmd(5))) // '.cf'
 call OfOpen (ierr,file,.true.,.false.,append,.false.,lufo,luMon)

 call wround (string,sb,3,.true.)
 if (sbhor) then
write(lufo,30) xsb,ysb+sb*0.12d0*scale, &
       string(1:lnblk(string)),' ',unit(1:lnblk(unit))
30        format(1x,'pt -l ',1PG12.4,1X,1PG12.4,1X,80A)
 else
write(lufo,30) xsb+sb*0.12d0*scale,ysb, &
       string(1:lnblk(string)),' ',unit(1:lnblk(unit))
 endif
 call OfClose (ierr,file,.false.,lufo,luMon)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xelmc (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use modeltopology,   only: NSD,NDOF,NEN
use iomodule
use postProcessData, only: deffac,deformed
use monitorVars,     only: elmon
use plotControl,     only: ielbas, ielstp
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "elmc" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)




!-locl
#include "lenchar.i"
 character(len=LENCHAR) file
 integer, external :: lnblk
!-init
 file = ' '

 if (.not.elmon .and. ncmd.lt.2) then
write(stderr,10)
10        format(1x,'Usage: elmc[enter] file')
ierr = SYNTAX
goto 1000
 endif
 if (ncmd.gt.1) file = cmd(2)

if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'elmcenter: first read FE-input file')
ierr = EXEC
goto 1000
 endif

 call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000
!
 call centrs (meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
 modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
   modeldatactx%SKEW,NEN,NSD, &
   NDOF,deformed,deffac,lufo,ierr,ielbas,ielstp,verbose,luMon)
!
 if (ierr.ne.OK) goto 1000

 call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xdt (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR, nsur
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed, cst
use monitorVars,     only: elmon
use plotControl,     only: iDtfl, npldof, nplt
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "dt" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
 integer, external :: lnblk
integer :: i
!-init
 ierr = OK
 i    = 2
 file = ' '

10   format(1x,'Usage: dt file',1x, &
 '[acc[uracy]|sta[bility]|osc[illation]]')
100    if (ncmd.ge.i) then
if (cmd(i)(1:3).eq.'sta') then
iDtfl = 0
elseif (cmd(i)(1:3).eq.'osc') then
iDtfl = 1
elseif (cmd(i)(1:3).eq.'acc') then
iDtfl = 2
else
if (lnblk(file).gt.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
file = cmd(i)
endif
i = i + 1
goto 100
 endif
 if (.not.elmon .and. lnblk(file).le.0) then
write(stderr,10)
ierr = SYNTAX
goto 1000
 endif

 if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'dt: read FE input and output files first')
ierr = EXEC
goto 1000
 endif

 if (ifeout.eq.0) then
if (MODE.ne.2 .and. MODE.ne.6) then
write(stderr,30)
30        format(1x,'dt: no mechanical calculations')
else
write(stderr,40)
40        format(1x,'dt: first read FE output file')
endif
ierr = EXEC
goto 1000
 endif

 if (MODE.ne.2 .and. MODE.ne.6) then
write(stderr,30)
ierr = EXEC
goto 1000
 endif

 if (NSTEP.le.0 .or. (.not.Mout)) then
write(stderr,50) NSTEP
50       format(1x,'dt: no mechanical output at itime = ', &
     I5)
ierr = EXEC
goto 1000
 endif

 if (triang) then
write(stderr,'("dt: triangulation not implemented")')
goto 1000
 endif

 if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
setria = .false.
if (NUMTR.le.0) then
write(stderr,60)
60        format(1x,'dt: no contour triangles generated')
goto 1000
endif
if (iecho.ge.4) then
  call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
endif
 endif

 call OfOpen (ierr,file,ascout,elmon,append,triang,lufo,luMon)
 if (ierr.ne.OK) goto 1000

 call DtChk (modeldatactx%STN,meshdatactx%MAT,matdatactx%prpMAT,plotdatactx%IJSUR,plotdatactx%XPT, &
   meshdatactx%IEN,NEN,NSTR,NSUR,NPLDOF,NPLT,verbose,iDtfl,luMon,CST)

 if (.not.elmon) then
 call writria (plotdatactx%XPT,NUMTR,lufo, 12,ierr)
 endif
 if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xedit (cmd,ncmd,ierr)

use iomodule
use postProcessData, only: opn
use errors,          only: ok, syntax, prgerr

! Routine that X-ecutes the "edit" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)





!-locl
 integer csystem
#include "lenchar.i"
 character(len=LENCHAR) file,syscmd
 integer, external :: lnblk, nextlu
integer :: lu
!-init
 ierr = OK
 file = ' '

!    call only allowed at lowest command stack level from stdin

 if (lucmd(ilu).ne.stdin) then
write(stderr,1)
 1        format(1x,'edit: command has no effect in command file')
ierr = SYNTAX
goto 1000
 endif
 if (ilu.gt.0) then
write(stderr,2)
 2        format(1x,'write: unresolvable command file stack')
ierr = PRGERR
goto 1000
 endif

 if (ncmd.le.1) then
file = fstack(0)
if (file.eq.'stdin') then
write(stderr,3)
 3        format(1x,'edit: no file')
ierr = SYNTAX
goto 1000
endif
 else
file = cmd(2)
 endif
 if (lucmd(ilu).ne.stdin) then
inquire(unit=lucmd(ilu),opened=opn)
if (opn) call closef(lucmd(ilu))
 endif
 lu = nextlu(0)

 if (bufedit) then
if (cmd(1).eq.'e!') then
call cptobuf (lu,file,ierr)
if (ierr.ne.OK) goto 1000
else
if (fstack(0).eq.'stdin') then
    call cptobuf (lu,file,ierr)
    if (ierr.ne.OK) goto 1000
else
    write(stderr,4)
 4               format(1x,'File has been changed: use "e!" to',1x, &
              'overwrite')
    ierr = SYNTAX
    goto 1000
endif
endif
 endif

 inquire(unit=lu,opened=opn)
 if (opn) call closef(lu)

 syscmd = 'vi ' // file(1:lnblk(file))
 ierr = csystem(syscmd)

 fcmd = fstack(0)
 if (fcmd.eq.'stdin') then
fcmd  = ' '
lucmd(ilu) = stdin
 else
lucmd(ilu) = 11
 endif
 bufedit = .false.

1000    return
 end
!-------------------------------------------------------------------------------
 subroutine xeceof (ierr)

use iomodule
use errors,          only: ok, stop, prgerr, eof

! Routine that X-ecutes the "eof" command.

 implicit none
!-pass
 integer ierr



!-locl
 integer, external :: lnblk
!-init

 if (ilu.gt.0) then
ierr = OK
!        'cf'-call ended
call closef(lucmd(ilu))
ilu = ilu - 1
fcmd = fstack(ilu)
if (fcmd.eq.'stdin') then
!        'cf' from stdin ended
fcmd = ' '
if (ilu.ne.0) then
    write(stderr,1)
 1            format(1x,'eof: unresolvable cf-stack')
    ierr = PRGERR
endif
lucmd(ilu) = stdin
goto 100
endif
 else
if (lucmd(ilu).eq.stdin) then
!        end of commands from standard input; stop program
    ierr = STOP
else
if (ierr.eq.EOF) then
!            eof on command file: divert command input to stdin
        call closef(lucmd(ilu))
        lucmd(ilu) = stdin
        ierr = EOF
else
!            quit on command file
    ierr = STOP
endif
endif
 endif

100    return
 end
!-------------------------------------------------------------------------------
 subroutine xecq (cmd,ncmd,ierr)

use iomodule
use errors,          only: eof, ok, stop

! Routine that X-ecutes "q" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-init
 ierr = OK

!    close current command file and divert to previous one
 call xeceof (ierr)
!    currently at lowest level of command files stack?
 if (ierr.ne.STOP) goto 10
!    lowest level.
!    check required whether buffer differs from command file?
 if (cmd(1).eq.'q!' .or. cmd(1).eq.'quit!') then
call bufcls()
bufedit = .false.
 else
if (bufedit) then
    if (fstack(0).eq.'stdin') then
    write(stderr,1)
 1                format(1x,'No write: use "q!" or "wq file"')
    else
    write(stderr,2)
 2                format(1x,'No write: use "q!" or "wq"')
    endif
    ierr = OK
else
call bufcls()
bufedit = .false.
endif
 endif

10    if (ierr.eq.EOF) ierr = OK
 return
 end
!-------------------------------------------------------------------------------
 subroutine xecw (cmd,ncmd,ierr)

use iomodule
use errors,          only: ok, syntax, prgerr

! Routine that X-ecutes "w" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
 character(len=255)file
 integer, external :: lnblk, nextlu
integer :: lu
!-init
 ierr = OK

!    call only allowed at lowest command stack level from stdin

 if (lucmd(ilu).ne.stdin) then
write(stderr,1)
 1        format(1x,'write: command has no effect in command file')
ierr = SYNTAX
goto 1000
 endif
 if (ilu.gt.0) then
write(stderr,2)
 2        format(1x,'write: unresolvable command file stack')
ierr = PRGERR
goto 1000
 endif

 if (ncmd.le.1) then
file = fstack(0)
if (file.eq.'stdin') then
write(stderr,3)
 3           format(1x,'write: no file')
ierr = SYNTAX
goto 1000
endif
 else
file = cmd(2)
 endif
 lu = nextlu(0)
!    write buffer contents to file
 call bufwrt (lu,file,ierr)
 bufedit = .false.

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xecr (cmd,ncmd,ierr)

USE ALGEBRA
use debugmodule,     only: iecho
use iomodule
use postProcessData, only: opn
use errors,          only: ok, syntax, prgerr

! Routine that X-ecutes the "restart" command.

implicit none

!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk, nextlu
integer :: lu
!-init
 ierr = OK
 file = ' '

!       call only allowed at lowest command stack level from stdin

if (lucmd(ilu).ne.stdin) then
    write(stderr,1)
 1       format(1x,'restart: command has no effect in command file')
    ierr = SYNTAX
    goto 1000
endif

if (ilu.gt.0) then
    write(stderr,2)
 2       format(1x,'restart: unresolvable command file stack')
    ierr = PRGERR
    goto 1000
endif

if (ncmd.le.1) then
    file = fstack(0)
    if (file.eq.'stdin') then
        write(stderr,3)
 3           format(1x,'restart: no file')
        ierr = SYNTAX
        goto 1000
    endif
else
    file = cmd(2)
endif
lu = nextlu(0)

if (bufedit) then
    if (cmd(1).eq.'r!') then
        call cptobuf (lu,file,ierr)
        if (ierr.ne.OK) goto 1000
    else
        write(stderr,4)
 4           format(1x,'File has been changed: use "r!" to',1x, &
         'overwrite')
        ierr = SYNTAX
        goto 1000
    endif
endif

inquire(unit=lu,opened=opn)
if (opn) then
    call closef(lu)
endif

fcmd = fstack(0)

if (fcmd.eq.'stdin') then
    fcmd  = ' '
    lucmd(ilu) = stdin
else
    lucmd(ilu) = nextlu(0)
endif

bufedit = .false.

if (iecho.eq.1) then
    write(stderr,5,advance='no')
 5    format(1x,'clearing arrays ... ')
endif

1000    return
end subroutine
!-------------------------------------------------------------------------------
 subroutine xech (cmd,ncmd,ierr)

use iomodule
use errors,          only: ok, syntax, prgerr

! Routine that X-ecutes the "history" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
logical, external ::  numeric
 integer, external :: lnblk
double precision, external :: chreal
integer :: nhist
!-init
 ierr = OK

!       call only allowed at lowest command stack level from stdin

if (lucmd(ilu).ne.stdin) then
    write(stderr,1)
 1       format(1x,'history: command has no effect in command file')
    ierr = SYNTAX
    goto 1000
endif
if (ilu.gt.0) then
    write(stderr,2)
 2       format(1x,'restart: unresolvable command file stack')
    ierr = PRGERR
    goto 1000
endif

 if (ncmd.gt.1) then
if (.not.numeric(cmd(2))) then
write(stderr,3)
 3        format(1x,'Usage: h[istory] [number]')
ierr = SYNTAX
goto 1000
endif
nhist = INT(chreal(cmd(2)))
 else
nhist = -1
 endif

 call bufcat (nhist,ierr)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xecu (cmd,ncmd,ierr)

use iomodule
use errors,          only: ok, syntax, prgerr

! Routine that X-ecutes the "undo" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
logical, external :: numeric
integer, external :: lnblk
double precision, external :: chreal
integer :: nhist
!-init
 ierr = OK

!       call only allowed at lowest command stack level from stdin

if (lucmd(ilu).ne.stdin) then
    write(stderr,1)
 1       format(1x,'undo: command has no effect in command file')
    ierr = SYNTAX
    goto 1000
endif
if (ilu.gt.0) then
    write(stderr,2)
 2       format(1x,'undo: unresolvable command file stack')
    ierr = PRGERR
    goto 1000
endif

 if (ncmd.gt.1) then
if (.not.numeric(cmd(2))) then
write(stderr,3)
 3        format(1x,'Usage: u[ndo] [number]')
ierr = SYNTAX
goto 1000
endif
nhist = INT(chreal(cmd(2)))
 else
nhist = 1
 endif

 call bufudo (nhist,ierr)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xyield (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "yield" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
 integer, external :: lnblk
!-init
 ierr = OK
 file = ' '

 if ((.not.elmon) .and. ncmd.lt.2) then
write(stderr,10)
10        format(1x,'Usage: yi[eld] file')
ierr = SYNTAX
goto 1000
 endif

 if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'yield: read FE input and output files first')
ierr = EXEC
goto 1000
 endif
 if (ifeout.eq.0) then
if (MODE.ne.2 .and. MODE.ne.6) then
write(stderr,30)
30        format(1x,'yield: no stresses calculated')
else
write(stderr,40)
40        format(1x,'yield: first read FE output file')
endif
ierr = EXEC
goto 1000
 endif
 if (MODE.ne.2 .and. MODE.ne.6) then
write(stderr,30)
ierr = EXEC
goto 1000
 endif
 if (.not.Mout) then
write(stderr,50) NSTEP
50       format(1x,'yield: no mechanical output at itime = ', &
     I5)
ierr = EXEC
goto 1000
 endif

 if (ncmd.ge.2) file = cmd(2)
 if (triang) then
write(stderr,51)
51        format(1x,'yield: no triangulation available')
goto 1000
 endif

 if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE,&
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
setria = .false.
if (NUMTR.le.0) then
write(stderr,60)
60        format(1x,'yield: no contour triangles generated')
goto 1000
endif
if (iecho.ge.4) then
     call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
endif
 endif

 call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000

 call effyld (plotdatactx%IJSUR,plotdatactx%XPT, verbose,luMon,CST)

 if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo,13,ierr)
 endif
 if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xTc (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "Tc" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
 integer, external :: lnblk
!-init
 ierr = OK
 file = ' '

10    format(1x,'Usage: tc file')
 if ((.not.elmon) .and. ncmd.lt.2) then
write(stderr,10)
ierr = SYNTAX
goto 1000
 endif
 if (ncmd.gt.1) file = cmd(2)

 if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'tc: read FE input and output files first')
ierr = EXEC
goto 1000
 endif
 if (ifeout.eq.0) then
if (MODE.ne.5 .and. MODE.ne.6) then
write(stderr,30)
30        format(1x,'tc: no temperatures calculated')
else
write(stderr,40)
40        format(1x,'tc: first read FE output file')
endif
ierr = EXEC
goto 1000
 endif
 if (MODE.ne.5 .and. MODE.ne.6) then
write(stderr,30)
ierr = EXEC
goto 1000
 endif
 if (.not.Tout) then
write(stderr,50) NSTEP
50        format(1x,'tc: no temperature output at itime = ',I5)
ierr = EXEC
goto 1000
 endif


 if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
setria = .false.
if (NUMTR.le.0) then
write(stderr,60)
60        format(1x,'tc: no contour triangles generated')
goto 1000
endif
if (iecho.ge.4) then
    call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
endif
 endif

 call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000

 call ctrtemp (meshdatactx%IEN, &
          meshdatactx%X,   &
         modeldatactx%D,   &
 modeldatactx%LMF,modeldatactx%tfault, &
   modeldatactx%T, &
   modeldatactx%TX,modeldatactx%LMTX,plotdatactx%IJSUR,plotdatactx%XPT,modeldatactx%SKEW,verbose,luMon, &
   CST)

 if (.not.elmon) then
call writria (plotdatactx%XPT,NUMTR,lufo, 14,ierr)
 endif
 if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xpr (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "pr" command. (pressure)

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
!-init
ierr = OK
file = ' '

 if ((.not.elmon) .and. ncmd.lt.2) then
write(stderr,10)
10        format(1x,'Usage: pr[essure] file')
ierr = SYNTAX
goto 1000
 endif

 if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'pressure: read FE input and output files first')
ierr = EXEC
goto 1000
 endif
 if (ifeout.eq.0) then
if (MODE.ne.2 .and. MODE.ne.6) then
write(stderr,30)
30        format(1x,'pressure: no stresses calculated')
else
write(stderr,40)
40        format(1x,'pressure: first read FE output file')
endif
ierr = EXEC
goto 1000
 endif
 if (MODE.ne.2 .and. MODE.ne.6) then
write(stderr,30)
ierr = EXEC
goto 1000
 endif
 if (.not.Mout) then
write(stderr,50) NSTEP
50       format(1x,'pressure: no mechanical output at itime = ', &
     I5)
ierr = EXEC
goto 1000
 endif

 if (ncmd.ge.2) file = cmd(2)
 if (triang) then
write(stderr,51)
51        format(1x,'pressure: no triangulation available')
goto 1000
 endif

 if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE,&
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
setria = .false.
if (NUMTR.le.0) then
write(stderr,60)
60        format(1x,'pressure: no contour triangles generated')
goto 1000
endif
if (iecho.ge.4) then
   call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
endif
 endif

 call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000

 call presur (meshdatactx%mat, matdatactx%prpmat, &
              meshdatactx%IEN, &
              modeldatactx%stn, & 
              plotdatactx%IJSUR, &
              plotdatactx%XPT,verbose,luMon,CST)

 if (.not.elmon) then
call writria (plotdatactx%XPT,NUMTR,lufo, 15,ierr)
 endif
 if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xVolumetricStrain (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule
use constants,       only: strain, strainRate
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "volume" command.
! Computes the volumetric strain

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
 integer, external :: lnblk
 logical vrate
integer :: i
!-init
 ierr = OK
 vrate = .false.
i    = 2
file = ' '
!

!write(*,*) "entered xVolumetricStrain"

100     if (i.le.ncmd) then
    if (cmd(i).eq.'rate') then
        vrate = .true.
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

 if ((.not.elmon) .and. ncmd.lt.2) then
write(stderr,10)
10        format(1x,'Usage: vol[ume_change] [rate] file')
ierr = SYNTAX
goto 1000
 endif

 if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'volume: read FE input and output files first')
ierr = EXEC
goto 1000
 endif

 if (ifeout.eq.0) then
if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30        format(1x,'volume: no strains calculated')
else
    write(stderr,40)
40        format(1x,'volume: first read FE output file')
endif
ierr = EXEC
goto 1000
 endif

 if (MODE.ne.2 .and. MODE.ne.6) then
write(stderr,30)
ierr = EXEC
goto 1000
 endif

 if (.not.Mout) then
write(stderr,50) NSTEP
50       format(1x,'volume: no mechanical output at itime = ', &
     I5)
ierr = EXEC
goto 1000
 endif

 if (vrate .and. NSTEP.eq.0) then
write(stderr,60)
60        format(1x,'volume: no rate calculation at time zero')
ierr = EXEC
goto 1000
 endif

 if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
setria = .false.

if (NUMTR.le.0) then
    write(stderr,70)
70        format(1x,'volume: no contour triangles generated')
    goto 1000
endif

if (iecho.ge.4) then
   call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
endif
 endif

 call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
 if (ierr.ne.OK) then
goto 1000
 endif
!

!write(*,*) "vrate: ", vrate

if (vrate) then
!    write(*,*) "modeldatactx%rate: ", modeldatactx%rate(1:3, 1:10)
    ! call with rate
    call VolumetricStrain (meshdatactx%IEN, &
             plotdatactx%IJSUR,plotdatactx%XPT,meshdatactx%MAT,matdatactx%PrPMAT, &
#ifndef EXPERIMENTAL_NORATESTORAGE
       modeldatactx%rate,&
#else
    modeldatactx%STR, &
#endif 
    verbose,luMon,CST,vrate,strainRate)

else
!    write(*,*) "modeldatactx%str: ", modeldatactx%str(1:3, 1:10)
    ! call with STR (strain)
    call VolumetricStrain (meshdatactx%IEN, &
             plotdatactx%IJSUR,plotdatactx%XPT,meshdatactx%MAT,matdatactx%PrPMAT, &
       modeldatactx%STR,verbose,luMon,CST,vrate, strain)
endif
!

!write(*,*) "B entries: "
!    call Printvec(plotdatactx%XPT,1,20)

 if (.not.elmon) then
call writria (plotdatactx%XPT,NUMTR,lufo, 16,ierr)
 endif

 if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
 end
!-------------------------------------------------------------------------------
subroutine xcp (cmd,ncmd,ierr)

use modeldefinition, only: MODE
use timestepmodule,  only: nmprt, ntprt, nintg
use debugmodule,     only: iecho
use iomodule,        only: verbose, lufo, ifein, ifeout, ascout, ascin, stderr, &
                           FEDSK, luin, stdout, append
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "cp" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) DskOut
character(len=1)bs
logical AscO,AscI
integer openform,ifopen,ifbopen
integer, external :: lnblk
integer :: i, iend, ifeold, ireq, lufold, luiold
!-init
bs = CHAR(8)

if (ifein.eq.0) then
    write(stderr,1)
 1        format(1x,'cp: first read FE-input file')
    ierr = EXEC
    goto 1000
endif

if (ncmd.lt.3) then
    write(stderr,2)
 2        format(1x,'Usage: cp FromFile ToFile')
    ierr = SYNTAX
    goto 1000
endif

ifeold = ifeout
luiold = luin
lufold = lufo
AscO   = ascout
AscI   = ascin

luin   = 19
lufo   = 20
ifeout = 0
FEDSK  = cmd(2)
DskOut = cmd(3)

call shvar (FEDSK)
call shvar (DskOut)

ierr = openform (luin,FEDSK,ascin)
if (ierr.ne.OK) then
    goto 900
endif

! Make sure the file is converted into another format (ascii/bin) than
! the original file (FEDSK.DAT)
ascout = (.not.ascin)
if (ascout) then
    ierr = ifopen  (lufo,DskOut,'unknown')
else
    ierr = ifbopen (lufo,DskOut,'unknown')
endif

!    Determine last output time
call LSTSTP (NINTG,MODE,NMPRT, NTPRT,Iend)

if (iecho.eq.1) then
    write(stdout,3,advance='no')
 3    format(1x,'Copying block     0')
endif

do I=0,Iend,1
    Ireq = I
    if (iecho.eq.1) then
        write(stdout,4,advance='no') bs,bs,bs,bs,bs,Ireq
 4        format(A,A,A,A,A,I5)
    endif
! MvdZ
! Third argument denotes call by xcp (0)
    call FOUTCP (Ireq,Iend,0,ierr)
    if (ierr.ne.OK) then
        if (iecho.eq.1) then
            write(stdout,'(/)')
            call flushf(stdout)
        endif
        write(stderr,5)
 5          format(1x,'cp: failed')
        call closef (lufo)
        ierr = EXEC
        goto 900
    endif
enddo

if (iecho.eq.1) then
    write(stdout,'(/)')
endif

call closef(lufo)

900    ascout = AscO
ascin  = AscI
ifeout = ifeold
luin   = luiold
lufo   = lufold

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xdu (cmd,ncmd,ierr)

use modeldefinition, only: mode
use timestepmodule,  only: nmprt, ntprt, nintg
use debugmodule,     only: iecho
use iomodule,        only: verbose, lufo, luMon, stdout, luin, append, &
                           ifein, ifeout, ascout, ascin, stderr, FEDSK
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "du" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) DskOut
character(len=1)bs
logical AscO,AscI
integer openform,ifopen,ifbopen
integer, external :: lnblk
integer :: i, iend, ifeold, ireq, lufold, luiold
!-init
bs = CHAR(8)

if (ifein.eq.0) then
    write(stderr,1)
 1        format(1x,'cp: first read FE-input file')
    ierr = EXEC
    goto 1000
endif

if (ncmd.lt.3) then
    write(stderr,2)
 2        format(1x,'Usage: du FromFile ToFile')
    ierr = SYNTAX
    goto 1000
endif

ifeold = ifeout
luiold = luin
lufold = lufo
AscO   = ascout
AscI   = ascin

luin   = 19
lufo   = 20
ifeout = 0
FEDSK  = cmd(2)
DskOut = cmd(3)

call shvar (FEDSK)
call shvar (DskOut)

ierr = openform(luin,FEDSK,ascin)
if (ierr.ne.OK) then
    goto 900
endif

! MvdZ Not used in du
!    ascout = (.not.ascin)
if (ascout) then
    ierr = ifopen  (lufo,DskOut,'unknown')
else
    ierr = ifbopen (lufo,DskOut,'unknown')
endif

!    Determine last output time
call LSTSTP (NINTG,MODE,NMPRT, NTPRT,Iend)

if (iecho.eq.1) then
    write(stdout,3,advance='no')
 3    format(1x,'Writing block     0')
endif

do I=0,Iend,1
    Ireq = I
    if (iecho.eq.1) then
        write(stdout,4,advance='no') bs,bs,bs,bs,bs,Ireq
 4        format(A,A,A,A,A,I5)
    endif

! MvdZ
! Third argument denotes call by xdu (1)
    call FOUTCP (Ireq,Iend,1,ierr)
    if (ierr.ne.OK) then
        if (iecho.eq.1) then
            write(stdout,'(/)')
            call flushf(stdout)
        endif
        write(stderr,5)
 5        format(1x,'du: failed')
        call closef (lufo)
        ierr = EXEC
        goto 900
    endif
enddo

if (iecho.eq.1) then
    write(stdout,'(/)')
endif

call closef(lufo)

900    ascout = AscO
ascin  = AscI
ifeout = ifeold
luin   = luiold
lufo   = lufold

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xcore (cmd,ncmd,ierr)

USE MODELDATAMODULE
use materialsmodule, only: matdatactx
use modeltopology,   only: NSD,NDOF,NSTR
use iomodule,        only: ifein, lufo, stderr
use errors,          only: exec


! Routine that X-ecutes the "core" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer ifopen
integer, external :: lnblk
!-init
file = 'plncore'

if (ifein.eq.0) then
    write(stderr,1)
 1        format(1x,'core: first read FE-input file')
    ierr = EXEC
    return
endif

if (ncmd.gt.1) then
    file = cmd(2)
endif

ierr = ifopen  (lufo,file,'old')

if (ierr.ne.0) then
    ierr = EXEC
    return
endif

call RdCore (modeldatactx%T,modeldatactx%Thermalflux, &
             matdatactx%prpmat,matdatactx%PRPLAS, &
             modeldatactx%D, &
             modeldatactx%DX,modeldatactx%STN,lufo,stderr,ierr)

call closef(lufo)

return
end subroutine
!-------------------------------------------------------------------------------
subroutine xheat (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, mout, append
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "heat" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
integer, external :: lnblk
!-init
ierr = OK
file = ' '

if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: heat file')
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'heat: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'heat: no stresses and strain rates',1x, &
         'calculated')
    else
        write(stderr,40)
40        format(1x,'heat: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'heat: no mechanical output at itime = ',I5)
    ierr = EXEC
    goto 1000
endif

if (NSTEP.le.0) then
    write(stderr,60) NSTEP
60        format(1x,'heat no strain rates at itime = ',I5)
    ierr = EXEC
    goto 1000
endif

if (ncmd.ge.2) then
    file = cmd(2)
endif

if (triang) then
    write(stderr,51)
51        format(1x,'heat: no triangulation available')
    goto 1000
endif

if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
    setria = .false.
    if (NUMTR.le.0) then
        write(stderr,70)
70        format(1x,'heat: no contour triangles generated')
        goto 1000
    endif

    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif

call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) then
    goto 1000
endif

call frheat (meshdatactx%IEN,modeldatactx%STN, &
#ifndef EXPERIMENTAL_NORATESTORAGE
  modeldatactx%rate,&
#endif
   plotdatactx%IJSUR,plotdatactx%XPT, &
   meshdatactx%mat,matdatactx%prpmat,verbose,luMon,CST)

if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 17,ierr)
endif

if (ierr.ne.OK) then
    goto 1000
endif

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xsrflx (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, mout, append
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "surflux" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file         
 logical, external :: numeric
 integer, external :: lnblk
double precision, external :: chreal
integer :: i, iside
!-init
ierr = OK

if (ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: surflux iside')
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'surflux: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'surflux: no thermal calculations done')
    else
        write(stderr,40)
40        format(1x,'surflux: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'surflux: no thermal output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif

if (.not.elmon) then
    write(stderr,60)
60        format(1x,'surflux: first set monitor elements')
    ierr = EXEC
    goto 1000
endif

i = 2
file = ' '
100    if (ncmd.ge.i) then
    if (numeric(cmd(i))) then
        iside = int(chreal(cmd(i)))
        if (iside.lt.1 .or. iside.gt.8) then
            write(stderr,70)
70            format(1x,'surflx: iside out of range [1,8]')
            ierr = EXEC
            goto 1000
        endif
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) then
    goto 1000
endif

call surflx (meshdatactx%X,   &
         meshdatactx%IEN, &
         meshdatactx%mat,modeldatactx%T,matdatactx%PRPTEM, &
        modeldatactx%D,   &
   modeldatactx%LMF,modeldatactx%tfault, &
   modeldatactx%TX,modeldatactx%LMTX,modeldatactx%SKEW,verbose,luMon,iside, &
   ierr)

if (ierr.ne.OK) then
    goto 1000
endif

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xbeta (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR, nType
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, mout, append
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt  
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "beta" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file         

integer, external :: lnblk
integer :: i, ix, nix

!-init
ierr = OK
file = ' '
i    = 2
ix   = 1
nix  = 0

10    format(1x,'Usage: be[ta] [x|y|z] file')

100    if (ncmd.ge.i) then
    if (cmd(i).eq.'x' .or. &
        cmd(i).eq.'y' .or. &
        cmd(i).eq.'z') then
        if (nix.gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif

        if (NTYPE.eq.1) then

            if (cmd(i).eq.'x') then
                ix = 1
            else if (cmd(i).eq.'y') then
                ix = 2
            else
                write(stderr,20)
20                   format(1x,'beta: plane stress/strain: only x-', &
                 1x,'and y-components supported')
                ierr = SYNTAX
                goto 1000
            endif
        else
            write(stderr,30) NTYPE
30            format(1x,'beta: PrgType=',I1,' not supported')
            ierr = EXEC
            goto 1000
        endif
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif

    i = i + 1
    goto 100
endif

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,40)
40        format(1x,'beta: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,50)
50        format(1x,'beta: no strains calculated')
    else
        write(stderr,60)
60        format(1x,'beta: first read FE output file')
    endif
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,50)
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,70) NSTEP
70       format(1x,'beta: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif
!
if (triang) then
!        invoke "triangulate" to generate smoothly interpolated triangles
    call trian (meshdatactx%IEN, &
                meshdatactx%X,   &
               modeldatactx%D,   &
               modeldatactx%STN, &
               modeldatactx%STR, &
#ifndef EXPERIMENTAL_NORATESTORAGE
           modeldatactx%rate,&
#endif
       modeldatactx%LMF,modeldatactx%LMX, &
       modeldatactx%tfault,meshdatactx%mat,matdatactx%prpmat, &
       matdatactx%PRPLAS,modeldatactx%DXE,modeldatactx%SKEW,1,deformed, &
       deffac,file)

    goto 1000
endif
!

if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)

    setria = .false.

    if (NUMTR.le.0) then
        write(stderr,80)
80        format(1x,'beta: no contour triangles generated')
        goto 1000
    endif

    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif
!
call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
if (ierr.ne.OK) then
    goto 1000
endif
!
call plbeta (meshdatactx%IEN,plotdatactx%IJSUR,plotdatactx%XPT,modeldatactx%STR,verbose,luMon, &
             ix,CST)
!
if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 18,ierr)
endif

if (ierr.ne.OK) then
    goto 1000
endif
!
call OfClose (ierr,file,elmon,lufo,luMon)
!
1000 return
end subroutine
!-------------------------------------------------------------------------------
subroutine xgref (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, ITMODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule,        only: ifein, ifeout, ascout, stderr
use plotControl,     only: g0, &
                           ixpnd, iypnd, &
                           nxpnd, nypnd
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "gref" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



!-locl
logical, external :: numeric
integer, external :: lnblk
double precision, external :: chreal
double precision :: x0
!-init
ierr = OK

10    format(1x,'Usage: gref x')

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'gref: read FE input file(s) first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    if (ITMODE.ne.1 .or. MODE.lt.3) then
        write(stderr,30)
30        format(1x,'gref: read FE output file first')
        ierr = EXEC
        goto 1000
    endif
endif

if (ncmd.lt.2) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (numeric(cmd(2))) then
    x0 = chreal(cmd(2))
else
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

call gref (meshdatactx%IEN, &
           meshdatactx%X,   &
          modeldatactx%D,   &
   meshdatactx%MAT,matdatactx%PrPMAT,matdatactx%PRPTEM, &
   modeldatactx%T,modeldatactx%LMF,modeldatactx%tfault,&
  modeldatactx%DXE,ixpnd,iypnd,modeldatactx%TX, &
 modeldatactx%LMTX,nxpnd,nypnd,x0)

if (iecho.eq.1) then
    write(stderr,40) g0
40    format(1x,'Gref = ',1PG12.4,' m/s^2')
endif

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xgrav (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx, grav
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, append
use monitorVars,     only: elmon
use errors,          only: ok, syntax
use plotControl,     only: crdmin, crdmax, dispmx, &
                           ixpnd, iypnd, &
                           nxpnd, nypnd

! Routine that X-ecutes the "grav" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)
!-locl
#include "lenchar.i"
character(len=LENCHAR) file
logical, external :: numeric
integer, external :: lnblk
double precision, external :: chreal
integer :: i, nnum
double precision :: xmin, xmax, xstp

!-init
ierr = OK
file = ' '

xmin = crdmin(1)
xmax = crdmax(1) + dispmx(1)
xstp = (xmax-xmin)/15d0

10    format(1x,'Usage: gra[vity] [xstep [xmin xmax]] file')

nnum = 0
do i=2,ncmd
    if (numeric(cmd(i))) then
        nnum = nnum + 1
        if (nnum.eq.1) then
            xstp = chreal(cmd(i))
        else if (nnum.eq.2) then
            xmin = chreal(cmd(i))
        else if (nnum.eq.3) then
            xmax = chreal(cmd(i))
        else
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
    endif
enddo

if (lnblk(file).le.0 .and. (.not.elmon)) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

call OfOpen (ierr,file,.true.,elmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) then
    goto 1000
endif

call gravty (meshdatactx%IEN, &
             meshdatactx%X,   &
            modeldatactx%D,   &
             meshdatactx%MAT, &
              matdatactx%PrPMAT, &
              matdatactx%PRPTEM, &
            modeldatactx%T, &
            modeldatactx%surf, &
           grav,ixpnd,iypnd, &
            modeldatactx%TX, &
            modeldatactx%LMTX, &
            modeldatactx%LMF, &
            modeldatactx%tfault, &
            modeldatactx%DXE, &
           nxpnd,nypnd,xmin, xmax,xstp,lufo,verbose,luMon)

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xcons (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
use materialsmodule, only: matdatactx
USE MESHDATAMODULE
use debugmodule,     only: iecho
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, mout, append
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax


! Routine that X-ecutes the "cons" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk
integer :: i

!-init
ierr = OK
file = ' '
i    = 2

10    format(1x,'Usage: cons[ervation] file')

100    if (ncmd.ge.i) then
    if (lnblk(file).gt.0) then
        write(stderr,10)
        ierr = SYNTAX
        goto 1000
    endif
    file = cmd(i)
    i = i + 1
    goto 100
endif

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'cons: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30        format(1x,'cons: no mechanical solution')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    write(stderr,40)
40        format(1x,'cons: first read FE output file')
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'cons: no mechanical output at itime = ',I5)
    ierr = EXEC
    goto 1000
endif

if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
            plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)

    setria = .false.
    if (NUMTR.le.0) then
        write(stderr,60)
60        format(1x,'cons: no contour triangles generated')
        goto 1000
    endif

    if (iecho.ge.4) then
        call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
    endif
endif

call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) then
    goto 1000
endif

call MassConserv (meshdatactx%X,   &
             modeldatactx%D,   &
              meshdatactx%IEN, &
   meshdatactx%MAT,matdatactx%PrPMAT, &
   matdatactx%PRPTEM,plotdatactx%XPT,plotdatactx%IJSUR,modeldatactx%T,&
   modeldatactx%T0,modeldatactx%TX, &
   modeldatactx%STN,modeldatactx%LMTX,verbose,luMon,CST)

if (.not.elmon) then
    call writria (plotdatactx%XPT,NUMTR,lufo, 19,ierr)
endif

if (ierr.ne.OK) then
    goto 1000
endif

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xsave (ierr)

USE MODELDEFINITION, only: MODE, NUMTR, nsurf
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: ifein, ifeout, ascout, stderr
use timestepmodule,  only: nstep
use plotControl,     only: isvsrf
use errors,          only: exec, ok

! Routine that X-ecutes the "save" command.

implicit none
!-pass
integer ierr




!-init
ierr = OK

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'save: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30        format(1x,'save: no mechanical solution')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    write(stderr,40)
40        format(1x,'save: first read FE output file')
    ierr = EXEC
    goto 1000
endif

if (NSURF.le.0) then
    write(stderr,50)
50        format(1x,'save: no surface points defined')
    ierr = EXEC
    goto 1000
endif

if (NSTEP.eq.0) then
!        copy into SURF array
    call ADDSRF (modeldatactx%ISURF, &
                 modeldatactx%SURF,  &
                  meshdatactx%X,     &
                 modeldatactx%D,     &
               NSD,NDOF, &
                NSURF,0)
endif

call SVSRF (modeldatactx%ISURF, &
            modeldatactx%SURF,  &
             plotdatactx%SRF1,  &
             meshdatactx%X,     &
            modeldatactx%D,     &
               NSURF)

isvsrf = NSTEP

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xsedi (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, mout, append
use timestepmodule,  only: nstep
use monitorVars,     only: elmon, npmon
use plotControl,     only: isvsrf
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "sedi" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
logical sedsav
integer, external :: lnblk
save sedsav
data sedsav/.false./
integer :: i

!-init
ierr = OK
file = ' '
i    = 2

10    format(1x,'Usage: sedi[ment] [no[save]] file')
100 if (ncmd.ge.i) then
    if (cmd(i)(1:4).eq.'nosa') then
        sedsav = .false.
    else if (cmd(i)(1:4).eq.'save') then
        sedsav = .true.
    else
        if (lnblk(file).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'sedi: read FE input and output files first')
    ierr = EXEC
    goto 1000
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30        format(1x,'sedi: no mechanical solution')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    write(stderr,40)
40        format(1x,'sedi: first read FE output file')
    ierr = EXEC
    goto 1000
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'sedi: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    goto 1000
endif

if (sedsav) then
    if (isvsrf.lt.0) then
        write(stderr,60)
60        format(1x,'sedi: first store a surface (see "save")')
        ierr = EXEC
        goto 1000
    else if (isvsrf.gt.NSTEP) then
        write(stderr,70) isvsrf,NSTEP
70           format(1x,'sedi: restoration? Saved erosion surface', &
         ' dates from itime = ',I5/1x,'sedi: Currently we''re', &
         ' living at itime = ',I5)
    endif
endif

call OfOpen (ierr,file,.true.,npmon,append,.false.,lufo,luMon)

if (ierr.ne.OK) then
    goto 1000
endif

call OutSrf (modeldatactx%ISURF, &
             modeldatactx%SURF,  &
              plotdatactx%SRF1,  &
              meshdatactx%X,     &
             modeldatactx%D,     &
       NSD,NDOF,lufo,luMon,verbose,sedsav)

call OfClose (ierr,file,npmon,lufo,luMon)

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xread (cmd,ncmd,ierr)

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: plotdatactx
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: ifein, stderr
use errors,          only: exec, ok, syntax
! Routine that X-ecutes the "read" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk
integer :: i
!-init
ierr = OK
file = ' '
i    = 2

10    format(1x,'Usage: read NpFile')

100    if (ncmd.ge.i) then
    if (lnblk(file).gt.0) then
        write(stderr,10)
        ierr = SYNTAX
        goto 1000
    endif
    file = cmd(i)
    i = i + 1
    goto 100
endif

if (lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'read: read FE input file first')
    ierr = EXEC
    goto 1000
endif

call RdXinp (file,plotdatactx%xinp,NSD,meshdatactx%nvglobal,ierr)

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xsample (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx, &
                           MaterialIDOfElement, &
                           lmat
use debugmodule,     only: AllocateError
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, stdout, mout, tout, append
use conversions,     only: oneCharIsNumeric, &
                           resize_dblarray
use constants,       only: SampleDisplacement, &
                           SampleVelocity, &
                           SampleViscosity, &
                           SampleTemperature, &
                           SampleStress, &
                           SampleStrain, &
                           SampleRotation, &
                           SampleStrainRate, &
                           SampleHeatFlux, &
                           SampleOnCubicGrid, &
                           SampleOnAnyLine, &
                           SampleFromFile
use timestepmodule,  only: nstep
use postProcessData, only: deformed
use errors,          only: exec, ok, syntax



! Routine that X-ecutes the "sample" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




! for stdout

!-locl
logical, external :: numeric
double precision :: S(5)
double precision, dimension(:), allocatable :: xlist, ylist
double precision :: x,y
integer  :: npoints, nx, ny, idx
integer  :: i,ii,j,jj, iostat
integer  :: xysize, oldsize, newsize
integer  :: error

character(len=255) :: WriteToFileName
integer, parameter :: writeToFileHandle=47
integer            :: writeTarget ! becomes either filehandle, or standard out
logical            :: writeToFile

character(len=255) :: argumentString

logical            :: samplePointFound

integer            :: lastSampledElement

logical :: mo, to, there
double precision :: xmin, xmax, xinc
double precision :: ymin, ymax, yinc
double precision :: xr, yr
integer :: idata, iElement
integer :: nCommas, dummyInt
double precision :: dummyDouble

integer          :: matNumber
double precision, allocatable :: viscosities(:)

!character(len=255) :: outputFile
!-init

dummyDouble = 0d0

mo = .false.
to = .false.

xysize = 10

WriteToFileName = " "
writeToFile     = .false.

argumentString = " "
lastSampledElement = 1

!
!        write(*,*) 'entered xsample'

10   format(1x,'Usage: sa[mple] d|v|T|s|e|r|f|w coordinate [range] [filename]')
!
! if the number of arguments is 4; then the argument is of the shape char, coor, xrange, yrange
! if the number of arguments is 3; then the argument is of the shape char, filename


!   d = displacement   --->   idata = 1
!   v = velocity       --->   idata = 2
!   T = temperature    --->   idata = 3
!   s = stress         --->   idata = 4
!   e = strain         --->   idata = 5
!   w = rotation       --->   idata = 8
!   r = strain rate    --->   idata = 6
!   f = heat flux      --->   idata = 7

if (ncmd.lt.3) then
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

!write(*,*) "cmd 2", trim(cmd(2))


! determine the quantity to sample
if (cmd(2).eq.'d') then
    idata = SampleDisplacement
    mo = .true.
else if (cmd(2).eq.'v') then
    idata = SampleVelocity
    mo = .true.
else if (cmd(2).eq.'c') then
    idata = SampleViscosity
    mo = .true.
else if (cmd(2).eq.'T') then
    idata = SampleTemperature
    to = .true.
else if (cmd(2).eq.'s') then
    idata = SampleStress
    mo = .true.
else if (cmd(2).eq.'e') then
    idata = SampleStrain
    mo = .true.
else if (cmd(2).eq.'r') then
    idata = SampleStrainRate
    mo = .true.
else if (cmd(2).eq.'f' .or. cmd(2).eq.'q') then
    ! geophysicists tend to use q for heat flux.
    ! oh, well, nobody is perfect. ;-)
    idata = SampleHeatFlux
    to = .true.
else if (cmd(2).eq.'w') then
    idata = SampleRotation
    mo = .true.
else
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif

!*** make a list of the coordinates to be samples.


argumentString = cmd(3)

!write(*,*) "cmd3, first char", trim(cmd(3))
!write(*,*) "argumentString", trim(argumentString(1:1))


if      (oneCharIsNumeric(argumentString(1:1)) .or. argumentString(1:1).eq."-") then

        ! there are a few possible syntaxes:
        ! sample d|v|T|s|e|w|r|f <xrange> <yrange>
        ! sample d|v|T|s|e|w|r|f <xrange> <yrange> filename


    ! an xrange and yrange has been given
    call parseSampleRange(cmd(3),xmin,xmax,xinc,ierr)
    if (ierr.ne.OK) then
        goto 1000
    endif

    call parseSampleRange(cmd(4),ymin,ymax,yinc,ierr)
    if (ierr.ne.OK) then
        goto 1000
    endif

!    write(*,*) "found x range: ", xmin,xmax,xinc
!    write(*,*) "found y range: ", ymin,ymax,yinc

    if (xmin.eq.xmax) then
        nx = 1
    else
        nx = floor((xmax-xmin)/xinc)+1
    endif

    if (ymin.eq.ymax) then
        ny = 1
    else
        ny = floor((ymax-ymin)/yinc)+1
    endif


    ! build a 2D list
    npoints = (nx * ny)

    allocate(xlist(npoints), stat=error)
    call AllocateError("xlist", error)

    allocate(ylist(npoints), stat=error)
    call AllocateError("ylist", error)

    idx = 1

    do jj=1,ny
        yr = ymin + DBLE(jj-1)*yinc
        do ii=1,nx
            xr = xmin + DBLE(ii-1)*xinc
            xlist(idx) = xr
            ylist(idx) = yr
            idx = idx + 1
        enddo
    enddo

    ! 2D list completed
    if (ncmd.eq.5) then
        ! write to file
        WriteToFileName = trim(cmd(5))
        writeToFile = .true.
    else if (ncmd.eq.4) then
        ! write to standard out
        writeToFile = .false.
    else
        write(*,*) "Syntax error; see: ? plnplt"
    endif


else if (argumentString(1:4) .eq. "line") then
    write(*,*) "sampling a line; not yet correctly implemented"

    ! there are a few possible syntaxes:
    ! sample d|v|T|s|e|w|r|f line <xstart,ystart> <xend,yend>
    ! sample d|v|T|s|e|w|r|f line <xstart,ystart> <xend,yend> [npoints] filename
    if (ncmd.eq.6) then

            call parseSampleRange(cmd(3),xmin,ymin,dummyDouble,ierr)
            if (ierr.ne.OK) goto 1000
            call parseSampleRange(cmd(4),xmax,ymax,dummyDouble,ierr)
            if (ierr.ne.OK) goto 1000

            read(cmd(5),*) nPoints


        ! build a 1D list
            allocate(xlist(npoints), stat=error)
            allocate(ylist(npoints), stat=error)


            xinc = dble(xmax - xmin) / dble(nPoints-1)
            yinc = dble(ymax - ymin) / dble(nPoints-1)

            do i=1,nPoints
                xlist(i) = xmin + (i-1) * xinc
                ylist(i) = ymin + (i-1) * yinc
            enddo
                ! 1D list completed 

            WriteToFileName = trim(cmd(6))
            writeToFile = .true.

    else if (ncmd.eq.5) then

        writeToFile = .false.

    else
        write(*,*) "Syntax error; see: ? plnplt"
    endif

else if (argumentString(1:4) .eq. "file") then

!    write(*,*) "thinks that argument ", trim(argumentString), " is a file"

!    write(*,*) "allocate size: ", xysize

        ! the third argument is a file name, containing coordinate pairs
        ! problem is that we do not know how many 
        ! coordinate pairs there are in this file.
        ! We will use an array, starting at size 10,
        ! that we will double in size every time we run out.

        ! After completely running through the file,
        ! we arrays will be scaled back to fit the
        ! data excatly.

        allocate(xlist(xysize), stat=error)
        call AllocateError("xlist", error)

        allocate(ylist(xysize), stat=error)
        call AllocateError("ylist", error)

!        write(*,*) "allocated stuff, checking existence of: ", trim(cmd(3))


        ! check whether the file exists.
        inquire(file=trim(cmd(4)), exist=there)

!        write(*,*) "finished", there

        if (there) then


            if (ncmd.eq.5) then
                ! write to file
                WriteToFileName = trim(cmd(5))
                writeToFile = .true.
            else if (ncmd.eq.4) then
                ! write to standard out
                writeToFile = .false.
            else
                write(*,*) "Syntax error; see ? sample"
            endif




!            write(*,*) "opening file"

            open(123, file=trim(cmd(4)))

!            write(*,*) "opened file"

            idx = 1

            do
                READ(123,*,end=13) x,y
        
!                write(*,*) "read coords: ", x, y

                ! check whether xsize and ysize are full.
                if (idx > xysize) then
                    write(*,*) "resizing"
                    oldsize = xysize
                    xysize  = xysize * 2
                    newsize = xysize
                    ! xlist and ylist will overflow it not enlarged. Double its size.
                    call resize_dblarray(xlist, oldsize, newsize)
                    call resize_dblarray(ylist, oldsize, newsize)
!                    write (*,*) 'resize completed. new size: ', xysize
                else
!                    write (*,*) 'array still suff. large: ', xysize
                endif

!                write(*,*) "adding value", idx, x, y

                xlist(idx) = x
                ylist(idx) = y

                idx = idx + 1

            enddo

13                  write(0,*) 'finished reading file'
            npoints = idx

            call resize_dblarray(xlist, xysize, npoints)
            call resize_dblarray(ylist, xysize, npoints)

        else
                
                write(*,*) "Sample: Could not find input file: ", trim(cmd(3))
                stop "Leaving GTecton..."
        endif


endif

!***** list completed *******



if (ifein.eq.0) then
    write(stderr,20)
20       format(1x,'sample: read FE input file first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    write(stderr,30)
30       format(1x,'sample: read FE output file first')
    ierr = EXEC
    goto 1000
endif

! do mechanical stuff
if (mo.and.(MODE.ne.2.and.MODE.ne.6)) then
    write(stderr,40)
40       format(1x,'sample: no mechanical solution calculated')
    ierr = EXEC
    goto 1000
else
    if (.not.Mout) then
        if (idata.eq.1) then
            write(stderr,50) 'displacements',NSTEP
50               format(1x,'sample: no ',A13,' at itime = ',I5)
            ierr = EXEC
            goto 1000
        else if (idata.eq.2) then
            write(stderr,50) 'velocities   ',NSTEP
            ierr = EXEC
            goto 1000
        endif
    endif
endif

! do thermal stuff
if (MODE.ne.5.and.MODE.ne.6) then
    if (to) then
        write(stderr,60)
60           format(1x,'sample: no thermal solution calculated')
        ierr = EXEC
        goto 1000
    endif
else
    if (.not.Tout) then
        if (idata.eq.3) then
            write(stderr,50) 'temperatures ',NSTEP
            ierr = EXEC
            goto 1000
        endif
    endif
endif
!


!        ni = NINT((xmax-xmin)/xinc)+1
!        nj = NINT((ymax-ymin)/yinc)+1
!        do jj=1,nj
!         yr = ymin + DBLE(jj-1)*yinc
!         do ii=1,ni
!          xr = xmin + DBLE(ii-1)*xinc


if (writeToFile) then
    if (append) then 
        open(unit=writeToFileHandle, file=trim(WriteToFileName), iostat=ierr, position="append")
    else 
        open(unit=writeToFileHandle, file=trim(WriteToFileName), iostat=ierr, position="rewind")
    endif

    if (ierr.ne.0) then
        write(*,*) "Could not open file ", trim(WriteToFileName), &
                   "to write sample output to."
        write(*,*) "Returning..."
        return
    else
        writeTarget = writeToFileHandle
    endif
else
    writeTarget = stdout
endif

if (npoints .gt. 1e6) then
    write(*,*) "Warning: sampling ", npoints, "points. This may take a while"
endif



! There are three different sample routines.
! One of them is called based on thje quantity that is desired.
! dsmpl : for displacement and velocities
!         (because it takes slippery and faulted displacement into account)
! nsmpl : for quantities that are only known on nodal points,
!         such as temperatures
! esmpl : for quantities that are only known on elements,
!         such as stress and strain


! viscosity sampling is a special case,
! because the viscosities have to be computed first.
! It would be wasteful to this again for each sample point,
! so we do this only once, before looping over the points.
if (idata.eq.SampleViscosity) then

    ! make an array with viscosity for the global 
    ! elements.
    allocate(viscosities(meshdatactx%neglobal))

    if (LMAT) then
        ! Each element has its own material.
        ! Used when viscosities change throughout the run
        viscosities = matdatactx%PRPMAT(3,:)
    else
        ! Viscosities are grouped by material
        ! Used when viscosities are static throughout the run
        do iElement = 1, meshdatactx%neglobal
            matNumber = MaterialIDOfElement(iElement)
            viscosities(iElement) = matdatactx%PRPMAT(3, matNumber)
        enddo
    endif

endif




do j=1,npoints  ! n point to be sampled

    xr = xlist(j)
    yr = ylist(j)

!**** do this for every combination of points ******************
    if (idata.eq.SampleDisplacement) then
        call dsmpl (meshdatactx%IEN, &
                    meshdatactx%X,   &
                   modeldatactx%D,   &
                   modeldatactx%DELD,&
     modeldatactx%LMF, &
     modeldatactx%tFault,modeldatactx%dFault,modeldatactx%SKEW,xr,yr,S,1, &
     deformed,lastSampledElement,samplePointFound)

        if (samplePointFound) then
            if (verbose) then
                write(writeTarget,70) xr,yr,(S(i),i=1,NDOF)
70           format(1x,'Displacement at ',2(1PG24.16),' = ',3(1PG14.6))
            else
                write(writeTarget,*) xr,yr,(S(i),i=1,NDOF)
            endif
        else
    !        write(FileHandle,65) xr,yr
            write(*,65) xr,yr
65           format(1x,'Sample: coordinate not found ',2(1PG14.6))
        endif

    else if (idata.eq.SampleVelocity) then

        call dsmpl (meshdatactx%IEN, &
                    meshdatactx%X,   &
                   modeldatactx%D,   &
                   modeldatactx%DELD,&
     modeldatactx%LMF, &
     modeldatactx%tFault,modeldatactx%dFault,modeldatactx%SKEW,xr,yr,S,2, &
     deformed,lastSampledElement,samplePointFound)

        if (samplePointFound) then
!        write(FileHandle,80) xr,yr,(S(i),i=1,NDOF)
            if (verbose) then
                write(writeTarget,80) xr,yr,(S(i),i=1,NDOF)
80               format(1x,'Velocity at ',2(1PG14.6),' = ',3(1PG14.6))
            else
                write(writeTarget,*) xr,yr,(S(i),i=1,NDOF)
            endif
        else
            write(*,65) xr,yr
        endif

    else if (idata.eq.SampleViscosity) then

    call esmpl (meshdatactx%IEN, &
                meshdatactx%X,   &
               modeldatactx%D,   &
               modeldatactx%LMF, &
             modeldatactx%tFault, &
     modeldatactx%SKEW,xr,yr, &
               viscosities, &
               1,S,deformed,.false.,lastSampledElement,samplePointFound)

    if (samplePointFound) then

        if (verbose) then
            write(writeTarget,81) xr,yr,(S(1))
81           format(1x,'Viscosity at ',2(1PG14.6),' = ',1PG14.6)
        else
            write(writeTarget,*) xr,yr,(S(1))
        endif
    else
!        write(FileHandle,65) xr,yr
        write(*,65) xr,yr
    endif




  else if (idata.eq.SampleTemperature) then
    call nsmpl (meshdatactx%IEN, &
                meshdatactx%X,   &
               modeldatactx%D,   &
     modeldatactx%LMF,modeldatactx%tFault, &
     modeldatactx%SKEW,xr,yr,modeldatactx%T,1,S,deformed,lastSampledElement,samplePointFound)

    if (samplePointFound) then

        if (verbose) then
            write(writeTarget,90) xr,yr,S(1)
90           format(1x,'Temperature at ',2(1PG14.6),' = ',1PG14.6)
        else
            write(writeTarget,*) xr,yr,S(1)
        endif

    else
!        write(FileHandle,65) xr,yr
        write(*,65) xr,yr
    endif

  else if (idata.eq.SampleStress) then
    call esmpl (meshdatactx%IEN, &
                meshdatactx%X,   &
               modeldatactx%D,   &
               modeldatactx%LMF, &
             modeldatactx%tFault, &
     modeldatactx%SKEW,xr,yr, &
               modeldatactx%STN, &
               NSTR,S,deformed,.false.,lastSampledElement,samplePointFound)


    if (samplePointFound) then
!        write(FileHandle,101) xr,yr,(S(i),i=1,NSTR)
        if (verbose) then
            write(writeTarget,101) xr,yr,(S(i),i=1,NSTR)
  101           format(1x,'Stress at ',2(1PG14.6),' = ',5(1PG14.6))
        else
            write(writeTarget,*) xr,yr,(S(i),i=1,NDOF)
        endif

    else

!        write(FileHandle,65) xr,yr
        write(*,65) xr,yr


    endif

  else if (idata.eq.SampleStrain) then
    call esmpl (meshdatactx%IEN, &
                meshdatactx%X,   &
               modeldatactx%D,   &
               modeldatactx%LMF, &
     modeldatactx%tFault, &
     modeldatactx%SKEW,xr,yr,&
               modeldatactx%STR, &
               NSTR,S,deformed,.false.,lastSampledElement,samplePointFound)


    if (samplePointFound) then

!        write(FileHandle,110) xr,yr,(S(i),i=1,NSTR)
        if (verbose) then
            write(writeTarget,110) xr,yr,(S(i),i=1,NSTR)

  110           format(1x,'Strain at ',2(1PG14.6),' = ',5(1PG14.6))
        else
            write(writeTarget,*) xr,yr,(S(i),i=1,NSTR)
        endif

    else
    !        write(FileHandle,65) xr,yr
        write(*,65) xr,yr

    endif

  else if (idata.eq.SampleStrainRate) then
    call esmpl (meshdatactx%IEN, &
                meshdatactx%X,   &
               modeldatactx%D,   &
     modeldatactx%LMF,modeldatactx%tFault, &
     modeldatactx%SKEW,xr,yr,&
#ifndef EXPERIMENTAL_NORATESTORAGE
           modeldatactx%rate,&
#else
     modeldatactx%str, &
#endif
     NSTR,S,&
     deformed,.false.,lastSampledElement,samplePointFound)


    if (samplePointFound) then

        if (verbose) then
!        write(FileHandle,120) xr,yr,(S(i),i=1,NSTR)
            write(writeTarget,120) xr,yr,(S(i),i=1,NSTR)

  120           format(1x,'Strain rate at ',2(1PG14.6),' = ',5(1PG14.6))
        else
            write(writeTarget,*) xr,yr,(S(i),i=1,NSTR)
        endif

    else
    !        write(FileHandle,65) xr,yr
        write(*,65) xr,yr

    endif

  else if (idata.eq.SampleHeatFlux) then

!    write(*,*) "samping heat flux"

    call esmpl (meshdatactx%IEN, &
                meshdatactx%X,   &
               modeldatactx%D,   &
     modeldatactx%LMF,modeldatactx%tFault, &
     modeldatactx%SKEW,xr,yr,modeldatactx%Thermalflux, &
     NSD,S,deformed,.false.,lastSampledElement,samplePointFound)



    if (samplePointFound) then

!        write(FileHandle,130) xr,yr,(S(i),i=1,NSD)
        if (verbose) then
            write(writeTarget,130) xr,yr,(S(i),i=1,NSD)

  130           format(1x,'Heat flow at ',2(1PG14.6),' = ',2(1PG14.6))
        else
            write(writeTarget,*) xr,yr,(S(i),i=1,NSD)
        endif

    else
    !        write(FileHandle,65) xr,yr
        write(*,65) xr,yr

    endif

  else if (idata.eq.SampleRotation) then

    call esmpl (meshdatactx%IEN, &
                meshdatactx%X,   &
               modeldatactx%D,   &
     modeldatactx%LMF,modeldatactx%tFault, &
     modeldatactx%SKEW,xr,yr,modeldatactx%AA,1,S,deformed,.true.,lastSampledElement,samplePointFound)


    if (samplePointFound) then

        if (verbose) then
!        write(FileHandle,140) xr,yr,S(1)
            write(writeTarget,140) xr,yr,S(1)
  140           format(1x,'x-y rotation at ',2(1PG14.6),' = ',2(1PG14.6))
        else
            write(writeTarget,*) xr,yr,S(1)
        endif
    else
    !        write(FileHandle,65) xr,yr
        write(*,65) xr,yr

    endif

  endif
!*** end of doing this for every point **************

enddo


if (idata.eq.SampleViscosity) then
    deallocate(viscosities)
endif



!
1000    if(allocated(xlist)) deallocate(xlist)
if(allocated(ylist)) deallocate(ylist)


if (writeToFile) then
    close(writeToFileHandle)
endif


return
end
!-------------------------------------------------------------------------------
subroutine xshsl (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use modeldefinition, only: NUMSLP
use iomodule,        only: ifein, stderr
use errors,          only: exec, syntax

! Routine that X-ecutes the "shsl" command.

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
#include "lenchar.i"
 character(len=LENCHAR) file
 integer, external :: lnblk, nextlu
 logical parallel
integer :: i
integer :: lu
!-init
 save parallel
 data parallel /.false./
 file = ' '
 i = 2

10    format(1x,'Usage: shsl [p|n] file')

100    if (i.le.ncmd) then
if (cmd(i).eq.'p') then
    parallel = .true.
else if (cmd(i).eq.'n') then
    parallel = .false.
else if (lnblk(file).le.0) then
    file = cmd(i)
else
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
i = i + 1
goto 100
 endif
!
if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'shsl: read FE input file first')
ierr = EXEC
goto 1000
 endif
!
 if (lnblk(file).ne.0) then
lu = nextlu(0)
call openf (lu,file,'unknown')
 endif

! TODO, dus numslp work correctly?
 call shslp (modeldatactx%NSLIP,5,  &
        meshdatactx%ien,  &
        meshdatactx%x,     &
    NUMSLP,lu,parallel)


 if (lnblk(file).ne.0) call closef (lu)
!
1000    return
 end
!-------------------------------------------------------------------------------
subroutine xecfo (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, mout, append
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "force" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
integer monside
 logical, external :: numeric
double precision, external :: chreal
 character(len=LENCHAR) file
double precision :: PCUR(3)
integer :: IDOFTYP
!-init
 file = ' '

 if (ncmd.lt.3) then
write(stderr,10)
10        format(1x,'Usage: fo[rce] monside [file]')
ierr = SYNTAX
goto 1000
 endif

if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'force: read FE input and output files first')
ierr = EXEC
goto 1000
 endif
 if (ifeout.eq.0) then
write(stderr,40)
40        format(1x,'force: first read FE-output file')
ierr = EXEC
goto 1000
 else
if (.not.Mout) then
write(stderr,50) NSTEP
50        format(1x,'force: no mechanical output at itime = ',I5)
ierr = EXEC
goto 1000
endif
 endif
 if (.not.elmon) then
write(stderr,60)
60        format(1x,'force: no monitor elements')
ierr = EXEC
goto 1000
 endif

if (numeric(cmd(2))) then
    monside = int(chreal(cmd(2)))
    if (monside.lt.0 .or. monside.gt.4) then
        write(stderr,70)
70           format(1x,'force: monside out of range [1,4]')
        ierr = EXEC
        goto 1000
    endif
else
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif


 call CURFOR (meshdatactx%mat, &
         meshdatactx%IEN, &
         meshdatactx%X,   &
 modeldatactx%LMF, &
   modeldatactx%tfault,matdatactx%prpmat, &
   modeldatactx%STN, &
 modeldatactx%D, &
 modeldatactx%SKEW,modeldatactx%DXE,monside,pcur)

IDOFTYP=1

 if (ncmd.gt.2) file = cmd(3)
 call OfOpen (ierr,file,.true.,.true.,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000

if (IDOFTYP.eq.1) then
80       format(1x,A1,'-direction: current force = ',1PG14.6)
    write(luMon,80) 'X',PCUR(1)
elseif (IDOFTYP.eq.2) then
write(stderr,90)
90        format(1x,'force: y-direction not implemented yet')
ierr = EXEC
goto 1000
 endif

1000    return
 end
!-------------------------------------------------------------------------------
! subroutine xregrid (cmd,ncmd,ierr)

! USE MODELDATAMODULE
! USE MESHDATAMODULE

! Routine that X-ecutes the "regrid" command.

! implicit none
!-pass
! integer ncmd,ierr
! character cmd*(*)
! dimension cmd(*)

!
!
!
!-locl
!#include "lenchar.i"
! logical, external :: numericic,there
 !character(len=255)command


! if (ncmd.lt.2) then
!write(stderr,10)
!10        format(1x,'Usage: regrid material_number[s]')
!ierr = SYNTAX
!goto 1000
! endif

!if (ifein.eq.0) then
!write(stderr,20)
!20        format(1x,'regrid: read FE input file first')
!ierr = EXEC
!goto 1000
! endif
!
! if (ISTART.eq.0) then
!write(stderr,25)
!25        format(1x,'regrid: ISTART<>0 error')
!ierr = EXEC
!goto 1000
! endif
! inquire(file='RESTART.DAT',exist=there)
! if (.not.there) then
!write(stderr,27)
!27        format(1x,'regrid: "RESTART.DAT" not in $cwd')
!ierr = EXEC
!goto 1000
! endif
! ipbtot = IPXPT
! ipnp   = ipbtot + NEQ*IWORD
! ipcon  = ipnp + NUMNP
! maxcon = MBTOT - ipcon
! if (maxcon.lt.1) then
!write(stderr,28)
!28        format(1x,'regrid: overflow on B array')
!ierr = EXEC
!goto 1000
! endif
! setria = .true.
!
! do nar=2,ncmd
!    if (numeric(cmd(nar))) then
!        matno = int(chreal(cmd(nar)))
!        if (matno.le.0 .or. matno.gt.NUMAT) then
!            write(stderr,30) NUMAT
!30               format(1x,'regrid: material number out of range',1x, &
!             '[1,',I3,']')
!            ierr = EXEC
!            goto 1000
!        endif
!    else
!        write(stderr,10)
!        ierr = SYNTAX
!        goto 1000
!    endif
!
!if (nar.eq.2) then
!    call RESTIO ( &
! modeldatactx%D, &
! modeldatactx%DXE,modeldatactx%STN,modeldatactx%STR, &
!       B(ipbtot),modeldatactx%T,A(IPFTOT),A(IPVPSTR),NDOF,NEN,NSTR,1)
!endif
!            write(stdout,*) " "
!    write(stdout,50) matno
!50   format(1x, 'regrid: regridding material',1x, I5)
!call regrid (meshdatactx%IEN, &
!             meshdatactx%X,   &
!            modeldatactx%D,   &
!     meshdatactx%mat,modeldatactx%NSLIP, &
!       matno,B(ipnp),B(ipcon),maxcon,ierr)

!if (ierr.ne.OK) goto 1000

!call mapelq (meshdatactx%IEN, &
!             meshdatactx%X,   &
!            modeldatactx%D,   &
!     modeldatactx%STN,modeldatactx%STR, &
!       modeldatactx%DXE,ierr)
!if (ierr.ne.OK) goto 1000
!call renew (meshdatactx%IEN, &
!     modeldatactx%STN,modeldatactx%STR,modeldatactx%DXE,modeldatactx%NSLIP, &
!       A(IPNSELSD),ierr)
!if (ierr.ne.OK) goto 1000
! enddo
!!
! call RESTIO (modeldatactx%D,   &!A(IPD), &
! modeldatactx%DXE,modeldatactx%STN,modeldatactx%STR,B(ipbtot), &
!   modeldatactx%T,A(IPFTOT),A(IPVPSTR),NDOF,NEN,NSTR,0)
! if (iecho.eq.1) write(stdout,60)
!60    format(1x,'regrid: succesfully written "RESTART.DAT"')
! call WRITEELM (meshdatactx%ien,meshdatactx%mat)
! if (iecho.eq.1) write(stdout,70)
!70    format(1x,'regrid: succesfully written "tecin.dat.elm_reg"')
! if (NUMSLP.ne.0) then
!call WRITESLP(modeldatactx%NSLIP,A(IPDIFF),A(IPNSELSD),A(IPXDAT), &
!       A(IPWINX))
!if (iecho.eq.1) write(stdout,80)
!80        format(1x,'regrid: succesfully written "tecin.dat.sns_reg"')
! endif
!
!
!1000    return
! end
!-------------------------------------------------------------------------------
 subroutine xconn(cmd,ncmd,ierr)

use debugmodule,     only: iecho, parallellog
use iomodule,        only: verbose, lufo, luMon, ifein, ascout, stderr, &
                           append
use postProcessData, only: GMT, connout
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "connection" command needed to plot pscontour files with gmt

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
#include "lenchar.i"
 character(len=LENCHAR) file
integer, external :: lnblk
integer :: i
!-init
 i     = 2
 file = ' '

if (iecho.eq.8) then
    call ParallelLog("xconn", "Entered")
endif


10    format(1x,'Usage: conn[ection] file ')

100    if (i.le.ncmd) then
    if (lnblk(file).le.0) then
 file = cmd(i)
    else
 write(stderr,10)
 ierr = SYNTAX
 goto 1000
    endif
    i = i + 1
    goto 100
 endif
 if (lnblk(file).le.0) then
write(stderr,10)
ierr = SYNTAX
goto 1000
 endif
!
if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'conn: first read FE-input file')
ierr = EXEC
goto 1000
 endif
 if (.not.GMT) write(stderr,30)
30    format(1x,'"P" does not need a connectivity file')
!
 call OfOpen (ierr,file,ascout,.false.,append,.false., &
 lufo,luMon)
 if (ierr.ne.OK) goto 1000
!
 call conn (lufo)
!
 call OfClose (ierr,file,.false.,lufo,luMon)

if (ierr.eq.OK) then
    connout = .true.
endif

!
1000  if (iecho.eq.8) then
    call ParallelLog("xfemo", "Finished")
endif
return
 end

!-------------------------------------------------------------------------------

subroutine xBoundaries (cmd,ncmd,ierr)
! routine to plot the partition boundaries

use modelctx, only: getsize
use meshdatamodule, only: meshdatactx, &
                          BuildElementNeighbourTable2D

implicit none


integer ncmd,ierr
character cmd*(*)
dimension cmd(*)

! boundaries will be stored as lines from node to node
! the storage will be chaotic and not much use for
! analysis. They are fine for plotting, though.
integer, parameter :: nodeSpace = 10000000  ! ugly... do more elegant
integer, allocatable :: fromNode(:), toNode(:)
integer :: iBoundary
integer :: boundaryFileID
integer, external :: nextlu
character(len=1024) :: boundaryFileName

allocate(fromNode(nodeSpace))
allocate(toNode(nodeSpace))

if (ncmd.eq.2) then
    boundaryFileName = cmd(2)
else
    write(*,*) "usage: bound[aries] file"
    return
endif

! initialize the neighbour list
call BuildElementNeighbourTable2D()

! the number of boundaries is not known,
! hence storage arrays will be allocated to
! worst case scenario, which is nvglobal.
! for very small meshes that may go awry,
! so we set a minimum size of 100

!nodeSpace = max(100,meshdatactx%nvglobal)
!allocate(fromNode(nodeSpace))
!allocate(  toNode(nodeSpace))

! from the neighbour list, plot the boundaries
call findBoundaries(nodeSpace, fromNode, toNode)

boundaryFileID = nextlu(0)

open(unit=boundaryFileID, file=trim(boundaryFileName))

! and write them to file.
do iBoundary =1,nodeSpace
    if (fromNode(iBoundary).eq. 0.0) then
        ! we are done
        return
    else
        ! write an entry to file
        write(boundaryFileID,*) ">"
        write(boundaryFileID,*) meshdatactx%X(1,fromNode(iBoundary)), meshdatactx%X(2,fromNode(iBoundary))
        write(boundaryFileID,*) meshdatactx%X(1,  toNode(iBoundary)), meshdatactx%X(2,  toNode(iBoundary))
    endif
enddo

close(boundaryFileID)

end subroutine xBoundaries


!-------------------------------------------------------------------------------
subroutine xdens (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, ITMODE, NUMTR
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use debugmodule,     only: iecho
use iomodule,        only: verbose, lufo, luMon, ifein, ifeout, ascout, &
                           stderr, tout, append
use postProcessData, only: GMT, connout, &
                           cst, deffac, deformed, setria
use monitorVars,     only: elmon
use plotControl,     only: nplt
use errors,          only: exec, ok, syntax

! Routine that X-ecutes the "density" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
logical thermal
character(len=LENCHAR) file
integer, external :: lnblk
integer :: i
!-init
 ierr = OK
i    = 2
file = ' '
!
100     if (i.le.ncmd) then
    if (lnblk(file).gt.0) then
        write(stderr,10)
        ierr = SYNTAX
        goto 1000
    endif
    file = cmd(i)
    i = i + 1
    goto 100
endif
 if ((.not.elmon) .and. ncmd.lt.2) then
write(stderr,10)
10        format(1x,'Usage: density file')
ierr = SYNTAX
goto 1000
 endif

 if (ifein.eq.0) then
write(stderr,20)
20        format(1x,'density: read FE input files first')
ierr = EXEC
goto 1000
 endif

 if (GMT .and. .not.connout .and. .not.elmon) then
write(stderr,30)
30        format(1x,'density: first write connectivity file')
ierr = EXEC
goto 1000
 endif

 if (setria .and. .not.elmon) then
!        calculate triangle coordinates in XPT
    call tricrd (meshdatactx%IEN, &
                 meshdatactx%X,   &
                modeldatactx%D,   &
     plotdatactx%IJSUR,plotdatactx%XPT, &
     modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
     modeldatactx%SKEW,deformed,deffac,NUMTR,CST)
setria = .false.

if (NUMTR.le.0) then
write(stderr,70)
70        format(1x,'density: no contour triangles generated')
goto 1000
endif

if (iecho.ge.4) then
    call tstria (plotdatactx%XPT,NPLT,NUMTR, 'trel.dat')
endif
 endif

 call OfOpen (ierr,file,ascout,elmon,append,.false.,lufo,luMon)
 if (ierr.ne.OK) goto 1000
!
 thermal = (Tout.or.((MODE.eq.5.or.MODE.eq.6).and.ITMODE.ne.0))
 call TDNSTY (meshdatactx%MAT,meshdatactx%IEN, &
                  meshdatactx%X,   &
                 modeldatactx%D,   &
    modeldatactx%LMF, &
   modeldatactx%tfault,modeldatactx%T,modeldatactx%TX,modeldatactx%LMTX,&
   matdatactx%PrPMAT,matdatactx%PRPTEM, &
   plotdatactx%IJSUR,plotdatactx%XPT,thermal,verbose,luMon,ierr)
 if (ierr.ne.OK) goto 1000
!
 if (.not.elmon) then
call writria (plotdatactx%XPT,NUMTR,lufo, 20,ierr)
 endif

 if (ierr.ne.OK) goto 1000

call OfClose (ierr,file,elmon,lufo,luMon)

1000    return
 end


!-------------------------------------------------------------------------------

subroutine xElementSize(writeToFile, WriteToFileName)

use modeldefinition, only: lgdef
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: nsd, nen, ndof
use constants,       only: useDXE


implicit none
!-pass
! integer ncmd,ierr
! character cmd*(*)
! dimension cmd(*)

integer            :: iElem, iPoint, iError
integer            :: lien(nen)
double precision   :: xl(nsd, nen), dl(nsd, nen)
character(len=255) :: WriteToFileName
logical            :: writeToFile
integer            :: fileID
double precision   :: area
double precision   :: a(2), b(2)

integer, external  :: nextlu


writeToFile = .false.



if (writeToFile) then
    fileID=nextlu(0) 
    open(unit=fileID, file=WriteToFileName, iostat=iError)
    if (iError.ne.0) then
        write(*,*) "Could not open file", trim(WriteToFileName), ". Error: ", iError
    endif
endif

do iElem=1,meshdatactx%neglobal

    ! for every element, compute the area
    lien = meshdatactx%ien(:,iElem)

    call LCOORD (meshdatactx%x,xl,meshdatactx%IEN(1,iElem),0)
    if (LGDEF.ne.0) then
        call LDISP  (DL,modeldatactx%D,meshdatactx%IEN(1,iElem),NDOF,NEN)
        call AddFaultDisplacement (DL,iElem,modeldatactx%TFAULT,NDOF,NEN,132)
        call ADDSNE (DL,iElem,NDOF,NEN,useDXE)
        call REZONE (XL,DL,1d0)
    endif

    a(1) = xl(1,2) - xl(1,1)
    a(2) = xl(2,2) - xl(2,1)
    b(1) = xl(1,3) - xl(1,2)
    b(2) = xl(2,3) - xl(2,2)
    AREA = 0.5d0*(a(1)*b(2)-a(2)*b(1))

    do iPoint=1,3
        if (writeToFile) then
            write(fileID,*) xl(1,iPoint), xl(2,iPoint), area
        else
            write(*,*) xl(1,iPoint), xl(2,iPoint), area
        endif
    enddo

enddo

if (writeToFile) then
    close(fileID)
endif

end subroutine

!-------------------------------------------------------------------------------

subroutine xElementRatio(writeToFile, WriteToFileName)

use modeldefinition, only: lgdef
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: nsd, nen, ndof
use constants,       only: pi, useDXE

implicit none


integer :: iElem, iPoint, iError
integer            :: lien(nen)
double precision   :: xl(nsd, nen),dl(nsd, nen)
character(len=255) :: WriteToFileName
logical            :: writeToFile
integer            :: fileID
double precision   :: dx, dy
double precision   :: side1, side2, side3
double precision   :: radius
double precision   :: areaOuterCircle, areaInnerCircle
double precision   :: quality
double precision   :: semiPerimeter

integer, external  :: nextlu

if (writeToFile) then
    fileID=nextlu(0) 
    open(unit=fileID, file=WriteToFileName, iostat=iError)
    if (iError.ne.0) then
        write(*,*) "Could not open file", trim(WriteToFileName), ". Error: ", iError
    endif
endif

do iElem=1,meshdatactx%neglobal

    ! for every element, compute the area of the circumscribed circle, 
    ! following: https://www.mathopenref.com/trianglecircumcircle.html
    lien = meshdatactx%ien(:,iElem)

    call LCOORD (meshdatactx%x,xl,meshdatactx%IEN(1,iElem),0)
    if (LGDEF.ne.0) then
        call LDISP  (DL,modeldatactx%D,meshdatactx%IEN(1,iElem),NDOF,NEN)
        call AddFaultDisplacement (DL,iElem,modeldatactx%TFAULT,NDOF,NEN,132)
        call ADDSNE (DL,iElem,NDOF,NEN,useDXE)
        call REZONE (XL,DL,1d0)
    endif


    ! compute the lengths of the sides
    ! side 1
    dx = xl(1,1) - xl(1,2)
    dy = xl(2,1) - xl(2,2)
    side1 = sqrt(dx**2 + dy**2)

    ! side 2
    dx = xl(1,3) - xl(1,2)
    dy = xl(2,3) - xl(2,2)
    side2 = sqrt(dx**2 + dy**2)

    ! side 3
    dx = xl(1,1) - xl(1,3)
    dy = xl(2,1) - xl(2,3)
    side3 = sqrt(dx**2 + dy**2)

    radius = side1 * side2 * side3 / &
             sqrt(( side1 + side2 + side3) * &
                  ( side1 + side2 - side3) * &
                  ( side1 - side2 + side3) * &
                  (-side1 + side2 + side3)) 

    areaOuterCircle = pi * radius**2


    semiPerimeter = 0.5d0 * (side1 + side2 + side3)

    radius = semiPerimeter * &
             sqrt((1d0 - side1/semiPerimeter) * &
                  (1d0 - side2/semiPerimeter) * &
                  (1d0 - side3/semiPerimeter))
    
    areaInnerCircle = pi * radius**2

    quality = areaInnerCircle / areaOuterCircle

    ! normalize the quality by multiplying with the ratio of the circles
    ! of an equilateral triangle, so that the maximum becomes 1.
    ! for an equilat. triiangle with all sides 1, we have
    ! outer circle, radius: 1 / sqrt(3)
    ! inner circle, radius: 0.5 / sqrt(3)
    ! ratio is the square of the circles: times 4

    quality = quality * 4d0


    do iPoint=1,3
        if (writeToFile) then
            write(fileID,*) xl(1,iPoint), xl(2,iPoint), quality
        else
            write(*,*) xl(1,iPoint), xl(2,iPoint), quality
        endif
    enddo

enddo

if (writeToFile) then
    close(fileID)
endif


end subroutine

!-------------------------------------------------------------------------------

subroutine xSideRatio(writeToFile, WriteToFileName)

use modeldefinition, only: lgdef
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: nsd, nen, ndof
use constants,       only: pi, useDXE

implicit none


integer :: iElem, iPoint, iError
integer            :: lien(nen)
double precision   :: xl(nsd, nen), dl(nsd, nen)
character(len=255) :: WriteToFileName
logical            :: writeToFile
integer            :: fileID
double precision   :: dx, dy
double precision   :: side1, side2, side3
double precision   :: quality

integer, external  :: nextlu

if (writeToFile) then
    fileID=nextlu(0) 
    open(unit=fileID, file=WriteToFileName, iostat=iError)
    if (iError.ne.0) then
        write(*,*) "Could not open file", trim(WriteToFileName), ". Error: ", iError
    endif
endif

do iElem=1,meshdatactx%neglobal

    ! for every element, compute the area of the circumscribed circle, 
    ! following: https://www.mathopenref.com/trianglecircumcircle.html
    lien = meshdatactx%ien(:,iElem)

    call LCOORD (meshdatactx%x,xl,meshdatactx%IEN(1,iElem),0)
    if (LGDEF.ne.0) then
        call LDISP  (DL,modeldatactx%D,meshdatactx%IEN(1,iElem),NDOF,NEN)
        call AddFaultDisplacement (DL,iElem,modeldatactx%TFAULT,NDOF,NEN,132)
        call ADDSNE (DL,iElem,NDOF,NEN,useDXE)
        call REZONE (XL,DL,1d0)
    endif


    ! compute the lengths of the sides
    ! side 1
    dx = xl(1,1) - xl(1,2)
    dy = xl(2,1) - xl(2,2)
    side1 = sqrt(dx**2 + dy**2)

    ! side 2
    dx = xl(1,3) - xl(1,2)
    dy = xl(2,3) - xl(2,2)
    side2 = sqrt(dx**2 + dy**2)

    ! side 3
    dx = xl(1,1) - xl(1,3)
    dy = xl(2,1) - xl(2,3)
    side3 = sqrt(dx**2 + dy**2)


    quality = min(side1, side2, side3) / max(side1, side2, side3)

    ! normalize the quality by multiplying with the ratio of the circles
    ! of an equilateral triangle, so that the maximum becomes 1.
    ! for an equilat. triiangle with all sides 1, we have
    ! outer circle, radius: 1 / sqrt(3)
    ! inner circle, radius: 0.5 / sqrt(3)
    ! ratio is the square of the circles: times 4

    do iPoint=1,3
        if (writeToFile) then
            write(fileID,*) xl(1,iPoint), xl(2,iPoint), quality
        else
            write(*,*) xl(1,iPoint), xl(2,iPoint), quality
        endif
    enddo

enddo

if (writeToFile) then
    close(fileID)
endif


end subroutine


!-------------------------------------------------------------------------------

subroutine xmesh(cmd,ncmd,iError)

implicit none

character cmd*(*)
dimension cmd(*)
integer :: iError, ncmd

character(len=255) :: WriteToFileName
logical            :: writeToFile
integer            :: fileID


!if (ncmd.eq.3) then
!    writeToFile = .true.
!    WriteToFileName = trim(cmd(3))
!else if (ncmd.eq.2) then
!    writeToFile = .false.
!else
!    write(*,*) "Usage: mesh informationType [filename]"
!endif


!if (writeToFile) then
!    fileID=nextlu(0)
!    open(unit=fileID, file=WriteToFileName, iostat=iError)
!    if (iError.ne.0) then
!        write(*,*) "Could not open file", trim(WriteToFileName), ". Error: ", iError
!    endif
!endif

if (ncmd.eq.2 .or. ncmd.eq.3) then

    if (ncmd.eq.3) then
        writeToFile = .true.
        WriteToFileName = trim(cmd(3))
    else if (ncmd.eq.2) then
        writeToFile = .false.
    else
        write(*,*) "Usage: mesh informationType [filename]"
    endif

    if       (cmd(2)(1:9) .eq. "arearatio") then
        call xElementRatio(writeToFile, WriteToFileName)
    else if  (cmd(2)(1:9) .eq. "sideratio") then
        call    xSideRatio(writeToFile, WriteToFileName)
    else if  (cmd(2)(1:4) .eq. "area") then
        call  xElementSize(writeToFile, WriteToFileName)
    else
        write(*,*) "did not recognize informationType:  ", trim(cmd(2))
        write(*,*) "Must be one of:"
        write(*,*) "  arearatio"
        write(*,*) "  sideratio"
        write(*,*) "  area"
    endif

else
    write(*,*) "Usage: mesh informationType [filename]"
endif


end subroutine

!-------------------------------------------------------------------------------


#ifdef EXPERIMENTAL_FRICTION

subroutine xfriction(cmd,ncmd,iError)

use frictionmodule,  only: frictionTraction, &
                           friction
use modeldefinition, only: nFaultSegmentEntriesGlobal,NTYPE
use modeldatamodule, only: modeldatactx
use modeltopology,   only: NDOF
use conversions,     only: oneCharIsNumeric
use debugmodule,     only: iecho

implicit none

character cmd*(*)
dimension cmd(*)
integer :: iError, ncmd

integer :: faultID
logical :: IDfound
integer :: iEntry, iFault

logical, allocatable :: IDsPresent(:)
integer              :: nFaults
integer, allocatable :: faultIDs(:)

character(len=255) :: WriteToFileName
logical            :: writeToFile
integer            :: fileID


if (.not. allocated(faultIDs)) then

    ! check which fault iDs are available
    ! do this only on the first fric call

    allocate(IDsPresent(maxval(friction(:)%faultNumber)))
    IDsPresent = .false.

    nFaults = 0

    do iEntry=1,nFaultSegmentEntriesGlobal
        if (.not. IDsPresent(friction(iEntry)%faultNumber) ) then
            nFaults = nFaults + 1
            IDsPresent(friction(iEntry)%faultNumber) = .true.
        endif
    enddo

    allocate(faultIDs(nFaults))
    faultIDs = 0

    iFault = 0
    do iEntry=1,maxval(friction(:)%faultNumber)
        if (IDsPresent(iEntry)) then
            iFault = iFault + 1
            faultIDs(iFault) = iEntry
        endif
    enddo

    deallocate(IDsPresent)

endif

! process command, and use available fault IDs when necessary

if (oneCharIsNumeric(cmd(2)(1:1))) then
    read(cmd(2),*,iostat=iError) faultID
    if (iError.ne.0) then
        write(*,*) "Usage: fric[tion] <friction ID>"
        write(*,*) "Could not convert argument", trim(cmd(2)), "to an integer"
        return
    endif
else if (cmd(2)(1:3).eq."ids") then
    write(*,*) "Available fault IDs: ", faultIDs
    return
else
    write(*,*) "Usage: fric[tion] <friction ID>"
    write(*,*) "Usage: fric[tion] ids for a list of available friction IDs"
    return
endif


if (ncmd.eq.3) then
    writeToFile = .true.
    WriteToFileName = trim(cmd(3))
else if (ncmd.eq.2) then
    writeToFile = .false.
else
    write(*,*) "Usage: fric[tion] <friction ID>"
endif


IDfound = .false.
do iEntry=1,nFaultSegmentEntriesGlobal
    if (friction(iEntry)%faultNumber.eq.faultID) then
        IDfound = .true.
    endif
enddo

if (.not.IDfound) then
    write(*,*) "fault ID", faultID, "unknown"
    write(*,*) "Available fault IDs: ", faultIDs
    return
endif

if (iecho.gt.0) then
    write(*,*) "calling with fault ID: ", faultID
endif

call frictionTraction(faultID, modeldatactx%DIFORC,nFaultSegmentEntriesGlobal, &
                      NTYPE,NDOF,writeToFile, WriteToFileName)

end subroutine

#endif

!-------------------------------------------------------------------------------


subroutine xsides (cmd,ncmd,ierror)
! prints sides of faults in contour format

use modeldatamodule, only: modeldatactx
use  meshdatamodule, only:  meshdatactx
use modeldefinition, only: NUMSLPglobal

implicit none

character cmd*(*)
dimension cmd(*)
integer :: iError, ncmd


character(len=255)   :: WriteToFileName
logical              :: writeToFile
integer              :: fileID

integer              :: iElem, iPoint, iSlipEntry
integer, allocatable :: elemLookup(:)

integer, external    :: nextlu

if (ncmd.eq.2) then
    writeToFile = .true.
    WriteToFileName = trim(cmd(2))
else if (ncmd.eq.2) then
    writeToFile = .false.
else
    write(*,*) "Usage: mesh informationType [filename]"
endif


! set up looking table for the elements
! for every element, set value to -1, 0, or 1, depending on the weights
! in nslipglobal
allocate(elemLookup(meshdatactx%neglobal))
elemLookup = 0


do iSlipEntry=1,NUMSLPglobal
    ! set weight to x
    elemLookup(modeldatactx%NSLIPglobal(1,iSlipEntry)) = &
               modeldatactx%NSLIPglobal(3,iSlipEntry)
    ! if x weight is 0, maybe y weight.
    if (elemLookup(modeldatactx%NSLIPglobal(1,iSlipEntry)).eq.0) then
        elemLookup(modeldatactx%NSLIPglobal(1,iSlipEntry)) = &
                   modeldatactx%NSLIPglobal(4,iSlipEntry)

    endif
enddo

if (writeToFile) then
    fileID=nextlu(0)
    open(unit=fileID, file=WriteToFileName, iostat=iError)
endif

do iElem=1,meshdatactx%neglobal
    ! and print it in contour format
    do iPoint=1,3
        if (writeToFile) then
            write(fileID,*) meshdatactx%X(1,meshdatactx%IEN(iPoint,iElem)), &
                            meshdatactx%X(2,meshdatactx%IEN(iPoint,iElem)), &
                            elemLookup(iElem)
        else
            write(*,*)      meshdatactx%X(1,meshdatactx%IEN(iPoint,iElem)), &                           
                            meshdatactx%X(2,meshdatactx%IEN(iPoint,iElem)), &
                            elemLookup(iElem)
        endif

    enddo

enddo


close(fileID)

deallocate(elemLookup)


end subroutine xsides
