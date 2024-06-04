subroutine xfemi (cmd,ncmd,ierr)

use algebra,         only: clear
use meshdatamodule,  only: coordsminmax,meshdatactx
use modeldefinition, only: mtot, mbtot, NUMFNGlobal, NUMSLPGlobal
use debugmodule,     only: iecho, parallellog,itest
use errors,          only: exec
use modeldefinition, only: MODE,NUMFN
use iomodule,        only: icoord, ifctr, ifelm, ifeout, &
                           ifein, fein, tempin, luin, stderr
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use modeldatamodule, only: modeldatactx
use filemodule,      only: partinfo
!use vtkmodule,       only: turnOffVelo
use timestepmodule,  only: timestepdatactx
use postProcessData, only: setria
use plotControl,     only: npbased, t0calc, CRDMIN,CRDMAX
!
! Routine that X-ecutes the "femi" command.
!
implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)



!-locl
integer fnblnk
integer, external :: lnblk
!
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
    call closef(luin)
endif

! establish the basic geometric dimensions and element type
call INTCON (NSD,NDOF,NSTR,NEN)


! read in system and mechanical data
if (iecho.eq.1 .and. MODE.gt.2) then
    write(stderr,1) 'reading "', &
     TEMPIN(fnblnk(TEMPIN):lnblk(TEMPIN)),'"'
endif

PARTINFO = "partition.info"

if (iecho.eq.8) then
    call parallelLog("xfemi","calling readfe")
endif

!------------------------------------------------
call READFE(NSD,NDOF,NSTR,NEN,ierr)
!------------------------------------------------

if (iecho.eq.8) then
    call parallelLog("xfemi","finished readfe")
endif
if (ierr.ne.0) then
    ierr = EXEC
    return
endif
if (itest.eq.1) then
    call ELMTEST (meshdatactx%X, meshdatactx%IEN)
    if (iecho.eq.1) write(*,*) 'elements tested OK'
endif


!write(*,*) "bla", timestepdatactx%MAXSTP(1), max(1,NINTG)
!if ((timestepdatactx%MAXSTP(1).eq.0) .and. max(1,NINTG).eq.1) then
!    call turnOffVelo()
!endif
!write(*,*) "velo turned off"

! recompute differential displacements and test input
if (NUMFN.gt.0) then
    call faultequalize (modeldatactx%NFAULT,modeldatactx%FAULT,NUMFN,NDOF)
endif

npbased = (NUMFNGlobal.eq.0 .and. NUMSLPGlobal.eq.0)

! set pointers for mechanical problem arrays
if (iecho.eq.8) then
    call parallelLog("xfemi","calling rdthrm")
endif

! read in thermal data
if (iecho.eq.2) then
    write(*,*) 'read thermal data'
endif

call RDTHRM (T0calc,ierr)

if (iecho.eq.8) then
    call parallelLog("xfemi","finished rdthrm")
endif

if (ierr.ne.0) then
    ierr = EXEC
    return
endif

! set pointers for thermal problem arrays
if (iecho.eq.2) then
    write(*,*) 'set thermal array pointers'
endif

if (iecho.eq.2) then
    write(*,*) 'Transfer material properties to individual elements'
endif
call ELMAT ()

ifein  = 1
ifeout = 0
icoord = 0
ifelm  = 0
ifctr  = 0
setria = .true.
 
! set coordinate extrema
if (iecho.eq.8) then
    call parallelLog("xfemi","minmaxing")
endif

call coordsMinMax (NSD,CRDMIN,CRDMAX)

if (iecho.eq.8) then
    call parallelLog("xfemi","done minmaxing")
endif

if (iecho.eq.8) then
    call parallelLog("xfemi","leaving sub")
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine xReadMesh(cmd,ncmd,ierr)

! read .msg file created by Gmesh

 implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)
!-local
character(len=255) :: meshFileName

if (ncmd.ne.2) then
    write(*,*) "Usage: readmesh [name of mesh file]"
    return
else
    meshFileName = cmd(2)
endif

call readGmsh(meshFileName)

return

end subroutine

!-------------------------------------------------------------------------------

subroutine xTraction(cmd,ncmd,ierr)

use algebra,         only: vmprd, arrayContains
use vtkmodule,       only: outputTags, nOutputTags, gmshRead
use meshdatamodule,  only: meshdatactx
use modeldefinition, only: lgdef
use modeldatamodule, only: modeldatactx
use iomodule,        only: ifein, ifeout, stderr
use errors,          only: exec

implicit none

!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-local
character(len=255) :: outFileName

integer            :: nTags
integer            :: elem1, elem2
integer            :: iTriangle, iTag

integer            :: pos1, pos2, pos3, pos4

double precision   :: fullStress(3,3)
double precision   :: xl(3,4), xn(3), traction(3)
double precision   :: area
integer            :: sideID
integer, external  :: sideIDofTriangle
integer            :: i, j, k
integer            :: lien(3)

if (nOutputTags.eq.0) then
    write(*,*) "traction only available when outputTag is set to non-zero)."
    write(*,*) "set tag [number]"
    return
endif

if (ncmd.ne.2) then
    write(*,*) "Usage: traction traction_file"
    return
else
    outFileName = cmd(2)
endif

if (ifein.eq.0) then
    write(stderr,*) "traction: first read FE-input file"
    ierr = 1
    return
endif

if (ifeout.eq.0) then
    write(stderr,*) "traction: first read FE-output file"
    ierr = EXEC
    return
endif

if (.not. gmshRead) then
    write(stderr,*) "First read mesh file"
!    write(stderr,*) "
    ierr = EXEC
    return
endif


! now we are sure that everything that we need has been read

! filter through the triangles that have the tag, and compute the traction for each
! of those triangles

open(unit=122, file=outFileName)

nTags = size(meshdatactx%triangleTags,1)

do iTriangle = 1,meshdatactx%nTriangles
    do iTag =1,nTags

                call arraycontains(outputTags, &
                           nOutputTags, &
                           meshdatactx%triangleTags(iTag, iTriangle), &
                           pos1)

                if(pos1.le.nOutputTags) then
!        if(meshdatactx%triangleTags(iTag, iTriangle) .eq. outputTag) then
            ! this triangle has the proper tag

                        ! abs to take the slippery weight out
            elem1 = abs(meshdatactx%triangleElements(1,iTriangle))
            elem2 = abs(meshdatactx%triangleElements(2,iTriangle))
!            write(*,*) "triangle ", iTriangle, &
!                       "has tag", outputTag, &!
!                       "and elems: ", elem1, elem2


            ! for now, only elem 1 is filled.

            ! fill the full stress tensor, following labels
            fullStress(1,1) = modeldatactx%STN(1,elem1)
            fullStress(1,2) = modeldatactx%STN(4,elem1)
            fullStress(1,3) = modeldatactx%STN(5,elem1)
            fullStress(2,1) = modeldatactx%STN(4,elem1)
            fullStress(2,2) = modeldatactx%STN(2,elem1)
            fullStress(2,3) = modeldatactx%STN(6,elem1)
            fullStress(3,1) = modeldatactx%STN(5,elem1)
            fullStress(3,2) = modeldatactx%STN(6,elem1)
            fullStress(3,3) = modeldatactx%STN(3,elem1)

            
            ! let  subroutine FCGEOM (XL,LIEN,AREA,XN)
            ! compute the normals
            call LCOORD (meshdatactx%X,xl,meshdatactx%IEN(1,elem1),0)

!            if (LGDEF.ne.0) then
!                call LDISP  (dl,D,IEN(1,k),NDOF,NEN)
!                call AddFaultDisplacement (dl,k,TFAULT,NDOF,NEN,36)
!                call ADDSNE (dl,k,NDOF,NEN,useDXE)
!                call REZONE (xl,dl,1d0)
!            endif


            sideID = sideIDofTriangle(iTriangle, elem1)

            ! compute lien
            call SIDENP(sideID,.true.,lien)

            call FCGEOM (XL,LIEN,AREA,XN)

            call VMPRD(fullStress,XN,traction,3,3)

!            write(122, "(i12,6e14.6)") iTriangle, xn, traction
            i = lien(1)
            j = lien(2)
            k = lien(3)
            ! write(122,*) sideId, lien ! Check sideId and lien output (was wrong previously)
            write(122, "(21e18.9e3)") fullStress(1,1), fullStress(2,2), fullStress(3,3), fullStress(1,2), fullStress(1,3), fullStress(2,3),  xl(1,i), xl(2,i), xl(3,i), xl(1,j), xl(2,j), xl(3,j), xl(1,k), xl(2,k), xl(3,k), xn, traction


            ! Do the same for elem 2 (if it exists)
            if (elem2>0) then
                ! fill the full stress tensor, following labels
                fullStress(1,1) = modeldatactx%STN(1,elem2)
                fullStress(1,2) = modeldatactx%STN(4,elem2)
                fullStress(1,3) = modeldatactx%STN(5,elem2)
                fullStress(2,1) = modeldatactx%STN(4,elem2)
                fullStress(2,2) = modeldatactx%STN(2,elem2)
                fullStress(2,3) = modeldatactx%STN(6,elem2)
                fullStress(3,1) = modeldatactx%STN(5,elem2)
                fullStress(3,2) = modeldatactx%STN(6,elem2)
                fullStress(3,3) = modeldatactx%STN(3,elem2)

                
                ! let  subroutine FCGEOM (XL,LIEN,AREA,XN)
                ! compute the normals
                call LCOORD (meshdatactx%X,xl,meshdatactx%IEN(1,elem2),0)

        !            if (LGDEF.ne.0) then
        !                call LDISP  (dl,D,IEN(1,k),NDOF,NEN)
        !                call AddFaultDisplacement (dl,k,TFAULT,NDOF,NEN,36)
        !                call ADDSNE (dl,k,NDOF,NEN,useDXE)
        !                call REZONE (xl,dl,1d0)
        !            endif


                sideID = sideIDofTriangle(iTriangle, elem2)

                ! compute lien
                call SIDENP(sideID,.true.,lien)

                call FCGEOM (XL,LIEN,AREA,XN)

                call VMPRD(fullStress,XN,traction,3,3)

        !            write(122, "(i12,6e14.6)") iTriangle, xn, traction
                i = lien(1)
                j = lien(2)
                k = lien(3)
                write(122, "(21e18.9e3)") fullStress(1,1), fullStress(2,2), fullStress(3,3), fullStress(1,2), fullStress(1,3), fullStress(2,3),  xl(1,i), xl(2,i), xl(3,i), xl(1,j), xl(2,j), xl(3,j), xl(1,k), xl(2,k), xl(3,k), xn, traction
            endif

        endif
    enddo
enddo

close(122)

end subroutine


!-------------------------------------------------------------------------------
subroutine xfemo (cmd,ncmd,ierr)

use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use filemodule,      only: BaseFedsk, &
                           index_filename
use timestepmodule,  only: maxitime, nstep
use debugmodule,     only: iecho, &
                           parallelLog
use errors,          only: exec, ok, syntax
use iomodule,        only: ifein, ifeout, ascin, fedsk, stderr, &
                           fedskrt, luin

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

integer :: i, ireq, j, k, n
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
subroutine xvtu(cmd,ncmd,ierr)

 ! subroutine to write output to vtk file. 
 ! This is likely not to be the definitive way, but a quick setup.
 ! Also, this version will only support a single piece-file.

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use vtkmodule,       only: WriteVTUfile, &
                           VTUOutputFileName
use iomodule,        only: ifein, stderr
use postProcessData, only: deformed

implicit none

! for stderr

#include "lenchar.i"

 
! for deformed

    
! for EXEC

integer ncmd,ierr
character cmd*(*)
dimension cmd(*)
character(len=LENCHAR) file

if (ifein.eq.0) then
    write(stderr,*) "VTU writer: first read FE-input file"
    ierr = 1
    return
endif

!if (ifeout.eq.0) then
!    write(stderr,*) "VTU writer: first read FE-output file"
!    ierr = EXEC
!    return
!endif


file = VTUOutputFileName
if (ncmd.gt.1) then
    file = cmd(2)
endif

call WriteVTUfile(file, deformed)

return

end subroutine

!-------------------------------------------------------------------------------

subroutine xvtutime(cmd,ncmd,ierr)
 ! subroutine to write timestamps to file
 ! still TODO
implicit none
integer   :: ncmd,ierr
character :: cmd*(*)
dimension :: cmd(*)

write(*,*) "VTU timestepping not yet implemented"

end subroutine

!-------------------------------------------------------------------------------
subroutine xedisp(cmd,ncmd,ierr)

! writes scalar data related to displacements
! needs to be rewritten particularly for slippery and faulted displacements (and
! similar subroutines below)

use modeldefinition, only: NUMFNGlobal, NUMSLPGlobal
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use vtkmodule,       only: outputTags, nOutputTags, connDone
use algebra,         only: arraycontains
use iomodule,        only: stderr,ifeout
use monitorVars,     only: elmon, npmon, monnp, nnp
use plotControl,     only: npbased, idisp
use errors,          only: syntax

implicit none
integer   :: ncmd,ierr
character :: cmd*(*)
dimension :: cmd(*)
integer   :: i, coord, iTriangle, iPoint, iTag, nTags, pos
logical   :: slippery,TagAndFault
double precision :: totalDisp





character(len=255) :: fileName
integer :: nodeID, iProgress
integer, allocatable :: outputTagLookup(:)
integer :: node1, node2, node3
double precision :: thisEdi

if (nOutputTags.ne.0 .and. (NUMSLPGlobal.gt.0 .or. NUMFNGlobal.gt.0)) then
    TagAndFault = .true.
    iDisp = 110 + MOD(iDisp,100)
else
    TagAndFault = .false.
    iDisp = 100 + MOD(iDisp,100)
endif

if (ncmd.lt.2) then
    write(stderr,10)
10  format(1x,'Usage: edi[sp] [total|normal|faulted|slippery] [x|y|z|e] file')
    ierr = SYNTAX
    return
endif

if (ifeout.eq.0) then
    write(stderr,*) "edi: first read FE output file"
    return
endif

i = 2
fileName = ' '

slippery = .false.
coord = -1

100 if (i.le.ncmd) then
    if (cmd(i).eq.'total') then
        if (TagAndFault) then
            write(stderr,101) 'total   '
101         format(1x,'edisp: writing ',A8, &
             ' displacements of tagged element sides not implemented')
            return
        endif
        coord = 0
        iDisp = 100 + MOD(iDisp,10)
    elseif (cmd(i).eq.'normal') then
        iDisp = 110 + MOD(iDisp,10)
    elseif (cmd(i).eq.'slip') then
        if (TagAndFault) then
            write(stderr,101) 'slippery'
            return
        endif
        slippery = .true.
        iDisp = 120 + MOD(iDisp,10)
    elseif (cmd(i).eq.'split') then
        if (TagAndFault) then
            write(stderr,101) 'split   '
            return
        endif
        iDisp = 130 + MOD(iDisp,10)
    elseif (cmd(i).eq.'e') then
        coord = 0
    elseif (cmd(i).eq.'x') then
        coord = 1
        iDisp = INT(iDisp/10)*10 + 1
    elseif (cmd(i).eq.'y') then
        coord = 2
        iDisp = INT(iDisp/10)*10 + 2
    elseif (cmd(i).eq.'z') then
        coord = 3
        iDisp = INT(iDisp/10)*10 + 3
    else
        if (len(trim(fileName)).gt.0) then
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
        fileName = cmd(i)
    endif
    i = i + 1
    goto 100
endif

open(unit=729, file=fileName)

if (nOutputTags.ne.0) then
    ! we need the connectivity to determine whether 
    ! points are needed or not. So check whether that has been done.
    if (.not. connDone) then
        write(stderr,*) "connectivity must be written first"
        ierr = SYNTAX
        close(729)
        return
    endif

    ! set a reference table for all the nodes, to filter out the ones that are part of the 
    ! connectivity. If a node is not part of the filtered connectivity, it is given
    ! coordinates (0,0,0)

!    allocate(outputTagLookup(meshdatactx%nvglobal))
!    outputTagLookup = 0
!    do iTriangle = 1, meshdatactx%nTriangles
!        do iPoint = 1,3
!            outputTagLookup(meshdatactx%triangleConnectivity(iPoint,iTriangle)) = &
!              meshdatactx%triangleConnectivity(iPoint,iTriangle)
!        enddo
!    enddo

    ! number of tags that can be associated with an element side ("triangle")
    nTags = size(meshdatactx%triangleTags,1)

    iProgress = 0

    do iTriangle = 1,meshdatactx%nTriangles
        do iTag=1,nTags

            ! search whether requested iTag appears in the gmsh triangleTags
            call arraycontains(outputTags, &
                               nOutputTags, &
                               meshdatactx%triangleTags(iTag, iTriangle), &
                               pos)

            if (pos.le.nOutputTags) then
                iProgress = iProgress + 1
                node1 = meshdatactx%triangleConnectivity(1,iTriangle)
                node2 = meshdatactx%triangleConnectivity(2,iTriangle)
                node3 = meshdatactx%triangleConnectivity(3,iTriangle)

                thisEdi = sqrt( modeldatactx%D(1,node1)**2 + &
                                modeldatactx%D(2,node1)**2 + &
                                modeldatactx%D(3,node1)**2)

                                 write(729,*) 3*iProgress-2, &
                     meshdatactx%X(:,meshdatactx%triangleConnectivity(1,iTriangle)), &
                                         thisEdi

                thisEdi = sqrt( modeldatactx%D(1,node2)**2 + &
                                modeldatactx%D(2,node2)**2 + &
                                modeldatactx%D(3,node2)**2)

                write(729,*) 3*iProgress-1, &
                     meshdatactx%X(:,meshdatactx%triangleConnectivity(2,iTriangle)), &
                     thisEdi

                thisEdi = sqrt( modeldatactx%D(1,node3)**2 + &
                                modeldatactx%D(2,node3)**2 + &
                                modeldatactx%D(3,node3)**2)

                write(729,*) 3*iProgress , &
                     meshdatactx%X(:,meshdatactx%triangleConnectivity(3,iTriangle)), &
                     thisEdi

            endif
        enddo
    enddo

    close(729)
    return

endif

if (.not.npbased) then
   ! faulted and/or slippery nodes in mesh
   if (npmon) then
        write(stderr,110)
110     format('edisp: element based writing not implemented')
        return
   endif
endif
if (elmon) then
    write(stderr,120)
120 format('edisp: element monitoring not implemented')
    return
endif

if (npmon) then 
    do i=1,nnp ! number of monitored nodes
        nodeID = monnp(i)
        
        if      (coord.eq.1) then
            write(729,*) nodeID, modeldatactx%D(1,nodeID)
        else if (coord.eq.2) then
            write(729,*) nodeID, modeldatactx%D(2,nodeID)
        else if (coord.eq.3) then
            write(729,*) nodeID, modeldatactx%D(3,nodeID)
        else
            totalDisp = sqrt(modeldatactx%D(1,nodeID)**2 + &
                             modeldatactx%D(2,nodeID)**2 + &
                             modeldatactx%D(3,nodeID)**2)
            write(729,*) nodeID, totalDisp
		endif

    enddo

else
    do i=1,meshdatactx%nvglobal


        if      (coord.eq.1) then
            write(729,*) i, modeldatactx%D(1,i)
        else if (coord.eq.2) then
            write(729,*) i, modeldatactx%D(2,i)
        else if (coord.eq.3) then
            write(729,*) i, modeldatactx%D(3,i)
		else
            totalDisp = sqrt(modeldatactx%D(1,i)**2 + &
                             modeldatactx%D(2,i)**2 + &
                             modeldatactx%D(3,i)**2)
            write(729,*) i, totalDisp
		endif

    enddo
endif




close(729)

return

end subroutine xedisp


!-------------------------------------------------------------------------------
subroutine xconn (cmd,ncmd,ierr)

use modeldatamodule, only: modeldatactx, plotdatactx
use meshdatamodule,  only: meshdatactx
use iomodule,        only: ifelm
use vtkmodule,       only: outputTags, nOutputTags, connDone
use modeldefinition, only: NUMSLPglobal
use algebra,         only: arraycontains
use iomodule,        only: ifein, ascout, stderr, elmfil
use monitorVars,     only: elmon, fmon
use errors,          only: exec, ok, syntax
!
! Routine X-excutes the "connections" command
!
implicit none
#include "lenchar.i"
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)




!-locl
 logical ExpandPath
 character(len=LENCHAR) fileName,dir,dirname
external          :: dirname
integer, external :: fnblnk,lnblk, nextlu
integer           :: n
integer           :: iTriangle, iTag, iSlipEntry
logical           :: includeThisOne
integer           :: fileID
integer, allocatable :: slipperyElemLookup(:,:)
integer           :: runningIndex, pos
!-init
 fileName = ' '
 ierr = OK
 ExpandPath = .false.

if (ncmd.lt.2 .and. (.not.elmon)) then
    write(stderr,10)
10        format(1x,'Usage: co[nnections] fileName')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'connections: first read FE-input fileName')
    ierr = EXEC
    return
endif

if (ncmd.gt.1) then
    fileName = cmd(2)
endif

!    Open monitor file if neccessary

! LWI: with the same name as the conn.dat ?????

call OfOpen(fileName,elmon.and.(.not.fmon),ierr)

if (nOutputTags.ne.0) then

    allocate(slipperyElemLookup(2,meshdatactx%neglobal))
    slipperyElemLookup = 0
    ! modeldatactx%NSLIPglobal contains the entries of the slippery nodes.
    ! We need to create a lookup table to match those to the elements.
    ! 1 contains the element, 3 contains the first weight.
    ! 4 has the second weight
    do iSlipEntry =1,NUMSLPglobal
        slipperyElemLookup(1,modeldatactx%NSLIPglobal(1,iSlipEntry)) = &
            modeldatactx%NSLIPglobal(3,iSlipEntry)
        slipperyElemLookup(2,modeldatactx%NSLIPglobal(1,iSlipEntry)) = &
            modeldatactx%NSLIPglobal(4,iSlipEntry)
    enddo

    ! user has requested only the information of a subsection
    ! write only the connectivity of elements in this region.
    fileID=481
    open(unit=fileID, file=fileName)

    runningIndex = 0

!    allocate(tagNodeLookup(3,meshdatactx%nTriangles))

    do iTriangle = 1,meshdatactx%nTriangles
        ! for every triangle...
        includeThisOne = .false.

        do iTag=1,size(meshdatactx%triangleTags,1)
            ! ... walk through the tags of this triangle, and see whether it is included.

            call arraycontains(outputTags, &
                               nOutputTags, &
                               meshdatactx%triangleTags(iTag, iTriangle), &
                               pos)

            if(pos.le.nOutputTags) then

!            if (meshdatactx%triangleTags(iTag, iTriangle) .eq. outputTag) then
                includeThisOne = .true.
            endif
        enddo

        if (includeThisOne) then
            
!            write(fileID,"(5i12)") meshdatactx%triangleConnectivity(:,iTriangle), &
!                                     slipperyElemLookup(:,iTriangle)
           write(fileID,"(5i12)") runningIndex+1, &
                                  runningIndex+2, &
                                  runningIndex+3, &
                                  slipperyElemLookup(:,iTriangle)

!            tagNodeLookup(:,runningIndex) = &
!                 meshdatactx%triangleConnectivity(:,iTriangle)

            runningIndex = runningIndex + 3
           


        endif

    enddo

    deallocate(slipperyElemLookup)
    close(fileID)
    connDone = .true.
    return
endif


! TODO pass of B mabe wrong
!call wrconn (fileName,meshdatactx%ien,B,ierr)
call wrconn (fileName,meshdatactx%ien,plotdatactx%IPT,ierr)

!
if (.not.(elmon.and.(.not.fmon))) then
    if (ierr.eq.OK) then

        if (ascout) then
            ifelm = 2
        else
            ifelm = 1
        endif

        if (ExpandPath) then
!            Expand file name to full path
            dir = dirname(fileName)
            n = fnblnk(dir)
            if (dir(n:n).eq.'/') then
                elmfil = fileName
            else
                call getwdir(dir)
                n = lnblk(dir)
                if (dir(n:n).ne.'/') then
                    dir(n+1:n+1) = '/'
                endif

                elmfil = dir(fnblnk(dir):lnblk(dir)) // &
               fileName(fnblnk(fileName):lnblk(fileName))
            endif
        else
            elmfil = fileName
        endif
    else
        ifelm  = 0
    endif
else
!        Close monitor file again
    call OfClose(fileName,elmon.and.(.not.fmon))
endif
!
connDone = .true.

 return
 end
!-------------------------------------------------------------------------------
subroutine xnp (cmd,ncmd,ierr)

USE MODELDATAMODULE, only: modeldatactx, plotdatactx
USE MESHDATaMODULE,  only: meshdatactx
use iomodule,        only: icoord, ifein, ascout, stderr, crdfil
use debugmodule,     only: iecho, parallelLog
use errors,          only: exec, OK, syntax
use monitorVars,     only: elmon, fmon, npmon
use plotControl,     only: npbased
!
! Routine that X-ecutes the "np" command.
!
implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)





!-locl
logical monitor,ExpandPath
integer, external :: fnblnk, lnblk
#include "lenchar.i"
 character(len=LENCHAR) file,dir,dirname
 external dirname
integer :: n
!-init
 file = ' '
 ierr = OK
 ExpandPath = .false.


if (iecho.eq.8) then
    call parallelLog("f3d; xnp", "entered sub")
endif

!
if (npbased) then
    monitor = npmon
else
    monitor=((.not.fmon).and.elmon)
endif
!
10    format(1x,'Usage: np file')
!
if (ncmd.gt.1) then
    file = cmd(2)
endif

if (lnblk(file).le.0 .and. (.not.monitor)) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'np: first read FE-input file')
    ierr = EXEC
    return
endif
!
!       Open monitor file if neccessary
call OfOpen(file,monitor,ierr)
!

if (iecho.eq.8) then
    call parallelLog("f3d; xnp", "calling nodals")
endif


!if (outputTag.ne.0) then
    ! user has requested only the information of a subsection
    ! write only the coordinates of elements in this region.
    

!    return
!endif


 call nodals (file, &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
 modeldatactx%LMF,modeldatactx%tfault, &
   modeldatactx%DXE,modeldatactx%SKEW, &
   plotdatactx%XPT,0,ierr)

if (iecho.eq.8) then
    call parallelLog("f3d; xnp", "finished nodals")
endif


!
if (.not.monitor) then
    if (ierr.eq.OK) then
if (ascout) then
    icoord = 2
else
    icoord = 1
endif

if (ExpandPath) then
!            Expand file name to full path
    dir = dirname(file)
    n = fnblnk(dir)
    if (dir(n:n).eq.'/') then
                crdfil = file
    else
        call getwdir(dir)
        n = lnblk(dir)
        if (dir(n:n).ne.'/') dir(n+1:n+1) = '/'
        crdfil = dir(fnblnk(dir):lnblk(dir)) // &
               file(fnblnk(file):lnblk(file))
    endif
else
    crdfil = file
endif
    else
        icoord = 0
    endif
else
!           Close monitor file again
    call OfClose(file,monitor)
endif
!
if (iecho.eq.8) then
    call parallelLog("f3d; xnp", "leaving sub")
endif


 return
 end
!-------------------------------------------------------------------------------
subroutine xelmc (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use debugmodule,     only: iecho, parallelLog
use errors,          only: exec, ok, syntax
use iomodule,        only: ifctr, ascout, ifein, stderr
use monitorVars,     only: elmon

!
! Routine that X-ecutes the "elmc" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) :: file, ctrfil
!-init
file = ' '
ierr = OK
!
if (iecho.eq.8) then
    call parallelLog("xelmc","entered")
endif


if (.not.elmon .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: elmc[enter] file')
    ierr = SYNTAX
    goto 1000
endif

if (ncmd.gt.1) then
    file = cmd(2)
endif


if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'elmcenter: first read FE-input file')
    ierr = EXEC
    goto 1000
endif
!
if (ncmd.gt.1) then
    file = cmd(2)
endif
!
if (iecho.eq.8) then
    call parallelLog("xelmc","opening output file")
endif

!       Open monitor file if neccessary
call OfOpen(file,elmon,ierr)
!
if (iecho.eq.8) then
    call parallelLog("xelmc","calling centrs")
endif


!    write center coordinates to file
call centrs (file, &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
        modeldatactx%LMF,modeldatactx%tfault, &
   modeldatactx%DXE,modeldatactx%SKEW, &
   plotdatactx%XPT,ierr)
!
if (.not.elmon) then
    if (ierr.eq.OK) then
        if (ascout) then
            ifctr = 2
        else
            ifctr = 1
        endif
        ctrfil = file
    else
        ifctr = 0
    endif
else
!           Close monitor file again
call OfClose(file,elmon)
endif

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xcf (cmd,ncmd,ierr)
!

use modeltopology, only: nen
use iomodule,      only: stderr, stdin, fcmd, ilu, fstack, maxStack, lucmd
use errors,        only: exec, ok, syntax


! Routine that X-ecutes the "cf" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
logical there
!    integer ftell
integer, external :: lnblk, nextlu, ifopen ! from tools
!-init
ierr = OK
!
if (ncmd.lt.2) then
    write(stderr,1)
 1        format(1x,'Usage: cf file')
    ierr = SYNTAX
    goto 1000
endif

file = cmd(2)
call shvar(file)
inquire (file=file,exist=there)

if (.not.there) then
    write(stderr,2) 'cf: file does not exist :"', &
     cmd(2)(1:lnblk(cmd(2))),'"'
 2        format(1x,80A)
    ierr = EXEC
    goto 1000
endif

if (ilu+1.gt.MaxStack) then
    write(stderr,3) MaxStack
3        format(1x,'cf: nesting depth limit = ',I2)
    ierr = EXEC
    goto 1000
endif

if (lucmd(ilu).ne.stdin) then
    fstack(ilu) = fcmd
!        lupos(ilu)  = ftell(lucmd(ilu))
else
    fstack(ilu) = 'stdin'
endif
!
ilu = ilu + 1
lucmd(ilu) = nextlu(0)
ierr = ifopen(lucmd(ilu),file,'old')

if (ierr.eq.0) then
    fcmd = file
else
    ierr = EXEC
    ilu  = ilu - 1
    fcmd = fstack(ilu)
    if (fcmd.eq.'stdin') then
        fcmd = ' '
    endif
endif
!
1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xmo (cmd,ncmd,ierr)

  USE MODELDATAMODULE
  USE MESHDATAMODULE
  use modeltopology,   only: nen
  use debugmodule,     only: iecho, parallellog
  use errors,          only: exec, syntax
  use iomodule,        only: ifein, stderr
  
  ! Routine that X-ecutes the "mo" command.
  !
  implicit none
  !-pass
  integer ncmd,ierr
  character cmd*(*)
  dimension cmd(*)
  
  if (ncmd.lt.3) then
      write(stderr,10)
10    format(1x,'Usage: mo[nitor] flag[s] number[s]')
      ierr = SYNTAX
      return
  endif

  if (ifein.eq.0) then
      write(stderr,20)
20    format(1x,'monitor: read FE input file first')
      ierr = EXEC
      return
  endif

  if (iecho.eq.8) then
      call parallelLog("xmo","calling setmon")
  endif

  call setmon (cmd, &
             ncmd, &
             meshdatactx%IEN, &
             meshdatactx%neglobal, &
             meshdatactx%nvglobal, &
             ierr)

  if (iecho.eq.8) then
      call parallelLog("xmo","finished setmon")
  endif

  return

end subroutine xmo
!-------------------------------------------------------------------------------
subroutine xnomo (cmd,ncmd,ierr)

  use iomodule,      only: icoord, stderr
  use monitorVars,   only: elmon, npmon, monnp, nnp, nnel, npmon, fmon
  use errors,        only: syntax
  use monitorVars,   only: dealloc_elem_monitor_data,dealloc_node_monitor_data
   
  ! Routine that X-ecutes the "nomo" command.
 
  implicit none
  !-pass
  integer ncmd,ierr
  character cmd*(*)
  dimension cmd(*)
  !-local
  integer :: ifelm

  if (ncmd.lt.2) then
      write(stderr,10)
10    format(1x,'Usage: nomo[nitor] flag'/5X,'with flag one of: n[ode], e[lement]')
      ierr = SYNTAX
      goto 1000
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

1000 return

end subroutine xnomo
!-------------------------------------------------------------------------------
subroutine xedit (cmd,ncmd,ierr)

use iomodule,      only: stderr, stdin, bufedit, fcmd, ilu, fstack, lucmd
use errors,        only: prgerr, ok, syntax


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
logical opned
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
    inquire(unit=lucmd(ilu),opened=opned)
    if (opned) then
        call closef(lucmd(ilu))
    endif
endif

lu = nextlu(0)

if (bufedit) then
    if (cmd(1).eq.'e!') then
        call cptobuf (lu,file,ierr)
        if (ierr.ne.OK) then
            goto 1000
        endif
    else
        if (fstack(0).eq.'stdin') then
            call cptobuf (lu,file,ierr)
            if (ierr.ne.OK) then
                goto 1000
            endif
        else
            write(stderr,4)
 4               format(1x,'File has been changed: use "e!" to',1x, &
             'overwrite')
            ierr = SYNTAX
            goto 1000
        endif
    endif
endif

inquire(unit=lu,opened=opned)
if (opned) then
    call closef(lu)
endif

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
end subroutine
!-------------------------------------------------------------------------------
subroutine xeceof (cmd,ncmd,ierr)

use iomodule,      only: stderr, stdin, fcmd, ilu, lucmd, fstack
use errors,        only: eof, ok, stop, prgerr

! Routine that X-ecutes the "eof" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)



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
end subroutine
!-------------------------------------------------------------------------------
subroutine xecq (cmd,ncmd,ierr)

use iomodule,      only: stderr, stdin, bufedit, fstack
use errors,        only: eof, ok, stop

! Routine that X-ecutes "q" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)


!-init
ierr = OK

!    close current command file and divert to previous one
call xeceof (cmd,ncmd,ierr)
!    currently at lowest level of command files stack?
if (ierr.ne.STOP) then
    goto 10
endif
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
end subroutine
!-------------------------------------------------------------------------------
subroutine xecw (cmd,ncmd,ierr)

use iomodule,       only: stderr, stdin, fstack, bufedit, lucmd, ilu
use errors,         only: ok, syntax, prgerr

! Routine that X-ecutes "w" command.

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
 3        format(1x,'write: no file')
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
end subroutine
!-------------------------------------------------------------------------------
subroutine xecr (cmd,ncmd,ierr)

USE ALGEBRA
USE MODELDEFINITION, only: MBTOT
use debugmodule,     only: iecho
use errors,          only: ok, prgerr, syntax
use iomodule,        only: stderr, stdin, fstack, bufedit, lucmd, ilu, fcmd

! Routine that X-ecutes the "restart" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)


!-locl
#include "lenchar.i"
character(len=LENCHAR) file
logical opned
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

inquire(unit=lu,opened=opned)
if (opned) then
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

use iomodule,      only: stderr, stdin, fstack, bufedit, lucmd, ilu
use errors,        only: ok, prgerr, syntax

! Routine that X-ecutes the "history" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)
!-locl
logical numeric
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
    nhist = NINT(chreal(cmd(2)))
else
    nhist = -1
endif

call bufcat (nhist,ierr)

1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xecu (cmd,ncmd,ierr)

use iomodule,      only: stderr, stdin, fstack, bufedit, lucmd, ilu
use errors,        only: ok, syntax, prgerr

! Routine that X-ecutes the "undo" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)


!-locl
 logical numeric
integer, external :: lnblk
real, external :: chreal
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
    nhist = NINT(chreal(cmd(2)))
else
    nhist = 1
endif

call bufudo (nhist,ierr)
1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xcp (cmd,ncmd,ierr)

use debugmodule,     only: iecho
use errors,          only: exec, OK, syntax
use modeldefinition, only: mode
use timestepmodule,  only: nmprt, ntprt, nintg
use iomodule,        only: ifein, luin, lufo, ifeout, ascin, ascout, &
                           stderr, fedsk, stdout

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
integer :: i, iend, ifeold, lufold, luiold, ireq
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

    call FOUTCP (Ireq,Iend,ierr)
    if (ierr.ne.OK) then
        if (iecho.eq.1) then
            write(stdout,'(/)')
            call flushf(stdout)
        endif
        write(stderr,5)
 5        format(1x,'cp: failed')
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
subroutine xgrid (cmd,ncmd,ierr)

use modeldefinition, only: numnp, numel, iword
use meshdatamodule,  only: meshdatactx
use formatsmodule
use iomodule,        only: icoord, ifelm, icoord, ifein, stderr
use modeltopology,   only: NEN
use monitorVars,     only: fmon, nnel
use plotControl,     only: npbased
use errors,          only: exec, syntax

!
! Writes a DataExplorer data file
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
integer, external :: lnblk, nextlu, ifopen
character(len=6) ftype
integer :: num
integer :: lu
!
if (ncmd.le.1) then
    write(stderr,10)
10        format(1x,'Usage: gr[id] file')
    ierr = SYNTAX
    return
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'grid: first read FE-input file')
    ierr = EXEC
    return
endif
!
if (icoord.eq.0) then
    write(stderr,30)
30        format(1x,'grid: first write nodal points file with "np"')
    ierr = EXEC
    return
endif
!
if (ifelm.eq.0) then
    write(stderr,40)
40        format(1x,'grid: first write connections file with "conn"')
    ierr = EXEC
    return
endif
!
lu = nextlu(0)
ierr = ifopen(lu,cmd(2),'unknown')

if (ierr.ne.0) then
    call closef(lu)
    ierr = EXEC
    return
endif

if (icoord.eq.1) then
    ftype='binary'
else
    ftype='ascii '
endif

if (fmon) then
    num = nnel*4
else
    if (npbased) then
        num = NUMNP
    else
        num = NUMEL*NEN
    endif
endif


!    write(lu,1) NDOF,num,ftype,crdfil(1:lnblk(crdfil))
!    1   format('# nodal points:'/'object 1 class array type float',
!     >   ' rank 1 shape ',I1,' items ',i8,1X,A6/'data file ',256A)

if (ifelm.eq.1) then
    ftype='binary'
else
    ftype='ascii '
endif

!    if (fmon) then
!        write(lu,2) 4,nnel,ftype,elmfil(1:lnblk(elmfil))
!    else
!        write(lu,2) NEN,NUMEL,ftype,elmfil(1:lnblk(elmfil))
!    endif

 2   format(/'# connectivity:'/'object 2 class array type int', &
 ' rank 1 shape ',I1,' items ',i8,1X,A6/'data file ',256A)

!    if (NTYPE.eq.5) then
!        if (fmon) then
!            write(lu,3) '"quads"'
!        else
!            write(lu,3) '"cubes"'
!        endif
!    else
!        write(lu,3) '"tetrahedra"'
!    endif

 3   format('attribute "element type" string ',A12/ &
 'attribute "ref" string "positions"')

!    if (fmon) then
!        if (icoord.eq.1) then
!               ftype='binary'
!          else
!                ftype='ascii '
!            endif
!            write(lu,4) NDOF,nnel*4,ftype,crdfil(1:lnblk(crdfil))
 4       format(/'# data:'/'object 3 class array type float', &
     ' rank 1 shape ',I1,' items ',i8,1X,A6/'data file ',256A)
!        write(lu,5)
 5       format('attribute "dep" string "positions"'// &
     'object', &
     ' "irregular positions irregular connections" class',1x, &
     'field'/ &
     'component "positions" value 1'/ &
     'component "connections" value 2'/ &
     'component "data" value 3'/ &
     'end')
!    else
!        write(lu,6)
 6       format(//'object', &
     ' "irregular positions irregular connections" class',1x, &
     'field'/ &
     'component "positions" value 1'/ &
     'component "connections" value 2'/ &
     'end')
!        endif


 call closef(lu)


 return
 end
!-------------------------------------------------------------------------------
subroutine xestress (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, ifein, ifeout, stderr, mout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax

!
! Routine that X-ecutes the "estre" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
character(len=2) flag
integer, external :: lnblk
save flag
integer :: i
!-init
data flag/'e'/
ierr = OK
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: estre[ss] [flag] file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'estre: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'estre: no stresses calculated')
    else
        write(stderr,40)
40        format(1x,'estre: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'estre: no mechanical output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'estre: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'estre: first write connections file with "conn"')
    ierr = EXEC
    return
endif
!
i=2

100    if (i.le.ncmd) then
    if (cmd(i).eq.'e '.or.cmd(i)(1:3).eq.'eff') then
        flag = 'e'
    else if (cmd(i).eq.'xx') then
        flag = 'xx'
    else if (cmd(i).eq.'yy') then
        flag = 'yy'
    else if (cmd(i).eq.'zz') then
        flag = 'zz'
    else if (cmd(i).eq.'xy'.or.cmd(i).eq.'yx') then
        flag = 'xy'
    else if (cmd(i).eq.'xz'.or.cmd(i).eq.'zx') then
        flag = 'xz'
    else if (cmd(i).eq.'yz'.or.cmd(i).eq.'zy') then
        flag = 'yz'
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif
    i = i + 1
    goto 100
endif
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

!    write effective stresses file(s)
call wreffs (file,modeldatactx%STN,plotdatactx%XPT,'stress',flag,ierr)
!
!    Close monitor file again
call OfClose(file,elmon)
!
return
end subroutine

!-------------------------------------------------------------------------------

subroutine xeffr (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: ifelm, icoord, ifein, ifeout, mout, stderr
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax

!
! Routine that X-ecutes the "erate" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file       
character(len=2) flag
integer, external :: lnblk
save flag
integer :: i
!-init
data flag/'e'/
ierr = OK
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: era[te] [flag] file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'erate: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'erate: no strain rates calculated')
    else
        write(stderr,40)
40        format(1x,'erate: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'erate: no mechanical output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'erate: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'erate: first write connections file with "conn"')
    ierr = EXEC
    return
endif
!
i=2
100     if (i.le.ncmd) then
    if (cmd(i).eq.'e '.or.cmd(i)(1:3).eq.'eff') then
        flag = 'e'
    else if (cmd(i).eq.'xx') then
        flag = 'xx'
    else if (cmd(i).eq.'yy') then
        flag = 'yy'
    else if (cmd(i).eq.'zz') then
        flag = 'zz'
    else if (cmd(i).eq.'xy'.or.cmd(i).eq.'yx') then
        flag = 'xy'
    else if (cmd(i).eq.'xz'.or.cmd(i).eq.'zx') then
        flag = 'xz'
    else if (cmd(i).eq.'yz'.or.cmd(i).eq.'zy') then
        flag = 'yz'
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif
    i = i + 1
    goto 100
endif
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

!    write effective strain rate file(s)
#ifndef EXPERIMENTAL_NORATESTORAGE
 call wreffs (file,modeldatactx%rate,plotdatactx%XPT,'strain rate',flag,ierr)
#endif
!
!    Close monitor file again
call OfClose(file,elmon)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xstress (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, stderr, ifeout, ifein, mout
use meshdatamodule,  only: meshdatactx
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use plotControl,     only: Sscale
use errors,          only: ok, exec, syntax        

!
! Routine that X-ecutes the "stress" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)


!-locl
#include "lenchar.i"
character(len=LENCHAR) fileName
logical tecton
integer, external :: lnblk
integer :: i

integer :: iTriangle

!-init
ierr = OK
fileName = ' '
tecton = .false.
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: stre[ss] file [tecton]')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'stress: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'stress: no stresses calculated')
    else
        write(stderr,40)
40        format(1x,'stress: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'stress: no mechanical output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'stress: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'stress: first write connections file',1x, &
   'with "conn"')
    ierr = EXEC
    return
endif
!
i = 2
do while (i.le.ncmd)
    if (cmd(i).eq.'tecton') then
        tecton = .true.
    else if (lnblk(fileName).le.0) then
        fileName = cmd(i)
    else
        write(stderr,10)
        ierr = SYNTAX
        return
    endif
    i = i + 1
enddo

if ((lnblk(fileName).le.0.and.(.not.elmon)).or. &
    (lnblk(fileName).le.0.and.tecton)) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif


!    Open monitor file if neccessary
call OfOpen(fileName,elmon.and.(.not.tecton),ierr)

!    write stresses file(s)
call wrtens(fileName,modeldatactx%STN,plotdatactx%XPT,'stress',Sscale,tecton,ierr)

!    Close monitor file again
call OfClose(fileName,elmon.or.tecton)
!
return
end subroutine

!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP
subroutine xtaylorstress (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, stderr, ifeout, ifein, mout
use meshdatamodule,  only: meshdatactx
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use plotControl,     only: Sscale
use errors,          only: ok, exec, syntax  
use modeltopology,   only: NSTRTAY      

!
! Routine that X-ecutes the "taylorstress" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)


!-locl
#include "lenchar.i"
character(len=LENCHAR) fileName
logical tecton
integer, external :: lnblk
integer :: i



integer :: iTriangle

!-init
ierr = OK
fileName = ' '
tecton = .false.
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: taylorstress file [tecton]')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'taylorstress: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'taylorstress: no stresses calculated')
    else
        write(stderr,40)
40        format(1x,'taylorstress: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'taylorstress: no mechanical output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'taylorstress: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'taylorstress: first write connections file',1x, &
   'with "conn"')
    ierr = EXEC
    return
endif
!
i = 2
do while (i.le.ncmd)
    if (cmd(i).eq.'tecton') then
        tecton = .true.
    else if (lnblk(fileName).le.0) then
        fileName = cmd(i)
    else
        write(stderr,10)
        ierr = SYNTAX
        return
    endif
    i = i + 1
enddo

if ((lnblk(fileName).le.0.and.(.not.elmon)).or. &
    (lnblk(fileName).le.0.and.tecton)) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif


!    Open monitor file if neccessary
 call OfOpen(fileName,elmon.and.(.not.tecton),ierr)

!    write taylorstresses file(s)
if (NSTRTAY.eq.1) then
! taylor stress is a scalar value

  !  do i=1,meshdatactx%neglobal
!       XPT(i) = REAL(modeldatactx%STNTAY(i))
 !   enddo
   ! write(stderr,*) modeldatactx%STNTAY
    call wrcval (fileName,(modeldatactx%STNTAY),'taylorstress',ierr)

elseif (NSTRTAY.eq.6) then
! full tensor
    write(stderr,*) 'no tensor output implemented yet, need a double precision version for writing tensors'
!   call wrtens(fileName,modeldatactx%STNTAY,plotdatactx%XPT,'taylorstress',Sscale,tecton,ierr)
else
    write(stderr,80)
80        format(1x,'taylorstress: only scalar and full tensor taylor stress supported, pleasexest contact model support')
    ierr = EXEC
    return
endif

!    Close monitor file again
call OfClose(fileName,elmon.or.tecton)
!
return
end subroutine
#endif
!-------------------------------------------------------------------------------
subroutine xestrain (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, stderr, ifein, ifeout, mout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon                   
use errors,          only: ok, exec, syntax        

!
! Routine that X-ecutes the "estra" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file       
character(len=2) flag
integer, external :: lnblk
save flag
integer :: i
!-init
data flag/'e'/
ierr = OK
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: estra[in] [flag] file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'estra: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'estra: no strains calculated')
    else
        write(stderr,40)
40        format(1x,'estra: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'estra: no mechanical output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'estra: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'estra: first write connections file with "conn"')
    ierr = EXEC
    return
endif
!
i=2
100     if (i.le.ncmd) then
    if (cmd(i).eq.'e '.or.cmd(i)(1:3).eq.'eff') then
        flag = 'e'
    else if (cmd(i).eq.'xx') then
        flag = 'xx'
    else if (cmd(i).eq.'yy') then
        flag = 'yy'
    else if (cmd(i).eq.'zz') then
        flag = 'zz'
    else if (cmd(i).eq.'xy'.or.cmd(i).eq.'yx') then
        flag = 'xy'
    else if (cmd(i).eq.'xz'.or.cmd(i).eq.'zx') then
        flag = 'xz'
    else if (cmd(i).eq.'yz'.or.cmd(i).eq.'zy') then
        flag = 'yz'
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif
    i = i + 1
    goto 100
endif
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

!    write effective strain file
call wreffs (file,modeldatactx%STR,plotdatactx%XPT,'strain',flag,ierr)
!
!    Close monitor file again
call OfClose(file,elmon)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xstrain (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, ifein, ifeout, stderr, mout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon                   
use plotControl,     only: Escale
use errors,          only: ok, exec, syntax        

!
! Routine that X-ecutes the "strain" command.
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
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: stra[in] file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'strain: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'strain: no strains calculated')
    else
        write(stderr,40)
40        format(1x,'strain: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'strain: no mechanical output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'strain: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'strain: first write connections file',1x, &
   'with "conn"')
    ierr = EXEC
    return
endif

if (ncmd.ge.2) then
    file = cmd(2)
endif
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

!    write strain file(s)
call wrtens (file,modeldatactx%STR,plotdatactx%XPT,'strain',Escale,.false.,ierr)
!
!    Close monitor file again
call OfClose(file,elmon)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xrate (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, ifeout, ifein, stderr, mout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon                   
use plotControl,     only: Escale
use errors,          only: ok, exec, syntax        

!
! Routine that X-ecutes the "rate" command.
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
!-init
ierr = OK
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: ra[te] file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'rate: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'rate: no strains calculated')
    else
        write(stderr,40)
40        format(1x,'rate: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'rate: no mechanical output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'rate: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'rate: first write connections file',1x, &
   'with "conn"')
    ierr = EXEC
    return
endif

if (ncmd.ge.2) then
    file = cmd(2)
endif
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

!    write rate file(s)
#ifndef EXPERIMENTAL_NORATESTORAGE
call wrtens (file,modeldatactx%rate,plotdatactx%XPT,'rate',Escale,.false.,ierr)
#endif
!
!    Close monitor file again
call OfClose(file,elmon)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xte (cmd,ncmd,ierr)
!
USE MODELDEFINITION, only: MODE, ITMODE
use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, tout, stderr, ifeout, ifein
use timestepmodule,  only: nstep
use monitorVars,     only: elmon, npmon
use plotControl,     only: itemp, npbased
use errors,          only: exec, syntax        

! Routine that X-ecutes the "te" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)


!-locl
logical tecton,monitor
#include "lenchar.i"
character(len=LENCHAR) file
integer, external :: lnblk, nextlu, ifopen
integer :: i, lu

!-init
file = ' '
tecton = .false.
i = 2

if (npbased) then
    monitor = npmon
else
    monitor = elmon
endif
!
10    format(1x,'Usage: te[mperature] [tecton] file')
!
100    if (i.le.ncmd) then
    if (cmd(i).eq.'tecton') then
        tecton = .true.
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif
    i = i + 1
    goto 100
endif

if ((.not.monitor) .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'temp: read FE input file first')
    ierr = EXEC
    return
endif
!
if (ifeout.eq.0) then
    if (ITMODE.ne.0 .and. MODE.ge.3) then
        write(stderr,30)
30        format(1x,'temperatures: write temperature INPUT data')
    else
        if (MODE.ne.5 .and. MODE.ne.6) then
            write(stderr,40)
40            format(1x,'temp: no temperatures calculated')
            ierr = EXEC
            return
        else
            write(stderr,50)
50            format(1x,'temp: first read FE output file')
            ierr = EXEC
            return
        endif
    endif
else
    if (MODE.ne.5 .and. MODE.ne.6) then
        write(stderr,40)
        ierr = EXEC
        return
    endif
    if (.not.Tout) then
        write(stderr,60) NSTEP
60        format(1x,'temp: no thermal output at itime = ',I5)
        ierr = EXEC
        return
    endif
endif
!
if (tecton) then
    lu = nextlu(0)
    ierr = ifopen(lu,file,'unknown')
    if (ierr.ne.0) then
        ierr = EXEC
        call closef(lu)
        return
    endif
    call tecout (modeldatactx%T,lu)
    call closef(lu)
    return
endif
!
if (icoord.eq.0.and.(.not.monitor)) then
    write(stderr,70)
70       format(1x,'temp: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0.and.(.not.monitor)) then
    write(stderr,80)
80       format(1x,'temp: first write connections file',1x, &
     'with "conn"')
    ierr = EXEC
    return
endif
!
!       Open monitor file if neccessary
call OfOpen(file,monitor,ierr)
!
if (itemp.eq.0) then
    call temprs (file,meshdatactx%IEN,modeldatactx%T,modeldatactx%TLST,modeldatactx%TX, &
       modeldatactx%TLSTX,modeldatactx%LMTX,plotdatactx%XPT,0,ierr)
else if (itemp.eq.1) then
    call temprs (file,meshdatactx%IEN,modeldatactx%T,modeldatactx%TLST,modeldatactx%TX, &
       modeldatactx%TLSTX,modeldatactx%LMTX,plotdatactx%XPT,1,ierr)
else if (itemp.eq.2) then
! todo check this
! NOT WORKING: TX0 array not defined
    call temprs (file,meshdatactx%IEN,modeldatactx%T,modeldatactx%T0,modeldatactx%TX, &
       modeldatactx%TLSTX,modeldatactx%LMTX,plotdatactx%XPT,2,ierr)
endif
!
call OfClose (file,monitor)
!
return
end
!-------------------------------------------------------------------------------
subroutine xdispl (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMFNGlobal, numel, NUMSLPGlobal
use modeldatamodule, only: modeldatactx, plotdatactx
use meshdatamodule,  only: meshdatactx
use iomodule,        only: icoord, ifelm, ifein, ifeout, stderr, mout
use modeltopology,   only: ndof, nen
use vtkmodule,       only: outputTags, nOutputTags, connDone
use algebra,         only: arraycontains,clear
use constants,       only: useDXE,useDELDXE
use timestepmodule,  only: nstep
use postProcessData, only: deffac,deformed
use monitorVars,     only: elmon, npmon
use plotControl,     only: idisp
use errors,          only: exec, syntax

! Routine that X-ecutes the "di" command.
 
implicit none

!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)





!-locl
#include "lenchar.i"
character(len=LENCHAR) fileName
logical monitor
integer, external :: lnblk,nextlu
integer :: i,j,iTag, nTags, iProgress, iTriangle, pos, iPoint
integer :: point1, point2, point3
integer :: nenPoint1, nenPoint2, nenPoint3
integer :: elem1, elem2
integer :: sideID, elemID, fileID,lien(3)
integer, external :: sideIDofTriangle
logical, parameter :: FaultWeightPositive=.false.
double precision :: XL(3,4),DL(3,4)

!-init
fileName = ' '
i = 2
iDisp = 10 + MOD(iDisp,10)

10   format(1x,'Usage: di[splacements] file',1x, &
 '[total|normal|slippery|faulted]')

do while (i.le.ncmd) 
    if (cmd(i).eq.'total') then
        iDisp = 10
    else if (cmd(i).eq.'normal') then
        iDisp = 11
    else if (cmd(i).eq.'slippery') then
        iDisp = 12
    else if (cmd(i).eq.'split') then
        iDisp = 13
    else if (cmd(i).eq.'faulted') then
        iDisp = 13
    else
        if (lnblk(fileName).le.0) then
            fileName = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif
    i = i + 1
enddo


if (nOutputTags.ne.0) then

    ! process outputs for user requested output for specific GMSH tags

    fileID = nextlu(0)
    open(unit=fileID, file=fileName)

    if (.not. connDone) then
        write(stderr,*) "connectivity must be written first"
        ierr = SYNTAX
        return
    endif

    ! set a reference table for all the nodes, to filter out the ones that are part of the 
    ! connectivity. If a node is not part of the filtered connectivity, it is given
    ! coordinates (0,0,0)

    nTags = size(meshdatactx%triangleTags,1)
    iProgress = 0
    do iTriangle = 1,meshdatactx%nTriangles
        do iTag=1,nTags

            call arraycontains(outputTags, &
                               nOutputTags, &
                               meshdatactx%triangleTags(iTag, iTriangle), &
                               pos)

            if (pos.le.nOutputTags) then

                ! This is a triangle (element side) with a requested tag.
                ! To write the slippery or faulted displacements, we need to know which element is on which side.
                point1 = meshdatactx%triangleConnectivity(1,iTriangle)  ! global node number
                point2 = meshdatactx%triangleConnectivity(2,iTriangle)
                point3 = meshdatactx%triangleConnectivity(3,iTriangle)

                elem1 = meshdatactx%triangleElements(1,iTriangle)       ! global element number
                elem2 = meshdatactx%triangleElements(2,iTriangle)

                iProgress = iProgress + 1
#ifdef OBSOLETE
                ! the conn runs from 1,2,3   4,5,6   etc.
                ! go back to displcaments from real global node IDs
                ! using the triangle connectivity

                ! triangle connectivity does not necessarily correspond to the
                ! the point sequence of NEN in modeldatactx%DELDXE(NDOF,NEN,neglobal)
#endif
                if (NUMSLPGlobal.gt.0 .or. NUMFNGlobal.gt.0) then

!                   Check whether to write elem1 or elem2, depending on minus determined in readmesh.f
                    if (FaultWeightPositive) then
                        ! Pick the element with positive weight.
                        if (elem2.eq.0 .or. elem1.gt.0) then
                            sideID = sideIDofTriangle(iTriangle, abs(elem1))
                            elemID = abs(elem1)
                        else
                            sideID = sideIDofTriangle(iTriangle, abs(elem2))
                            elemID = abs(elem2)
                        endif
                    else
                        ! Pick the element with negative weight.
                        if (elem2.eq.0 .or. elem1.lt.0) then
                            ! absolute values because minus signs are carried 
                            ! along to indicate slippery weights, to save memory.
                            sideID = sideIDofTriangle(iTriangle, abs(elem1))
                            elemID = abs(elem1)
                        else
                            sideID = sideIDofTriangle(iTriangle, abs(elem2))
                            elemID = abs(elem2)
                        endif
                    endif

                    ! get element coordinates into XL
                    call LCOORD (meshdatactx%X,XL,meshdatactx%IEN(1,elemID),0)
                    if (deformed) then
                        call LDISP (DL,modeldatactx%D,meshdatactx%IEN(1,elemID),NDOF,NEN)
                        call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
                        call AddFaultDisplacement (DL,elemID,modeldatactx%TFAULT,NDOF,NEN,36)
                        call REZONE (XL,DL,deffac)
                    endif

                    ! get element displacements into DL
                    if (iDisp.eq.10 .or. iDisp.eq.11) then
                        ! normal or total displacement
                        call LDISP (DL,modeldatactx%D,meshdatactx%IEN(1,elemID),NDOF,NEN)
                    else
                        call CLEAR (DL,NDOF*NEN,"DL")
                    endif
                    if (iDisp.eq.10 .or. iDisp.eq.12) then
                        ! slippery or total displacement
                        call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
                    endif
                    if (iDisp.eq.10 .or. iDisp.eq.13) then
                        ! faulted or total displacement
                        call AddFaultDisplacement (DL,elemID,modeldatactx%TFAULT,NDOF,NEN,36)
                    endif

                    ! Get element-local point numbers of the element side that makes up this triangle
                    call SIDENP(SideID,.true.,lien)
#ifdef OBSOLETE
!                   ! triangle connectivity does not necessarily correspond to the
!                   ! the point sequence of NEN in modeldatactx%DELDXE(NDOF,NEN,neglobal)
!                   ! We have to make sure that the proper differential displacement
!                   ! is added to the proper node.
!                   !  | \            face 1: nodes 1 2 4
!                   !  |__\2_____X         2:       1 3 2
!                   !  1\ /                3:       1 4 3
!                   !    \4                4:       2 3 4
!                   do iPoint=1,4
!                       if (meshdatactx%ien(iPoint,elemID).eq.point1) then
!                           nenPoint1 = iPoint
!                       endif
!                       if (meshdatactx%ien(iPoint,elemID).eq.point2) then
!                           nenPoint2 = iPoint
!                       endif
!                       if (meshdatactx%ien(iPoint,elemID).eq.point3) then
!                           nenPoint3 = iPoint
!                       endif
!                   enddo
!
!                   write(fileID,*) 3*iProgress-2,meshdatactx%X(:,point1), &
!                    modeldatactx%D(:,point1), &
!                    modeldatactx%DXE(:,nenPoint1,elemID)
!
!                   write(fileID,*) 3*iProgress-1, &
!                    meshdatactx%X(:,point2), &
!                    modeldatactx%D(:,point2), &
!                    modeldatactx%DXE(:,nenPoint2,elemID)

!                   write(fileID,*) 3*iProgress , &
!                    meshdatactx%X(:,point3), &
!                    modeldatactx%D(:,point3), &
!                    modeldatactx%DXE(:,nenPoint3,elemID)
#endif
                    write(fileID,*) 3*iProgress-2,XL(:,lien(1)),DL(:,lien(1))  
                    write(fileID,*) 3*iProgress-1,XL(:,lien(2)),DL(:,lien(2))  
                    write(fileID,*) 3*iProgress-0,XL(:,lien(3)),DL(:,lien(3))  
                else
                    ! no slippery and faulted nodes in the model: nodal point-centered data model.
                    do j=1,3
                        do i=1,NDOF
                            XL(i,j) = meshdatactx%X(i,meshdatactx%triangleConnectivity(j,iTriangle))
                            DL(i,j) = modeldatactx%D(i,meshdatactx%triangleConnectivity(j,iTriangle))
                        enddo
                    enddo
                    if (deformed) then
                        call REZONE (XL,DL,deffac)
                    endif
                    write(fileID,*) 3*iProgress-2,XL(:,1),DL(:,1)
                    write(fileID,*) 3*iProgress-1,XL(:,1),DL(:,1)
                    write(fileID,*) 3*iProgress-0,XL(:,1),DL(:,1)
                endif

            endif
        enddo
    enddo

    close(629)
    return

endif

monitor = (npmon.or.elmon)

if (lnblk(fileName).le.0 .and. (.not.monitor)) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20  format(1x,'displ: read FE input and output files first')
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30  format(1x,'displ: no displacements calculated')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    write(stderr,40)
40  format(1x,'displ: first read FE output file')
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50  format(1x,'displ: no mechanical output at itime = ',I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0.and.(.not.monitor)) then
    write(stderr,60)
60  format(1x,'displ: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0.and.(.not.monitor)) then
    write(stderr,70)
70  format(1x,'displ: first write connections file',1x, &
     'with "conn"')
    ierr = EXEC
    return
endif

!   Open monitor file if neccessary
call OfOpen(fileName,monitor,ierr)

! RG 2021 if (NUMFNGlobal.gt.0) then
! RG 2021     call LOADFT (meshdatactx%ien,&
! RG 2021                  modeldatactx%NFAULT,&
! RG 2021                  modeldatactx%FAULT,&
! RG 2021                  modeldatactx%TFAULT,&
! RG 2021                  NUMEL)
! RG 2021 endif

call wrdisp (fileName,meshdatactx%IEN,modeldatactx%D,modeldatactx%LMF,&
             modeldatactx%tfault,modeldatactx%DXE, &
             plotdatactx%XPT,ierr)

call OfClose (fileName,monitor)
 
return
end
!-------------------------------------------------------------------------------
 subroutine xvelo (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMFNGlobal, numel, NUMSLPGlobal
use modeldatamodule, only: modeldatactx, plotdatactx
use meshdatamodule,  only: meshdatactx
use iomodule,        only: icoord, ifelm, ifein, ifeout, stderr, mout
use modeltopology,   only: ndof, nen
use vtkmodule,       only: outputTags, nOutputTags, connDone
use algebra,         only: arraycontains,clear
use constants,       only: useDXE,useDELDXE
use timestepmodule,  only: nstep
use postProcessData, only: deffac,deformed
use monitorVars,     only: elmon, npmon
use plotControl,     only: idisp
use errors,          only: exec, syntax

! Routine that X-ecutes the "ve" command.

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)





!-locl
#include "lenchar.i"
character(len=LENCHAR) fileName
logical monitor
integer, external :: lnblk,nextlu
integer :: i,j,iTag, nTags, iProgress, iTriangle, pos, iPoint
integer :: point1, point2, point3
integer :: nenPoint1, nenPoint2, nenPoint3
integer :: elem1, elem2
integer :: sideID, elemID, fileID,lien(3)
integer, external :: sideIDofTriangle
logical, parameter :: FaultWeightPositive=.false.
double precision :: XL(3,4),DL(3,4),VL(3,4)

!-init
fileName = ' '
i = 2
monitor = (npmon.or.elmon)
iDisp = 20 + MOD(iDisp ,10)
 
10   format(1x,'Usage: ve[locity] fileName',1x, &
 '[total|normal|slippery|faulted]')
 
do while (i.le.ncmd) 
     if (cmd(i).eq.'total') then
          iDisp = 20
     else if (cmd(i).eq.'normal') then
          iDisp = 21
     else if (cmd(i).eq.'slippery') then
          iDisp = 22
     else if (cmd(i).eq.'faulted') then
          iDisp = 23
     else if (cmd(i).eq.'split') then
          iDisp = 23
     else
          if (lnblk(fileName).le.0) then
              fileName = cmd(i)
         else
              write(stderr,10)
              ierr = SYNTAX
              return
         endif
     endif
     i = i + 1
enddo

if (nOutputTags.ne.0) then

    fileID = nextlu(0)
    open(unit=fileID, file = fileName)

    ! we need the connectivity to determine whether 
    ! points are needed or not. So check whether that has been done.
    if (.not. connDone) then
        write(stderr,*) "connectivity must be written first"
        ierr = SYNTAX
        close(fileID)
        return
    endif

    nTags = size(meshdatactx%triangleTags,1)

    iProgress = 0

    do iTriangle = 1,meshdatactx%nTriangles
        do iTag=1,nTags

            call arraycontains(outputTags,nOutputTags, &
             meshdatactx%triangleTags(iTag, iTriangle),pos)

            if (pos.le.nOutputTags) then

                ! get global node and element number belonging to this element face
                point1 = meshdatactx%triangleConnectivity(1,iTriangle)
                point2 = meshdatactx%triangleConnectivity(2,iTriangle)
                point3 = meshdatactx%triangleConnectivity(3,iTriangle)

                elem1 = meshdatactx%triangleElements(1,iTriangle)
                elem2 = meshdatactx%triangleElements(2,iTriangle)

                iProgress = iProgress + 1

                if (NUMSLPGlobal.gt.0 .or. NUMFNGlobal.gt.0) then

!                   Check whether to write elem1 or elem2, depending on minus determined in readmesh.f
                    if (FaultWeightPositive) then
                        ! Pick the element with positive weight.
                        if (elem2.eq.0 .or. elem1.gt.0) then
                            sideID = sideIDofTriangle(iTriangle, abs(elem1))
                            elemID = abs(elem1)
                        else
                            sideID = sideIDofTriangle(iTriangle, abs(elem2))
                            elemID = abs(elem2)
                        endif
                    else
                        ! Pick the element with negative weight.
                        if (elem2.eq.0 .or. elem1.lt.0) then
                            ! absolute values because minus signs are carried 
                            ! along to indicate slippery weights, to save
                            ! memory.
                            sideID = sideIDofTriangle(iTriangle, abs(elem1))
                            elemID = abs(elem1)
                        else
                            sideID = sideIDofTriangle(iTriangle, abs(elem2))
                            elemID = abs(elem2)
                        endif
                    endif

                    ! get element coordinates into XL
                    call LCOORD (meshdatactx%X,XL,meshdatactx%IEN(1,elemID),0)
                    if (deformed) then
                        call LDISP (DL,modeldatactx%D,meshdatactx%IEN(1,elemID),NDOF,NEN)
                        call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
                        call AddFaultDisplacement (DL,elemID,modeldatactx%TFAULT,NDOF,NEN,36)
                        call REZONE (XL,DL,deffac)
                    endif

                    ! get element velocities into VL
                    if (iDisp.eq.20 .or. iDisp.eq.21) then
                        ! normal or total velocity
                        call LDISP (VL,modeldatactx%DELD,meshdatactx%IEN(1,elemID),NDOF,NEN)
                    else
                        call CLEAR (VL,NDOF*NEN,"VL")
                    endif
                    if (iDisp.eq.20 .or. iDisp.eq.22) then
                        ! slippery or total velocity
                        call ADDSNE (VL,elemID,NDOF,NEN,useDELDXE)
                    endif
                    if (iDisp.eq.20 .or. iDisp.eq.23) then
                        ! faulted or total velocity
                        call AddFaultDisplacement (VL,elemID,modeldatactx%DFAULT,NDOF,NEN,36)
                    endif

                    ! Get element-local point numbers of the element side that
                    ! makes up this triangle
                    call SIDENP(SideID,.true.,lien)

                    write(fileID,*) 3*iProgress-2,XL(:,lien(1)),VL(:,lien(1))
                    write(fileID,*) 3*iProgress-1,XL(:,lien(2)),VL(:,lien(2))
                    write(fileID,*) 3*iProgress-0,XL(:,lien(3)),VL(:,lien(3))
                else
                    ! no slippery and faulted nodes in the model: nodal point-centered data model.
                    do j=1,3
                        do i=1,NDOF
                            XL(i,j) = meshdatactx%X(i,meshdatactx%triangleConnectivity(j,iTriangle))
                            if (deformed) then
                                DL(i,j) = modeldatactx%D(i,meshdatactx%triangleConnectivity(j,iTriangle))
                            endif
                            VL(i,j) = modeldatactx%DELD(i,meshdatactx%triangleConnectivity(j,iTriangle))
                        enddo
                    enddo
                    if (deformed) then
                        call REZONE (XL,DL,deffac)
                    endif
                    write(fileID,*) 3*iProgress-2,XL(:,1),VL(:,1)
                    write(fileID,*) 3*iProgress-1,XL(:,1),VL(:,1)
                    write(fileID,*) 3*iProgress-0,XL(:,1),VL(:,1)
                endif
            endif
        enddo
    enddo
    close(fileID)
    return

endif

if (lnblk(fileName).le.0 .and. (.not.monitor)) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20  format(1x,'velo: read FE input and output files first')
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
30  format(1x,'velo: no velocities calculated')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    write(stderr,40)
40  format(1x,'velo: first read FE output fileName')
    ierr = EXEC
    return
 endif

 if (.not.Mout) then
    write(stderr,50) NSTEP
50  format(1x,'velo: no mechanical output at itime = ',I5)
    ierr = EXEC
    return
 endif
 
 if (icoord.eq.0.and.(.not.monitor)) then
    write(stderr,60)
60  format(1x,'velo: first write nodal points file with "np"')
    ierr = EXEC
    return
 endif

 if (ifelm.eq.0.and.(.not.monitor)) then
    write(stderr,70)
70  format(1x,'velo: first write connections file',1x, &
     'with "conn"')
    ierr = EXEC
    return
 endif
 
 if (NSTEP.eq.0) then
    write(stderr,80)
80  format(1x,'velo: no velocities available at time step 0')
    ierr = EXEC
    return
endif
 
! Open monitor file if neccessary
call OfOpen(fileName,monitor,ierr)

call wrdisp (fileName,meshdatactx%IEN,modeldatactx%DELD,modeldatactx%LMF,&
   modeldatactx%DFAULT,modeldatactx%DELDXE,plotdatactx%XPT,ierr)

call OfClose (fileName,monitor)
 
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xtc (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, ITMODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule,        only: icoord, ifelm, ifeout, ifein, stderr, tout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax
!
! Routine that X-ecutes the "tc" command.
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
!-init
ierr = OK
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: tc file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'tc: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (ITMODE.ne.0 .and. MODE.ge.3) then
        write(stderr,30)
30           format(1x,'tc: write temperature INPUT data')
    else
        if (MODE.ne.5 .and. MODE.ne.6) then
            write(stderr,40)
40               format(1x,'tc: no temperatures calculated')
            ierr = EXEC
            return
        else
            write(stderr,50)
50               format(1x,'tc: first read FE output file')
            ierr = EXEC
            return
        endif
    endif
else
    if (MODE.ne.5 .and. MODE.ne.6) then
        write(stderr,40)
        ierr = EXEC
        return
    endif
    if (.not.Tout) then
        write(stderr,60) NSTEP
60           format(1x,'tc: no thermal output at itime = ',I5)
        ierr = EXEC
        return
    endif
endif
!
if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'tc: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,80)
80       format(1x,'tc: first write connections file',1x, &
   'with "conn"')
    ierr = EXEC
    return
endif
!
!       Compute center temperatures into XPT
call tempc( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
 modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
   modeldatactx%T,modeldatactx%LMTX,modeldatactx%TX,plotdatactx%XPT,ierr)

if (ierr.ne.OK) then
    return
endif
!
if (ncmd.ge.2) then
    file = cmd(2)
endif
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

!    write file
call wrcval (file,plotdatactx%XPT,'tc',ierr)
!
!    Close monitor file again
call OfClose(file,elmon)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xmat (cmd,ncmd,ierr)

use modeldatamodule, only: plotdatactx
use meshdatamodule,  only: meshdatactx
use iomodule,        only: icoord, ifelm, ifein, stderr
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax

! Routine that X-ecutes the "mat" command.
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
!-init
ierr = OK
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: mat file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'mat: read FE input file first')
    ierr = EXEC
    return
endif
!
if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'mat: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,80)
80       format(1x,'mat: first write connections file',1x, &
   'with "conn"')
    ierr = EXEC
    return
endif
!
!       Compute center temperatures into XPT
call matc(meshdatactx%mat,plotdatactx%XPT)
!
if (ncmd.ge.2) then
    file = cmd(2)
endif
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

!    write file
call wrcval (file,plotdatactx%XPT,'mat',ierr)
!
!    Close monitor file again
call OfClose(file,elmon)
!
return
end
!-----------------------------------------------------------------------
subroutine xflux (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, ifein, ifeout, stderr, tout
use modeltopology,   only: nsd
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax

! Routine that X-ecutes the "flux" command.
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
!-init
ierr = OK
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: fl[ux] file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'flux: read FE input and output files first')
    ierr = EXEC
    return
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
    return
endif

if (MODE.ne.5 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Tout) then
    write(stderr,50) NSTEP
50       format(1x,'flux: no thermal output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'flux: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'flux: first write connections file with "conn"')
    ierr = EXEC
    return
endif
!
 if (ncmd.ge.2) file = cmd(2)
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

!    transfer heat flow data to XPT array, or echo monitored data
call setflux (modeldatactx%Thermalflux,plotdatactx%XPT,ierr)

if (ierr.ne.OK) return
if (.not.elmon) call wrvect (file,NSD,plotdatactx%XPT,ierr)
if (ierr.ne.OK) return
!
!    Close monitor file again
call OfClose(file,elmon)
!
return
end

!-------------------------------------------------------------------------------

subroutine xError(cmd,ncmd,ierr)

! Compute the norm of the a posteriori error vector following Verfurth J. Comp. Appl. Math. 1994:

use algebra,         only: vectorLength
use  meshdatamodule, only:  meshdatactx
use modeldatamodule, only: modeldatactx
use errors,          only: ok, exec, syntax

implicit none

integer ncmd,ierr
character cmd*(*)
dimension cmd(*)

integer :: iElem, fileID


! get filename from  command line args
if (ncmd.ne.2) then
    write(*,*) "usage: error [filename]"
else
    fileID = 382
    open(unit=fileID, file=cmd(2))
endif

allocate(modeldatactx%elemError(meshdatactx%neglobal,3))
call sumFnError3D(modeldatactx%elemError)

do iElem = 1, meshdatactx%neglobal
        write(*,*) iElem, vectorLength(modeldatactx%elemError(iElem,:), 3)
enddo

deallocate(modeldatactx%elemError)

close(382)

end subroutine

!-------------------------------------------------------------------------------
subroutine setflux (FLUX,XPT,ierr)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: icoord, ifelm, verbose, lumon, stderr
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: Fscale,Fsb
use errors,          only: exec
!
! Transfers heat flow vector data to XPT array, or echoes monitor data
!
implicit none
!-pass
integer :: ierr
double precision :: FLUX, XPT
dimension FLUX(NSD,*),XPT(NSD,*)
!-local
integer :: i, n, nn
integer :: ios


if (.not.elmon) then
    do n=1,NUMEL
        XPT(1,n) = REAL(FLUX(1,n)*Fscale)
        XPT(2,n) = REAL(FLUX(2,n)*Fscale)
        XPT(3,n) = REAL(FLUX(3,n)*Fscale)
    enddo
else
    do nn=1,nnel
        n = monel(nn)

        if (VERBOSE) then
            write(luMon,10,err=1100,iostat=ios) n,(FLUX(i,n), &
           i=1,NSD)
10               format(1x,'Heat flow in element ',I5,':', &
             3(1X,1PG15.6))
        else
            write(luMon,*,err=1100,iostat=ios) (FLUX(i,n), &
           i=1,NSD)
        endif
    enddo
endif
!
return
!
1100    write(stderr,12) ios,luMon
12   format(1x,'flux: write err ',I3,' on monitor logical unit ', &
    I3)
call prerr('flux')
ierr = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xvi (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule,        only: icoord, ifelm
use modeldefinition, only: ivlim, icvis
use modeltopology,   only: nen, nstr, ndof, nsd
use algebra,         only: idig
use iomodule,        only: ifein, ifeout, stderr, mout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use plotControl,     only: viscli, viscnu, vislog, visraw
use errors,          only: ok, exec, syntax

!
! Routine that X-ecutes the "viscosity" command.
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
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,1)
 1       format(1x,'Usage: vi[scosity] [[no]log] [no[nu[mbers]]]',1x, &
     '[no[li[near]]] file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,2)
 2        format(1x,'viscosity: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    write(stderr,3)
 3        format(1x,'viscosity: first read FE output file')
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,4) NSTEP
 4       format(1x,'no viscosities output at itime = ', &
     I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,5)
 5       format(1x,'viscosities: first write nodal points file',1x, &
   'with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,6)
 6       format(1x,'viscosities: first write connections file',1x, &
   'with "conn"')
    ierr = EXEC
    return
endif
!

i = 2
file = ' '
100     if (ncmd.ge.i) then
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
            write(stderr,1)
            ierr = SYNTAX
            goto 1000
        endif
        file = cmd(i)
    endif
    i = i + 1
    goto 100
endif

if (.not.elmon .and. lnblk(file).le.0) then
    write(stderr,1)
    ierr = SYNTAX
    goto 1000
endif
!
!    Recompute viscosities without upper and lower limits
if (visraw.and.IVLIM.ne.0.and.IDIG(ICVIS,4,1).ne.0) then
    call setvisc ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
     modeldatactx%T,modeldatactx%TX, &
     modeldatactx%STN,modeldatactx%STR, &
       modeldatactx%LMF,modeldatactx%tfault,modeldatactx%LMX,modeldatactx%DXE, &
       modeldatactx%LMTX,modeldatactx%SKEW,NSD,NSTR,NDOF, &
       NEN)
endif

if (visraw.and.IVLIM.ne.0.and.IDIG(ICVIS,4,3).ne.0) then
    call vweak (modeldatactx%STR,NSTR)
endif
!
!       Transfer viscosities into B
call viscosityTransfer(plotdatactx%XPT)
!
!       Open monitor file if neccessary
call OfOpen(file,elmon,ierr)
if (ierr.ne.OK) then
    goto 1000
endif

!       Write file
call wrcval (file,plotdatactx%XPT,'visc',ierr)
if (ierr.ne.OK) then
    goto 1000
endif
!
!       Close monitor file again
call OfClose(file,elmon)
!
1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine viscosityTransfer(XPT)


USE ALGEBRA
use materialsmodule, only: MaterialIDOfElement, &
                           matdatactx
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: nstr
use timestepmodule,  only: NINTG,MAXIT
use plotControl,     only: Vscale, vislog, viscli
!
!    transfers viscosities from PRPMAT to XPT
!
implicit none
!-pass
double precision XPT(*)
!-local
integer          :: m, n
double precision :: emhu, anpwr, es
double precision, external :: effs, etaLin
!
do n=1,meshdatactx%neglobal

    m = MaterialIDOfElement(n)
    EMHU  = matdatactx%PRPMAT(3,m)
    ANPWR = matdatactx%PRPMAT(4,m)

    if (viscli) then
        ES = effs (modeldatactx%STN(1,n))
        EMHU = etaLin (ES,EMHU,ANPWR)
    endif

    if (Vscale.gt.0.) then
        EMHU = Vscale*EMHU
    endif

    if (vislog) then
        EMHU = TENLOG(EMHU)
    endif

    ! save to output variable
    XPT(n) = REAL(EMHU)
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine xnbc (cmd,ncmd,ierr)
!
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifein, stderr
use monitorVars,     only: npmon
use plotControl,     only: npbased
use errors,          only: exec, syntax

! Routine that X-ecutes the "nbc" command.
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
integer, save :: iflg
integer :: i
!-init
data iflg/6/
file = ' '
i = 2
!
10    format(1x,'Usage: nbc [a|d|v|f|i|r] file')
!
100 if (i.le.ncmd) then
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
            ierr = SYNTAX
            return
        endif
    endif

    i = i + 1
    goto 100
endif

if ((.not.npmon) .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'nbc: read FE input file first')
    ierr = EXEC
    return
endif
!
if (icoord.eq.0.and.(.not.npmon)) then
    write(stderr,30)
30       format(1x,'nbc: first write nodal points file with "np"')
    ierr = EXEC
    return
endif
!
if (.not.npbased) then
    write(stderr,31)
31       format(1x,'nbc: element based node definitions not ', &
     'implemented')
    ierr = EXEC
    return
endif
!
!       Open monitor file if neccessary
call OfOpen(file,npmon,ierr)
!
call wrnbc (file,modeldatactx%ibond, modeldatactx%bond, &
                  modeldatactx%SKEW,plotdatactx%XPT,iflg,ierr)
!
call OfClose (file,npmon)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xebc (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use modeldefinition, only: iword, numel
use modeltopology,   only: nsd
use iomodule,        only: icoord, ifein, stderr
use errors,          only: exec, syntax
!
! Routine that X-ecutes the "ebc" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer :: i
integer, external :: lnblk
integer, save :: iflg
!-init
file = ' '
iflg = -1
10    format(1x,'Usage: ebc w|s file')
!
i = 2
do while (i.le.ncmd)
    if (cmd(i).eq.'w') then
        iflg = 1
    else if (cmd(i).eq.'s') then
        iflg = 2
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif
    i = i + 1
enddo

if (lnblk(file).le.0 .or. iflg.eq.-1) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'ebc: read FE input file first')
    ierr = EXEC
    return
endif
!
if (iflg.eq.1) then
    call BCWNKP (modeldatactx%IWELM,modeldatactx%IWSIDE,&
        modeldatactx%IWTIME,modeldatactx%WPRES, &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
       modeldatactx%LMF,modeldatactx%tfault, &
       modeldatactx%SKEW, plotdatactx%xpt,file)
else if (iflg.eq.2) then

    call BCSTRS (modeldatactx%IELSTR, modeldatactx%ISSIDE, &
             modeldatactx%ISTR, modeldatactx%STRS, &   ! does this not conflict with modeldatactx%STR ?
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
       modeldatactx%LMF,modeldatactx%tfault, &
       modeldatactx%SKEW,file)
endif
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xreslv (cmd,ncmd,ierr)
!
USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule,        only: ifeout, ifein, stderr, mout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax

! Routine that X-ecutes the "rt" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
integer           :: i, j
double precision,    external :: chreal
logical, external :: numeric
integer, external :: lnblk
integer, save     :: iflg,iside
!-init
data iflg,iside/0,1/
i = 2
file = ' '
!
10    format(1x,'Usage: rt [stress|strain|rate] [facenumber]')
!
100    if (i.le.ncmd) then
    if (cmd(i).eq.'stress') then
        iflg = 0
    else if (cmd(i).eq.'strain') then
        iflg = 1
    else if (cmd(i).eq.'rate') then
        iflg = 2
    else if (numeric(cmd(i))) then
        j = NINT(chreal(cmd(i)))
        if (j.ge.1 .and. j.le.6) then
            iside = j
        else
            write(stderr,20)
20               format(1x,'rt: face numbers should be in the',1x, &
             'range [1,6]')
            ierr = SYNTAX
            return
        endif
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif

    i = i + 1
    goto 100

endif
!
if (ifein.eq.0) then
    write(stderr,30)
30        format(1x,'rt: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,40)
40        format(1x,'rt: no mechanical solution calculated')
    else
        write(stderr,50)
50        format(1x,'rt: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,40)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,60) NSTEP
60       format(1x,'rt: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    return
endif

if (iflg.eq.2.and.NSTEP.eq.0) then
    write(stderr,70)
70        format(1x,'rt: no elastic strain rates')
    ierr = EXEC
    return
endif
!
if (.not.elmon) then
    return
endif
!
!       Open monitor file if neccessary
call OfOpen(file,elmon,ierr)
if (ierr.ne.OK) then
    return
endif
!
if (iflg.eq.0) then
    call wreslv ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
     modeldatactx%LMF,modeldatactx%tfault, &
       modeldatactx%DXE,modeldatactx%STN,'  stress   ',iside,ierr)
else if (iflg.eq.1) then
    call wreslv ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
     modeldatactx%LMF,modeldatactx%tfault, &
       modeldatactx%DXE,modeldatactx%STR,'  strain   ',iside,ierr)
#ifndef EXPERIMENTAL_NORATESTORAGE
else
    call wreslv ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
     modeldatactx%LMF,modeldatactx%tfault, &
       modeldatactx%DXE,modeldatactx%rate,'strain rate',iside,ierr)
#endif
endif
!
call OfClose (file,elmon)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xdebrah (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use meshdatamodule,  only: meshdatactx
use materialsmodule, only: matdatactx
use iomodule,        only: icoord, ifelm, ifeout, ifein, stderr, mout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax
!
! Routine that X-ecutes the "deborah" command.
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
!-init
ierr = OK
file = ' '
!
if ((.not.elmon) .and. ncmd.lt.2) then
    write(stderr,10)
10        format(1x,'Usage: deb[orah] file')
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'deborah: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'deborah: no strain rates calculated')
    else
        write(stderr,40)
40        format(1x,'deborah: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'deborah: no mechanical output at itime = ', I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'deborah: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'deborah: first write connections file with "conn"')
    ierr = EXEC
    return
endif
!
if (ncmd.ge.2) file = cmd(2)
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

write(*,*) "Deborah has been temporarily disabled"
! call Deborah (modeldatactx%rate,modeldatactx%STN,meshdatactx%mat,matdatactx%prpmat,B,ierr)
 if (ierr.ne.OK) goto 1000
!
!       Write file
call wrcval (file,plotdatactx%XPT,'deborah',ierr)
 if (ierr.ne.OK) goto 1000
!
!    Close monitor file again
 call OfClose(file,elmon)
!
1000    return
 end
!-------------------------------------------------------------------------------
subroutine xcoulomb (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
use modeldatamodule, only: modeldatactx, plotdatactx
use iomodule,        only: icoord, ifelm, ifeout, ifein, stderr, mout
use modeltopology,   only: nsd
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax
!
! Routine that X-ecutes the "coulomb" command.
!
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
double precision,    external :: chreal
double precision, save :: frmhu
real :: effe
integer :: i
!-init
data frmhu /0.4d0/
ierr = OK
i = 2
file = ' '
!
10    format(1x,'Usage: coulomb [mhu] file')
!
100    if (i.le.ncmd) then
    if (numeric(cmd(i))) then
        effe = chreal(cmd(i))
        if (effe.ge.0d0 .and. effe.le.1d0) then
            frmhu = effe
        else
            write(stderr,15)
15               format(1x,'coulomb: friction coefficient should', &
             'be in the range [0,1]')
            ierr = SYNTAX
            return
        endif
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif
    i = i + 1
    goto 100
endif
!
if ((.not.elmon) .and. lnblk(file).le.0) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif

if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'coulomb: read FE input and output files first')
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30        format(1x,'coulomb: no stresses calculated')
    else
        write(stderr,40)
40        format(1x,'coulomb: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'coulomb: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0 .and. (.not.elmon)) then
    write(stderr,60)
60       format(1x,'coulomb: first write nodal points file with', &
     ' "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0 .and. (.not.elmon)) then
    write(stderr,70)
70       format(1x,'coulomb: first write connections file with', &
   ' "conn"')
    ierr = EXEC
    return
endif
!
!    Open monitor file if neccessary
call OfOpen(file,elmon,ierr)

call coulomb (modeldatactx%STN,frmhu,plotdatactx%XPT,ierr)
if (ierr.ne.OK) goto 1000
!
if (.not.elmon) call wrvect (file,NSD,plotdatactx%XPT,ierr)
if (ierr.ne.OK) goto 1000
!
!    Close monitor file again
 call OfClose(file,elmon)
!
1000    return
end subroutine
!-------------------------------------------------------------------------------
subroutine xrot (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE
USE MODELDATAMODULE
USE MESHDATAMODULE
use iomodule,        only: icoord, ifelm, ifein, ifeout, stderr, mout
use timestepmodule,  only: nstep
use monitorVars,     only: elmon
use errors,          only: exec, syntax

!
! Routine that X-ecutes the "rot" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
#include "lenchar.i"
character(len=LENCHAR) file
logical RATE
integer, external :: lnblk
integer :: i
save RATE
!-init
data RATE /.false./
file = ' '
i = 2
!
10    format(1x,'Usage: rot[ation] [no[rate]] file')
!
100    if (i.le.ncmd) then
    if (cmd(i).eq.'norate') then
        RATE = .false.
    else if (cmd(i).eq.'rate') then
        RATE = .true.
    else
        if (lnblk(file).le.0) then
            file = cmd(i)
        else
            write(stderr,10)
            ierr = SYNTAX
            return
        endif
    endif
    i = i + 1
    goto 100
endif

if (lnblk(file).le.0 .and. (.not.elmon)) then
    write(stderr,10)
    ierr = SYNTAX
    return
endif
!
if (ifein.eq.0) then
    write(stderr,20)
20        format(1x,'rot: read FE input and output files first')
    ierr = EXEC
    return
endif
!
if (ifeout.eq.0) then
    if (MODE.ne.2 .and. MODE.ne.6) then
        write(stderr,30)
30           format(1x,'rot: no mechanical solution computed')
    else
        write(stderr,40)
40           format(1x,'rot: first read FE output file')
    endif
    ierr = EXEC
    return
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    write(stderr,30)
    ierr = EXEC
    return
endif

if (.not.Mout) then
    write(stderr,50) NSTEP
50       format(1x,'rot: no mechanical output at itime = ', &
     I5)
    ierr = EXEC
    return
endif

if (icoord.eq.0.and.(.not.elmon)) then
    write(stderr,60)
60       format(1x,'rot: first write nodal points file with "np"')
    ierr = EXEC
    return
endif

if (ifelm.eq.0.and.(.not.elmon)) then
    write(stderr,70)
70       format(1x,'rot: first write connections file',1x, &
     'with "conn"')
    ierr = EXEC
    return
endif

if (NSTEP.eq.0 .and. RATE) then
    write(stderr,80)
80        format(1x,'rot: no rotation rates available at time step 0')
    ierr = EXEC
    return
endif
!
!       Open monitor file if neccessary
call OfOpen(file,elmon,ierr)
!
call INFROT (file,modeldatactx%B, &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
        modeldatactx%DELD,&
     modeldatactx%LMF, &
   modeldatactx%tFault, &
   modeldatactx%SKEW,modeldatactx%DELDXE,&
#ifndef EXPERIMENTAL_NORATESTORAGE
         modeldatactx%rate,&
#endif
ierr)

!
call OfClose (file,elmon)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine xsample (cmd,ncmd,ierr)

USE MODELDEFINITION, only: MODE, NUMFNGlobal, numel
USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx, &
                           MaterialIDOfElement, &
                           lmat
use debugmodule,     only: AllocateError
use iomodule,        only: icoord, ifelm, ifeout, stderr, tout, mout, ifein, verbose
use modeltopology,   only: nen, ndof, nstr, nsd
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
use conversions,     only: oneCharIsNumeric, &
                           resize_dblarray
use debugmodule,     only: parallelLog, iecho
use timestepmodule,  only: nstep
use postProcessData, only: deformed
use errors,          only: ok, exec, syntax

! Routine that X-ecutes the "sample" command.

implicit none
!-pass
integer          :: ncmd,ierr
character        :: cmd*(*)
dimension        :: cmd(*)
integer          :: xi, yi, zi  ! loop indices
double precision :: xmin, xmax, xinc
double precision :: ymin, ymax, yinc
double precision :: zmin, zmax, zinc




!-locl
logical            :: numeric,mo,to
double precision   :: S(6)
character(len=255) :: WriteToFileName, readFromFileName
logical            :: writeToFile, there
character(len=255) :: argumentString
integer, parameter :: FileHandle=53
integer            :: idata,lastSampledElement
integer            :: sampleType
integer            :: xyzsize
double precision   :: x, y, z
double precision   :: xr, yr, zr
integer            :: nx, ny, nz
integer            :: i,ii,j,jj,k,kk, iostatm, iElement
double precision, dimension(:), allocatable :: xlist, ylist, zlist
integer            :: npoints, idx
integer            :: oldsize, newsize
double precision   :: dummyDouble
logical            :: samplePointFound
integer            :: matNumber
double precision, allocatable :: viscosities(:)



!-init
mo = .false.
to = .false.
xyzsize = 10
samplePointFound = .false.

!
10    format(1x,'Usage: sa[mple] d|v|T|s|e|r|f coordinate [range] [filename]')
!
!if (ncmd.lt.6) then
!    write(stderr,10)
!    ierr = SYNTAX
!    goto 1000
!endif
!
if (iecho.eq.8) then
    call parallelLog("xsample, 3D","entered sub")
endif



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
else if (cmd(2).eq.'w') then
    idata = SampleRotation
    mo = .true.
else if (cmd(2).eq.'f') then
    idata = SampleHeatFlux
    mo = .true.
else
    write(stderr,10)
    ierr = SYNTAX
    goto 1000
endif
 

argumentString = cmd(3)

if    (oneCharIsNumeric(argumentString(1:1)) .or. argumentString(1:1).eq."-") then
    ! sample by range
    sampleType = SampleOnCubicGrid

    call parseSampleRange(cmd(3),xmin,xmax,xinc,ierr)

    if (ierr.ne.OK) goto 1000

    call parseSampleRange(cmd(4),ymin,ymax,yinc,ierr)

    if (ierr.ne.OK) goto 1000

    call parseSampleRange(cmd(5),zmin,zmax,zinc,ierr)

    if (ierr.ne.OK) goto 1000


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

    if (zmin.eq.zmax) then
        nz = 1
    else
        nz = floor((zmax-zmin)/zinc)+1
    endif

    ! build a 3D list
    npoints = (nx * ny * nz)

    allocate(xlist(npoints), stat=ierr)
    call AllocateError("xlist", ierr)

    allocate(ylist(npoints), stat=ierr)
    call AllocateError("ylist", ierr)

    allocate(zlist(npoints), stat=ierr)
    call AllocateError("zlist", ierr)

    idx = 1

    do kk=1,nz
        zr = zmin + DBLE(kk-1)*zinc
        do jj=1,ny
            yr = ymin + DBLE(jj-1)*yinc
            do ii=1,nx
                xr = xmin + DBLE(ii-1)*xinc
                xlist(idx) = xr
                ylist(idx) = yr
                zlist(idx) = zr
                idx = idx + 1
            enddo
        enddo
    enddo

    ! 3D list completed
    if (ncmd.eq.6) then
        ! write to file
        WriteToFileName = trim(cmd(6))
        writeToFile = .true.
        open(unit=FileHandle, file=WriteToFileName)

    else if (ncmd.eq.5) then
        ! write to standard out
        writeToFile = .false.
    else
        write(*,*) "Syntax error; see: ? plnplt"
    endif

else if (argumentString(1:4) .eq. "line") then

    write(*,*) "sampling a line in 3D; not yet correctly implemented"


    sampleType = SampleOnAnyLine

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
            allocate(xlist(npoints), stat=ierr)
            allocate(ylist(npoints), stat=ierr)
            allocate(zlist(npoints), stat=ierr)

            xinc = dble(xmax - xmin) / dble(nPoints-1)
            yinc = dble(ymax - ymin) / dble(nPoints-1)
            zinc = dble(zmax - zmin) / dble(nPoints-1)


            do i=1,nPoints
                xlist(i) = xmin + (i-1) * xinc
                ylist(i) = ymin + (i-1) * yinc
                zlist(i) = zmin + (i-1) * zinc
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

    sampleType = SampleFromFile

    readFromFileName = cmd(4)

!    write(*,*) "cmd1: ", trim(cmd(1))
!   write(*,*) "cmd2: ", trim(cmd(2))
!   write(*,*) "cmd3: ", trim(cmd(3))
!   write(*,*) "cmd4: ", trim(cmd(4))




    ! the third argument is a file name, containing coordinate pairs
    ! problem is that we do not know how many 
    ! coordinate pairs there are in this file.
    ! We will use an array, starting at size 10,
    ! that we will double in size every time we run out.

    ! After completely running through the file,
    ! we arrays will be scaled back to fit the
    ! data excatly.

    allocate(xlist(xyzsize), stat=ierr)
    call AllocateError("xlist", ierr)

    allocate(ylist(xyzsize), stat=ierr)
    call AllocateError("ylist", ierr)

    allocate(zlist(xyzsize), stat=ierr)
    call AllocateError("zlist", ierr)

    ! check whether the file exists.
    inquire(file=trim(cmd(4)), exist=there)

    if (there) then


        if (ncmd.eq.5) then
            ! write to file
            WriteToFileName = trim(cmd(5))
            writeToFile = .true.
            open(unit=FileHandle, file=WriteToFileName)
        else if (ncmd.eq.4) then
            ! write to standard out
            writeToFile = .false.
        else
            write(*,*) "Syntax error; see ? sample"
        endif

        open(123, file=trim(cmd(4)))

        idx = 1

        do
            READ(123,*,end=13) x,y,z

            ! check whether xsize and ysize are full.
            if (idx > xyzsize) then
                oldsize = xyzsize
                xyzsize  = xyzsize * 2
                newsize = xyzsize
                ! xlist and ylist will overflow it not enlarged. Double its size.
                call resize_dblarray(xlist, oldsize, newsize)
                call resize_dblarray(ylist, oldsize, newsize)
                call resize_dblarray(zlist, oldsize, newsize)
            else
                ! the array is still large anough. Do nothing
            endif

            xlist(idx) = x
            ylist(idx) = y
            zlist(idx) = z

            idx = idx + 1
        enddo

13      write(0,*) 'finished reading file'
        npoints = idx

        call resize_dblarray(xlist, xyzsize, npoints)
        call resize_dblarray(ylist, xyzsize, npoints)
        call resize_dblarray(zlist, xyzsize, npoints)

    else

        write(*,*) "Sample: Could not find input file: ", trim(cmd(3))
        stop "Leaving GTecton..."
    endif

else
    write(stderr, *) "argument ", trim(argumentString), " not recognized."
    write(stderr, *) "must be either: "
    write(stderr, *) "numeric, to indicate a range,"
    write(stderr, *) "'line',  to indicate sampling over an arbitrary line,"
    write(stderr, *) "'file',  to read sample points from a file."
endif


if (ifein.eq.0) then
    write(stderr,20)
20   format(1x,'sample: read FE input file first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    write(stderr,30)
30        format(1x,'sample: read FE output file first')
    ierr = EXEC
    goto 1000
endif

if (mo.and.(MODE.ne.2.and.MODE.ne.6)) then
    write(stderr,40)
40        format(1x,'sample: no mechanical solution calculated')
    ierr = EXEC
    goto 1000
else
    if (.not.Mout) then
        if (idata.eq.1) then
            write(stderr,50) 'displacements',NSTEP
50                format(1x,'sample: no ',A13,' at itime = ',I5)
            ierr = EXEC
            goto 1000
        else if (idata.eq.2) then
            write(stderr,50) 'velocities   ',NSTEP
            ierr = EXEC
            goto 1000
        endif
    endif
endif

if (MODE.ne.5.and.MODE.ne.6) then
    if (to) then
        write(stderr,60)
60        format(1x,'sample: no thermal solution calculated')
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

 
lastSampledElement = 1
! find some way to to this cleverly... Will save sample times
!lastSampledElementY = 1
!lastSampledElementZ = 1

if (iecho.eq.8) then
    call parallelLog("xsample, 3D","enterin sample loop")
endif


if (npoints .gt. 1e6) then
    write(*,*) "Warning: sampling ", npoints, "points. This may take a while"
endif


if (writeToFile) then
    open(unit=FileHandle, file=trim(writeToFileName))
endif


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
    zr = zlist(j)

!      ========================= displacement =======
    if (idata.eq.SampleDisplacement) then


! RG 2021-01 call LOADFT (meshdatactx%ien,modeldatactx%NFAULT,modeldatactx%FAULT,&
! RG 2021-01             modeldatactx%TFAULT,NUMEL)


!     do i=1,NUMFNglobal
!        write(stderr,*) i,modeldatactx%TFAULT(1,i),modeldatactx%TFAULT(2,i),modeldatactx%TFAULT(3,i)
!     enddo

!    write(*,*) "calling nsmpl for displacement"

   call nsmpl ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
       modeldatactx%LMF,modeldatactx%tfault, &
       modeldatactx%SKEW,xr,yr,zr,modeldatactx%D,   &
       NDOF,S,deformed, idata,lastSampledElement, samplePointFound)

!    write(*,*) "finished displacement"


    if (samplePointFound) then

       if (verbose) then
           write(FileHandle,70) xr,yr,zr,(S(i),i=1,NDOF)
70            format(1x,'Displacement at ',3(1PG24.16),' = ',3(1PG14.6))
       else
           write(FileHandle,*) xr,yr,zr,(S(i),i=1,NDOF)
       endif
    else
       write(FileHandle,65) xr,yr,zr
65          format(1x,'Sample: coordinate not found ',3(1PG14.6))

   endif

!      =========================== velocity ============
    else if (idata.eq.SampleVelocity) then

   call nsmpl ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
         modeldatactx%D,   &
         modeldatactx%LMF,modeldatactx%tfault, &
         modeldatactx%SKEW,xr,yr,zr, &
         modeldatactx%DELD,&
         NDOF,S,deformed, idata,lastSampledElement, samplePointFound)


    if (samplePointFound) then

       if (verbose) then
           write(FileHandle,80) xr,yr,zr,(S(i),i=1,NDOF)
80            format(1x,'Velocity at ',3(1PG14.6),' = ',3(1PG14.6))
       else
           write(FileHandle,*) xr,yr,zr,(S(i),i=1,NDOF)
       endif
    else
       write(FileHandle,65) xr,yr,zr
   endif


!      =========================== viscosity ============
    else if (idata.eq.SampleViscosity) then

   call nsmpl ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
         modeldatactx%D,   &
         modeldatactx%LMF,modeldatactx%tfault, &
         modeldatactx%SKEW,xr,yr,zr, &
         viscosities,&
         NDOF,S,deformed, idata,lastSampledElement, samplePointFound)


    if (samplePointFound) then

       if (verbose) then
           write(FileHandle,81) xr,yr,zr,(S(i),i=1,NDOF)
81            format(1x,'Viscosity at ',3(1PG14.6),' = ',3(1PG14.6))
       else
           write(FileHandle,*) xr,yr,zr,(S(i),i=1,NDOF)
       endif
    else
       write(FileHandle,65) xr,yr,zr
   endif





!      =========================== temperature ============
    else if (idata.eq.SampleTemperature) then
   call nsmpl ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
         modeldatactx%D,   &
         modeldatactx%LMF,modeldatactx%tfault, &
         modeldatactx%SKEW,xr,yr,zr,modeldatactx%T,1,S,&
         deformed, idata,lastSampledElement, samplePointFound)

    if (samplePointFound) then

       if (verbose) then
           write(FileHandle,90) xr,yr,zr,S(1)
90            format(1x,'Temperature at ',3(1PG14.6),' = ',1PG14.6)
       else
           write(FileHandle,*) xr,yr,zr,S(1)
       endif

    else
       write(FileHandle,65) xr,yr,zr



   endif

    else if (idata.eq.SampleStress) then
   call esmpl ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
         modeldatactx%D,   &
         modeldatactx%LMF,modeldatactx%tfault, &
         modeldatactx%SKEW,xr,yr,zr,modeldatactx%STN,NSTR,S,&
         deformed,lastSampledElement, samplePointFound)


    if (samplePointFound) then

       if (verbose) then
           write(FileHandle,101) xr,yr,zr,(S(i),i=1,NSTR)
  101            format(1x,'Stress at ',3(1PG14.6),' = ',6(1PG14.6))
       else
           write(FileHandle,*) xr,yr,zr,(S(i),i=1,NSTR)
       endif

    else
       write(FileHandle,65) xr,yr,zr

   endif

    else if (idata.eq.SampleStrain) then
   call esmpl ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
        modeldatactx%D,   &
       modeldatactx%LMF,modeldatactx%tfault, &
       modeldatactx%SKEW,xr,yr,zr,modeldatactx%STR,NSTR,S,&
       deformed,lastSampledElement, samplePointFound)

    if (samplePointFound) then

           if (verbose) then
               write(FileHandle,110) xr,yr,zr,(S(i),i=1,NSTR)
  110            format(1x,'Strain at ',3(1PG14.6),' = ',6(1PG14.6))
           else
              write(FileHandle,*) xr,yr,zr,(S(i),i=1,NSTR)
           endif

        else
           write(FileHandle,65) xr,yr,zr


       endif

    else if (idata.eq.SampleStrainRate) then
#ifndef EXPERIMENTAL_NORATESTORAGE
    call esmpl ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
         modeldatactx%D,   &
         modeldatactx%LMF,modeldatactx%tfault, &
         modeldatactx%SKEW,xr,yr,zr,modeldatactx%rate,NSTR,S,&
         deformed,lastSampledElement, samplePointFound)
#endif
    if (samplePointFound) then

       if (verbose) then
           write(FileHandle,120) xr,yr,zr,(S(i),i=1,NSTR)
  120            format(1x,'Strain rate at ',3(1PG14.6),' = ',6(1PG14.6))
       else
           write(FileHandle,*) xr,yr,zr,(S(i),i=1,NSTR)
       endif
    else
        write(FileHandle,65) xr,yr,zr



    endif
    else if (idata.eq.SampleHeatFlux) then
   call esmpl ( &
         meshdatactx%IEN, &
         meshdatactx%X,   &
         modeldatactx%D,   &
         modeldatactx%LMF,modeldatactx%tfault, &
         modeldatactx%SKEW,xr,yr,zr,modeldatactx%ThermalFlux,NSD,S,&
         deformed,lastSampledElement, samplePointFound)

    if (samplePointFound) then

           if (verbose) then
               write(FileHandle,130) xr,yr,zr,(S(i),i=1,NSD)
  130              format(1x,'Heat flow at ',3(1PG14.6),' = ',3(1PG14.6))

           else
               write(FileHandle,*) xr,yr,zr,(S(i),i=1,NSD)
           endif

       endif
    else
           write(FileHandle,65) xr,yr,zr
    endif


enddo ! enddo loop over de points from the sample list


if (idata.eq.SampleViscosity) then
    deallocate(viscosities)
endif


!
if (iecho.eq.8) then
    call parallelLog("xsample, 3D","leaving sub")
endif


  
1000 close(FileHandle)

 return
 end

!-------------------------------------------------------------------------------

subroutine writeSampledData(fileHandle, quantity, xr, yr, zr, nValues, values)
! TODO, actually use this
use iomodule, only: verbose

implicit none

integer          :: fileHandle
character(len=*) :: quantity
double precision :: xr, yr, zr
integer          :: nValues
double precision :: values(nValues)

integer          :: i

if (values(1).eq.-999.0 .and. &
    values(2).eq.-999.0 .and. &
    values(3).eq.-999.0) then
    write(FileHandle,*) xr,yr,zr
else
    if (verbose) then
        write(FileHandle,*) trim(quantity) // " at", xr,yr,zr,(values(i),i=1,nValues)
    else
        write(FileHandle,*) xr,yr,zr,(values(i),i=1,nValues)
    endif
endif


end subroutine



!-------------------------------------------------------------------------------
subroutine xgrav (cmd,ncmd,ierr)

USE MODELDATAMODULE
USE MODELDEFINITION, only: NUMFNglobal, numel
USE MESHDATAMODULE
USE MATERIALSMODULE, only: matdatactx
use modeltopology,   only: nen, ndof
use iomodule,        only: ifein, stderr
use monitorVars,     only: elmon
use errors,          only: ok, exec, syntax
!
! Routine that X-ecutes the "grav" command.
!
implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
integer arglog
#include "lenchar.i"
character(len=LENCHAR) file
logical numeric
double precision lonstep,lonmin,lonmax,latstep,latmin,latmax,height
integer i,nnum
! logical tecton
integer, external :: lnblk
double precision,    external :: chreal
!-init
ierr = OK
file = ' '

 10 format(1x,'Usage: gra[vity] lonstep lonmin lonmax latstep latmin latmax heigth file')
 ! read arguments
! nnum = 2
! read(cmd(nnum),'(i1)') arglog
! write(stderr,*) '2',arglog
! if (arglog.eq.0) then
!  detpot = .false.
! else if (arglog.eq.1) then
!  detpot = .true.
! else
!  write(stderr,10)
!  ierr = SYNTAX
!  goto 1000
! endif
!   write(stderr,*) 'potential',detpot
! nnum = 3
! read(cmd(nnum),'(i1)') arglog
! write(stderr,*) '3',arglog
! if (arglog.eq.0) then
!  detgrav = .false.
! else if (arglog.eq.1) then
!  detgrav = .true.
! else
!  write(stderr,10)
!  ierr = SYNTAX
!  goto 1000
! endif
! write(stderr,*) 'gravity',detgrav
! nnum = 4
! read(cmd(nnum),'(i1)') arglog
! write(stderr,*) '4',arglog
! if (arglog.eq.0) then
!  detgrad = .false.
! else if (arglog.eq.1) then
!  detgrad = .true.
!  write(stderr,*) 'gradient',detgrad
! else
!  write(stderr,10)
!  ierr = SYNTAX
!  goto 1000
! endif
! write(stderr,*) 'gradient',detgrad
! nnum = 5
! coordtype = cmd(nnum)(1:4)
nnum = 1
do i=2,ncmd
    if (numeric(cmd(i))) then
        nnum = nnum + 1
        if (nnum.eq.2) then
            lonstep = chreal(cmd(i))
!           write(stderr,*) 'lonstep',lonstep
        else if (nnum.eq.3) then
            lonmin = chreal(cmd(i))
!           write(stderr,*) 'lonmin',lonmin
        else if (nnum.eq.4) then
            lonmax = chreal(cmd(i))
!           write(stderr,*) 'lonmax',lonmax
        else if (nnum.eq.5) then
            latstep = chreal(cmd(i))
!           write(stderr,*) 'latstep',latstep
        else if (nnum.eq.6) then
            latmin = chreal(cmd(i))
!           write(stderr,*) 'latmin',latmin
        else if (nnum.eq.7) then
            latmax = chreal(cmd(i))
!           write(stderr,*) 'latmax',latmax
        else if (nnum.eq.8) then
            height = chreal(cmd(i))
!           write(stderr,*) 'heigth',height
        else
            write(stderr,*) 'unrecognized option ',i,cmd(i)
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

write(stderr,*) 'output file',file

if (lnblk(file).le.0 .and. (.not.elmon)) then
   write(stderr,10)
   ierr = SYNTAX
   goto 1000
endif

! check if input and output of FE has already been read
if (ifein.eq.0) then
   write(stderr,20) 
   20        format(1x,'displ: read FE input and output files first')
   ierr = EXEC
   return
endif

 ! check on sample points
if (lonmin.gt.lonmax) then
    write(stderr,*) 'longitude min larger than max',lonmin,lonmax
    return
endif

if (latmin.gt.latmax) then
   write(stderr,*) 'latitude min larger than max',latmin,latmax
   return
endif 
 
 ! update faulted displacements
 
! RG 2021-01 call LOADFT (meshdatactx%ien,modeldatactx%NFAULT,modeldatactx%FAULT,&
! RG 2021-01               modeldatactx%TFAULT,NUMEL)

 call OfOpen (file,elmon,ierr)
if (ierr.ne.OK) then
    goto 1000
endif

call gravitytetrahedron (meshdatactx%IEN,meshdatactx%X,&
 modeldatactx%D,meshdatactx%mat, matdatactx%PRPMAT,&
 modeldatactx%STR,modeldatactx%TFAULT,modeldatactx%DXE,&
 lonmin,lonmax,lonstep,latmin,latmax,latstep,height,&
 file)

 call OfClose (file,elmon)

1000    return
end subroutine
!--------------------------------------------------------------------------------------------------!

subroutine xfullgrav (cmd,ncmd,ierr)

use iomodule,        only: ifein, ifeout, stderr
use errors,          only: exec

! computes the gravity of the domain into a list of points, given in a file.
! This can for example be a spherical around an Earth model.
! This list should contain X,Y,Z pairs

implicit none

!-pass
integer ncmd, ierr
character cmd*(*)
dimension cmd(*)

!local
integer :: nSamplePoints, iSamplePoint
double precision :: samplePoint(3)

double precision, allocatable :: samplePoints(:,:)
double precision, allocatable :: gravityVec(:,:)
double precision, allocatable :: Potential(:)
double precision, allocatable :: GravGradient(:,:,:)





if (ncmd.ne.3) then
    write(*,*) "Usage: fullgrav inputfile outputfil"
    return
endif

if (ifein.eq.0) then
    write(stderr,*) 'xsec: read FE input file first'
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    write(stderr,*) 'xsec: read FE output file first'
    ierr = EXEC
    return
endif

!--- read sample points

nSamplePoints = 0
open(unit=68, file=trim(cmd(2)))
do while (.true.)
    read(68,*,end=10) samplePoint(1), samplePoint(2), samplePoint(3)
    nSamplePoints = nSamplePoints + 1
enddo

10 rewind(68)

allocate (samplePoints(nSamplePoints,3  ))
allocate (  gravityVec(nSamplePoints,3  ))
allocate (   Potential(nSamplePoints    ))
allocate (GravGradient(nSamplePoints,3,3))     

if (nSamplePoints.eq.0) then
    write(*,*) "Plt3d was not able to read points from file ",trim(cmd(2))
    return
else
!   write(*,*) "found ", nSamplePoints, "nSamplePoints"
endif

open(unit=69, file=trim(cmd(3)))


do iSamplePoint=1,nSamplePoints
    read(68,*) samplePoints(iSamplePoint,:)
enddo

!$OMP PARALLEL
!$OMP DO
do iSamplePoint=1,nSamplePoints
    call gravityOfDomain(samplePoints(iSamplePoint,:), &
                         gravityVec(iSamplePoint,:), &
                         Potential(iSamplePoint), &
                         GravGradient(iSamplePoint,:,:))


enddo
!$OMP END DO
!$OMP END PARALLEL

do iSamplePoint=1,nSamplePoints
    ! write vector and magnitude
    write(69,*) gravityVec(iSamplePoint,:), &
                sqrt(sum(gravityVec(iSamplePoint,:)**2)), &
                Potential(iSamplePoint), &
                GravGradient(iSamplePoint,:,:)

enddo
close(69)
close(68)

deallocate(samplePoints)
deallocate(gravityVec)
deallocate(Potential)
deallocate(GravGradient)

end subroutine


!-------------------------------------------------------------------------------



subroutine xcrosssec(cmd,ncmd,ierr)
!----
! Plot nodal or elemental quantities (interpolated) on a planar cross-section
!----

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use vtkmodule,       only: gmshRead
use constants,       only: SampleDisplacement, &
                           SampleVelocity, &
                           SampleTemperature, &
                           SampleStress, &
                           SampleStrain, &
                           SampleRotation, &
                           SampleStrainRate, &
                           SampleHeatFlux
use iomodule,        only: ifein, ifeout, stderr
use errors,          only: exec

implicit none


integer :: ierr

integer :: ncmd
character :: cmd*(*)
dimension :: cmd(*)
double precision :: xsec_origin(3), xsec_normal(3), p(3), q(3), nh

integer :: idata, nQty

integer, allocatable :: doesThisSegmentIntersectXsec(:)
double precision, allocatable :: interpolationFactor(:)
double precision, allocatable :: interpolatedCoord(:,:)

integer :: iFile, iNode, iLine, iElem
integer :: i, j, k, n
integer :: np1, np2
double precision :: u(3), v(3), uproj, vproj
integer :: nIntersect
integer :: intersectSum
integer :: quadOrder(4)
double precision :: quadPt(3,4)
double precision :: qty(6)
double precision :: factor
integer :: whichLineInElem(6)


! Make sure the tecin.dat files and mesh file are read
if (ifein.eq.0) then
    write(stderr,*) 'xsec: read FE input file first'
    ierr = EXEC
    return
endif

if (ifeout.eq.0) then
    write(stderr,*) 'xsec: read FE output file first'
    ierr = EXEC
    return
endif

if (.not. gmshRead) then
    write(stderr,*) 'xsec: read mesh file first'
    ierr = EXEC
    return
endif

! Read the cross-section parameters (format: xsec d|v|T|s|e|r|f x0,y0,z0 nx,ny,nz [file])
if (ncmd.lt.4) then
    write(0,*) 'xsec: missing arguments'
    write(0,*) 'Usage: xsec d|v|s|e x0,y0,z0 nx,ny,nz [file]'
    return
endif

if (cmd(2).eq.'d') then
    idata = SampleDisplacement
elseif (cmd(2).eq.'v') then
    idata = SampleVelocity
elseif (cmd(2).eq.'s') then
    idata = SampleStress
elseif (cmd(2).eq.'e') then
    idata = SampleStrain
else
    write(0,*) 'xsec: sorry, no option to get cross section for quantity ',trim(cmd(2)),' yet'
    return
endif

i = index(cmd(3),',')
cmd(3)(i:i) = ' '
i = index(cmd(3),',')
cmd(3)(i:i) = ' '
read(cmd(3),*,err=9001,iostat=ierr) (xsec_origin(i),i=1,3)
9001 if (ierr.ne.0) then
    write(0,*) 'xsec: error parsing cross-section origin coordinates ',trim(cmd(3))
    return
endif

i = index(cmd(4),',')
cmd(4)(i:i) = ' '
i = index(cmd(4),',')
cmd(4)(i:i) = ' '
read(cmd(4),*,err=9002,iostat=ierr) (xsec_normal(i),i=1,3)
9002 if (ierr.ne.0) then
    write(0,*) 'xsec: error parsing cross-section normal vector ',trim(cmd(3))
    return
endif

! Normalize the normal vector to cross-section plane
if (abs(xsec_normal(1)).lt.1.0d-6 .and. &
    abs(xsec_normal(2)).lt.1.0d-6 .and. &
    abs(xsec_normal(3)).lt.1.0d-6) then
    write(0,*) 'xsec: normal vector must be non-zero'
    return
endif
factor = xsec_normal(1)*xsec_normal(1)+xsec_normal(2)*xsec_normal(2)+xsec_normal(3)*xsec_normal(3)
factor = sqrt(factor)
xsec_normal(1) = xsec_normal(1)/factor
xsec_normal(2) = xsec_normal(2)/factor
xsec_normal(3) = xsec_normal(3)/factor

! Projected coordinates (p in strike-parallel direction, q in up-dip direction)
nh = sqrt(xsec_normal(1)*xsec_normal(1)+xsec_normal(2)*xsec_normal(2))
if (nh.lt.1e-6) then
    p(1) = 1.0
    p(2) = 0.0
    p(3) = 0.0
    q(1) = 0.0
    q(2) = 1.0
    q(3) = 0.0
else
    p(1) = -xsec_normal(2)
    p(2) = xsec_normal(1)
    p(3) = 0.0
    q(1) = -xsec_normal(3)*xsec_normal(1)/nh
    q(2) = -xsec_normal(3)*xsec_normal(2)/nh
    q(3) = nh
endif


if (ncmd.eq.5) then
    iFile = 57
    open(unit=iFile,file=cmd(5),status='unknown')
else
    iFile = 6
endif


! Specify the number of values in each quantity
nQty = 0
if (idata.eq.SampleDisplacement) then
    nQty = NSD
elseif (idata.eq.SampleVelocity) then
    nQty = NSD
elseif (idata.eq.SampleStress) then
    nQty = NSTR
elseif (idata.eq.SampleStrain) then
    nQty = NSTR
endif
                           !SampleTemperature, &
                           !SampleRotation, &
                           !SampleStrainRate, &
                           !SampleHeatFlux


! Determine which line segments intersect the cross-section
! Store the coordinate of the intersection and the interpolation factor
allocate(doesThisSegmentIntersectXsec(meshdatactx%nLineSegments))
allocate(interpolatedCoord(5,meshdatactx%nLineSegments))
allocate(interpolationFactor(meshdatactx%nLineSegments))
doesThisSegmentIntersectXsec = 0
interpolatedCoord = 0.0d0
interpolationFactor = 0.0d0

! Possible cases for intersections:
!     0: Line does not intersect plane at all
!     1: Line intersects cross-section
!     10: One point lies in cross-section
!     100: Both points lie in cross-section
nIntersect = 0
do iLine = 1,meshdatactx%nLineSegments
    np1 = meshdatactx%lineSegments(1,iLine)
    np2 = meshdatactx%lineSegments(2,iLine)

    ! Vectors from cross-section origin to points of line segment
    do j = 1,NSD
        u(j) = meshdatactx%X(j,np1) - xsec_origin(j)
        v(j) = meshdatactx%X(j,np2) - xsec_origin(j)
    enddo

    ! Project these vectors onto the cross-section normal
    uproj = u(1)*xsec_normal(1) + u(2)*xsec_normal(2) + u(3)*xsec_normal(3)
    vproj = v(1)*xsec_normal(1) + v(2)*xsec_normal(2) + v(3)*xsec_normal(3)

    ! Both points lie in the plane
    if (abs(uproj).lt.1e-3.and.abs(vproj).lt.1e-3) then
        nIntersect = nIntersect+1
        doesThisSegmentIntersectXsec(iLine) = 100
        factor = 99

    ! One point lies in the plane (keep track of which point with sign)
    elseif (abs(uproj).lt.1e-3) then
        nIntersect = nIntersect+1
        doesThisSegmentIntersectXsec(iLine) = -10
        factor = (-uproj)/(vproj-uproj)
    elseif (abs(vproj).lt.1e-3) then
        nIntersect = nIntersect+1
        doesThisSegmentIntersectXsec(iLine) = 10
        factor = (-uproj)/(vproj-uproj)

    ! Line segment crosses plane
    elseif (uproj*vproj.lt.0) then
        nIntersect = nIntersect+1
        doesThisSegmentIntersectXsec(iLine) = 1
        factor = (-uproj)/(vproj-uproj)

    else
        ! Does not intersect
    endif

    ! Store some things for line segments that intersect plane
    if (abs(doesThisSegmentIntersectXsec(iLine)).ne.0) then
        ! Calculate and store the coordinate of the intersection
        do j = 1,NSD
            interpolatedCoord(j,iLine) = &
                meshdatactx%X(j,np1) + (meshdatactx%X(j,np2)-meshdatactx%X(j,np1))*factor
        enddo

        ! Also calculate the coordinate on the projected coordinate system (p,q)
        ! Vector from origin to projected coordinate
        do j = 1,3
            u(j) = interpolatedCoord(j,iLine) - xsec_origin(j)
        enddo
        ! Project this vector onto p, q axes
        interpolatedCoord(4,iLine) = u(1)*p(1) + u(2)*p(2) + u(3)*p(3)
        interpolatedCoord(5,iLine) = u(1)*q(1) + u(2)*q(2) + u(3)*q(3)

        ! Store the interpolation factor for this line segment
        interpolationFactor(iLine) = factor
    endif
enddo
if (nIntersect.eq.0) then
    write(0,*) 'xsec: found no line segments that intersect plane'
    if (iFile.eq.57) then
        close(iFile)
    endif
    return
endif
!do i = 1,meshdatactx%nLineSegments
!    write(0,*) i,doesThisSegmentIntersectXsec(i)
!enddo


! NOTE: COULD PROBABLY DO THIS WITHOUT STORING LINE SEGMENT INFORMATION BY STARTING
! AT THIS POINT AND COMPUTING INTERSECTIONS FOR EACH NODE PAIR IN EACH ELEMENT. THE
! UPSIDE OF STORING THE LINE SEGMENT INFORMATION IS THAT THERE ARE FEWER LINE SEGMENTS
! THAN ELEMENTS, SO IT SHOULD SAVE A FEW PROCESSING LOOPS (BUT I AM NOT SURE OF THIS).
! THE DOWNSIDE IS USING MORE MEMORY AND THE OVERHEAD OF LOADING THESE ARRAYS.


! Check which elements appear in the cross-section and print the desired quantity for them
! They can intersect the plane in different ways, which we need to handle differently
! Possible cases:
!     - Three lines intersect plane (forms a triangle): intersectSum=3
!     - Four lines intersect plane (quadrilateral):     intersectSum=4
!     - One point lies in plane:                        intersectSum=30 (DO NOT PLOT!)
!     - One point lies in plane, two lines intersect:   intersectSum=32
!     - One line lies in plane:                         intersectSum=140 (DO NOT PLOT!)
!     - One line lies in plane, one line intersects:    intersectSum=141
!     - Three lines lie in plane:                       intersectSum=330
do iElem = 1,meshdatactx%neglobal

    ! The intersectSum variable encodes the geometry of intersection
    intersectSum = 0
    do j = 1,6
        iLine = meshdatactx%linesInElements(j,iElem)
        intersectSum = intersectSum + abs(doesThisSegmentIntersectXsec(iLine))
    enddo

    whichLineInElem = 0
    n = 0

    ! Three lines intersect the cross-section plane to form a triangle
    if (intersectSum.eq.3) then
        ! Get the lines that intersect
        do j = 1,6
            iLine = meshdatactx%linesInElements(j,iElem)
            if (doesThisSegmentIntersectXsec(iLine).eq.1) then
                n = n + 1
                whichLineInElem(n) = j
            endif
            if (n.eq.3) then
                exit
            endif
        enddo
        ! Write coordinates and quantities on cross-section
        do n = 1,3
            j = whichLineInElem(n)
            iLine = meshdatactx%linesInElements(j,iElem)
            call getInterpolatedQuantity(qty,idata,iLine,iElem,interpolationFactor(iLine))
            write(iFile,*) (interpolatedCoord(i,iLine),i=1,3),(qty(i),i=1,nQty),&
                           (interpolatedCoord(i,iLine),i=4,5),&
                           qty(1)*p(1)+qty(2)*p(2)+qty(3)*p(3),&
                           qty(1)*q(1)+qty(2)*q(2)+qty(3)*q(3),&
                           qty(1)*xsec_normal(1)+qty(2)*xsec_normal(2)+qty(3)*xsec_normal(3)
        enddo

    ! Four lines intersect the cross-section plane to form a quadrilateral
    ! Divide it into two triangles for plotting
    elseif (intersectSum.eq.4) then
        ! Get the four lines that intersect
        do j = 1,6
            iLine = meshdatactx%linesInElements(j,iElem)
            if (doesThisSegmentIntersectXsec(iLine).eq.1) then
                n = n + 1
                whichLineInElem(n) = j
                quadPt(1,n) = interpolatedCoord(1,iLine)
                quadPt(2,n) = interpolatedCoord(2,iLine)
                quadPt(3,n) = interpolatedCoord(3,iLine)
            endif
            if (n.eq.4) then
                exit
            endif
        enddo
        ! Put the points in order
        call getQuadOrder(quadOrder,quadPt)
        ! Write coordinates and quantities on cross-section for triangle points 123
        do n = 1,3
            j = whichLineInElem(quadOrder(n))
            iLine = meshdatactx%linesInElements(j,iElem)
            call getInterpolatedQuantity(qty,idata,iLine,iElem,interpolationFactor(iLine))
            write(iFile,*) (interpolatedCoord(i,iLine),i=1,3),(qty(i),i=1,nQty),&
                           (interpolatedCoord(i,iLine),i=4,5),&
                           qty(1)*p(1)+qty(2)*p(2)+qty(3)*p(3),&
                           qty(1)*q(1)+qty(2)*q(2)+qty(3)*q(3),&
                           qty(1)*xsec_normal(1)+qty(2)*xsec_normal(2)+qty(3)*xsec_normal(3)
        enddo
        ! Write coordinates and quantities on cross-section for triangle points 134
        do n = 1,4
            if (n.eq.2) then
                cycle
            endif
            j = whichLineInElem(quadOrder(n))
            iLine = meshdatactx%linesInElements(j,iElem)
            call getInterpolatedQuantity(qty,idata,iLine,iElem,interpolationFactor(iLine))
            write(iFile,*) (interpolatedCoord(i,iLine),i=1,3),(qty(i),i=1,nQty),&
                           (interpolatedCoord(i,iLine),i=4,5),&
                           qty(1)*p(1)+qty(2)*p(2)+qty(3)*p(3),&
                           qty(1)*q(1)+qty(2)*q(2)+qty(3)*q(3),&
                           qty(1)*xsec_normal(1)+qty(2)*xsec_normal(2)+qty(3)*xsec_normal(3)
        enddo

    ! One point intersects the plane
    elseif (intersectSum.eq.30) then
        ! do not plot

    ! One point and two lines intersect the plane, forming a triangle
    elseif (intersectSum.eq.32) then
        ! Get the two lines that intersect
        do j = 1,6
            iLine = meshdatactx%linesInElements(j,iElem)
            if (doesThisSegmentIntersectXsec(iLine).eq.1) then
                n = n + 1
                whichLineInElem(n) = j
            endif
            if (n.eq.2) then
                exit
            endif
        enddo
        ! Get the point that intersects
        ! Once we find it, we can exit the loop
        do j = 1,6
            iLine = meshdatactx%linesInElements(j,iElem)
            if (abs(doesThisSegmentIntersectXsec(iLine)).eq.10) then
                n = n + 1
                whichLineInElem(n) = j
                exit
            endif
        enddo
        ! Write coordinates and quantities on cross-section
        do n = 1,3
            j = whichLineInElem(n)
            iLine = meshdatactx%linesInElements(j,iElem)
            call getInterpolatedQuantity(qty,idata,iLine,iElem,interpolationFactor(iLine))
            write(iFile,*) (interpolatedCoord(i,iLine),i=1,3),(qty(i),i=1,nQty),&
                           (interpolatedCoord(i,iLine),i=4,5),&
                           qty(1)*p(1)+qty(2)*p(2)+qty(3)*p(3),&
                           qty(1)*q(1)+qty(2)*q(2)+qty(3)*q(3),&
                           qty(1)*xsec_normal(1)+qty(2)*xsec_normal(2)+qty(3)*xsec_normal(3)
        enddo

    ! One line lies in the plane
    elseif (intersectSum.eq.140) then
        ! do not plot

    ! One line lies in the plane and another intersects, forming a triangle
    elseif (intersectSum.eq.141) then
        ! Get the two points that lie in the plane
        ! We should only get each node once, so keep track of which one we grabbed first
        iNode = 0
        do j = 1,6
            iLine = meshdatactx%linesInElements(j,iElem)
            if (doesThisSegmentIntersectXsec(iLine).eq.-10) then
                if (n.eq.0) then
                    n = n + 1
                    iNode = meshdatactx%lineSegments(1,iLine)
                    whichLineInElem(n) = j
                else
                    if (iNode.ne.meshdatactx%lineSegments(1,iLine)) then
                        n = n + 1
                        whichLineInElem(n) = j
                    endif
                endif
            elseif (doesThisSegmentIntersectXsec(iLine).eq.10) then
                if (n.eq.0) then
                    n = n + 1
                    iNode = meshdatactx%lineSegments(2,iLine)
                    whichLineInElem(n) = j
                else
                    if (iNode.ne.meshdatactx%lineSegments(2,iLine)) then
                        n = n + 1
                        whichLineInElem(n) = j
                    endif
                endif
            endif
            if (n.eq.2) then
                exit
            endif
        enddo
        ! Get the line segment intersection
        ! There is only one of these
        do j = 1,6
            iLine = meshdatactx%linesInElements(j,iElem)
            if (doesThisSegmentIntersectXsec(iLine).eq.1) then
                n = n + 1
                whichLineInElem(n) = j
                exit
            endif
        enddo
        ! Write coordinates and quantities on cross-section
        do n = 1,3
            j = whichLineInElem(n)
            iLine = meshdatactx%linesInElements(j,iElem)
            call getInterpolatedQuantity(qty,idata,iLine,iElem,interpolationFactor(iLine))
            write(iFile,*) (interpolatedCoord(i,iLine),i=1,3),(qty(i),i=1,nQty),&
                           (interpolatedCoord(i,iLine),i=4,5),&
                           qty(1)*p(1)+qty(2)*p(2)+qty(3)*p(3),&
                           qty(1)*q(1)+qty(2)*q(2)+qty(3)*q(3),&
                           qty(1)*xsec_normal(1)+qty(2)*xsec_normal(2)+qty(3)*xsec_normal(3)
        enddo

    ! Three lines lie in the plane, forming a triangle
    elseif (intersectSum.eq.330) then
        ! Get the three points that lie in the plane
        ! This is nice and easy because we can just get the endpoints of line segments that stick out of the plane
        do j = 1,6
            iLine = meshdatactx%linesInElements(j,iElem)
            if (abs(doesThisSegmentIntersectXsec(iLine)).eq.10) then
                n = n + 1
                whichLineInElem(n) = j
            endif
            if (n.eq.3) then
                exit
            endif
        enddo
        ! Write coordinates and quantities on cross-section
        do n = 1,3
            j = whichLineInElem(n)
            iLine = meshdatactx%linesInElements(j,iElem)
            call getInterpolatedQuantity(qty,idata,iLine,iElem,interpolationFactor(iLine))
            write(iFile,*) (interpolatedCoord(i,iLine),i=1,3),(qty(i),i=1,nQty),&
                           (interpolatedCoord(i,iLine),i=4,5),&
                           qty(1)*p(1)+qty(2)*p(2)+qty(3)*p(3),&
                           qty(1)*q(1)+qty(2)*q(2)+qty(3)*q(3),&
                           qty(1)*xsec_normal(1)+qty(2)*xsec_normal(2)+qty(3)*xsec_normal(3)
        enddo

    elseif (intersectSum.ne.0) then
        write(0,*) 'xsec: I am not sure how you got intersectSum=',intersectSum,' in element',iElem
        write(0,*) '      Moving on to the next element'

    else
    endif
enddo

! Free memory back up
deallocate(doesThisSegmentIntersectXsec)
deallocate(interpolatedCoord)
deallocate(interpolationFactor)

if (iFile.eq.57) then
    close(iFile)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine getInterpolatedQuantity(qty,idata,iLine,iElem,factor)
!----
! Calculate the desired quantity (specified by idata) at the point where a line segment (iLine)
! intersects a plane. The quantity may also depend on the element (iElem). Some of the work has
! already been done, in calculating the interpolation factor (factor).
!----

use meshdatamodule, only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMSLPglobal
use modeltopology, only: NSD,NDOF,NSTR,NEN
use constants, only: SampleDisplacement, &
                     SampleVelocity, &
                     SampleTemperature, &
                     SampleStress, &
                     SampleStrain, &
                     SampleRotation, &
                     SampleStrainRate, &
                     SampleHeatFlux
implicit none
! pass
integer :: idata, iLine, iElem
double precision :: factor, qty(6)
! local
integer :: np1, np2, nen1, nen2, j, n

! Calculate the quantity interpolated onto the cross-section plane
qty = 0.0d0

! Nodes making up line segment
np1 = meshdatactx%lineSegments(1,iLine)
np2 = meshdatactx%lineSegments(2,iLine)

! Displacements
if (idata.eq.SampleDisplacement) then
    do j = 1,NSD
        qty(j) = &
            modeldatactx%D(j,np1) + (modeldatactx%D(j,np2)-modeldatactx%D(j,np1))*factor
    enddo

    if (NUMSLPglobal.gt.0) then
        do n = 1,NEN
            if (meshdatactx%ien(n,iElem).eq.np1) then
                nen1 = n
            endif
            if (meshdatactx%ien(n,iElem).eq.np2) then
                nen2 = n
            endif
        enddo
        do j = 1,NSD
            qty(j) = qty(j) + &
                     modeldatactx%DXE(j,nen1,iElem) + &
                     (modeldatactx%DXE(j,nen2,iElem)-modeldatactx%DXE(j,nen1,iElem))*factor
        enddo
    endif

! Velocities
elseif (idata.eq.SampleVelocity) then
    do j = 1,NSD
        qty(j) = &
            modeldatactx%DELD(j,np1) + (modeldatactx%DELD(j,np2)-modeldatactx%DELD(j,np1))*factor
    enddo

    if (NUMSLPglobal.gt.0) then
        do n = 1,NEN
            if (meshdatactx%ien(n,iElem).eq.np1) then
                nen1 = n
            endif
            if (meshdatactx%ien(n,iElem).eq.np2) then
                nen2 = n
            endif
        enddo
        do j = 1,NSD
            qty(j) = qty(j) + &
                     modeldatactx%DELDXE(j,nen1,iElem) + &
                     (modeldatactx%DELDXE(j,nen2,iElem)-modeldatactx%DELDXE(j,nen1,iElem))*factor
        enddo
    endif

! Stresses
elseif (idata.eq.SampleStress) then
    do j = 1,NSTR
        qty(j) = modeldatactx%STN(j,iElem)
    enddo

! Strains
elseif (idata.eq.SampleStrain) then
    do j = 1,NSTR
        qty(j) = modeldatactx%STR(j,iElem)
    enddo

! Sorry, Dave. I'm afraid I can't do that.
else
    write(0,*) 'This quantity has not been implemented in xsec yet...stopping'
    stop
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine getQuadOrder(order,points)
!----
! Arrange four planar points in spatial order
! This workflow assumes the points form a convex quadrilateral
!----
implicit none
! pass
integer :: order(4)
double precision :: points(3,4)
! local
double precision :: x12,y12,z12,x23,y23,z23,x34,y34,z34,x41,y41,z41
double precision :: cpx1,cpy1,cpz1,cpx2,cpy2,cpz2,cpx3,cpy3,cpz3
double precision :: dp1,dp2

order = 0

! Vector 1->2
x12 = points(1,2)-points(1,1)
y12 = points(2,2)-points(2,1)
z12 = points(3,2)-points(3,1)
! Vector 2->3
x23 = points(1,3)-points(1,2)
y23 = points(2,3)-points(2,2)
z23 = points(3,3)-points(3,2)
! Vector 3->4
x34 = points(1,4)-points(1,3)
y34 = points(2,4)-points(2,3)
z34 = points(3,4)-points(3,3)
! Vector 4->1
x41 = points(1,1)-points(1,4)
y41 = points(2,1)-points(2,4)
z41 = points(3,1)-points(3,4)

! Cross-product of 1->2 and 2->3
cpx1 = y12*z23 - z12*y23
cpy1 = z12*x23 - x12*z23
cpz1 = x12*y23 - y12*x23
! Cross-product of 2->3 and 3->4
cpx2 = y23*z34 - z23*y34
cpy2 = z23*x34 - x23*z34
cpz2 = x23*y34 - y23*x34
! Cross-product of 3->4 and 4->1
cpx3 = y34*z41 - z34*y41
cpy3 = z34*x41 - x34*z41
cpz3 = x34*y41 - y34*x41

! Dot product of cross-product 1 and 2
dp1 = cpx1*cpx2 + cpy1*cpy2 + cpz1*cpz2
! Dot product of cross-product 2 and 3
dp2 = cpx2*cpx3 + cpy2*cpy3 + cpz2*cpz3

! Dot products determine ordering of nodes
if (dp1.gt.0.and.dp2.gt.0) then
    order(1) = 1
    order(2) = 2
    order(3) = 3
    order(4) = 4
elseif (dp1.lt.0.and.dp2.gt.0) then
    order(1) = 1
    order(2) = 2
    order(3) = 4
    order(4) = 3
elseif (dp1.gt.0.and.dp2.lt.0) then
    order(1) = 1
    order(2) = 3
    order(3) = 2
    order(4) = 4
else
    write(0,*) '!! Problem in getQuadOrder. Stopping'
    stop
endif

return
end subroutine







!subroutine getXsecTaggedSurfaces()
!use meshdatamodule, only: meshdatactx
!use vtkmodule, only: outputTags, nOutputTags, gmshRead
!implicit none
!integer :: iTriangle, iTag, pos
!logical :: includeThisOne
!
!do iTriangle = 1,meshdatactx%nTriangles
!    ! for every triangle...
!    includeThisOne = .false.
!    do iTag=1,size(meshdatactx%triangleTags,1)
!        ! ... walk through the tags of this triangle, and see whether it is included.
!        call arraycontains(outputTags, &
!                           nOutputTags, &
!                           meshdatactx%triangleTags(iTag, iTriangle), &
!                           pos)
!        if(pos.le.nOutputTags) then
!            includeThisOne = .true.
!        endif
!    enddo
!    if (includeThisOne) then
!    endif
!enddo
!
!return
!end subroutine
