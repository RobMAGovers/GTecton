! this subroutine will read a Gmsh 4.0 file and put its output in a data structure
! which can be accessed for further postprocessing

! The connectivity and the nodal point coordinates have already been read
! by readfe. So all we need to do is to take out to lines and the triangles.

! we have, in meshdatamodule: 
!         for each line, its start and end point
!  PetscInt, ALLOCATABLE :: lineEndPointIDs(2,:)
!         for each line, its marker    
!  PetscInt, ALLOCATABLE :: lineMarker(:)

!         for each triangle, its marker
!  PetscInt, ALLOCATABLE :: triangleMarker(:)
!         for each triangle, its boundaring elements
!  PetscInt, ALLOCATABLE :: triangleElements(2,:)





subroutine readGmsh(meshFileName)

use modeldefinition, only: mode

use constants,  only: elemIsPoint, &
                      elemIsLine, &
                      elemIsTriangle, &
                      elemIsTetrahedron
use filemodule, only: meshFileID
use vtkmodule,  only: gmshRead
use iomodule,      only: stderr

implicit none

integer, external :: nextlu ! from tools
character(len=*)  :: meshFileName
character(len=255):: string

integer           :: gmshVersion
character(len=3)  :: gmshVersionString
integer           :: datasize, filetype
integer           :: ios

meshFileID = nextlu(0)

open(unit=meshFileID, file=meshFileName,iostat=ios)

if (ios.ne.0) then
    write(stderr,*) "could not open meshfile. ", trim(meshFileName)
    write(stderr,*) "Error: ", ios
    return
else
!    write(stderr,*) "opened ", meshFileName(1:11), "with ID", meshFileID
endif


! read the header: 
! read first line
read(meshFileID,*) string(1:11)

if (ios.ne.0) then
    write(stderr,*) "could not read meshfile. Error: ", ios
    return
endif

if (string(1:11).eq.'$MeshFormat') then
    ! read second line
    read(meshFileID,*,err=1020,end=1030,iostat=ios) gmshVersionString,&
                                             filetype,   &
                                             datasize
    if (gmshVersionString(1:3).eq."2.2") then
        gmshVersion = 2
    endif

    if (filetype.ne.0) then
        write(stderr,1)
1       format(1x,'gmsh2tecton: non-ASCII GMSH mesh file')
        call exitp(1)
    endif

    if (datasize.ne.8) then
        write(stderr,2)
2       format(1x,'gmsh2tecton: only 8-sized supported')
        call exitp(1)
    endif

    ! read third line
    read(meshFileID,'(a14)',err=1040,end=1050,iostat=ios) string
    if (string(1:14).ne."$EndMeshFormat") then
        write(stderr,3)
3       format(1x,'gmsh2tecton: "$MeshFormat" block not properly closed')
        call exitp(1)
    endif

    if (INT(gmshVersion).eq.2 .or. INT(gmshVersion).eq.3) then
        if (MODE.eq.5) then
            write(*,10) 2.2
10      format('gmsh file format:',F3.1)
        endif

        ! read the actual file:
        call rdmsh2 ()


    else if (INT(gmshVersion).eq.4) then
        write(*,12)
12      format(1x,'gmsh2tecton: support for gmsh file format 4',&
         ' is under development'/1x,&
         'For now, you can save in format 2.2 by ',&
         'adding "Mesh.MshFileVersion = 2.2;" to ',&
         'your "~/.gmsh-options" file')
         call exitp(1)
    else
        write(stderr,4) INT(gmshVersion)
4       format(1x,'gmsh2tecton: unknown mesh file format ',I2)
        call exitp(1)
    endif
else
    write(stderr,5)
5   format(1x,'gmsh2tecton: "$MeshFormat" header missing in mesh file')
    call exitp(1)
endif

call lineSegmentsInTetrahedra()

gmshRead = .true.

return

1000 write(stderr,*) "could not read meshfile 1. ", ios
return

1010 write(stderr,*) "prematurely reached end of mesh file 1", ios
return

1020 write(stderr,*) "could not read meshfile 2. ", ios
return

1030 write(stderr,*) "prematurely reached end of mesh file 2", ios
return

1040 write(stderr,*) "could not read meshfile 3. ", ios
return

1050 write(stderr,*) "prematurely reached end of mesh file 3", ios
return



end subroutine


subroutine rdmsh2 ()

!NODE,X,IEL,ITYPE,IPHYS,IEN,NUMNP,MAXNP,&
!                  NNEL,MAXEL,NUMAT,vertexTag,verbose,logfileID,MODE)

!
! read GMSH mesh data from Gmsh version 2.2.
! May also read other versions, but definitely not 4 or later.
!
use meshdatamodule, only: meshdatactx
use constants, only: elemIsPoint, &
                     elemIsLine, &
                     elemIsTriangle, &
                     elemIsTetrahedron
use filemodule, only: meshFileID
use algebra,    only: threeInFour
use iomodule,   only: stderr

    implicit none
!-pass
    integer NODE,IEL,ITYPE,IPHYS,IEN,NUMNP,MAXNP,NNEL,MAXEL,&
            NUMAT,vertexTag
    real X
!    dimension NODE(MAXNP),X(3,MAXNP),IEL(MAXEL),ITYPE(MAXEL),&
!         IPHYS(MAXEL),IEN(8,MAXEL),vertexTag(MaxNP)
    logical :: verbose
    integer :: logfileID,MODE

    
!-locl
    character(len=256) string
    character(len=256) elementString
    integer nsz(19),ios,i,n,j,numel,filetype,datasize,itag,lien(27),&
     numline,numtria
    real meshformat
    integer :: iPoint, iElem, iStep, iIndex, iTriangle
    integer :: ntags, tags(10)

    integer :: nPoints       
    integer :: nLines        
    integer :: nTetrahedrons 
    integer :: nTagsPoints       
    integer :: nTagsLines        
    integer :: nTagsTriangles    
    integer :: nTagsTetrahedrons
    

    integer :: dummyInt, dummyInt2, elemType
    logical :: markerFOund


    integer, allocatable :: sortedElemsIndex(:), newIndex(:)
    integer, allocatable :: sortedIen(:,:)
    integer :: highestLowestIndex
    integer, allocatable :: elemsPerLowestIndex(:)

    integer :: lowestTriangleNode
    logical :: gotya
    integer :: offset, startoffset, endoffset

    integer :: nElemsFound

    
    ! number of vertices per element.
        ! element 1 is a line -> 2 vertices
        ! element 2 is a triangle -> 3 vertices
        ! etc.
        data nsz/2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13/

    if (verbose) then
        write(0,*) "reading format version 2 mesh file"
    endif

    NUMNP = -1
    NNEL = -1

! we read the file line by line, and when we run across 
! header lines, we know what is being read, so that it
! may be put in the appropriate structure.


100 read(meshFileID,'(a256)',err=1000,end=200,iostat=ios) string

    ! nodes
    if (string(1:7).eq.'$Nodes ') then
        ! read number of GMSH nodes.
        read(meshFileID,*,err=1020,end=1030,iostat=ios) NUMNP

        

!        if (NUMNP.le.0) then
!            write(stderr,10)
!  10        format(1x,'gmsh2tecton: no nodal coordinates in mesh file')
!            call exitp(1)
!        else if (NUMNP.gt.MAXNP) then
!            write(stderr,20)
!  20        format(1x,'gmsh2tecton: MAXNP dimensioned too small')
!            call exitp(1)
!        endif

        !*** read the actual nodes
        do i=1,NUMNP
!            read(meshFileID,*,err=1040,end=1050,iostat=ios) &
!            NODE(i),X(1,i),X(2,i),X(3,i)
            read(meshFileID,*,err=1040,end=1050,iostat=ios) string
        enddo
        !***

        ! make sure that the node block is properly closed.
        read(meshFileID,'(a256)',err=1060,end=1070,iostat=ios) string

        if (string(1:10).ne.'$EndNodes ') then
            write(stderr,40)
  40        format(1x,'gmsh2tecton: mesh file "$EndNodes" missing')
            call exitp(1)
        endif

        if (verbose) then
            write(stderr,42) NUMNP
  42    format('gmsh2tecton read ',I7,' GMSH nodes from mesh file')
        endif

        if (MODE.eq.5) then
            write(*,44) NUMNP
  44    format('Number of GMSH nodes = ',I7)
        endif

    ! elements
    elseif (string(1:10).eq.'$Elements ') then
        read(meshFileID,*,err=1100,end=1110,iostat=ios) NNEL   ! number of GMSH elements

!        if (NNEL.le.0) then
!            write(stderr,50)
!  50        format(1x,'gmsh2tecton: no GMSH elements in mesh file')
!            call exitp(1)
!        else if (NNEL.gt.MAXEL) then
!            write(stderr,60)
!  60        format(1x,'gmsh2tecton: MAXEL dimensioned too small')
!            call exitp(1)
!        endif

        ! GMSH calls any (1D,2D,3D) connectivity an element. GTECTON is only
        ! interested in the 3D tetrahedry.

        !*** read the actual elements
        ! element notation is as follows:
        ! elm-number elm-type number-of-tags < tags > node-number-list

        ! we do this in two passes.

        ! On the first pass, we count the numbers of points and lines and triangles and elements
        ! For each type, we also monitor the largest numbers of tags.
        ! Then we allocate the proper size of memory

        nPoints       = 0
        nLines        = 0
        meshdatactx%nTriangles    = 0
        nTetrahedrons = 0
        nTagsPoints       = 0
        nTagsLines        = 0
        nTagsTriangles    = 0
        nTagsTetrahedrons = 0

        do i=1,NNEL
!!!         Potential out-of-array boundaries nsz problem here for new GMSH element types
            read(meshFileID,*,err=1120,end=1130,iostat=ios) dummyInt,&
                                                     elemType,&
                                                     ntags,  &
                                                     (tags(j),j=1,ntags), &
                                                     (lien(j),j=1,nsz(elemType))

            if (elemType.eq.elemIsPoint) then
                nPoints = nPoints + 1
                if (ntags.gt.nTagsPoints) then
                    nTagsPoints = ntags
                endif
            else if (elemType.eq.elemIsLine) then
                nLines = nLines + 1
                 if (ntags.gt.nTagsLines) then
                    nTagsLines = ntags
                endif
            else if (elemType.eq.elemIsTriangle) then
                meshdatactx%nTriangles = meshdatactx%nTriangles + 1
                if (ntags.gt.nTagsTriangles) then
                    nTagsTriangles = ntags
                endif
            else if (elemType.eq.elemIsTetrahedron) then
                nTetrahedrons = nTetrahedrons + 1
                if (ntags.gt.nTagsTetrahedrons) then
                    nTagsTetrahedrons = ntags
                endif
            else
                write(stderr,*) "gmshReader: element type ", elemType, &
                                "not supported by GTECTON"
                call exitp(1)
            endif
        enddo

!        write(*,*) "counted nPoints", nPoints
!        write(*,*) "counted nTagsPoints: ", nTagsPoints

!        write(*,*) "counted nLines", nLines
!        write(*,*) "counted nTagsLines: ", nTagsLines

!        write(*,*) "counted nTriangles", meshdatactx%nTriangles
!        write(*,*) "counted nTagsTriangles: ", nTagsTriangles

!        write(*,*) "counted nTetrahedrons", nTetrahedrons
!        write(*,*) "counted nTagsTetraherons: ", nTagsTetrahedrons


        ! allocate memory for the lines and the triangles and the tags
        allocate(meshdatactx%pointTags(nTagsPoints, nPoints))
        meshdatactx%pointTags = 0

        allocate(meshdatactx%lineEndPointIDs(2,nLines))
        allocate(meshdatactx%lineTags(nTagsLines, nLines))
        meshdatactx%lineEndPointIDs = 0        
        meshdatactx%lineTags = 0

        allocate(meshdatactx%triangleConnectivity(3,meshdatactx%nTriangles))
        allocate(meshdatactx%triangleElements(2,meshdatactx%nTriangles))
!        allocate(meshdatactx%triangleElementSide(2,meshdatactx%nTriangles))
        allocate(meshdatactx%triangleTags(nTagsTriangles,meshdatactx%nTriangles))
        meshdatactx%triangleConnectivity = 0
        meshdatactx%triangleElements = 0
        meshdatactx%triangleTags = 0

        allocate(meshdatactx%tetrahedronTags(nTagsTetrahedrons,nTetrahedrons))
        meshdatactx%tetrahedronTags = 0


!************************************************************
        rewind(meshFileID)
!************************************************************


        ! read until we reach the elements marker in the mesh file again
        markerFound = .false.
        do while (.not. markerFound)
            read(meshFileID,'(a256)',err=1000,end=200,iostat=ios) string
            if (string(1:10).eq.'$Elements ') then
                markerFound = .true.
            endif
        enddo

        ! reread number of elements.
        read(meshFileID,*,err=1120,end=1130,iostat=ios) elementString


        ! we are now back where we were.
        ! read the elements again.

        ! In the second pass here we actually fill the array with the values we read.

        nPoints = 0
        nLines = 0
        meshdatactx%nTriangles = 0
        nTetrahedrons = 0

        do i=1,NNEL

            ! set them to 0, becaus the next read may be of a smaller
            ! element, so that not all entries in tags and lien will
            ! will be overwritten.
            tags = 0
            lien = 0

            read(meshFileID,*,err=1120,end=1130,iostat=ios) dummyInt, &
                                  elemType, &
                                  ntags, &
                                 (tags(j),j=1,ntags), &
                                 (lien(j),j=1,nsz(elemType))


            if (elemType.eq.elemIsPoint) then
                ! only administrate the tags, because the nodes are already 
                ! in meshdatactx%x
                nPoints = nPoints + 1
                do iTag=1,ntags
                    meshdatactx%pointTags(iTag, nPoints) = tags(iTag)
                enddo
            else if (elemType.eq.elemIsLine) then
                nLines = nLines + 1
                meshdatactx%lineEndPointIDs(:,nLines) = lien(1:2)
                do iTag=1,ntags
                    meshdatactx%lineTags(iTag, nLines) = tags(iTag)
                enddo
            else if (elemType.eq.elemIsTriangle) then
                meshdatactx%nTriangles = meshdatactx%nTriangles + 1
                meshdatactx%triangleConnectivity(:,meshdatactx%nTriangles) = lien(1:3)
                do iTag=1,ntags
                    meshdatactx%triangleTags(iTag, meshdatactx%nTriangles) = tags(iTag)
                enddo
            else if (elemType.eq.elemIsTetrahedron) then
                ! only administrate the tags, because the elements are already
                ! in meshdatactx%ien
                nTetrahedrons = nTetrahedrons + 1
                do iTag=1,ntags
                    meshdatactx%tetrahedronTags(iTag, nTetrahedrons) = tags(iTag)
                enddo
            else
                 write(stderr,*) "gmshReader: element type ", elemType, &
                                 "not supported by GTECTON"
                call exitp(1)
            endif
        enddo

        ! find the element side (triangle) corresponding to an element (tet)
        call TriangleTetrahedronRelation()


!       Set vertex tags first by surfaces, then by lines, then by points
!        do i=1,NNEL
!            if (ITYPE(i).eq.2 .and. IPHYS(i).ne.0) then
!                do iPoint=1,nsz(ITYPE(i))
!                    vertexTag(IEN(iPoint,i)) = IPHYS(i)
!                enddo
!            endif
!        enddo

!        do i=1,NNEL
!            if (ITYPE(i).eq.1 .and. IPHYS(i).ne.0) then
!                do iPoint=1,nsz(ITYPE(i))
!                    vertexTag(IEN(iPoint,i)) = IPHYS(i)
!                enddo
!            endif
!        enddo

!        do i=1,NNEL
!            if (ITYPE(i).eq.15 .and. IPHYS(i).ne.0) then
!                do iPoint=1,nsz(ITYPE(i))
!                    vertexTag(IEN(iPoint,i)) = IPHYS(i)
!                enddo
!            endif
!        enddo

!       read data block footer
        read(meshFileID,'(a12)',err=1140,end=1150,iostat=ios) string

        if (string(1:13).ne.'$EndElements ') then
            write(stderr,80)
  80        format(1x,'gmsh reader: "$Elements" block not properly closed')
            call exitp(1)
        endif

        if (nTetrahedrons.eq.0) then
            write(stderr,90)
  90        format('gmsh reader: no tetrahedral elements in mesh file')
            call exitp(1)
        endif

        if (verbose) then
            write(stderr,92) nTetrahedrons
  92    format('gmsh reader read ',I12,' tetrahedral elements from mesh file')
        endif

        if (MODE.eq.5) then
            write(*,94) NNEL,numline,numtria,NUMEL,NUMAT
  94        format('Number of GMSH elements = ',I12/&
                   'Number of line elements = ',I12/&
                   'Number of triangle elements = ',I12/&
                   'Number of tetrahedral elements = ',I12/&
                   'Number of materials = ',I12)
        endif

        if (NUMAT.eq.0) then
            write(stderr,95)
  95    format('gmsh2tecton: NUMAT=0, so you either did not define ',&
         'Physical Volumes yet, or you wrote your mesh file while ',&
         'Mesh.SaveAll=1')
        endif

    elseif (string(1:15).eq.'$PhysicalNames ') then
        ! names given to tags
        read(meshFileID,*,err=1140,end=1150,iostat=ios) j
        if (j.gt.0 .and. MODE.eq.5) then
            write(*,96) j
        endif
  96    format('Number of physical names = ',I12)
        string = '$End' // string(2:15)

        if (MODE.eq.5) then
            call skipblock(trim(string),.true.)
        else
            call skipblock(trim(string),.false.)
        endif

    elseif (string(1:10).eq.'$Periodic ') then

        read(meshFileID,*,err=1140,end=1150,iostat=ios) j

        if (j.gt.0 .and. MODE.eq.5) then
            write(*,98) j
        endif

  98    format('Number of periodic entities = ',I3)
        string = '$End' // string(2:252)
        call skipblock(trim(string),.false.)
    elseif (string(1:1).eq.'$') then

        if (verbose) then
            write(stderr,*) 'skipping msh block with ',trim(string)
        endif

        string = '$End' // string(2:252)
        call skipblock(trim(string),.true.)
    endif
    goto 100
!
200 if (NUMNP.eq.-1) then
        write(stderr,*) 'gmsh2tecton: no nodal points block in mesh file'
        call exitp(1)
    endif
    if (NNEL.eq.-1) then
        write(stderr,*) 'gmsh2tecton: no elements block in mesh file'
        call exitp(1)
    endif
    return
!
1000 write(stderr,1001)
1001 format(1x,'gmsh2tecton: cannot read block header information')
     call exitp(1)

1020 write(stderr,1021) ios
1021 format(1x,'gmsh2tecton: read error ',I5,' on NUMNP')
     call exitp(1)

1030 write(stderr,1031)
1031 format(1x,'gmsh2tecton: premature EOF while reading NUMNP')
     call exitp(1)

1040 write(stderr,1041) ios
1041 format(1x,'gmsh2tecton: read error ',I5,' on node coordinates')
     call exitp(1)

1050 write(stderr,1051)
1051 format(1x,'gmsh2tecton: premature EOF while reading node coordinates')
     call exitp(1)

1060 write(stderr,1061) ios
1061 format(1x,'gmsh2tecton: read error ',I5,' on "$EndNodes"')
     call exitp(1)

1070 write(stderr,1071)
1071 format(1x,'gmsh2tecton: premature EOF while reading "$EndNodes"')
     call exitp(1)

1100 write(stderr,1101) ios
1101 format(1x,'gmsh2tecton: read error ',I5,' on number of elements')
     call exitp(1)

1110 write(stderr,1111)
1111 format(1x,'gmsh2tecton: premature EOF while reading number of elements')
     call exitp(1)

1120 write(stderr,1121) ios
1121 format(1x,'gmsh2tecton: read error ',I5,' on element data')
     call exitp(1)

1130 write(stderr,1131)
1131 format(1x,'gmsh2tecton: premature EOF while reading element data')
     call exitp(1)

1140 write(stderr,1141) ios
1141 format(1x,'gmsh2tecton: read error ',I5,' on "$EndElements"')
     call exitp(1)

1150 write(stderr,1151)
1151 format(1x,'gmsh2tecton: premature EOF while reading "$EndElements"')
     call exitp(1)
!
    end subroutine

!-------------------------------------------------------------------------------
    subroutine skipblock(string,echo)

    use filemodule, only: meshFileID

    implicit none
!-pass
    logical echo
    character(len=*) string

!    integer stdin,stdout,stderr,lui,luo
!-locl
    integer ios
    character(len=256) input

100 read(meshFileID,'(a256)',err=300,end=200,iostat=ios) input
    if (trim(input).eq.string) return
    if (echo) write(*,10) trim(input)
10  format(256a)
    goto 100

200 write(0,1) trim(string)
  1 format(1x,'gmsh2tecton: section end not found in mesh file: ',a256)
    call exitp(1)
300 write(0,2) trim(string)
  2 format(1x,'gmsh2tecton: mesh read error in data section: ',a256)
    call exitp(1)
    end

!-------------------------------------------------------------------------------
subroutine TriangleTetrahedronRelation()

use meshdatamodule, only: meshdatactx
use algebra,        only: threeInFour, sort4vec
use modeldefinition, only: NUMSLPglobal
use modeldatamodule, only: modeldatactx

implicit none

integer :: iElem, iStep, iIndex, iTriangle, iPos, iMatch, iNode, iSlip, iSide
integer :: dummyInt, dummyInt2


    integer, allocatable :: newIndex(:)
    integer, allocatable :: sortedIen(:,:)
    integer :: highestLowestIndex, lowestHighestIndex
    integer, allocatable :: elemsPerLowestIndex(:)
    integer, allocatable :: elemsPerHighestIndex(:)
    integer, allocatable :: elementLookup(:)

    integer :: lowestTriangleNode, highestTriangleNode
    logical :: gotya
    integer :: offset, startoffset, endoffset

    integer :: nElemsFound
    integer :: elemsFound(4)

    integer, allocatable :: sortedLowElemsIndex(:)
    integer, allocatable :: sortedHighElemsIndex(:)

integer :: elemToCheck
integer :: mySign
logical, parameter :: debugThisSub = .false.


        !********* Now we need to fill meshdatactx%triangleElements
        ! which is a bit of a pain, as we have no correlation
        ! between the triangles and the tetrahedrons, except their
        ! points. This information is required for traction and such
        ! on triangles when selecting by tag.

        ! Sort the elements based on the lowest node ID
        ! After that, search for the lowest node ID of the triangle
        ! into a smll subsection of elements.
        allocate(sortedLowElemsIndex(meshdatactx%neglobal))
        allocate(sortedHighElemsIndex(meshdatactx%neglobal))

        allocate(sortedIen(4,meshdatactx%neglobal))

        sortedLowElemsIndex = 0
        sortedHighElemsIndex = 0

        sortedIen = meshdatactx%ien

        ! of every element, set the lowest index to the front
        ! the orientation is no ong relevant
        do iElem = 1,meshdatactx%neglobal
            call sort4vec(sortedIen(:,iElem)) ! routine from algebra module
        enddo

        ! now sort the actual elements bij lowest vertex
        ! and by highest vertex.



        ! First, count the number of tetrahedrons who have each vertex as lowest vertex
        highestLowestIndex = maxval(sortedIen(1,:))
        lowestHighestIndex = minval(sortedIen(4,:))

if (debugThisSub) then
        write(*,*) "highestLowestIndex", highestLowestIndex
        write(*,*) "lowestHighestIndex", lowestHighestIndex
endif

        allocate(elemsPerLowestIndex(highestLowestIndex))
        allocate(elemsPerHighestIndex(meshdatactx%nvglobal - lowestHighestIndex+1))

if (debugThisSub) then
        write(*,*) "size of lowest index array", highestLowestIndex
        write(*,*) "size of highest index array", meshdatactx%nvglobal - lowestHighestIndex+1
endif

        elemsPerLowestIndex = 0
        elemsPerHighestIndex = 0

        
        do iElem = 1,meshdatactx%neglobal
            elemsPerLowestIndex(sortedIen(1,iElem)) = &
            elemsPerLowestIndex(sortedIen(1,iElem)) + 1
            ! maybe   elemsPerLowestIndex(sortedIen(1,:)) without the loop? test this.

            elemsPerHighestIndex(sortedIen(4,iElem)-lowestHighestIndex+1) = &
            elemsPerHighestIndex(sortedIen(4,iElem)-lowestHighestIndex+1) + 1
        enddo

if (debugThisSub) then
        write(*,*) "elemsPerLowestIndex", elemsPerLowestIndex(1:100)
        write(*,*) "elemsPerHighestIndex", elemsPerHighestIndex(1:100)
endif

        ! compute the cumulative value to determine starting positions.
        ! Here there is not yet any sorting by tag, so all the entries should
        ! be in the array, although the highest numbers have zeroes,
        ! as those elements are counted by lower indices.


        ! turn lowest into cumulative
        dummyInt = elemsPerLowestIndex(1)
        elemsPerLowestIndex(1) = 0
        do iIndex = 2,highestLowestIndex
            dummyInt2 = elemsPerLowestIndex(iIndex)
            elemsPerLowestIndex(iIndex) = dummyInt
            dummyInt = dummyInt + dummyInt2
        enddo


        ! turn highest into cumulative
        dummyInt = elemsPerHighestIndex(1)
        elemsPerHighestIndex(1) = 0
        do iIndex = 2,meshdatactx%nvglobal - lowestHighestIndex+1
            dummyInt2 = elemsPerHIghestIndex(iIndex)
            elemsPerHighestIndex(iIndex) = dummyInt
            dummyInt = dummyInt + dummyInt2
        enddo

        
! tested and OK
if (debugThisSub) then
       write(*,*) "elemsPerLowestIndex sorted", elemsPerLowestIndex(1:100)
       write(*,*) "elemsPerHighestIndex sorted", elemsPerHighestIndex(1:100)
endif

        ! note that the tetrahedron node indices in the elemsPerLowestIndex
        ! correspond to the positional index in this array. They both start at 0.
        ! However, this is NOT the case for the elemsPerHighestIndex, where the 
        ! array index of 1 corresponds the lowest highest index, which does not
        ! correspond to 1. The difference between the two is hidden in the
        ! the variable lowestHighestIndex.
        do iElem = 1,meshdatactx%neglobal
            ! this elements has lowest index: sortedIen(1,iElem)
            ! and as such, its place in the sorted list is determined
            ! by the offset we just computed in elemsPerLowestIndex

            


            sortedLowElemsIndex(elemsPerLowestIndex(sortedIen(1,iElem)) + 1) = iElem
            sortedHighElemsIndex(elemsPerHighestIndex(sortedIen(4,iElem)-lowestHighestIndex+1)+ 1) = iElem

            ! update the arrays that now actually hold the offsets,
            ! so that th next element with the same highest or lowest point ID
            ! will be entered one firther in the array
            elemsPerLowestIndex(sortedIen(1,iElem)) = &
            elemsPerLowestIndex(sortedIen(1,iElem)) + 1

            elemsPerHighestIndex(sortedIen(4,iElem)-lowestHighestIndex+1) = &
            elemsPerHighestIndex(sortedIen(4,iElem)-lowestHighestIndex+1) + 1
        enddo

if (debugThisSub) then
       write(*,*) "sortedLowElemsIndex" , sortedLowElemsIndex(1:100)
       write(*,*) "sortedHighElemsIndex", sortedHighElemsIndex(1:100)
endif

        ! there are cases when we have an element where the triangle belongs
        ! to the highest three nodes. They are not found when selecting 
        ! by the lowest node, our current search strategy.

        ! There are also triangles that will be found by both the upper and
        ! the lower index. We must make sure to count these only once.

        ! A maximum of four triangles will be found, if both are double found.
        




        ! now the first elements indices in sortedElemsIndex are the 
        ! elements that are adjacent to node 1, then follow the elements that are not
        ! adjacent to 1, but to 2, etc. 

        ! to link the triangles to the tetrahedrons, we have:

        !  meshdatactx%ien                   : four points in tetrahedron we need to match
        !  meshdatactx%triangleConnectivity  : three points we need to match
        !  elemsPerLowestIndex               : offset of where we need to search
        !  sortedElemsIndex                  : the global elements indices at those offsets

        ! and we fill: meshdatactx%triangleElements, for every triangle it will contain  
        ! the one or two tetrahedrons to which it is adjacent

!       write(*,*) "first 10 offsets" ,elemsPerLowestIndex(1:10)

        do iTriangle = 1,meshdatactx%nTriangles
!        do iTriangle = 1,10

if (debugThisSub) then
            write(*,*) "*************************************************************"
            write(*,*) "triangle", &
                       iTriangle, &
                       meshdatactx%triangleConnectivity(:,iTriangle)
endif

            nElemsFound = 0
            elemsFound = 0

            ! check whether there are elements whose lowest index
            ! matches the lowest node of the triangle
            lowestTriangleNode  = minval(meshdatactx%triangleConnectivity(:, iTriangle))

            if (lowestTriangleNode.le.HighestLowestIndex) then

                if (lowestTriangleNode.eq.1) then
                    startoffset = 0
                else
                    startoffset = elemsPerLowestIndex(lowestTriangleNode-1)
                endif

                endoffset = elemsPerLowestIndex(lowestTriangleNode)

if (debugThisSub) then
             write(*,*) "search low from pos", startoffset, &
                           "to", endoffset, &
                           "all nodes: ", meshdatactx%triangleConnectivity(:, iTriangle)
endif

                ! walk through the array to check for matches on the lower end
                do iPos = startoffset+1, endoffset
                    elemToCheck = sortedLowElemsIndex(iPos)

if (debugThisSub) then
                    write(*,*) "Check elem", iPos, "which is", elemToCheck, &
                               "with connectivity", meshdatactx%ien(:,elemToCheck)
endif

                    if (threeInFour(meshdatactx%triangleConnectivity(:, iTriangle), &
                                    meshdatactx%ien(:,elemToCheck))) then
if (debugThisSub) then
                        write(*,*) "GOT ONE!"
endif
                        nElemsFound = nElemsFound + 1
                        elemsFound(nElemsFound) = elemToCheck
                    endif
                enddo

            else
if (debugThisSub) then
                write(*,*) "triangle ", iTriangle, &
                           "cannot overlap with a low point, because lowest point is higher  than lowest tetrahedron point"
                write(*,*) "lowestTriangleNode", lowestTriangleNode
                write(*,*) "HighestLowestIndex", HighestLowestIndex

endif
                ! the lowestTriangle node is higher than the HighestLowestIndex of the 
                ! tetrahedrons
                ! No match will be found here.
            endif


            ! check whether there are elements whose highest index
            ! matches the highest node of the triangle
            highestTriangleNode = maxval(meshdatactx%triangleConnectivity(:, iTriangle))

            ! the highestTriangle node can be lower than the lowestHighestIndex of all the
            ! elements in the domain!
            ! In that case, there is no chance the the element will match with the highest index.


            if (highestTriangleNode.ge.lowestHighestIndex) then

                if (highestTriangleNode.eq.meshdatactx%nvglobal - lowestHighestIndex+1) then
                    startoffset = 1
                elseif (highestTriangleNode.eq.lowestHighestIndex) then
                    startoffset = 1
                else
!                    write(*,*) "two vals", highestTriangleNode, lowestHighestIndex
                    startoffset = elemsPerHighestIndex(highestTriangleNode - lowestHighestIndex)
                endif

                endoffset = elemsPerHighestIndex(highestTriangleNode - lowestHighestIndex+1)

if (debugThisSub) then
                write(*,*) "search high from pos", startoffset, &
                           "to", endoffset, &
                           "all nodes: ", meshdatactx%triangleConnectivity(:, iTriangle)
endif

                do iPos = startoffset, endoffset
                    elemToCheck = sortedHighElemsIndex(iPos)

if (debugThisSub) then
                    write(*,*) "Check elem", iPos, "which is", elemToCheck, &
                               "with connectivity", meshdatactx%ien(:,elemToCheck)
endif

                    if (threeInFour(meshdatactx%triangleConnectivity(:, iTriangle), &
                                    meshdatactx%ien(:,elemToCheck))) then
if (debugThisSub) then
                        write(*,*) "GOT HIGH ONE!"
endif
                        nElemsFound = nElemsFound + 1
                        elemsFound(nElemsFound) = elemToCheck
                    endif
                enddo
            else
if (debugThisSub) then
                write(*,*) "triangle ", iTriangle, &
                           "cannot overlap with a high point, because highest point is lower than lowest high tetrahedron point "
                write(*,*) "highestTriangleNode", highestTriangleNode
                write(*,*) "lowestHighestIndex", lowestHighestIndex
endif
                ! the highestTriangle is lower than the lowestHighestIndex of the tetrahedrons
                ! No match will be found here.
            endif

            !make sure whether we did not mess up:
            if (nElemsFound.eq.0) then
                write(*,*) "Found no element adjacent to triangle", iTriangle
                write(*,*) "This is not supposed to happen"
                write(*,*) "Please contact model support"
                stop "Leaving Plt3d"
            endif


            ! Now that we have the elements matching this triangle,
            ! get the unique ones out, and stuff them where they 
            ! belong in meshdatactx%triangleElements

            do iMatch=1,nElemsFound
                if      (meshdatactx%triangleElements(1,iTriangle) .eq.0) then
                         meshdatactx%triangleElements(1,iTriangle) = elemsFound(iMatch)
                else if (meshdatactx%triangleElements(2,iTriangle) .eq.0) then
                    if  (meshdatactx%triangleElements(1,iTriangle) .ne. elemsFound(iMatch)) then
                         meshdatactx%triangleElements(2,iTriangle) = elemsFound(iMatch)
if (debugThisSub) then
                        write(*,*) "***********************************************"
                        write(*,*) "triangle", iTriangle, &
                                   "has two elems adjacent"
                        do iNode=1,3
                        write(*,*) "triangle node", iNode, &
                                   "is",  meshdatactx%triangleConnectivity(iNode,iTriangle), &
                                   "with coords", meshdatactx%X(:,meshdatactx%triangleConnectivity(iNode,iTriangle))
                        enddo

                        do iElem=1,2
                        write(*,*) "element ", iElem, &
                                   "is", meshdatactx%triangleElements(iElem,iTriangle), &
                                   "with nodes",  meshdatactx%ien(:, meshdatactx%triangleElements(iElem,iTriangle))

                        enddo
endif                    
                    endif
                else
                    ! we have found two different elements, all the rest must be doubles
                endif
            enddo


        enddo

        deallocate(elemsPerLowestIndex)
        deallocate(elemsPerHighestIndex)

        deallocate(sortedLowElemsIndex)
        deallocate(sortedHighElemsIndex)
        deallocate(sortedIen)

                ! now that memory has been cleaned, we can determine the slippery sides
                ! when this is the case.
                ! We will solve this by adding a minus sign to the
                ! meshdatactx%triangleElements
                ! for the elements that have a minus sign to the slippery weights.

                if (NUMSLPglobal.gt.0) then
                        allocate(elementLookup(meshdatactx%neglobal))
                        elementLookup = 0

                        ! walk through the slippery entries to set the lookup table
                        do iSlip = 1,NUMSLPglobal
                                ! 1 is elemID, 3 = first weight
                                elementLookup(modeldatactx%NSLIPglobal(1,iSlip)) = &
                                                modeldatactx%NSLIPglobal(3,iSlip)
                        enddo

                        ! use this info to adapt meshdatactx%triangleElements
                        do iTriangle=1,meshdatactx%nTriangles
                                do iSide=1,2
                                        if (meshdatactx%triangleElements(iSide,iTriangle).ne.0) then
                                                
                                                mySign = elementLookup(meshdatactx%triangleElements(iSide,iTriangle))

                                                if (mySign.eq.-1) then
                                                        meshdatactx%triangleElements(iSide,iTriangle) = &
                                   -meshdatactx%triangleElements(iSide,iTriangle)
                                                endif
                                        else
                                                ! can be 0 when triangle is on edge of domain, where it only has 
                                                ! a single element adjacent to it.
                                        endif
                                enddo
                        enddo
                        deallocate(elementLookup)
                endif

end subroutine



!--------------------------------------------------------------------------------------------------!

subroutine lineSegmentsInTetrahedra()
!----
! This subroutine defines line segments as pairs of nodes (lineSegments) and stores the four line
! segments that make up each tetrahedral element (linesInElements)
!----

use meshdatamodule, only: meshdatactx
use algebra, only: sort4vec
implicit none
integer :: nLineSegments
integer :: i, j, k, n, iLineSeg
integer :: np(4)
integer :: alreadyThere
integer, allocatable :: sortedIen(:,:)
integer, allocatable :: nNodeNeighbors(:)
integer, allocatable :: NodeNeighbors(:,:)
integer, allocatable :: cumulativeLineSegments(:)
integer :: maxNeighbors

! Allocate memory for some local arrays and initialize them
maxNeighbors = 25
allocate(sortedIen(4,meshdatactx%neglobal))
allocate(nNodeNeighbors(meshdatactx%nvglobal))
allocate(nodeNeighbors(maxNeighbors,meshdatactx%nvglobal))
nNodeNeighbors = 0
nodeNeighbors = 0


! Sort the nodes within each element in ascending order
sortedIen = meshdatactx%ien
do i = 1,meshdatactx%neglobal
    call sort4vec(sortedIen(:,i)) ! routine from algebra module
enddo


! Count the total number of line segments making up the mesh
meshdatactx%nLineSegments = 0

! We do not a priori know how many line segments (node pairs) there are
! First, go through the nodes in all tetrahedra and build node neighbor array
! Keep track of the total number to allocate memory for the lineSegment array
do i = 1,meshdatactx%neglobal

    do j = 1,4
        np(j) = sortedIen(j,i)
    enddo

    ! For each line segment in the tetrahedron, check if the node's neighbors should be added 
    ! to the neighbor array
    do j = 1,3
        do k = j+1,4

            ! Check if the neighbor is already in the array
            alreadyThere = 0
            do n = 1,nNodeNeighbors(np(j))
                if (nodeNeighbors(n,np(j)).eq.np(k)) then
                    alreadyThere = 1
                    exit
                endif
            enddo

            ! Did not find the neighbor... add it to the list
            if (alreadyThere.eq.0) then
                nNodeNeighbors(np(j)) = nNodeNeighbors(np(j)) + 1
                if (nNodeNeighbors(np(j)).gt.maxNeighbors) then
                    write(0,*) '!! Error: number of nodes neighboring node ',np(j), &
                               ' exceeds maxNeighbors=',maxNeighbors
                    stop
                endif
                nodeNeighbors(nNodeNeighbors(np(j)),np(j)) = np(k)
                
                ! Update line segment count
                meshdatactx%nLineSegments = meshdatactx%nLineSegments + 1
            endif
        enddo
    enddo

enddo
!do i = 1,meshdatactx%nvglobal
!    write(0,*) i,nNodeNeighbors(i),(nodeNeighbors(j,i),j=1,nNodeNeighbors(i))
!enddo


! We will need to remember where nodes are in the neighbor array for linking them to elements,
! so compute the cumulative number of line segments at each node
allocate(cumulativeLineSegments(meshdatactx%nvglobal))
cumulativeLineSegments = 0
do i = 2,meshdatactx%nvglobal
    cumulativeLineSegments(i) = nNodeNeighbors(i-1) + cumulativeLineSegments(i-1)
enddo
!do i = 1,meshdatactx%nvglobal
!    write(0,*) i,cumulativeLineSegments(i)
!enddo


! Load the node pairs into the line segment array in order of their appearance in nodeNeighbors
allocate(meshdatactx%lineSegments(2,meshdatactx%nLineSegments))
n = 0
do i = 1,meshdatactx%nvglobal
    if (nNodeNeighbors(i).gt.0) then
        do j = 1,nNodeNeighbors(i)
            n = n + 1
            meshdatactx%lineSegments(1,n) = i
            meshdatactx%lineSegments(2,n) = nodeNeighbors(j,i)
        enddo
    endif
enddo
!do i = 1,meshdatactx%nLineSegments
!    write(0,*) i,meshdatactx%lineSegments(1,i),meshdatactx%lineSegments(2,i)
!enddo


! Go back through the elements and their nodes one more time
! This time, keep track of which line segments are in each element
allocate(meshdatactx%linesInElements(6,meshdatactx%neglobal))
meshdatactx%linesInElements = 0

do i = 1,meshdatactx%neglobal

    do j = 1,4
        np(j) = sortedIen(j,i)
    enddo

    ! For each node pair (line segment) in the tetrahedron, figure out its line segment number
    iLineSeg = 0
    do j = 1,3
        do k = j+1,4
            do n = 1,nNodeNeighbors(np(j))
                if (nodeNeighbors(n,np(j)).eq.np(k)) then
                    iLineSeg = iLineSeg + 1
                    meshdatactx%linesInElements(iLineSeg,i) = cumulativeLineSegments(np(j)) + n
                    exit
                endif
            enddo
        enddo
    enddo
enddo
!do i = 1,meshdatactx%neglobal
!    write(0,*) i,(meshdatactx%linesInElements(j,i),j=1,6)
!enddo

! Finished with local arrays
deallocate(sortedIen)
deallocate(nNodeNeighbors)
deallocate(nodeNeighbors)
deallocate(cumulativeLineSegments)

return
end subroutine
