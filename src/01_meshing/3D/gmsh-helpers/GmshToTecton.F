module meshdata

implicit none

real,    allocatable    :: X(:,:), XN(:)

integer, allocatable :: NODE(:)
integer, allocatable :: IEL(:),ITYPE(:)
integer, allocatable :: IPHYS(:),IEN(:,:)
integer, allocatable :: NE(:), NP(:)
integer, allocatable :: vertexTag(:)

integer, allocatable :: iSpec(:)

integer, allocatable :: elementTag(:)
integer, allocatable :: isVertexSpec(:)

end module

module io

implicit none

integer :: stdin,stdout,stderr,lui,luo

end module



program main
!
! reads GMSH mesh file and outputs GTECTON data
! MODE = 0    generate node coordinate list (command line argument "-c")
! MODE = 1    generate element numbers list ("-e")
! MODE = 2    select node numbers with label iphys ("-n iphys")
! MODE = 3    select side numbers with label iphys ("-s iphys")
! MODE = 4    select fault element/node numbers ("-f iphys")
! MODE = 5    print major info on mesh file ("-I")
! MODE = 6    output VTK file with line segments

! NODE(MAXNP)    GMSH node number
! X(3,MAXNP)    node coordinate
! IEL(MAXEL)    GMSH element number
! ITYPE(MAXEL)    element type:   (by gmsh convention)
!         1 Line (2 nodes, 1 edge).
!         2 Triangle (3 nodes, 3 edges).
!         3 Quadrangle (4 nodes, 4 edges). 
!         4 Tetrahedron (4 nodes, 6 edges, 4 facets). 
!        15 Point (1 node)
!         .
!         .
! IPHYS(MAXEL)    physical number
! IEN(8,MAXEL)    GMSH node numbers list
!

use meshdata
use io,       only: stdin,stdout,stderr,lui,luo

implicit none

integer, parameter :: logfileID = 42
character(len=15), parameter :: logFileName="Gmsh2tecton.log"

!-locl
character(len=256) fin,fout,frmt
logical  :: verbose
logical :: MattCoords
logical :: TacoMaterials

!    integer iflu,mode,nextlu,NUMNP,NNEL,NODE,IEL,ITYPE,IPHYS,IREG,
!     >   IEN,ispec,i,NP,j,NUMEL,NUMAT,k,nsh,nfc,NE,n,ierr

integer(kind=8) :: NUMNP,NNEL,numel
integer :: iflu,mode,nextlu,IREG, &
    i,j,NUMAT,k,nsh,nfc,n,ierr

integer :: nspec, jspec
logical :: foundAMatch


!    real X,XN,AREA,det
real AREA,det


!    dimension :: NODE(MAXNP),X(3,MAXNP),
!     >   IEL(MAXEL),ITYPE(MAXEL),
!     >   IPHYS(MAXEL),IEN(8,MAXEL),NP(MXINP),
!     >   XN(3),NE(MXINE)

integer :: nrelevant

external iflu,nextlu

integer :: triVertex,tetVertex

allocate(iSpec(1000)) ! number of possible different tags

!-init
nrelevant = 0
stdin  = iflu('stdin')
stdout = iflu('stdout')
stderr = iflu('stderr')
lui = stdin
luo = stdout
!    
!     set defaults
mode = 0
ispec = 1
nspec = 1
fin  = '<stdin>'
fout = '<stdout>'
verbose = .false.

!     process command line arguments
call shell (mode, ispec, nspec, fin, fout, verbose, MattCoords, TacoMaterials)

if (verbose) then
    open(unit=logfileID, file=logFileName)
    rewind(unit=logfileID)
    write(stderr,*) 'Detailed information in file "Gmsh2tecton.log"'
    write(stderr,*) "Allocating arrays..."
endif

! Nobody I know ever reads from standard input in gmsh2tecton...
if (fin.ne.'<stdin>') then
!       open GMSH mesh file
    lui = nextlu(stdout)
    call openf(lui,fin,'old')
endif

if (fout.ne.'<stdout>') then
!       open output file
    luo = nextlu(stdout)
    call openf(luo,fout,'unknown')
endif

if (verbose) then
    write(stderr,*) "Reading mesh..."
endif

!   read GMSH mesh data
call readMeshFile (NUMNP,NNEL,NUMAT,verbose,logfileID,MODE)

if (verbose) then
    write(stderr,*) "finished reading mesh of size:"
    write(stderr,*) "numnp: ", numnp
    write(stderr,*) "nnel:  ", nnel
    write(stderr,*) "numat: ", numat
endif

if (verbose) then
    write(stderr,*) "found mode: ", MODE
endif



if (MODE.eq.5) then
    call exitp(0)
endif

!   generate lists of new node and element numbers (NP and NE)
if (verbose) then
    write(stderr,*) "Generating inverted node list...", NUMNP,NNEL
endif
call INVLST (NUMNP,NNEL)
!
!   Count number of elements
if (verbose) write(stderr,*) "Counting elements ..."
NUMEL = 0
do j=1,NNEL
    if (MODE.eq.6) then
        ! count number of line elements
        if (itype(j).eq.1) NUMEL = NUMEL + 1
    else
        ! count number of tetrahedral elements
        if (itype(j).eq.4) NUMEL = NUMEL + 1
    endif
enddo

if (verbose) write(stderr,*) "Checking material properties of elements..."

if (MODE.ne.6 .and. NUMAT.le.99999) then
    do 10 i=1,NUMAT
        do j=1,NNEL
            if (itype(j).eq.4 .and. IPHYS(j).eq.i) then
                goto 10
            endif
        enddo
        write(stderr,9) i
9        format(1x,'gmsh2tecton: no tetrahedral elements with ', &
       'material number ',I5,' in output')
10      continue
endif



if (MODE.eq.0) then

    !******** generate tecin.dat.nps *************

    if (verbose) then
        write(stderr,*) "mode 0; generating tecin.dat.nps"
    endif
    frmt = '(2I12,3E25.17))' ! very large format to take up all the precision created by gmsh
    do i=1,NUMNP
        write(luo,frmt,err=1000) i,vertexTag(i),X(1,i),X(2,i),X(3,i)
    enddo
    write(luo,20,err=1000)
   20   format('end nodal coordinate data')


else if (MODE.eq.1) then

    !******** generate tecin.dat.elm *************
    ! ****** write tetrahedral elements while simultaneously redefining node (?element?) numbers

    if (verbose) then
        write(stderr,*) "mode 1; generating tecin.dat.elm"
    endif
    frmt = '(6I12)'

    if (NUMAT.gt.99999) then
        write(stderr,30)
    endif
   30   format(14x,'ignoring physical geometry information')

    n = 0 ! counter for number of tetrahedra
    do j=1,NNEL
        if (itype(j).eq.4) then
            if (NUMAT.gt.99999) then
                IPHYS(j) = 1
            endif
            n = n + 1

            call elmtest (j,X,IEN(1:8,j),NP,NUMNP,det,ierr)
            if (ierr.ne.0) then
                write(stderr,40) n,IPHYS(j),det
            endif
   40           format(1x,'element = ',I12,', material number = ',I12, &
             ', volume = ',1PG14.6)

            if (IPHYS(j)>0) then
                write(luo,frmt,err=1000) n,IPHYS(j),&
                (NP(IEN(i,j)),i=1,4)
            else
                if (verbose) then
                    write(*,*) "No physical volumes detected in the mesh file."
                    write(*,*) "Setting material code 1 for every element."
                endif
                write(luo,frmt,err=1000) n,1,(NP(IEN(i,j)),i=1,4)
            endif
        endif
    enddo
    write(luo,60,err=1000)
   60   format('end element data')

    if (NUMAT.gt.99999) then
        write(stderr,70) 1
    else
        write(stderr,70) NUMAT
    endif
   70   format(1x,'NUMAT = ',I5)

else if (MODE.eq.2) then

    ! *****  select nodes with a specific tag

    if (verbose) then
        write(stderr,*) "mode 2; nodal boundary conditions"
    endif
    do i=1,NUMNP
        do j = 1,nspec
            if (VertexTag(i).eq.ISPEC(j)) then
                if (MattCoords) then
                    write(luo,'(I12,3F20.6)',err=1000) &
                    NP(i),X(1,NP(i)),X(2,NP(i)),X(3,NP(i))
                else
                    write(luo,'(I12)',err=1000) NP(i)
                endif
                exit ! break j loop over tags
            endif
        enddo
    enddo

else if (MODE.eq.3) then

    ! ****** MODE 3; select by side tag number

    if (verbose) then
        write(stderr,*) "mode 3; element boundary conditions"
    endif

    ! Element side boundary conditions: output element/surface numbers with a specific iphys value
    do i=1,NNEL
        do jspec = 1,nspec
            if (ITYPE(i).eq.2 .and. IPHYS(i).eq.ISPEC(jspec)) then
                ! Search for tetrahedra that have this face and establish the face number NFC
                do j=1,NNEL
                    nsh = 0
                    nfc = 0
                    if (ITYPE(j).eq.4) then

                        if (IEN(1,i).eq.IEN(1,j) .or. &
                            IEN(2,i).eq.IEN(1,j) .or. &
                            IEN(3,i).eq.IEN(1,j)) then
                            nsh = nsh + 1 ! # shared nodes
                            nfc = nfc + 1
                        endif

                        if (IEN(1,i).eq.IEN(2,j) .or.&
                            IEN(2,i).eq.IEN(2,j).or.&
                            IEN(3,i).eq.IEN(2,j)) then
                            nsh = nsh + 1
                            nfc = nfc + 2
                        endif

                        if (IEN(1,i).eq.IEN(3,j) .or.&
                            IEN(2,i).eq.IEN(3,j).or.&
                            IEN(3,i).eq.IEN(3,j)) then
                            nsh = nsh + 1
                            nfc = nfc + 3
                        endif

                        if (IEN(1,i).eq.IEN(4,j) .or.&
                            IEN(2,i).eq.IEN(4,j).or.&
                            IEN(3,i).eq.IEN(4,j)) then
                            nsh = nsh + 1
                            nfc = nfc + 4
                        endif

                        if (nsh.eq.3) then    
                            if (nfc.eq.6) then
                                nfc = 2
                            else if (nfc.eq.7) then
                                nfc = 1
                            else if (nfc.eq.8) then
                                nfc = 3
                            else if (nfc.eq.9) then
                                nfc = 4
                            else
                                write(stderr,90) (IEN(k,i),k=1,3),(IEN(k,j),k=1,4)
   90                               format(1x,'gmsh2tecton: can''t',1x,&
                                    'establish face number:'/14x,&
                                    'face node numbers: ',3i12/14x,&
                                    'element node numbers: ',4i12)
                                call exitp(1)
                            endif
                            call FCGEOM (X,NP,IEN(1,j),nfc,XN,AREA)
            
                            write(luo,100) NE(IEL(j)),nfc,&
                                            XN(1),XN(2),XN(3),AREA
  100                           format(2I12,4(1X,1PG14.6))
                        endif
                    endif
                enddo ! end loop over j (second loop over tetrahedra)
                exit ! exit jspec loop since we found the matching element
            endif
        enddo ! end loop over jspec (list of tags)
    enddo ! end loop over i (first loop over tetrahedra)


else if (MODE.eq.4) then

    !********* MODE 4; select fault element/nodes combinations

    if (verbose) then
        write(stderr,*) "mode 4, faulted node boundary conditions"
    endif

    allocate(isVertexSpec(NUMNP), stat=iErr)
    if (iErr .ne. 0) then
        write(*,*) "could not allocate isVertexSpec array. Error: ", iErr
        stop "Leaving Gmsh2tecton"
    endif


    ! Initialize each node as not specified
    do i = 1,NUMNP
        isVertexSpec(i) = 0
    enddo

    ! Loop through all elements to load tagged node lookup table
    do i = 1,NNEL
        ! pick only triangle elements
        if (ITYPE(i).ne.2) then
            cycle ! not a triangle, cycle loop through elements
        endif

        ! pick only elements that match a tag specified in ISPEC
        foundAMatch = .false.
        do j = 1,nspec
            if (IPHYS(i).eq.ISPEC(j)) then
                foundAMatch = .true.
                exit ! found an element with the correct label, exit loop over tags
            endif
        enddo
        if (.not.foundAMatch) then
            cycle ! no triangular element with that tag, cycle loop through elements
        endif

        ! fill the nodes of this triangle into the lookup table
        do j = 1,3
            triVertex = IEN(j,i)        ! original number of element node
            isVertexSpec(triVertex) = i ! designate this node as specified with triangular element number
        enddo
    enddo

    ! Use lookup array (isVertexSpec) to determine if specified node is in tetrahedral element
    do j = 1,NNEL

        ! pick only tetrahedral elements
        if (ITYPE(j).ne.4) then
            cycle ! not a tetrahedron, cycle loop through elements
        endif

        ! loop through nodes in tetrahedral element
        do k = 1,4
            tetVertex = IEN(k,j)        ! original number of element node
            if (isVertexSpec(tetVertex).ge.1) then
                call vector(X,NP,IEN(1,j),tetVertex,IEN(1,isVertexSpec(tetVertex)),XN, ierr)
                if (ierr.ne.0) then
                    write(stderr,101) NE(IEL(j)),NP(tetVertex)
  101                   format(1x,'element = ',I12,', node = ',I12)
                    call exitp(1)
                endif
                if (MattCoords) then
                    if (TacoMaterials) then
                        ! new element number, new node number, tag, normal vector, and vector to first node
                        write(luo,110) NE(IEL(j)),NP(tetVertex),IPHYS(j),XN(1),XN(2),XN(3),&
                         X(1,tetVertex),X(2,tetVertex),X(3,tetVertex)
  110                       format(3I12,3(1PG14.6),3(1PG14.6))
                    else
                        ! new element number, new node number, normal vector, and vector to first node
                        write(luo,111) NE(IEL(j)),NP(tetVertex),XN(1),XN(2),XN(3),&
                         X(1,tetVertex),X(2,tetVertex),X(3,tetVertex)
  111                       format(2I12,3(1PG14.6),3(1PG14.6))
                    endif
                else
                    if (TacoMaterials) then
                        ! new element number, new node number, tag, normal vector
                        write(luo,112) NE(IEL(j)),NP(tetVertex),IPHYS(j),XN(1),XN(2),XN(3)
  112                       format(3I12,3(1PG14.6))
                    else
                        ! new element number, new node number, normal vector
                        write(luo,113) NE(IEL(j)),NP(tetVertex),XN(1),XN(2),XN(3)
  113                       format(2I12,3(1PG14.6))
                    endif
                endif
            endif
        enddo
    enddo
    deallocate(isVertexSpec)


else if (MODE.eq.6) then

    !****** write line VTK file

    write(luo,210) 2*NUMEL
210     format('# vtk DataFile Version 3.0'/'Generated by gmsh2tecton'/&
     'ASCII'/'DATASET POLYDATA'/'POINTS ',I6,' float')
     do j = 1,NNEL
        if (ITYPE(j).eq.1) then
            write(luo,220) X(1,IEN(1,j)),X(2,IEN(1,j)),X(3,IEN(1,j)),&
             X(1,IEN(2,j)),X(2,IEN(2,j)),X(3,IEN(2,j))
220             format(6(1PG16.6,1X))
        endif
     enddo
     write(luo,230) NUMEL,NUMEL*3
230      format('LINES ',I6,1X,I6)
     i = 0
     do j = 1,NNEL
        if (ITYPE(j).eq.1) then
            write(luo,240) i,i+1
240             format('2 ',I6,1X,I6)
            i = i + 2
        endif
     enddo

else
    write(stderr,*) 'gmsh2tecton: no MODE number ',MODE
    call exitp(1)
endif

if (verbose) close(logfileID)

call exitp(0)
!
1000    write(stderr,1001)
 1001    format(1x,'gmsh2tecton: cannot write output')
call exitp(1)

end program
!-------------------------------------------------------------------------------
subroutine readMeshFile (NUMNP,NNEL,NUMAT,verbose,logfileID,MODE)

use io, only: lui, stderr

!
! read GMSH mesh data
!
implicit none
!-pass
integer(kind=8) :: NUMNP,NNEL
integer :: numat, mode
logical :: verbose
integer :: logfileID
!-locl
character(len=14) string
integer(kind=8) :: numel
integer nsz,ios,i,n,j,filetype,datasize,itag,ntag,lien
real meshformat
integer :: iPoint
dimension nsz(19),lien(27),itag(10)
    ! number of vertices per element.
    ! element 1 is a line -> 2 vertices
    ! element 2 is a triangle -> 3 vertices
    ! etc.
    data nsz/2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13/
!-init
numel = 0
!
! read '$MeshFormat'
read(lui,'(a14)',err=1000,end=1010,iostat=ios) string

! tead content of mesh format
if (string(1:11).eq.'$MeshFormat') then
    read(lui,*,err=1020,end=1030,iostat=ios) meshformat,&
    filetype,datasize
    if (filetype.ne.0) then
        write(stderr,1)
1           format(1x,'gmsh2tecton: non-ASCII GMSH mesh file')
        call exitp(1)
    endif

    if (datasize.ne.8) then
        write(stderr,2)
2       format(1x,'gmsh2tecton: only 8-sized supported')
        call exitp(1)
    endif

! read '$EndMeshFormat'

    read(lui,'(a14)',err=1040,end=1050,iostat=ios) string
    if (string(1:14).ne.'$EndMeshFormat') then
        write(stderr,3)
3       format(1x,'gmsh2tecton: "$MeshFormat" block not properly closed')
        call exitp(1)
    endif

    if (INT(meshformat).eq.2 .or. INT(meshformat).eq.3) then
        if (MODE.eq.5) then
            write(*,10) 2.2
10      format('gmsh file format:',F3.1)
        endif

        call readMeshData (NUMNP,NNEL,NUMAT,verbose,logfileID,MODE)

    else if (INT(meshformat).eq.4) then
        write(*,12)
12      format(1x,'gmsh2tecton: support for gmsh file format 4',&
     ' is under development'/1x,&
     'For now, you can save in format 2.2 by ',&
     'adding "Mesh.MshFileVersion = 2.2;" to ',&
     'your "~/.gmsh-options" file')
         call exitp(1)
    else
        write(stderr,4) INT(meshformat)
4       format(1x,'gmsh2tecton: unknown mesh file format ',I2)
        call exitp(1)
    endif
else
    write(stderr,5)
5   format(1x,'gmsh2tecton: "$MeshFormat" header missing in mesh file')
    call exitp(1)
endif

return

1000 write(stderr,1001)
1001 format(1x,'gmsh2tecton: empty file')
     call exitp(1)

1010 write(stderr,1011)
1011 format(1x,'gmsh2tecton: read error on 1st line of mesh file')
     call exitp(1)

1020 write(stderr,1021) ios
1021 format(1x,'gmsh2tecton: read error ',I5,' on 2nd line of mesh file')
     call exitp(1)

1030 write(stderr,1031)
1031 format(1x,'gmsh2tecton: premature EOF on 2nd line of mesh file ')
     call exitp(1)

1040 write(stderr,1041) ios
1041 format(1x,'gmsh2tecton: read error ',I5,' on 3rd line of mesh file')
     call exitp(1)

1050 write(stderr,1051)
1051 format(1x,'gmsh2tecton: premature EOF on 3rd line of mesh file')
     call exitp(1)

end subroutine
!-------------------------------------------------------------------------------
subroutine readMeshData (NUMNP,NNEL,NUMAT,verbose,logfileID,MODE)
!
! read GMSH mesh data: data format 2
!
use constants, only: elemIsPoint, &
                     elemIsLine, &
                     elemIsTriangle, &
                     elemIsTetrahedron

use meshData
use io,        only: lui, stderr

implicit none
!-pass
integer(kind=8) :: NUMNP,NNEL
integer :: NUMAT
logical :: verbose
integer :: logfileID,MODE
!-locl
character(len=256) string
integer nsz,ios,i,n,j,numel,filetype,datasize,itag,ntag,lien,&
 numline,numtria, iError
real meshformat
integer :: iPoint
dimension nsz(19),lien(27),itag(10)
    ! number of vertices per element.
    ! element 1 is a line -> 2 vertices
    ! element 2 is a triangle -> 3 vertices
    ! etc.
    data nsz/2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13/

if (verbose) then
    write(0,*) 'reading format version 2 mesh file'
endif

NUMNP = -1
NNEL = -1

100 read(lui,'(a256)',err=1000,end=200,iostat=ios) string

if (string(1:7).eq.'$Nodes ') then
    ! read number of GMSH nodes.
    read(lui,*,err=1020,end=1030,iostat=ios) NUMNP

    if (NUMNP.le.0) then
        write(stderr,10)
  10        format(1x,'gmsh2tecton: no nodal coordinates in mesh file')
        call exitp(1)
    else 

        if (verbose) then
            write(*,*) "Allocating nodal point data for ",numnp,"points"
        endif

        allocate(NODE(NUMNP), stat = iError)
        if (iError .ne. 0) then
            write(*,*) "could not allocate node array. Error: ", iError
            stop "Leaving Gmsh2tecton"
        endif
        allocate(X(3,NUMNP))
        if (iError .ne. 0) then
                  write(*,*) "could not allocate X array. Error: ", iError
                  stop "Leaving Gmsh2tecton"
        endif
        allocate(NP(10*NUMNP))
        if (iError .ne. 0) then
                  write(*,*) "could not allocate NP array. Error: ", iError
                  stop "Leaving Gmsh2tecton"
        endif
        allocate(XN(3))
        if (iError .ne. 0) then
                  write(*,*) "could not allocate XN array. Error: ", iError
                  stop "Leaving Gmsh2tecton"
        endif
        allocate(vertexTag(NUMNP))
        if (iError .ne. 0) then
                  write(*,*) "could not allocate vertexTag array. Error: ", iError
                  stop "Leaving Gmsh2tecton"
        endif

        NODE = 0
        X = 0d0
        NP = 0
        XN = 0d0
        vertexTag = 0

    endif

    do i=1,NUMNP
        read(lui,*,err=1040,end=1050,iostat=ios) &
        NODE(i),X(1,i),X(2,i),X(3,i)
    enddo

    read(lui,'(a256)',err=1060,end=1070,iostat=ios) string
    if (string(1:10).ne.'$EndNodes ') then
        write(stderr,40) 
  40        format(1x,'gmsh2tecton: mesh file "$EndNodes" missing')
        call exitp(1)
    endif

    if (verbose) write(stderr,42) NUMNP
  42    format('gmsh2tecton read ',I7,' GMSH nodes from mesh file')

    if (MODE.eq.5) write(*,44) NUMNP
  44    format('Number of GMSH nodes = ',I7)

elseif (string(1:10).eq.'$Elements ') then
    read(lui,*,err=1100,end=1110,iostat=ios) NNEL   ! number of GMSH elements
    if (NNEL.le.0) then
        write(stderr,50)
  50        format(1x,'gmsh2tecton: no GMSH elements in mesh file')
        call exitp(1)
    else

        if (verbose) then
            write(*,*) "Allocating elements data for ",nnel, "elements"
            write(*,*) "Allocating NE array for ", 10*NNEL, "elements"
        endif

        allocate(IEL(NNEL))
        if (iError .ne. 0) then
                  write(*,*) "could not allocate IEL array. Error: ", iError
                  stop "Leaving Gmsh2tecton"
        endif
        allocate(ITYPE(NNEL))
        if (iError .ne. 0) then
                  write(*,*) "could not allocate ITYPE array. Error: ", iError
                  stop "Leaving Gmsh2tecton"
        endif
        allocate(IPHYS(NNEL))
        if (iError .ne. 0) then
                  write(*,*) "could not allocate IPHYS array. Error: ", iError
                  stop "Leaving Gmsh2tecton"
        endif
        allocate(IEN(8,NNEL))
        if (iError .ne. 0) then
                  write(*,*) "could not allocate IEN array. Error: ", iError
                  stop "Leaving Gmsh2tecton"
        endif
        allocate(elementTag(NNEL))
        if (iError .ne. 0) then
            write(*,*) "could not allocate elementTag array. Error: ", iError
            stop "Leaving Gmsh2tecton"
        endif
        allocate(NE(10*NNEL))
        if (iError .ne. 0) then
            write(*,*) "could not allocate NE array. Error: ", iError
            stop "Leaving Gmsh2tecton"
        endif


        IEL = 0
        ITYPE = 0
        IPHYS = 0
        IEN = 0
        elementTag = 0

    endif
    ! GMSH calls any (1D,2D,3D) connectivity an element. GTECTON is only
    ! interested in the 3D tetrahedra.
    NUMEL = 0
    NUMAT = 0
    numline = 0
    numtria = 0
    do i=1,NNEL
!!!         Potential out-of-array boundaries nsz problem here for new GMSH element types
        read(lui,*,err=1120,end=1130,iostat=ios) &
         IEL(i),ITYPE(i),ntag,(itag(j),j=1,ntag),&
         (lien(j),j=1,nsz(ITYPE(i)))

! the Gmsh constants
! (taken from: http://gmsh.info/doc/texinfo/gmsh.html)
!integer,          parameter :: elemIsPoint          = 15
!integer,          parameter :: elemIsLine           = 1
!integer,          parameter :: elemIsTriangle       = 2
!integer,          parameter :: elemIsTetrahedron    = 4



        if (ITYPE(i).ne.elemIsPoint     .and. &   
            ITYPE(i).ne.elemIsLine      .and. &      
            ITYPE(i).ne.elemIstriangle  .and. & 
            ITYPE(i).ne.elemIsTetrahedron) then 
            write(stderr,70) ITYPE(i)
  70            format(1x,'gmsh2tecton: element type ',I2, &
             ' not supported by GTECTON')
            call exitp(1)
        endif
        IPHYS(i) = 0
        ! detect whether the element has a tag, and if so, possibly assign it to the nodes.
        ! make sure not to overwrite existing tags.
        ! ONLY looks at the first tag when there are multiple
        if (ntag.gt.0) then
            IPHYS(i) = itag(1)
        endif

        do j=1,nsz(ITYPE(i))
            IEN(j,i) = lien(j)
        enddo

        if (ITYPE(i).eq.4) then
            NUMEL = NUMEL + 1
            NUMAT = MAX(NUMAT,IPHYS(i))
        else if (ITYPE(i).eq.2) then
            numtria = numtria + 1
        else if (ITYPE(i).eq.1) then
            numline = numline + 1
        endif
    enddo

!       Set vertex tags first by surfaces, then by lines, then by points
    do i=1,NNEL
        if (ITYPE(i).eq.2 .and. IPHYS(i).ne.0) then
            do iPoint=1,nsz(ITYPE(i))
                vertexTag(IEN(iPoint,i)) = IPHYS(i)
            enddo
        endif
    enddo
    do i=1,NNEL
        if (ITYPE(i).eq.1 .and. IPHYS(i).ne.0) then
            do iPoint=1,nsz(ITYPE(i))
                vertexTag(IEN(iPoint,i)) = IPHYS(i)
            enddo
        endif
    enddo
    do i=1,NNEL
        if (ITYPE(i).eq.15 .and. IPHYS(i).ne.0) then
            do iPoint=1,nsz(ITYPE(i))
                vertexTag(IEN(iPoint,i)) = IPHYS(i)
            enddo
        endif
    enddo

!       read data block footer
    read(lui,'(a12)',err=1140,end=1150,iostat=ios) string
    if (string(1:13).ne.'$EndElements ') then
        write(stderr,80) 
  80        format(1x,'gmsh2tecton: "$Elements" block not properly closed')
        call exitp(1)
    endif

    if (NUMEL.eq.0) then
        write(stderr,90) 
  90        format('gmsh2tecton: no tetrahedral elements in mesh file')
        call exitp(1)
    endif

    if (verbose) write(stderr,92) NUMEL
  92    format('gmsh2tecton read ',I7,' tetrahedral elements from mesh file')
    if (MODE.eq.5) then
        write(*,94) NNEL,numline,numtria,NUMEL,NUMAT
  94        format('Number of GMSH elements = ',I7/&
               'Number of line elements = ',I7/&
               'Number of triangle elements = ',I7/&
               'Number of tetrahedral elements = ',I7/&
               'Number of materials = ',I7)
    endif

!        if (NUMAT.eq.0) write(stderr,95)
!  95    format('gmsh2tecton: NUMAT=0, so you either did not define ',&
!         'Physical Volumes yet, or you wrote your mesh file while ',&
!         'Mesh.SaveAll=1')

elseif (string(1:15).eq.'$PhysicalNames ') then
    read(lui,*,err=1140,end=1150,iostat=ios) j
    if (j.gt.0 .and. MODE.eq.5) write(*,96) j
  96    format('Number of physical names = ',I3)
    string = '$End' // string(2:15)
    if (MODE.eq.5) then
        call skipblock(trim(string),.true.)
    else
        call skipblock(trim(string),.false.)
    endif
elseif (string(1:10).eq.'$Periodic ') then
    read(lui,*,err=1140,end=1150,iostat=ios) j
    if (j.gt.0 .and. MODE.eq.5) write(*,98) j
  98    format('Number of periodic entities = ',I3)
    string = '$End' // string(2:252)
    call skipblock(trim(string),.false.)
elseif (string(1:1).eq.'$') then
    if (verbose) write(stderr,*) 'skipping msh block with ',trim(string)
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
subroutine shell (MODE,ISPEC,nspec,FIN,FOUT, verbose, MattCoords, TacoMaterials)

use io, only: stderr

!
! process command line arguments (see man page)
!
implicit none
!-pass
integer :: MODE, ISPEC(1000), nspec
character*(*) FIN,FOUT
logical :: verbose
!-locl
logical numeric
character(len=256) string,lfin,lfout
real chreal
integer narg,iargc,i,lmode
external numeric,chreal
logical :: MattCoords
logical :: TacoMaterials
!-init
verbose = .false.
MattCoords = .false.
TacoMaterials = .false.
lmode   = -999
lfin    = '<stdin>'
lfout   = '<stdout>'
nspec = 0

narg = COMMAND_ARGUMENT_COUNT()
i = 1

do while (i.le.narg) 
    call get_command_argument(i,string)
    if (string(1:3).eq.'-c ' .or. &
        string(1:3).eq.'-e ' .or. & 
        string(1:3).eq.'-l ' .or. &
        string(1:3).eq.'-n ' .or. &
        string(1:3).eq.'-I ' .or. &
        string(1:3).eq.'-s ') then
        if (lmode.ne.-999) then
            if (string(1:3).eq.'-n '.and.lmode.eq.2) then
                ! Okay to add another node tag
            elseif (string(1:3).eq.'-s '.and.lmode.eq.3) then
                ! Okay to add another label for element sides
            elseif (string(1:3).eq.'-f '.and.lmode.eq.4) then
                ! Okay to add another label for fault element/node pairs
            else
                write(stderr,10)
10                  format(1x,'gmsh2tecton: conflicting command line arguments.')
                call exitp(1)
            endif
        endif
    endif

    if (string(1:3).eq.'-c ') then
        lmode = 0

    else if (string(1:3).eq.'-e ') then
        lmode = 1

    else if (string(1:3).eq.'-n ') then

        if (i.eq.narg) then
            write(stderr,20)
20              format(1x,'Usage: gmsh2tecton -n number')
            call exitp(1)
        endif
        i = i + 1
        call get_command_argument(i,string)

        if (numeric(string)) then
            nspec = nspec + 1
!                if (nspec.gt.maxspec) then
!                    write(stderr,*) 'gmsh2tecton: too many physical labels specified'
!                    call exitp(1)
!                endif
            read( string, '(i10)' ) ISPEC(nspec)
        else
            write(stderr,20)
            call exitp(1)
        endif
        lmode = 2

    else if (string(1:3).eq.'-s ') then

        if (i.eq.narg) then
            write(stderr,30)
30              format(1x,'Usage: gmsh2tecton -s number')
            call exitp(1)
        endif

        i = i + 1
        call get_command_argument(i,string)

        if (numeric(string)) then
            nspec = nspec + 1
!                if (nspec.gt.maxspec) then
!                    write(stderr,*) 'gmsh2tecton: too many physical labels specified'
!                    call exitp(1)
!                endif
            read( string, '(i10)' ) ISPEC(nspec)
!                ISPEC = INT(chreal(string))
        else
            write(stderr,30)
            call exitp(1)
        endif

        lmode = 3

    else if (string(1:3).eq.'-f ') then

        if (i.eq.narg) then
            write(stderr,35)
35              format(1x,'Usage: gmsh2tecton -f number')
            call exitp(1)
        endif

        i = i + 1
        call get_command_argument(i,string)

        if (numeric(string)) then
            nspec = nspec + 1
!                if (nspec.gt.maxspec) then
!                    write(stderr,*) 'gmsh2tecton: too many physical labels specified'
!                    call exitp(1)
!                endif
            read( string, '(i10)' ) ISPEC(nspec)
!                ISPEC = INT(chreal(string))
        else
            write(stderr,35)
            call exitp(1)
        endif

        lmode = 4

    else if (string(1:3).eq.'-i ') then

        if (i.eq.narg) then
            write(stderr,40)
40              format(1x,'Usage: gmsh2tecton -i meshfile')
            call exitp(1)
        endif

        i = i + 1
        call get_command_argument(i,string)

        if (lfin.ne.'<stdin>') then
            write(stderr,10)
            call exitp(1)
        endif

        lfin = string

    else if (string(1:3).eq.'-o ') then

        if (i.eq.narg) then
            write(stderr,50)
50              format(1x,'Usage: gmsh2tecton -o outputfile')
            call exitp(1)
        endif
        i = i + 1
        call get_command_argument(i,string)

        if (lfout.ne.'<stdout>') then
            write(stderr,10)
            call exitp(1)
        endif

        lfout = string

    else if (string(1:3).eq.'-v ') then
        verbose = .true.

    else if (string(1:3).eq.'-m ') then
        MattCoords = .true.
!            write(stderr,52)
!52          format('gmsh2tecton: the "-m" option is deprecated')
    else if (string(1:3).eq.'-t ') then
        ! material numbers are printed as well.
        ! particulary useful for determining on which side of a fault
        ! an element is.
        TacoMaterials = .true.

    else if (string(1:3).eq.'-I ') then
        lmode = 5

    else if (string(1:3).eq.'-l ') then
        lmode = 6

    else

        write(stderr,*) ' no option ',trim(string)
        write(stderr,*)
        write(stderr,60)
60          format(1x,'Usage: gmsh2tecton [-c|-e|-n iphys|-s', &
         ' iphys|-f iphys] [-l] [-i meshfile] [-o tectonfile] [-I] [-m] [-t]')
            call exitp(1)
    endif

    i = i + 1

enddo

if (narg.eq.0) then
    write(stderr,60)
    call exitp(1)
endif

!
if (lmode.ne.-999) MODE = lmode
if (lfin.ne.'<stdin>') FIN = lfin
if (lfout.ne.'<stdout>') FOUT = lfout
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine INVLST (NUMNP,NNEL)
!
use meshData

implicit none
!-pass
integer(kind=8) :: NUMNP,NNEL
!-locl
integer i,j,n
logical test
data test/.false./
!
! * subroutine input:
!
!   NUMNP       number of (gmsh and gtecton) nodes
!   NODE(i)     node number list
!   X(3,i)      node coordinate list
!   
!   NNEL        number of gmsh elements
!   IEL(i)      gmsh element number list
!   ITYPE(i)    gmsh element type list
!   IPHYS(i)    element tag list
!   IEN(.,.)    list of nodes in gmsh element IEL
!
! * subroutine output:
!   NP(i)       list of new node numbers
!   NE(i)       list of new gmsh element numbers
!
!    Redefine node numbers to create a continuous and monotonous series.
!    To this end, create an "inverted" node number list; example:
!        i  NODE    NP
!        1     5    -
!        2     8    -
!        3     6    - 
!        4     4    4
!        5     -    1
!        6     -    3
!        7     -    -
!        8     -    2
j = NODE(1)
do i=2,NUMNP
    j = MAX(j,NODE(i))
enddo

!    if (j.gt.MXINP) then
!        write(stderr,10) j
!   10   format(1x,'gmsh2tecton: MXINP < ',I6)
!        call exitp(1)
!    endif

if (test) then
   open(unit=39,file='invnp.dat')
   rewind(39)
   write(39,11)
   11  format('    GMSH  TECTON NODES')
endif

do i=1,NUMNP
   NP(NODE(i)) = i
   if (test) write(39,12) NODE(i),i
   12  format(2I12)
enddo

if (test) close(39)

!   create inverted element number list with property NE(IEL(i))=i
j = 0
do i=1,NNEL
    if (ITYPE(i).eq.4) j = MAX(j,IEL(i))
enddo

!    if (j.gt.MXINE) then
!        write(stderr,20) j
!   20   format(1x,'gmsh2tecton: MXINE < ',I6)
!        call exitp(1)
!    endif

if (test) then
   open(unit=39,file='invelm.dat')
   rewind(39)
   write(39,21)
   21  format('    GMSH  TECTON ELEMENTS')
endif

n = 0
do i=1,NNEL
     if (ITYPE(i).eq.4) then
         n = n + 1
         NE(IEL(i)) = n
     endif
     if (test) write(39,22) IEL(i),NE(i)
   22    format(2I12)
enddo

if (test) then
    close(39)
endif

end subroutine
!-------------------------------------------------------------------------------
subroutine condadd (N,IA,NUM)
!
! conditional add of N to array IA: if not already in the array
!
implicit none
!-pass
integer N,IA,NUM
dimension IA(*)
!-locl
integer i
!
if (NUM.gt.0) then
    do i=1,NUM
        if (IA(i).eq.N) then
            return
        endif
    enddo
endif

NUM = NUM + 1
IA(NUM) = N
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine FCGEOM (X,NP,IEN,IFC,XN,AREA)

use io, only: stderr

!
! Calculates area of an element face, perpendicular to face normal,
! and the unit normal vector.
!
implicit none
!-pass
integer NP,IEN,IFC
real X,XN,AREA
dimension NP(*),IEN(8),X(3,*),XN(3)
!-locl
integer lien,i
real xl,xll
dimension xl(3,4),xll(3,3),lien(3)
!
!    localize coordinates
i = NP(IEN(1)) ! index number of coordinates
xl(1,1) = X(1,i)
xl(2,1) = X(2,i)
xl(3,1) = X(3,i)
i = NP(IEN(2))
xl(1,2) = X(1,i)
xl(2,2) = X(2,i)
xl(3,2) = X(3,i)
i = NP(IEN(3))
xl(1,3) = X(1,i)
xl(2,3) = X(2,i)
xl(3,3) = X(3,i)
i = NP(IEN(4))
xl(1,4) = X(1,i)
xl(2,4) = X(2,i)
xl(3,4) = X(3,i)

if (IFC.eq.1) then
    lien(1) = 1
    lien(2) = 2
    lien(3) = 4
elseif (IFC.eq.2) then
    lien(1) = 1
    lien(2) = 3
    lien(3) = 2
elseif (IFC.eq.3) then
    lien(1) = 1
    lien(2) = 4
    lien(3) = 3
elseif (IFC.eq.4) then
    lien(1) = 2
    lien(2) = 3
    lien(3) = 4
endif

!    establish face coordinates in xll
do i=1,3
    xll(1,i) = xl(1,lien(i))
    xll(2,i) = xl(2,lien(i))
    xll(3,i) = xl(3,lien(i))
enddo
!
XN(1)=(xll(2,2)-xll(2,1))*(xll(3,3)-xll(3,1))&
     -(xll(3,2)-xll(3,1))*(xll(2,3)-xll(2,1))
XN(2)=(xll(3,2)-xll(3,1))*(xll(1,3)-xll(1,1))&
     -(xll(1,2)-xll(1,1))*(xll(3,3)-xll(3,1))
XN(3)=(xll(1,2)-xll(1,1))*(xll(2,3)-xll(2,1))&
     -(xll(2,2)-xll(2,1))*(xll(1,3)-xll(1,1))

AREA = SQRT(XN(1)*XN(1)+XN(2)*XN(2)+XN(3)*XN(3))

if (AREA.le.0.0) then
    write(stderr,20) (xll(1,i),xll(2,i),xll(3,i),i=1,3)
   20        format(1x,'gmsh2tecton: zero AREA of element face. Face ',&
         'coordinates:',3(1x,/3(1PG14.6)))
    write(stderr,30) (xl(1,i),xl(2,i),xl(3,i),i=1,4)
   30        format(1x,'element coordinates:',4(1x,/3(1PG14.6)))
    write(stderr,40) lien(1),lien(2),lien(3)
   40        format(1x,'local node numbers: ',3i12)
       call exitp(1)
else
    XN(1)=XN(1)/AREA
    XN(2)=XN(2)/AREA
    XN(3)=XN(3)/AREA
    AREA = AREA*5E-1
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine vector (X,NP,IEN,NODE,IEM,V,ierr)
!
use io, only: stderr


implicit none
!-pass
integer NP,IEN,NODE,IEM,ierr
real X,V
dimension NP(*),X(3,*),IEN(4),IEM(3),V(3)
!-locl
integer k
double precision xce,yce,zce,p(3),xl(3,3),inprod,mag
!-init
ierr = 0
!
!    define vector from element center to slippery node
k = NP(IEN(1))

xce = DBLE(X(1,k))
yce = DBLE(X(2,k))
zce = DBLE(X(3,k))

k = NP(IEN(2))

xce = xce + DBLE(X(1,k))
yce = yce + DBLE(X(2,k))
zce = zce + DBLE(X(3,k))

k = NP(IEN(3))

xce = xce + DBLE(X(1,k))
yce = yce + DBLE(X(2,k))
zce = zce + DBLE(X(3,k))

k = NP(IEN(4))

xce = (xce + DBLE(X(1,k)))*0.25d0
yce = (yce + DBLE(X(2,k)))*0.25d0
zce = (zce + DBLE(X(3,k)))*0.25d0

k = NP(NODE)

p(1) = DBLE(X(1,k)) - xce
p(2) = DBLE(X(2,k)) - yce
p(3) = DBLE(X(3,k)) - zce

mag = p(1)*p(1)+p(2)*p(2)+p(3)*p(3)

if (mag.le.1e-19) then
    write(stderr,10)
   10        format(1x,'gmsh2tecton: zero direction vector error')
    ierr = 1
    return
endif
!
!    define normal vector to the fault face. Note that this is not necessarily the
!    face of the element that is being considered, so you cannot compute face area
!    from it

k = NP(IEM(1))

xl(1,1) = DBLE(X(1,k))
xl(2,1) = DBLE(X(2,k))
xl(3,1) = DBLE(X(3,k))

k = NP(IEM(2))

xl(1,2) = DBLE(X(1,k))
xl(2,2) = DBLE(X(2,k))
xl(3,2) = DBLE(X(3,k))

k = NP(IEM(3))

xl(1,3) = DBLE(X(1,k))
xl(2,3) = DBLE(X(2,k))
xl(3,3) = DBLE(X(3,k))

V(1)=(xl(2,2)-xl(2,1))*(xl(3,3)-xl(3,1))&
    -(xl(3,2)-xl(3,1))*(xl(2,3)-xl(2,1))
V(2)=(xl(3,2)-xl(3,1))*(xl(1,3)-xl(1,1))&
    -(xl(1,2)-xl(1,1))*(xl(3,3)-xl(3,1))
V(3)=(xl(1,2)-xl(1,1))*(xl(2,3)-xl(2,1))&
    -(xl(2,2)-xl(2,1))*(xl(1,3)-xl(1,1))

mag = DBLE(V(1)*V(1)+V(2)*V(2)+V(3)*V(3))

if (mag.le.1e-19) then
    write(stderr,20)
   20        format(1x,'gmsh2tecton: zero normal vector error')
       ierr = 1
    return
endif

inprod = p(1)*DBLE(V(1))+p(2)*DBLE(V(2))+p(3)*DBLE(V(3))

if (ABS(inprod).le.1e-19) then
    write(stderr,30)
   30        format(1x,'gmsh2tecton: direction vector perpendicular to ',&
           'fault normal vector')
    ierr = 1
    return
else if (inprod.lt.0.0) then
    V(1) = -V(1)
    V(2) = -V(2)
    V(3) = -V(3)
endif

mag = 1d0/SQRT(mag)

V(1) = SNGL(DBLE(V(1)*mag))
V(2) = SNGL(DBLE(V(2)*mag))
V(3) = SNGL(DBLE(V(3)*mag))

return
end subroutine
!-------------------------------------------------------------------------------
subroutine elmtest (elemID,X,IEN,NP,NUMNP,DET,IERR)

use io, only: stderr

!
! tests elements
!
implicit none

integer NEN,NSD
parameter (NEN=4,NSD=3)
real eps
parameter (EPS=1.0E-20)
!-pass
integer (kind=8) :: NUMNP
integer IEN(8),NP,IERR
real X,DET
dimension X(NSD,*),NP(*)
!-locl
integer i,j
integer p(4) ! indicent of the points
double precision :: xl(NSD,NEN)
integer :: tempnode ! for exchanging nodes when negative element
integer :: elemID   ! global(? really?) index of the element
!-init
IERR = 0
!

! let T be tetrahedron with 4 vertices with coordinates

! x1, y1, z1
! x2, y2, z2
! x3, y3, z3
! x4, y4, z4

! of which is moved with the fourth points to be on the origin,
! so that we get coordinates:

! x1-x4, y1-y4, z1-z4
! x2-x4, y2-y4, z2-z4
! x3-x4, y3-y4, z3-z4
!     0,     0,     0


! volume of element is determinant of matrix

! x1-x4, y1-y4, z1-z4
! x2-x4, y2-y4, z2-z4
! x3-x4, y3-y4, z3-z4




! determine global indices of nodal points
do i=1,4
   p(i) = NP(IEN(i))
enddo

! makes sure that p1...4 are sensible numbers
do i=1,4
    if (p(i).lt.1 ) then
        write(stderr,*) 'Error: node number ', p(i), 'less than 1'
        IERR = 1
        return
    else if (p(i).gt.NUMNP) then
        write(stderr,*) 'Error: node number ', p(i), &
  'greater than ', NUMNP
        IERR = 1
        return  
    endif
enddo

! fill the matrix of which we want the determinant

do i=1,3 ! over the points
    do j=1,3 ! over x,y,z
        ! X(coordinate, vertexID)
        xl(i,j) = X(j,p(i)) - X(j,p(4))
    enddo
enddo


det = xl(1,1) * ( xl(2,2)*xl(3,3) - xl(2,3)*xl(3,2) ) - &
      xl(2,1) * ( xl(1,2)*xl(3,3) - xl(1,3)*xl(3,2) ) + &
      xl(3,1) * ( xl(1,2)*xl(2,3) - xl(2,2)*xl(1,3) )

! somehow the determinant as computed here is the negative determinant
! This is an ugly solution...
det = -det

if (ABS(DET).le.eps) then
    write(*,*) "Error: element ", elemID, " is flat. Zero volume"
    write(*,*) "All the nodes are in a single plane."
    write(*,*) "xyz coords of node 1", ien(1), x(1,ien(1)), x(2,ien(1)), x(3,ien(1))
    write(*,*) "xyz coords of node 2", ien(2), x(1,ien(2)), x(2,ien(2)), x(3,ien(2)) 
    write(*,*) "xyz coords of node 3", ien(3), x(1,ien(3)), x(2,ien(3)), x(3,ien(3)) 
    write(*,*) "xyz coords of node 4", ien(4), x(1,ien(4)), x(2,ien(4)), x(3,ien(4)) 
    write(*,*) "the .msh file contains an error"
    ierr = 1
    return
endif

if (DET.lt.0.0) then
    ! element has the wrong order of nodes, so that elements becomes negative.
    ! exchange two nodes of the element, and all should be well.
    write(*,*) 'Warning: element ', elemID, ' is inverted.'
    write(*,*) 'Switching vertex 1 and 2 to fix'

    tempnode = IEN(1)
    IEN(1) = IEN(2)        
    IEN(2) = tempnode

    return
endif
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine skipblock(string,echo)

use io, only: lui

implicit none
!-pass
logical echo
character(len=*) string
!-locl
integer ios
character(len=256) input

100 read(lui,'(a256)',err=300,end=200,iostat=ios) input

if (trim(input).eq.string) then
    return
endif

if (echo) then
    write(*,10) trim(input)
10  format(256a)
endif

goto 100

200 write(0,1) trim(string)
  1 format(1x,'gmsh2tecton: section end not found in mesh file: ',a256)
call exitp(1)

300 write(0,2) trim(string)
  2 format(1x,'gmsh2tecton: mesh read error in data section: ',a256)
call exitp(1)

end subroutine
