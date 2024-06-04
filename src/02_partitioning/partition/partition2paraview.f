program partition2paraview

! this program reads in a partitioned mesh,
! assumed to be written to generated files
! tecin.dat.part.nps/elm

! It writes this to a form paraview which contains
! the coordinates, the connectivity and the partitioning.


! input
character(len=19), parameter   :: vertexFile    = "tecin.dat.partf.nps"
character(len=19), parameter   :: elementFile   = "tecin.dat.partf.elm"

character(len=44), parameter   :: vertexFormat  = "(I5,I12,1x,I5,1x,1e25.17,1x,1e25.17,I6,6I13)"
character(len=9),  parameter   :: elementFormat = "(I5,6I12)"

! output
character(len=8),  parameter   :: outputfile    = "mesh.vtu"



double precision,  allocatable :: x(:), y(:)
integer,           allocatable :: v1(:), v2(:), v3(:), v4(:)
integer,           allocatable :: vertexPartition(:)
integer,           allocatable :: elementPartition(:)

integer                        :: nVertices, nElements
integer                        :: iVertex,   iElem

integer                        :: dummyInt

character(len=12) nvert_char, nelem_char


! ************** read vertices ************

! first count the vertices
nVertices = 0
open(unit=42, file=vertexFile)
do
    read(42,*,end=10)
    nVertices = nVertices + 1
enddo
10 close(42)
nVertices = nVertices - 1

! make space
allocate(x(nVertices))
allocate(y(nVertices))
allocate(vertexPartition(nVertices))

! and read them
open(unit=44, file=vertexFile)
do iVertex=1,nVertices
    ! dummies are index and marker
    read(44,*) vertexPartition(iVertex), dummyInt, dummyInt, x(iVertex), y(iVertex)
enddo

! ************** read elements ************

! count the elements
nElements = 0
open(unit=45, file=elementFile)
do
    read(45,*,end=20)
    nElements = nElements + 1
enddo
20 close(45)
nElements = nElements - 1

! make space 
allocate(v1(nElements))
allocate(v2(nElements))
allocate(v3(nElements))
allocate(v4(nElements))
allocate(elementPartition(nElements))

! and read them
open(unit=46, file=elementFile)
do iElem=1,nElements
    read(46,*) elementPartition(iElem), &
               dummyInt, &
               dummyInt, &
               v1(iElem), &
               v2(iElem), &
               v3(iElem), &
               v4(iElem)
enddo

close(46)

! ******************* write vtk file **********

write (nvert_char, *) nVertices
write (nelem_char, *) nElements


open (unit = 47, file = outputfile)

write(47,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(47,*) '   <UnstructuredGrid>'
write(47,*) '      <Piece NumberOfPoints="'//nvert_char//'" NumberOfCells="'//nelem_char//' ">'

!***************************************************************

write(47,*) '         <PointData Scalars="scalars">'

write(47,*) '             <DataArray type="UInt32" NumberOfComponents="1" Name="vertex partitioning" Format="ascii">'
do iVertex=1,nvertices
    write(47,*) vertexPartition(iVertex)
enddo
write(47,*) '             </DataArray>'

write(47,*) '         </PointData>'

!***************************************************************

write(47,*) '         <CellData Scalars="scalars">'
write(47,*) '             <DataArray type="UInt32" NumberOfComponents="1" Name="Element partitioning" Format="ascii">'
do iElem=1,nElements
    write(47,*) (elementPartition(iElem))
enddo

    write(47,*) '             </DataArray>'

write(47,*) '         </CellData>'

!***************************************************************

write(47,*) '         <Points>'
write(47,*) '            <DataArray type="Float32" NumberOfComponents="3" Name="Vertex coordinates" Format="ascii">'
do iVertex=1,nvertices
    write(47,*) x(iVertex),y(ivertex),0
enddo
write(47,*) '            </DataArray>'

write(47,*) '         </Points>'

!***************************************************************

write(47,*) '         <Cells>'
write(47,*) '            <DataArray type="UInt64" Name="connectivity" Format="ascii">'
do iElem=1,nElements
    write(47,*) v1(iElem)-1, v2(iElem)-1, v3(iElem)-1
enddo
write(47,*) '            </DataArray>'


write(47,*) '            <DataArray type="UInt32" Name="offsets" Format="ascii">'
do iElem=1,nElements
    write(47,*) 3*iElem
enddo
write(47,*) '            </DataArray>'


write(47,*) '            <DataArray type="UInt32" Name="types" Format="ascii">'
do iElem=1,nElements-1
! 10 is paraviewinese for 'tetrahedron'.
! 5 is triangle
    write(47,'(a3)',advance="no") '5 ' ! write them all on one line to keep it compact
enddo
write(47,'(a2)') '5' ! to force the linebreak
write(47,*) '            </DataArray>'



write(47,*) '         </Cells>'
write(47,*) '      </Piece>'

write(47,*) '   </UnstructuredGrid>'
write(47,*) '</VTKFile>'


close(47)





end program
