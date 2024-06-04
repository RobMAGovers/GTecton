program useEncode

use base64Encode, only: EncodeDoublesToFile

implicit none

double precision, dimension(4)   :: doublesPyra = (/1.0d0, -2.0d0, 3.0d0, 4.0d0/)


open(unit=62, file = "oneTetrahedron.vtu")


write(62,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
write(62,*)'  <UnstructuredGrid>'
write(62,*)'    <Piece NumberOfPoints="4" NumberOfCells="1">'

write(62,*)'      <PointData Scalars="scalars">'
write(62,*)'        <DataArray type="Float32" Name="scalars" format="binary">'
!************************* lo and behold *******************
! call with filehandle, array of doubles, size of array, precision (1 = single; 2 = double).
call EncodeDoublesToFile(62, doublesPyra, 4, 1)
!***********************************************************
write(62,*)'        </DataArray>'
write(62,*)'      </PointData>'

write(62,*)'      <Points>'

write(62,*)'        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
write(62,*)'          0 0 0 0 0 1 0 1 0 1 0 0'
write(62,*)'        </DataArray>'

write(62,*)'      </Points>'

write(62,*)'      <Cells>'

write(62,*)'        <DataArray type="Int32" Name="connectivity" format="ascii">'
write(62,*)'          0 1 2 3'
write(62,*)'        </DataArray>'

write(62,*)'        <DataArray type="Int32" Name="offsets" format="ascii">'
write(62,*)'          4'
write(62,*)'        </DataArray>'

write(62,*)'        <DataArray type="UInt8" Name="types" format="ascii">'
write(62,*)'          10'
write(62,*)'        </DataArray>'

write(62,*)'      </Cells>'

write(62,*)'    </Piece>'

write(62,*)'  </UnstructuredGrid>'

write(62,*)'</VTKFile>'



close(62)


end program
