module VTKmodule

 implicit none

 integer :: nprocs  ! to determine nr of piece files in the VTK file

 integer fp

 character(len=8), parameter :: mainfile = 'main.vtk'
 character(len=8), parameter :: piecefile = 'data.vtp'

 character(len=8), parameter :: VTKOutputFileName = 'data.vtu'
 integer,          parameter :: vtufile = 833


 

 !TODO, generate a way to add a counter to data.vtp files.

 public :: VTKwrite
 public :: VTKInititialisePlotoptions

 type plotoptions
logical prestress
logical materials
logical velocity
logical displacement
logical stress
logical strain
logical elempartition
logical vertexpartition
logical writeASCII
 end type plotoptions

 type(plotoptions) plotthis

contains 

subroutine VTKInititialisePlotoptions
 implicit none
!    type(plotoptions plotthis

 plotthis%prestress       = .false.
 plotthis%materials       = .false.
 plotthis%velocity        = .false.
 plotthis%displacement    = .true.
 plotthis%stress          = .false.
 plotthis%strain          = .false.
 plotthis%elempartition   = .false.
 plotthis%vertexpartition = .false.
 plotthis%writeASCII      = .true.

end subroutine


subroutine VTKwrite(filename)

#include "lenchar.i"

 character(len=LENCHAR) :: filename

 open (unit = vtufile, file = filename)

 call WriteVTUfile(vtufile)

 close (vtufile)


end subroutine


subroutine WriteVTUfile(fp)

 use  meshdatamodule, only:  meshdatactx
 use modeldatamodule, only: modeldatactx

 implicit none

 integer :: fp, i, n

 ! write number of vertices and elements to characters, to ouput them in the file.
 character(len=12) nvert_char, nelem_char
 write (nvert_char, *) meshdatactx%nvglobal
 write (nelem_char, *) meshdatactx%neglobal


 write(fp,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
 write(fp,*) '   <UnstructuredGrid>'
 write(fp,*) '      <Piece NumberOfPoints="'//nvert_char//'" NumberOfCells="'//nelem_char//' ">'

!*********** point data ******

 ! write displacement
 write(fp,*) '         <PointData Scalars="scalars">'


!    write(fp,*) '             <DataArray type="Float32" NumberOfComponents="3" Name="Displacement" Format="ascii">'
!    do i=1,meshdatactx%nvglobal
!        write(fp,*) modeldatactx%D(1,i), modeldatactx%D(2,i), modeldatactx%D(3,i)
!    enddo
!    write(fp,*) '             </DataArray>'

 write(fp,*) '         </PointData>'

!*******  Cell data
! LWITODO; stress and strain tensors are being written here. This will have to be changed to good visualisation.

 write(fp,*) '         <CellData Scalars="scalars">'

  if (plotthis%stress) then
 write(fp,*) '             <DataArray type="Float32" NumberOfComponents="1" Name="Stress" Format="ascii">'
 do i=1,meshdatactx%neglobal
write(fp,*) (modeldatactx%STR(n,i),n=1,6)
 enddo
 write(fp,*) '             </DataArray>'
  endif

!    write(fp,*) '             <DataArray type="Float32" NumberOfComponents="1" Name="Strain" Format="ascii">'
!    do i=1,meshdatactx%neglobal
!        write(fp,*) (modeldatactx%STN(n,i),n=1,6)
!    enddo
!    write(fp,*) '             </DataArray>'

 write(fp,*) '         </CellData>'


!****** vertices **********

 write(fp,*) '         <Points>'

 ! write coordinate points
 write(fp,*) '            <DataArray type="Float32" NumberOfComponents="3" Name="coordinates" Format="ascii">'
 do i=1,meshdatactx%nvglobal
write(fp,*) meshdatactx%X(1,i), meshdatactx%X(2,i), meshdatactx%X(3,i)
 enddo
 write(fp,*) '            </DataArray>'

 write(fp,*) '         </Points>'

! ********** Cells & Connectivity *************

 write(fp,*) '         <Cells>'

 ! write connectivity
 write(fp,*) '            <DataArray type="Int32" Name="connectivity" Format="ascii">'
 do i=1,meshdatactx%neglobal
write(fp,*) meshdatactx%ien(1,i)-1, meshdatactx%ien(2,i)-1, meshdatactx%ien(3,i)-1, meshdatactx%ien(4,i)-1
 enddo
 write(fp,*) '            </DataArray>'

 write(fp,*) '            <DataArray type="Int32" Name="offsets" Format="ascii">'
 do i=1,meshdatactx%neglobal
write(fp,*) 4*i
 enddo
 write(fp,*) '            </DataArray>'

 write(fp,*) '            <DataArray type="Int32" Name="types" Format="ascii">'
 do i=1,meshdatactx%neglobal
write(fp,*) '10'
 enddo
 write(fp,*) '            </DataArray>'

 write(fp,*) '         </Cells>'

!************
 write(fp,*) '      </Piece>'
 write(fp,*) '   </UnstructuredGrid>'
 write(fp,*) '</VTKFile>'




end subroutine






subroutine WriteVTKHeader(fp)

 integer :: fp

 ! lwitodo; for niceties we should check en verify whether fp is opened
 ! and handle accordingly when it is not.

end subroutine


end module
