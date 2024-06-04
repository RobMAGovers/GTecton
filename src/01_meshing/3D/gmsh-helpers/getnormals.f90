module cli

    implicit none
    ! Inputs from command line
    character(len=200) :: nodelistfile      ! File with list of node numbers for computing normal vectors
    character(len=200) :: meshfile          ! Gmsh output mesh file
    integer,dimension(20) :: labels         ! Labels of surfaces to include in normal vector calculation
    integer :: nlabels                      ! Number of surface labels
end module

module arrays
    implicit none
    ! Arrays
    integer, allocatable :: nodelist(:)     ! array with list of specified nodes
    real, allocatable    :: coords(:,:)
    integer, allocatable :: triangles(:,:)
    integer, allocatable :: trianglesatnodes(:,:)
    real, allocatable :: normals(:,:)
    ! Array dimensions
    integer :: nnodelist                    ! number of specified nodes
    integer :: ncoords
    integer :: ntriangles
end module

!--------------------------------------------------------------------------------------------------!

program main

use cli
use arrays

implicit none
integer :: i

! Parse command line and check inputs (variables in module cli)
write(0,*) 'getnormals: getting arguments from command line'
call gcmdln()
write(0,*) '    Node list file: ',trim(nodelistfile)
write(0,*) '    Mesh file:      ',trim(meshfile)
write(0,*) '    Surface labels: ',(labels(i),i=1,nlabels)
call chkin()

! Read node list into array nodelist(nnodelist)
call readnodelist()

! Read arrays coords(ncoords,3), triangles(ntriangles,4) from Gmsh file
call readmesh()

! Compute normal vector for each triangle and put into array normals(ntriangles,3)
call trianglenormals()

! Make an array linking all nodes to triangles with specified labels named trianglesatnodes(ncoords,21)
call linknodestriangles()
deallocate(triangles)

! Compute mean normal at each node in list and print results to standard output
call normalatnode()

deallocate(nodelist)
deallocate(coords)
deallocate(normals)
deallocate(trianglesatnodes)
end

!--------------------------------------------------------------------------------------------------!

subroutine chkin()

use cli

implicit none

logical :: ex
integer :: i

if (nodelistfile.eq.'') then
    call usage('!! Error: Node list file is not defined')
else
    inquire(file=nodelistfile,exist=ex)
    if (.not.ex) then
        call usage('!! Error: No node list file found named '//trim(nodelistfile))
    endif
endif

if (meshfile.eq.'') then
    call usage('!! Error: Mesh file is not defined')
else
    inquire(file=meshfile,exist=ex)
    if (.not.ex) then
        call usage('!! Error: No mesh file found named '//trim(meshfile))
    endif
endif

if (nlabels.eq.0) then
    call usage('!! Error: No surface labels defined')
endif

do i = 1,nlabels
    if (labels(i).lt.0) then
        call usage('!! Error: Found negative label in input list')
    endif
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine readnodelist()

use cli, only : nodelistfile
use arrays, only : nodelist, nnodelist

implicit none

integer :: i

! count number of nodes in list
nnodelist = 0
open(unit=11,file=nodelistfile,status='old')

do
    read(11,*,end=111) 
    nnodelist = nnodelist + 1
enddo

111 rewind(11)
! allocate memory for the array
allocate(nodelist(nnodelist))

! read the list
do i = 1,nnodelist
    read(11,*,end=111) nodelist(i)
enddo

close(11)
return
end

!--------------------------------------------------------------------------------------------------!

subroutine readmesh()

use cli, only : meshfile
use arrays, only: coords, ncoords, triangles, ntriangles

implicit none

character(len=1000) :: text, input_line
integer :: nelements
integer :: i, j
integer, allocatable :: triarraytmp(:,:)

ntriangles = 0

open(unit=12,file=meshfile,status='old')

do
    read(12,*,end=121) text

    ! Locate node block
    if (trim(text).eq.'$Nodes') then
        write(0,*) 'getnormals: reading nodal coordinates'
        read(12,*) ncoords
        ! Allocate memory for coordinates
        allocate(coords(ncoords,3))
        ! Read nodal coordinates
        do i = 1,ncoords
            read(12,*) coords(i,1),(coords(i,j),j=1,3)
        enddo
        write(0,*) 'getnormals: finished reading nodal coordinates'
    endif

    ! Locate element block
    if (trim(text).eq.'$Elements') then
        write(0,*) 'getnormals: reading elements'

        ! Total number of elements
        read(12,*) nelements

        ! Allocate array for all elements (too big for now)
        allocate(triarraytmp(nelements,4))

        ! Read only 2D triangular elements into array
        do i = 1,nelements
            read(12,'(A)') input_line
            read(input_line,*) j,j
            if (j.eq.2) then
                ntriangles = ntriangles + 1
                read(input_line,*) j,j,j,triarraytmp(ntriangles,4),j,&
                                   (triarraytmp(ntriangles,j),j=1,3)
            endif
        enddo
        write(0,*) 'getnormals: finished reading triangles; found ntriangles=',ntriangles

        ! Allocate memory for triangle coordinates (1-3) and labels (4)
        allocate(triangles(ntriangles,4))

        ! Save values to new array, deallocate old big array
        triangles = triarraytmp(1:ntriangles,:)
        deallocate(triarraytmp)

    endif
enddo
121 close(12)
return
end

!--------------------------------------------------------------------------------------------------!

subroutine linknodestriangles()

use cli,    only : labels, nlabels
use arrays, only : ncoords, triangles, ntriangles, trianglesatnodes

implicit none

integer :: i, j, k
integer :: nadjacent

write(0,*) 'getnormals: making array linking nodes to triangles'
! allocate linking array
allocate(trianglesatnodes(ncoords,21))
trianglesatnodes = 0

do i = 1,ntriangles
    do j = 1,nlabels
        if (triangles(i,4).eq.labels(j)) then
            do k = 1,3
                ! First column is number of adjacent triangles
                nadjacent = trianglesatnodes(triangles(i,k),1)
                if (nadjacent.ge.20) then
                    call usage('!! Error: Node has more than 20 adjacent triangles')
                endif
                ! Load new entry
                nadjacent = nadjacent + 1
                trianglesatnodes(triangles(i,k),nadjacent+1) = i
                ! Update number of adjacent triangles in array
                trianglesatnodes(triangles(i,k),1) = nadjacent
            enddo
            exit
        endif
    enddo
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine trianglenormals()

use arrays, only : coords, triangles, ntriangles, normals

implicit none

real    :: v1(3), v2(3)
integer :: i

! allocate normal vector
allocate(normals(ntriangles,3))
do i = 1,ntriangles
    ! Make two vectors from point 1 to points 2 and 3
    v1(1) = coords(triangles(i,2),1) - coords(triangles(i,1),1)
    v1(2) = coords(triangles(i,2),2) - coords(triangles(i,1),2)
    v1(3) = coords(triangles(i,2),3) - coords(triangles(i,1),3)
    v2(1) = coords(triangles(i,3),1) - coords(triangles(i,1),1)
    v2(2) = coords(triangles(i,3),2) - coords(triangles(i,1),2)
    v2(3) = coords(triangles(i,3),3) - coords(triangles(i,1),3)
    ! Cross product is normal vector
    normals(i,1) = v1(2)*v2(3) - v1(3)*v2(2)
    normals(i,2) = v1(3)*v2(1) - v1(1)*v2(3)
    normals(i,3) = v1(1)*v2(2) - v1(2)*v2(1)
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine normalatnode()
use arrays, only : nodelist, nnodelist, coords, normals, trianglesatnodes
implicit none
integer :: i, j, k
integer :: nadjacent
integer :: itriangle
real :: totalnormal(3)
real :: dotproduct
real :: length
do i = 1,nnodelist
    ! Exit with error if node has no specified neighbors
    nadjacent = trianglesatnodes(nodelist(i),1)
    if (nadjacent.eq.0) then
        write(0,*) '!! Error: Node ',nodelist(i),' has no neighboring triangles with specified labels'
        stop
    endif
    ! Add normal vector from each neighboring triangle
    totalnormal = 0.0e0
    do j = 1,nadjacent
        itriangle = trianglesatnodes(nodelist(i),j+1)
        if (j.gt.1) then
            ! Make sure normals point in about the same direction
            dotproduct = totalnormal(1)*normals(itriangle,1) &
                             + totalnormal(2)*normals(itriangle,2) &
                             + totalnormal(3)*normals(itriangle,3)
            if (dotproduct.lt.0.0) then
                do k = 1,3
                    normals(itriangle,k) = -normals(itriangle,k)
                enddo
            endif
        endif
        do k = 1,3
            totalnormal(k) = totalnormal(k) + normals(itriangle,k)
        enddo
    enddo
    ! Normalize total normal vector to unit length
    length = sqrt(totalnormal(1)*totalnormal(1) &
                      + totalnormal(2)*totalnormal(2) &
                      + totalnormal(3)*totalnormal(3))
    do k = 1,3
        totalnormal(k) = totalnormal(k)/length
    enddo
    write(*,*) (coords(nodelist(i),k),k=1,3),(totalnormal(k),k=1,3)
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use cli

implicit none

integer :: narg
integer :: i
character(len=200) :: arg

! Initialize variables
nodelistfile = ''
meshfile = ''
do i = 1,20
    labels(i) = 0
enddo
nlabels = 0

! Parse command line
narg = command_argument_count()
i = 1
do while (i.le.narg)
    call get_command_argument(i,arg)
    if (trim(arg).eq.'-n'.or.trim(arg).eq.'--nodelist') then
        i = i + 1
        call get_command_argument(i,nodelistfile)
    elseif (trim(arg).eq.'-m'.or.trim(arg).eq.'--meshfile') then
        i = i + 1
        call get_command_argument(i,meshfile)
    elseif (trim(arg).eq.'-l'.or.trim(arg).eq.'--label') then
        i = i + 1
        nlabels = nlabels + 1
        call get_command_argument(i,arg)
        read(arg,*) labels(nlabels)
    else
        call usage('!! Error: No option '//trim(arg))
    endif
    i = i + 1
enddo

return
end

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)
implicit none
character(len=*) :: str
if (str.ne.''.and.str.ne.'long') then
    write(0,*) str
    write(0,*)
endif
write(0,*) 'Usage: getnormals -n|--nodelist NODELIST -m|--meshfile MESHFILE -l|--label LABEL [-l LABEL...]'
write(0,*) '    -n NODELIST   List of nodes where normals are desired'
if (str.eq.'long') then
    write(0,*) 
endif
write(0,*) '    -m MESHFILE   GMSH output mesh file'
if (str.eq.'long') then
    write(0,*) 
endif
write(0,*) '    -l LABEL      Surface labels to include in normal vector calculation'
if (str.eq.'long') then
    write(0,*) 
endif
stop
return
end

