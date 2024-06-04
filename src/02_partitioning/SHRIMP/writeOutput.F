subroutine writeOutput(outputFormat, nSteps)

use enumerates

implicit none

integer :: outputFormat, nSteps

if (outputFormat .eq. outputFormatSimple) then
    call writeSimple()
else if (outputFormat .eq. outputFormatGMT) then
    call writeGMT()
else if (outputFormat .eq. outputFormatVTK) then
    call writeVTK()
else if (outputFormat .eq. outputFormatGTECTON) then
    call writeGTECTON(nSteps)
else
     write (0,*) "output format:", outputFormat
      stop "outputformat not recognized. Leaving SHRIMP..."
endif


end subroutine


subroutine writeSimple()

use meshData
use moduleMPI, only: MPI_COMM_WORLD, &
                     myProcessorID, &
                     nProcessors

implicit none

integer :: iPass, iPoint, iElement, iErr
integer :: outputFileUnit

do iPass = 1, nProcessors

    outputFileUnit = 200 + iPass

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (myProcessorID .eq. iPass-1) then
        if (myProcessorID.eq.0) then
            ! open file to create a new file, or overwrite existing file.
            open(unit=outputFileUnit, file="partitionedElements.txt")
        else
            ! all higher thread append to the same file.
            open(unit=outputFileUnit, file="partitionedElements.txt", position='append', status='old')
        endif

   	    do iElement = 1,nElementsLocal
!            if (outputFormat .eq. outputFormatSimple) then
                write(outputFileUnit,*) partitionOfElements(iElement)
!            else if (outputFormat .eq. outputFormatGMT) then

                do iPoint = 1,nPointsPerElement
                    write(outputFileUnit,*) coordsForMyElems(:, iPoint, iElement), &
                                            partitionOfElements(iElement)
                enddo
!            else if (outputFormat .eq. outputFormatVTK) then
                write(outputFileUnit,*) iElement, partitionOfElements(iElement)
!            else
!                stop "outputformat not recognized. Leaving SHRIMP..."
!            endif

        enddo

	close(outputFileUnit)

    endif


    call MPI_Barrier(MPI_COMM_WORLD, ierr)

enddo

end subroutine


subroutine writeGMT()

use meshData
use moduleMPI, only: MPI_COMM_WORLD, &
       	       	     myProcessorID, &
       	       	     nProcessors
implicit none

integer :: iPass, iPoint, iElement, iErr
integer	:: outputFileUnit

outputFileUnit = 300 + iPass

call MPI_Barrier(MPI_COMM_WORLD, ierr)

if (myProcessorID .eq. 0) then
    open(unit=outputFileUnit, file="connectivityGMT.dat")


    do iElement = 1,nElementsGlobal
        write(outputFileUnit,*)  3*iElement-3, 3*iElement-2, 3*iElement-1
    enddo

close(outputFileUnit)

endif

call MPI_Barrier(MPI_COMM_WORLD, ierr)



end subroutine


subroutine writeVTK()

use meshData
use moduleMPI, only: MPI_COMM_WORLD, &
                     myProcessorID, &
                     nProcessors

implicit none

integer, parameter :: VTKfile = 72
integer :: iPass, iPoint, iElement, ierr
integer	:: outputFileUnit
character(len=12) nPoint_char, nElem_char




! write header 
if (myProcessorID.eq.0) then
!    write (nPoint_char, *) pointsPerElem * nElemsGlobal
    write (nPoint_char, *) 4 * nElementsGlobal
    write (nElem_char, *) nElementsGlobal

    open(unit=VTKfile, file="partitioning.vtu")
    write(VTKfile,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
    write(VTKfile,*) '   <UnstructuredGrid>'
    write(VTKfile,*) '	 <Piece NumberOfPoints="'//nPoint_char//'" NumberOfCells="'//nElem_char//' ">'
   !############# write the partitioning of the points
    write(VTKfile,*) '         <PointData Scalars="scalars">'
    write(VTKfile,*) '            <DataArray type="UInt32" NumberOfComponents="1" Name="vertex partitioning" Format="ascii">'

    close(VTKfile)

endif

call MPI_Barrier(MPI_COMM_WORLD, ierr)

! each thread writes it vertex partitioning in turn
do iPass = 1, nProcessors

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (iPass-1 .eq. myProcessorID) then

	    open(unit=VTKfile, file="partitioning.vtu", position='append', status='old')

        ! print the nodes per element.
        do iElement = 1, nElementsLocal
            do iPoint = 1, nPointsPerElement
                write(VTKfile,*) partitionOfElementPoints(iPoint,iElement)
            enddo
        enddo

	    close(VTKfile)

	endif

	call MPI_Barrier(MPI_COMM_WORLD, ierr)


enddo


! header for the element partitioning
if (myProcessorID.eq.0) then
    open(unit=VTKfile, file="partitioning.vtu", position='append', status='old')
    write(VTKfile,*) '            </DataArray>'
    write(VTKfile,*) '         </PointData>'
    write(VTKfile,*) '         <CellData Scalars="scalars">'
    write(VTKfile,*) '            <DataArray type="UInt32" NumberOfComponents="1" Name="Element partitioning" Format="ascii">'
    close(VTKfile)
endif


! each thread writes its own actual element partitioning

do iPass = 1, nProcessors

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (iPass-1 .eq. myProcessorID) then

        open(unit=VTKfile, file="partitioning.vtu", position='append', status='old')

        do iElement = 1, nElementsLocal
            write(VTKfile,*) partitionOfElements(iElement)
        enddo


        close(VTKfile)

    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)


enddo

! header for point coordinates

if (myProcessorID.eq.0) then
 open(unit=VTKfile, file="partitioning.vtu", position='append', status='old')
 write(VTKfile,*) '              </DataArray>'
 write(VTKfile,*) '          </CellData>'
 write(VTKfile,*) '          <Points>'
 write(VTKfile,*) '             <DataArray type="Float32" NumberOfComponents="3" Name="Vertex coordinates" Format="ascii">'
 close(VTKfile)

endif

! write the actual vertex coordinates

do iPass = 1, nProcessors

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (iPass-1 .eq. myProcessorID) then

        open(unit=VTKfile, file="partitioning.vtu", position='append', status='old')

		! print the nodes per element.
		do iElement = 1, nElementsLocal
            do iPoint = 1, nPointsPerElement
                write(VTKfile,*) coordsForMyElems(:,iPoint,iElement)
            enddo
        enddo

        close(VTKfile)

    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)


enddo


! write all the rest of the VTU file

if (myProcessorID.eq.0) then

    open(unit=VTKfile, file="partitioning.vtu", position='append', status='old')


    write(VTKfile,*) '              </DataArray>'
    write(VTKfile,*) '           </Points>'
    write(VTKfile,*) '           <Cells>'
    write(VTKfile,*) '              <DataArray type="UInt64" Name="connectivity" Format="ascii">'
    ! connectivity is not the actual connectivity of the mesh,
    ! but the separate connectivity for each element,
    ! in order to plot a sharp boundary between the different elements.
    do iElement = 1, nElementsGlobal
        write(VTKfile,*) 4 * iElement - 4, &
                         4 * iElement - 3, &
                         4 * iElement - 2, &
                         4 * iElement - 1
    enddo

    write(VTKfile,*) '              </DataArray>'
    write(VTKfile,*) '              <DataArray type="UInt32" Name="offsets" Format="ascii">'

    do iElement = 1, nElementsGlobal
        write(VTKfile,*) 4 * iElement
    enddo

    write(VTKfile,*) '             </DataArray>'
    write(VTKfile,*) '             <DataArray type="Int32" Name="types" Format="ascii">'
! element type... 10 is tetrahedron; see https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf, page 9
    do iElement = 1, nElementsGlobal
        write(VTKfile,*) 10
    enddo


    write(VTKfile,*) '             </DataArray>'
    write(VTKfile,*) '          </Cells>'
    write(VTKfile,*) '       </Piece>'
    write(VTKfile,*) '    </UnstructuredGrid>'
    write(VTKfile,*) ' </VTKFile>'

    close(VTKfile)


endif

end subroutine


subroutine writeGtecton(nSteps)

use meshData
use moduleMPI, only: MPI_COMM_WORLD, &
       	       	     myProcessorID, &
       	       	     nProcessors, &
                     MPI_INT, MPI_SUM
implicit none

character(len=12) :: outputTag_char

integer :: iPass, iPoint, iElement, iNeighbor, ierr
integer :: outputFileUnit, tagSize

integer :: nSteps

! the partition points have zigzagged around,
! and the last operation could van been either
! a zig or zag. Make sure we pic the right one.
if (modulo(nProcessors,2) .eq. 1) then
    partitionOfPoints => partitionOfPointsB
else
    partitionOfPoints => partitionOfPointsA
endif

if (modulo(nProcessors,2) .eq. 1) then
    pointCoords => pointCoordsB
else
    pointCoords => pointCoordsA
endif



do iPass = 1, nProcessors

    if (iPass-1 .eq. myProcessorID) then

        call setIndexTag(myProcessorID+1, nProcessors, outputTag_char, tagSize)

        outputFileUnit = 200 + iPass

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

!       -----------------------------------------------------------------------------
        open(unit=outputFileUnit, file='tecin.dat.partf.nps'//(outputTag_char(13-tagSize:12)))
!       -----------------------------------------------------------------------------


        do iPoint = 1, nPointsLocal

            write(outputFileUnit,'(I5,2I13,2x,1e24.17,2x,1e24.17,2x,1e24.17,I5)', advance='no') &
                  partitionOfPoints(iPoint) - 1, & ! count partitions from 0
                  iPoint + minNodeThisPartition - 1, &
                  pointMarkers(iPoint), &
                  pointCoords(:, iPoint), &
                  neighborCount(iPoint)

             do iNeighbor = 1, neighborCount(iPoint)-1
                    write(outputFileUnit,'(I13)',advance='no') neighborlist(iPoint)%neighbors(iNeighbor)

             enddo
             write(outputFileUnit,'(I13)') neighborlist(iPoint)%neighbors(neighborCount(iPoint))

        enddo

        close(outputFileUnit)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

    endif



enddo


do iPass = 1, nProcessors

    if (iPass-1 .eq. myProcessorID) then

        call setIndexTag(myProcessorID+1, nProcessors, outputTag_char, tagSize)

        outputFileUnit = 200 + iPass

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

!       -----------------------------------------------------------------------------
        open(unit=outputFileUnit, file='tecin.dat.partf.elm'//(outputTag_char(13-tagSize:12)))
!       -----------------------------------------------------------------------------

        do iElement = 1, nElementsLocal

            write(outputFileUnit, '(i5,6i13)') partitionOfElements(iElement) -1, & ! count partitions from 0
                                               iElement + elementOffset(myProcessorID+1), &
                                               elementMarkers(iElement), &
                                               connectivity(:,iElement)

        enddo

        close(outputFileUnit)

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

    endif
enddo



!------------------ create partition.info, containing the number of elements 
! and vertices for each partition, used by GTECTON to allocate local space

if (myProcessorID .eq. 0) then
	allocate(pointCountPerPartitionGlobal(2**nSteps))
	allocate(elementCountPerPartitionGlobal(2**nSteps))
	PointCountPerPartitionGlobal = 0
	elementCountPerPartitionGlobal = 0
endif

PointCountPerPartitionThread = 0
ElementCountPerPartitionThread = 0

do iPoint = 1, nPointsLocal
	PointCountPerPartitionThread(partitionOfPoints(iPoint)) = &
    PointCountPerPartitionThread(partitionOfPoints(iPoint)) + 1
enddo

do iElement = 1, nElementsLocal
    elementCountPerPartitionThread(partitionOfElements(iElement)) = &
    elementCountPerPartitionThread(partitionOfElements(iElement)) + 1
enddo

call MPI_Barrier(MPI_COMM_WORLD, ierr)

call MPI_Reduce( &
       PointCountPerPartitionThread, &  ! sending from
       PointCountPerPartitionGlobal, &  ! sending to
       2**nSteps, &    ! this many numbers in the array
       MPI_INT, &        ! data type
       MPI_SUM, &        ! Summation of weighted gives the average
       0, &              ! send to thread 0                
       MPI_COMM_WORLD, iErr)

call MPI_Reduce( &
       ElementCountPerPartitionThread, &  ! sending from
       ElementCountPerPartitionGlobal, &  ! sending to
       2**nSteps, &    ! this many numbers in the array
       MPI_INT, &        ! data type
       MPI_SUM, &        ! Summation of weighted gives the average
       0, &    	       	 ! send to thread 0                
       MPI_COMM_WORLD, iErr)

call MPI_Barrier(MPI_COMM_WORLD, ierr)


if (myProcessorID .eq. 0) then

	open(unit=401, file="partition.info")

	write(401,"(i10)") 2**nSteps

	do iPass = 1, 2**nSteps
		write(401,*) iPass-1, PointCountPerPartitionGlobal(iPass), ElementCountPerPartitionGlobal(iPass)
	enddo
	
	close(401)

	deallocate(ElementCountPerPartitionGlobal)
	deallocate(PointCountPerPartitionGlobal)


endif

deallocate(PointCountPerPartitionThread)
deallocate(ElementCountPerPartitionThread)

!-------------------- create a mergescript that will merge the output
! files into input files run by gtecton.

if (myProcessorID .eq. 0) then

	open(unit=402, file='mergeScript.sh')

	write(402,"(a)") "#!/bin/bash"

	call setIndexTag(1, nProcessors, outputTag_char, tagSize)

	write(402,"(a)") "cat tecin.dat.partf.nps"//(outputTag_char(13-tagSize:12))//" >| tecin.dat.partf.nps"
    write(402,"(a)") "cat tecin.dat.partf.elm"//(outputTag_char(13-tagSize:12))//" >| tecin.dat.partf.elm"

    do iPass=2,nProcessors

	    call setIndexTag(iPass, nProcessors, outputTag_char, tagSize)

        write(402,"(a)") "cat tecin.dat.partf.nps"//(outputTag_char(13-tagSize:12))//" >> tecin.dat.partf.nps"
        write(402,"(a)") "cat tecin.dat.partf.elm"//(outputTag_char(13-tagSize:12))//" >> tecin.dat.partf.elm"

	enddo

	write(402,"(a)") "echo 'end partitioned nodes' >> tecin.dat.partf.nps"
    write(402,"(a)") "echo 'end partitioned elements' >> tecin.dat.partf.elm"

	close(402)

endif

end subroutine

!-----------------------------------------------------------------

subroutine setIndexTag(idx, total, tag, tagSize)

implicit none

integer idx, total, tagsize, iZero
character(len=12) tag

write (tag, *) idx

tagSize = log10(real(total-1)) + 1

do iZero = 1, 12 - tagSize
	tag(iZero:iZero) = "0"
enddo

end subroutine
