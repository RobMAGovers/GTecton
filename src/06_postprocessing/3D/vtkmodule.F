module VTKmodule

! after all the options has been set in setcmd3D, this module the writing of the VTU file.

!*************************************************************************
!  A mesh contains every nodal point only once.
!  For a mesh without slippery or faulted nodes this is fine.
!  However, if we want to draw a mesh with a slippery surface or a fault,
!  the nodal points on the slippery/faulted surfaces are in multiple places
!  at once. Also, Paraview only wants one quantity per node.
!
!  The simple way to resolve this would be to plot every element separately,
!  and to clone each nodal point for every element to which it is attached.
!  This will multiply the total number of nodal points by about 13, making
!  the files much bigger, because the nodal quantities section will also
!  increase in size. In most models with slippery/faulted nodes,
!  almost all the of the nodes will be regular nodes, so this cloning 
!  is very inefficient.
!  
!  A more intelligent way to handle this is to clone only the 
!  slippery/faulted nodes, and leave all the other nodes as they were.
!  This is done below.
!
!  To do this, we have to modify the connectivity of the mesh.
!  We first count the extra space needed, then allocate this space,
!  then fill it with the modified mesh. That way we do not have to reallocate
!  the arrays for every node clone.


!  The slippery quantities will be defined based on which element it is, which can
!  be extracted easily from the DELDXE and DX arrays.

!  redefining the mesh will also changed the nodes on which the boundary conditions work!

!  It will also increase the number of nodal points. The number of elements will remain unchanged


implicit none

type plotoptions
logical :: prestress
logical :: materials
logical :: slipperyMarkers
logical :: slipperyVectors
logical :: slipperyWeights
logical :: faultedDisplacement
logical :: faultedVelocity
logical :: slipperyDisplacement
logical :: slipperyVelocity
logical :: velocity
logical :: velocityBC
logical :: displacement
logical :: displacementBC
logical :: forceBC
logical :: slideInPlaneBC ! special displacement BC
logical :: slideOnLineBC ! special displacement BC
logical :: fixedPointBC ! special displacement BC
logical :: EulerAnglesBC
logical :: stress
logical :: strain
logical :: strainrate
logical :: viscosity
logical :: volume
logical :: elempartition
logical :: vertexpartition
logical :: vertexmarker
logical :: writeASCII
integer :: precision
logical :: elemQualitySide
logical :: elemQualitySurface
logical :: elemQualityVolume
logical :: error
#ifdef EXPERIMENTALTRANSIENTDISCREEP
logical :: taylorstress
#endif
end type plotoptions

type(plotoptions) plotthis

logical :: connDone
integer :: nOutputTags
integer :: outputTags(100)
logical :: gmshRead



integer :: nExtraNodes ! the number of nodes that is created by duplication to 
                       ! accommodate the slippery and split nodes. 
integer :: nAllNodes   ! the number of nodes after the slippery and split nodes
                       ! have been duplicated.

! array for binary printing of offsets
! They come by multiples of fours 4,8,12,16 etc.
! and as such are easily generated for ascii writing
! However, the binary writing demands an array, which we have to 
! fill first.
integer, allocatable :: offsets(:)

integer, allocatable :: referenceNode(:)   ! newly created notes have the value of this one for disp and such
integer, allocatable :: referenceElem(:)   ! element nr of newly created entry to which the new belongs. 
                                           ! Also contains the element of the first node, to which all the order nodes refer.
                                           ! Because of this, the referenceElem array is longer than the referenceNode array.
                                           ! For ease, we allocate this to size nAllNodes.
integer, allocatable :: firstReferenceElem(:) ! to contain the reference element of the original node,
                                              ! when the total is not yet known.

character(len=8), parameter :: VTUOutputFileName = 'data.vtu' ! the default, that can be changed in subroutine xvtu in xeclib
integer, parameter          :: fp = 833 ! file pointer to the VTU file

integer, allocatable :: ienDuplicate(:,:)   ! a new connectivity array to put the modified connectivity
                                            ! in, when the split and slippery nodes have been duplicated.
double precision, allocatable :: allNodeCoords(:,:) ! coordinates including the duplicate nodes
integer, allocatable :: faultEntryLookup(:) ! ... to describe precisely

! to scale the default arrows for nodal boundary conditions
! so that paraview does not crash by havign to draw huuuuuuge arrows
double precision, parameter :: scaleFactor = 1.0e-0

logical :: meshAdjusted

interface nodalQtyToAllNodesQty
    module procedure nodalDoubleQtyToAllNodesQty, &
                     nodalIntegerQtyToAllNodesQty
end interface

contains


!-----------------------------------------------------------------------------


subroutine WriteVTUfile(fileName, deformed)

use modeldefinition, only: NUMFNglobal, &
                           NUMSLPglobal
use materialsmodule, only: matdatactx, &
                           MaterialIDOfElement
use meshdatamodule,  only: meshdatactx
use constants,       only: elemQualitySide, &
                           elemQualitySurface, &
                           elemQualityVolume, &
                           dispBC, &
                           veloBC, &
                           EulerAnglesBC, &
                           fixedPointBC, &
                           forceBC, &
                           slideInPlaneBC, &
                           slideOnLineBC
use     debugmodule, only: iecho
#ifdef EXPERIMENTALTRANSIENTDISCREEP
use modeltopology, only: NSTRTAY
#endif

implicit none

#include "lenchar.i"

character(len=LENCHAR) :: filename
logical                :: deformed

character(len=12)      :: nvert_char, nelem_char

open (unit = fp, file = fileName)

if (NUMSLPglobal.gt.0 .or. NUMFNglobal.gt.0) then
    if (.not. meshAdjusted) then
        ! it is possible to run multiple write commands in one 
        ! plt3d script. The mesh needs to be adapted only the first time.
	    call adjustMeshForSlipAndFaults()
        meshAdjusted = .true.
	endif
    write (nvert_char, *) nAllNodes
    ! the new connectivity will be written in IENduplicate
else
    write (nvert_char, *) meshdatactx%nvglobal
endif
write (nelem_char, *) meshdatactx%neglobal

! write the header of the VTU file
write(fp,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(fp,*) '   <UnstructuredGrid>'
write(fp,*) '      <Piece NumberOfPoints="'//nvert_char//'" NumberOfCells="'//nelem_char//' ">'

! the data order of the XML is determiend by the VTU standard:
!<PointData>...</PointData>   ! simulation result data per vertex,    optional
!<CellData>...</CellData>     ! simulation result data per element,   optional
!<Points>...</Points>         ! nodal coordinates of the vertices,    mandatory
!<Cells>...</Cells>           ! connectivity, element types,          mandatory


! <PointData> block
write(fp,*) '         <PointData Scalars="scalars">'
! displacement and velocity cannot be turned off.
call writeDisplacement()
call writeVelocity()
if (plotthis%vertexpartition) then
    call writeVertexPartition()
endif
if (plotthis%vertexMarker) then
    call writeVertexMarker()
endif
if (plotthis%slipperyDisplacement) then
    call writeSlipperyDisplacement()
endif
if (plotthis%slipperyVelocity) then
    call writeSlipperyVelocity()
endif
if (plotthis%slipperyWeights) then
    call writeSlipperyWeights()
endif
if (plotthis%faultedDisplacement) then
    call writeFaultedDisplacement()
endif
if (plotthis%faultedVelocity) then
    call writeFaultedVelocity()
endif
!-- boundary conditions
if (plotthis%velocityBC) then
    call writeNodalBC(veloBC)
endif
if (plotthis%displacementBC) then
    call writeNodalBC(dispBC)
endif
if (plotthis%forceBC) then
    call writeNodalBC(forceBC)
endif
if (plotthis%slideOnLineBC) then
    call writeNodalBC(slideOnLineBC)
endif
if (plotthis%slideInPlaneBC) then
    call writeNodalBC(slideInPlaneBC)
endif
if (plotthis%fixedPointBC) then
    call writeNodalBC(fixedPointBC)
endif
if (plotthis%EulerAnglesBC) then
    call writeNodalBC(EulerAnglesBC)
endif
write(fp,*) '         </PointData>'


! <CellData> block
write(fp,*) '         <CellData Scalars="scalars">'
if (plotthis%prestress) then
    call writePrestress()
endif
if (plotthis%stress) then
    call writeStress()
endif

#ifdef EXPERIMENTALTRANSIENTDISCREEP
if (plotthis%taylorstress) then
    call writeTaylorStress()
endif
#endif
if (plotthis%viscosity) then
    call writeViscosity()
endif
if (plotthis%strain) then
    call writeStrain()
endif
if (plotthis%strainrate) then
    call writeStrainRate()
endif
if (plotthis%elempartition) then
    call writeElementPartition()
endif
if (plotthis%volume) then
    call writeElementVolume()
endif
if (plotthis%materials) then
    call writeMaterials()
endif
if (plotthis%elemQualitySide) then
    call writeElemQuality(elemQualitySide)
endif
if (plotthis%elemQualitySurface) then
    call writeElemQuality(elemQualitySurface)
endif
if (plotthis%elemQualityVolume) then
    call writeElemQuality(elemQualityVolume)
endif
write(fp,*) '         </CellData>'


! <points> block
write(fp,*) '         <Points>'
call writePointCoordsToVTU(deformed)
write(fp,*) '         </Points>'

! <cells> block
write(fp,*) '         <Cells>'
call writeConnectivityToVTU()
write(fp,*) '         </Cells>'

! write the tail of the VTU file
write(fp,*) '       </Piece>'
write(fp,*) '   </UnstructuredGrid>'
write(fp,*) '</VTKFile>'

close (fp)

end subroutine

!--------------------------------------------------------------------------------

subroutine writePrestress()

! prestress is a little bit of work, because the STN0 array does not contain prestresses for every element,
! but only for elements for which a prestress has been set. This is memory efficient
! However, the VTU file wants the prestress for every element, so that it will get a lot zeroes.
! To do this, we need the element numbers for which the prestresses have been set, which are in modeldatactx%ISELM

! Moreover, the prestress array does not have to be in numerical order, but they have
! the order that is set by the boundary condition file, which may be any.

! we create a zero array and copy the prestress value in it.

use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use modeldefinition, only: NPRE

implicit none

double precision, allocatable :: prestress(:,:)
integer   :: i,j

allocate(prestress(6, meshdatactx%neglobal))
prestress(i,j) = 0d0

do i=1,NPRE
    do j=1,6
        prestress(j, modeldatactx%ISELM(i)) = modeldatactx%STN0(j,i)
    enddo
enddo

call  writeVTUFloatBlock(fp, "PreStress Boundary condition", 28, 6, prestress, meshdatactx%neglobal, .true.)

deallocate(prestress)

end subroutine

!--------------------------------------------------------------------------------

subroutine writeStress()

use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx

implicit none

call writeVTUFloatBlock(fp, "Stress", 6, 6, modeldatactx%STN, meshdatactx%neglobal, .true.)

end subroutine

!--------------------------------------------------------------------------------

#ifdef EXPERIMENTALTRANSIENTDISCREEP
subroutine writeTaylorStress()

use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: NSTRTAY

implicit none

call  writeVTUFloatBlock(fp, "Taylorstress", 12, NSTRTAY, &
                         modeldatactx%STNTAY, meshdatactx%neglobal, .true.)

end subroutine
#endif


!--------------------------------------------------------------------------------

subroutine writeViscosity()

use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use materialsmodule, only: matdatactx, MaterialIDOfElement

implicit none

integer :: iElem, iMat
double precision, allocatable :: RawViscosity(:)

! allocate raw viscosity
allocate(RawViscosity(meshdatactx%neglobal))
! this is writing out viscosity as defined in matdatactx%PRPMAT(3,m), and depending on settings this may not be the actual viscosity
! check viscositytransfer function

do iElem=1,meshdatactx%neglobal
    iMat = MaterialIDOfElement(iElem)
    RawViscosity(iElem)  = matdatactx%PRPMAT(3,iMat)
enddo

call  writeVTUFloatBlock(fp, "Viscosity", 9, 1, RawViscosity, meshdatactx%neglobal, .true.)

deallocate(RawViscosity)

end subroutine


!--------------------------------------------------------------------------------

subroutine writeStrain()

use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx

implicit none 

call writeVTUFloatBlock(fp, "Strain", 6, 6, modeldatactx%STR, meshdatactx%neglobal, .true.)

end subroutine

!--------------------------------------------------------------------------------

subroutine writeStrainRate()

use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx

implicit none

integer :: iElem

if (plotthis%writeASCII) then
    if (plotthis%precision.eq.1) then
        write(fp,*) '            <DataArray type="Float32" NumberOfComponents="6" Name="Strainrate" Format="ascii">'
    else
        write(fp,*) '            <DataArray type="Float64" NumberOfComponents="6" Name="Strainrate" Format="ascii">'
    endif

    do iElem=1,meshdatactx%neglobal
#ifdef EXPERIMENTAL_NORATESTORAGE

        call StrainrateOfElem (meshdatactx%X,&
                               modeldatactx%D,&
                              modeldatactx%DELD, &
                               strainRate,&
                               meshdatactx%IEN,&
                              modeldatactx%LMF, &
                              modeldatactx%TFAULT,&
                              modeldatactx%DFAULT,& ! faulted node velocity
                              modeldatactx%LMX, &
                              modeldatactx%DX, &
                              modeldatactx%SKEW,&
                              modeldatactx%DXE, &
                                   iElem)

        write(fp,*) (strainRate(1:6))  ! str = strain, not stress.
#else
        write(fp,*) (modeldatactx%rate(1:6,iElem))  ! str = strain, not stress.
#endif

      enddo
else
    write(*,*) "Writing strainrate in binary not yet supported"
        ! TODO, have to find out this concatenation...
!        if (plotthis%precision.eq.1) then
!            write(fp,*) '            <DataArray type="Float32" NumberOfComponents="6" Name="Strainrate" 
!Format="binary">'
!        else
!            write(fp,*) '            <DataArray type="Float64" NumberOfComponents="6" Name="Strainrate" 
!Format="binary">'
!        endif
!        call EncodeArrayToFile(fp, modeldatactx%rate, plotthis%precision)
endif

write(fp,*) '             </DataArray>'

end subroutine

!--------------------------------------------------------------------------------

subroutine writeDisplacement()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMFNglobal, &
                           NUMSLPglobal

implicit none

integer :: iElem, iPoint, nodeID
integer :: faultEntry

double precision, allocatable :: nodalDoubleQuantity(:,:)

if (NUMSLPglobal.eq.0 .and. NUMFNglobal.eq.0) then
    call  writeVTUFloatBlock(fp, "Displacement", 12, 3, modeldatactx%D, meshdatactx%nvglobal, .true.)
else
    allocate(nodalDoubleQuantity(3, nAllNodes))
    nodalDoubleQuantity = 0d0
    call nodalQtyToAllNodesQty(modeldatactx%D, nodalDoubleQuantity)

      ! Add slippery and split displacement if necessary.
    ! This gives a sharp transition across the interface in the plot


    if (NUMSLPglobal .gt. 0) then
        ! and then add the slippery displacements. Get this from DXE
        do iElem = 1, meshdatactx%neglobal
            do iPoint = 1, 4
                nodeID = ienDuplicate(iPoint, iElem)
                nodalDoubleQuantity(:, nodeID) = nodalDoubleQuantity(:, nodeID) + &
                modeldatactx%DXE(:, iPoint, iElem)
            enddo
        enddo
    endif

    if (NUMFNglobal .gt. 0) then
        ! and then add the faulted displacements.
        ! They still need to be rotated according to Euler angles
        do iPoint = 1, nAllNodes
            faultEntry = faultEntryLookup(iPoint)
            if (faultEntry .gt. 0) then
                ! this node is faulted.
                
                nodalDoubleQuantity(:, iPoint) = nodalDoubleQuantity(:, iPoint) + &
                    modeldatactx%tfault(:, faultEntry)
            endif
        enddo
    endif


    call  writeVTUFloatBlock(fp, "Displacement", 12, 3, nodalDoubleQuantity, nAllNodes, .true.)
    deallocate(nodalDoubleQuantity)
endif

end subroutine

!--------------------------------------------------------------------------------

subroutine writeSlipperyDisplacement()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMSLPglobal

implicit none

integer :: iElem, iPoint, nodeID

double precision, allocatable :: nodalDoubleQuantity(:,:)

if (NUMSLPglobal .gt. 0) then

    allocate(nodalDoubleQuantity(3, nAllNodes))
    nodalDoubleQuantity = 0d0

    do iElem = 1, meshdatactx%neglobal
        do iPoint = 1, 4
            ! nodeID is unique thanks to the duplication
            nodeID = ienDuplicate(iPoint, iElem)
            nodalDoubleQuantity(:, nodeID) = modeldatactx%DXE(:, iPoint, iElem)
        enddo
    enddo

    call  writeVTUFloatBlock(fp, "Slippery displacement", 21, 3, nodalDoubleQuantity, nAllNodes, .true.)
    deallocate(nodalDoubleQuantity)

endif

end subroutine


!--------------------------------------------------------------------------------

subroutine writeSlipperyVelocity()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMSLPglobal

implicit none

integer :: iElem, iPoint, nodeID

double precision, allocatable :: nodalDoubleQuantity(:,:)

if (NUMSLPglobal .gt. 0) then

    allocate(nodalDoubleQuantity(3, nAllNodes))
    nodalDoubleQuantity = 0d0

    do iElem = 1, meshdatactx%neglobal
        do iPoint = 1, 4
            nodeID = ienDuplicate(iPoint, iElem)
            nodalDoubleQuantity(:, nodeID) = modeldatactx%DELDXE(:, iPoint, iElem)
        enddo
    enddo

    call  writeVTUFloatBlock(fp, "Slippery velocity", 17, 3, nodalDoubleQuantity, nAllNodes, .true.)
    deallocate(nodalDoubleQuantity)

endif

end subroutine

!--------------------------------------------------------------------------------


subroutine writeFaultedDisplacement()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMFNglobal

implicit none

integer :: iElem, iPoint, nodeID

double precision, allocatable :: nodalDoubleQuantity(:,:)
integer ::  faultEntry, originalID
double precision :: rotatedVector(3)

if (NUMFNglobal .gt. 0) then

    allocate(nodalDoubleQuantity(3, nAllNodes))
    nodalDoubleQuantity = 0d0

    ! They still need to be rotated according to Euler angles
    do iPoint = 1, nAllNodes
        faultEntry = faultEntryLookup(iPoint)
        if (faultEntry .gt. 0) then
            ! this node is faulted.

            call originalNode(iPoint, originalID)
            ! take possible Euler angles into account
            call rotateVector(modeldatactx%tfault(:, faultEntry), &
                              originalID, &
                              rotatedVector)

            nodalDoubleQuantity(:, iPoint) = nodalDoubleQuantity(:, iPoint) + &
                   rotatedVector

        endif
    enddo

    call  writeVTUFloatBlock(fp, "Faulted displacement", 20, 3, nodalDoubleQuantity, nAllNodes, .true.)
    deallocate(nodalDoubleQuantity)

endif

end subroutine

!--------------------------------------------------------------------------------


subroutine writeFaultedVelocity()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMFNglobal

implicit none

integer :: iElem, iPoint, nodeID

double precision, allocatable :: nodalDoubleQuantity(:,:)
integer ::  faultEntry, originalID
double precision :: rotatedVector(3)

if (NUMFNglobal .gt. 0) then

    allocate(nodalDoubleQuantity(3, nAllNodes))
    nodalDoubleQuantity = 0d0

    ! They still need to be rotated according to Euler angles
    do iPoint = 1, nAllNodes
        faultEntry = faultEntryLookup(iPoint)
        if (faultEntry .gt. 0) then
            ! this node is faulted.

            call originalNode(iPoint, originalID)
            ! take possible Euler angles into account
            call rotateVector(modeldatactx%dfault(:, faultEntry), &
                              originalID, &
                              rotatedVector)

            nodalDoubleQuantity(:, iPoint) = nodalDoubleQuantity(:, iPoint) + &
                   rotatedVector

        endif
    enddo

    call  writeVTUFloatBlock(fp, "Faulted velocity", 16, 3, nodalDoubleQuantity, nAllNodes, .true.)
    deallocate(nodalDoubleQuantity)

endif

end subroutine


!--------------------------------------------------------------------------------

subroutine writeVelocity()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMFNglobal, &
                           NUMSLPglobal

implicit none

integer :: iElem, iPoint, nodeID
integer :: faultEntry

double precision, allocatable :: nodalDoubleQuantity(:,:)

if (NUMSLPglobal.eq.0 .and. NUMFNglobal.eq.0) then
    call  writeVTUFloatBlock(fp, "Velocity", 8, 3, modeldatactx%DELD, meshdatactx%nvglobal, .true.)
else
    allocate(nodalDoubleQuantity(3, nAllNodes))
    nodalDoubleQuantity = 0d0

    call nodalQtyToAllNodesQty(modeldatactx%DELD, nodalDoubleQuantity)

    if (NUMSLPglobal .gt. 0) then
        ! and then add the slippery displacements. Get this from DELD
        do iElem = 1, meshdatactx%neglobal
            do iPoint = 1, 4
                nodeID = ienDuplicate(iPoint, iElem)
                nodalDoubleQuantity(:, nodeID) = nodalDoubleQuantity(:, nodeID) + &
                  modeldatactx%DELDXE(:, iPoint, iElem)
            enddo
        enddo
    endif

    if (NUMFNglobal .gt. 0) then
        ! and then add the regular displacements
        do iPoint = 1, nAllNodes
            faultEntry = faultEntryLookup(iPoint)
            if (faultEntry .gt. 0) then
                ! this node is faulted.
                nodalDoubleQuantity(:, iPoint) = nodalDoubleQuantity(:, iPoint) + &
                    modeldatactx%dfault(:, faultEntry)
            endif
        enddo
    endif


    call  writeVTUFloatBlock(fp, "Velocity", 8, 3, nodalDoubleQuantity, nAllNodes, .true.)
    deallocate(nodalDoubleQuantity)
endif

end subroutine

!--------------------------------------------------------------------------------

subroutine rotateVector(vector, nodeID, rotatedVector)
! Some input variables are vectors that are rotated over Euler Angles,
! Such as displacements and velocities and fault displacements and velocities.
! This generic function is used to rotate all those vectors

use algebra,         only: formrt, equate, vmprd
use modeltopology,   only: ndof
use modeldatamodule, only: modeldatactx

implicit none

double precision :: vector(3), rotatedVector(3)
integer          :: nodeID

double precision :: ROT(3,3), thisFault(3)

if (ABS(modeldatactx%SKEW(1,nodeID)).gt.1e-8 .or. &
    ABS(modeldatactx%SKEW(2,nodeID)).gt.1e-8) then
    ! the is rotated. Rotate fault velocities as well

    ! 1) compute rotation matrix from Euler angles
    call FORMRT (modeldatactx%SKEW(:,nodeID),ROT,NDOF)

    ! 2) use rotation matrix to rotate the vector
    call VMPRD (ROT, vector, rotatedVector ,NDOF,3)
else
    ! the node is not rotated, the vector is dandy.
    rotatedVector = vector
endif

end subroutine

!--------------------------------------------------------------------------------

subroutine writeVertexPartition()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMFNglobal, &
                           NUMSLPglobal

implicit none

integer, allocatable :: nodalIntegerQuantity(:)

if (NUMSLPglobal.eq.0 .and. NUMFNglobal.eq.0) then
    call  writeVTUIntegerBlock(fp, "Vertex Partition", 16, 1, meshdatactx%v2p, meshdatactx%nvglobal, .true.)
else
    allocate(nodalIntegerQuantity(nAllNodes))
    nodalIntegerQuantity = 0
    call nodalQtyToAllNodesQty(meshdatactx%v2p, nodalIntegerQuantity)
    call  writeVTUIntegerBlock(fp, "Vertex Partition", 16, 1, nodalIntegerQuantity, nAllNodes, .true.)
    deallocate(nodalIntegerQuantity)
endif

end subroutine

!--------------------------------------------------------------------------------

subroutine writeVertexMarker()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMFNglobal, &
                           NUMSLPglobal

implicit none

integer, allocatable :: nodalIntegerQuantity(:)

if (NUMSLPglobal.eq.0 .and. NUMFNglobal.eq.0) then
    call  writeVTUIntegerBlock(fp, "Vertex marker", 13, 1, meshdatactx%globalMarker, meshdatactx%nvglobal, .true.)
else
    allocate(nodalIntegerQuantity(nAllNodes))
    nodalIntegerQuantity = 0
    call nodalQtyToAllNodesQty(meshdatactx%globalMarker, nodalIntegerQuantity)
    call  writeVTUIntegerBlock(fp, "Vertex marker", 13, 1, nodalIntegerQuantity, nAllNodes, .true.)
    deallocate(nodalIntegerQuantity)
endif

end subroutine


!--------------------------------------------------------------------------------

subroutine writeElementPartition()

use meshdatamodule,  only: meshdatactx

implicit none

call writeVTUIntegerBlock(fp, "Element partitioning", 20, 1, meshdatactx%e2p, meshdatactx%neglobal, .true.)

end subroutine

!--------------------------------------------------------------------------------

subroutine writeElementVolume()

use meshdatamodule,  only: meshdatactx
use algebra,         only: volumeOfTetrahedron

implicit none

integer :: iElem
double precision, allocatable :: nodalDoubleQuantity(:)

double precision :: x1, y1, z1
double precision :: x2, y2, z2
double precision :: x3, y3, z3
double precision :: x4, y4, z4

allocate(nodalDoubleQuantity(meshdatactx%neglobal))
nodalDoubleQuantity = 0d0

do iElem = 1,meshdatactx%neglobal
    call getLocalElementCoordinates(iElem, .false., x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
    nodalDoubleQuantity(iElem) = VolumeOfTetrahedron(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
enddo

call  writeVTUFloatBlock(fp, "Volume", 6, 1, nodalDoubleQuantity, nAllNodes, .true.)

deallocate(nodalDoubleQuantity)

end subroutine

!--------------------------------------------------------------------------------

subroutine writeElemQuality(qualityType)

use meshdatamodule,  only: meshdatactx
use algebra,         only: volumeOfTetrahedron, &
                           areaOfTriangleIn3D, &
                           radiusOfSphereInsideTetrahedron, &
       	       	       	   radiusOfSphereAroundTetrahedron
use constants,       only: elemQualitySide, &
                           elemQualitySurface, &
                           elemQualityVolume

implicit none

integer :: qualityType

double precision, allocatable :: nodalDoubleQuantity(:)

double precision :: x1, y1, z1
double precision :: x2, y2, z2
double precision :: x3, y3, z3
double precision :: x4, y4, z4

integer :: iElem

double precision :: sideSquare(6)
double precision :: face(4)
double precision :: innerradius, outerRadius

allocate(nodalDoubleQuantity(meshdatactx%neglobal))
nodalDoubleQuantity = 0d0

do iElem = 1, meshdatactx%neglobal
    ! deformed = .false.
    call getLocalElementCoordinates(iElem, .false., x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

    if (qualityType .eq. elemQualitySide) then
        sideSquare(1) = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        sideSquare(2) = (x1-x3)**2 + (y1-y3)**2 + (z1-z3)**2
        sideSquare(3) = (x1-x4)**2 + (y1-y4)**2 + (z1-z4)**2
        sideSquare(4) = (x2-x3)**2 + (y2-y3)**2 + (z2-z3)**2
        sideSquare(5) = (x2-x4)**2 + (y2-y4)**2 + (z2-z4)**2
        sideSquare(6) = (x3-x4)**2 + (y3-y4)**2 + (z3-z4)**2
        nodalDoubleQuantity(iElem) = minval(sideSquare) / maxval(sideSquare)

    else if (qualityType .eq. elemQualitySurface) then
        face(1) = areaOfTriangleIn3D(x1, y1, z1, x2, y2, z2, x3, y3, z3)
        face(2) = areaOfTriangleIn3D(x1, y1, z1, x2, y2, z2, x4, y4, z4)
        face(3) = areaOfTriangleIn3D(x1, y1, z1, x3, y3, z3, x4, y4, z4)
        face(4) = areaOfTriangleIn3D(x2, y2, z2, x3, y3, z3, x4, y4, z4)
        nodalDoubleQuantity(iElem) = minval(face) / maxval(face)

    else if (qualityType .eq. elemQualityVolume) then
        innerradius = radiusOfSphereInsideTetrahedron(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
        outerRadius = radiusOfSphereAroundTetrahedron(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
        nodalDoubleQuantity(iElem) = 9d0 * innerradius**2 / outerRadius**2
        ! times 9 to normalizer quality to the [0,1] iterval

    endif
enddo

if      (qualityType .eq. elemQualitySide) then
    call  writeVTUFloatBlock(fp, "Side quality", 12, 1, nodalDoubleQuantity, meshdatactx%neglobal, .true.)
else if (qualityType .eq. elemQualitySurface) then
    call  writeVTUFloatBlock(fp, "Surface quality", 15, 1, nodalDoubleQuantity, meshdatactx%neglobal, .true.)
else if (qualityType .eq. elemQualityVolume) then
    call  writeVTUFloatBlock(fp, "Volume quality", 14, 1, nodalDoubleQuantity, meshdatactx%neglobal, .true.)
endif    

deallocate(nodalDoubleQuantity)

end subroutine

!--------------------------------------------------------------------------------

subroutine writeElementError()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use algebra,         only: vectorLength

implicit none

integer :: iElem
double precision, allocatable :: nodalDoubleQuantity(:)

allocate(nodalDoubleQuantity(meshdatactx%neglobal))
nodalDoubleQuantity = 0d0

allocate(modeldatactx%elemError(meshdatactx%neglobal,3))
call sumFNerror3D()
! writes error in modeldatactx%elemError

do iElem = 1, meshdatactx%neglobal
    nodalDoubleQuantity(iElem) = vectorLength(modeldatactx%elemError(iElem,:),3)
enddo
call  writeVTUFloatBlock(fp, "Element Error", 13, 1, nodalDoubleQuantity, meshdatactx%neglobal, .true.)

deallocate(nodalDoubleQuantity)

end subroutine

!--------------------------------------------------------------------------------

subroutine writeSlipperyWeights()

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition, only: NUMSLPglobal

implicit none


integer              :: iVertex, iSlipEntry
integer	       	     ::	placeInElem
integer              :: thisElement, thisVertex

integer, allocatable :: nodalIntegerQuantity(:)
allocate(nodalIntegerQuantity(nAllNodes))
nodalIntegerQuantity = 0


do iSlipEntry = 1,NUMSLPglobal

    thisElement = modeldatactx%NSLIPglobal(1,iSlipEntry)
    thisVertex  = modeldatactx%NSLIPglobal(2,iSlipEntry)

    ! find the vertex place in the original IEN, without the duplicates,
    ! because the NSLIPglobal array does not have duplicates.
    do iVertex = 1,4
        if (meshdatactx%ien(iVertex,thisElement) .eq. thisVertex) then
            placeInElem = iVertex
        endif
    enddo

    ! Now get the duplicated vertex ID corresponding to this one.
    ! If the note has no duplicates, thevertex ID will remain the same.
    thisVertex = ienDuplicate(placeInElem,thisElement)

    nodalIntegerQuantity(thisVertex) = modeldatactx%NSLIPglobal(3,iSlipEntry)

enddo

call writeVTUIntegerBlock(fp, "Slippery weights", 16, 1, nodalIntegerQuantity, &
                          nAllNodes, .true.)

deallocate(nodalIntegerQuantity)

end subroutine

!--------------------------------------------------------------------------------

subroutine writeMaterials()

use meshdatamodule,  only: meshdatactx

implicit none

call writeVTUIntegerBlock(fp, "Materials", 9, 1, meshdatactx%mat, meshdatactx%neglobal, .false.)

end subroutine

!--------------------------------------------------------------------------------

subroutine writeNodalBC(BCType)

use modeldefinition, only: NUMSLPglobal, &
                           NUMFNglobal
use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use constants,       only: veloBC, &
                           dispBC, &
                           forceBC, &
                           slideInPlaneBC, &
                           slideOnLineBC, &
                           fixedPointBC, &
                           EulerAnglesBC
use algebra,         only: rotatePoint

implicit none

integer                       :: BCType

double precision, allocatable :: BCarray(:,:)
integer                       :: iNode, iDof
double precision              :: rotated(3)


! we must make an array containing the boundary conditions
! me do not need to take duplicate nodes into account,
! because nodal bcs do not work on slippery nodes.
! we must adjust the size of the array, though
if (NUMSLPglobal.gt.0 .or. NUMFNglobal.gt.0) then
    allocate(BCarray(3,nAllNodes))
else
    allocate(BCarray(3,meshdatactx%nvglobal))
endif

BCarray = 0d0

if (BCType.eq.veloBC.or.BCType.eq.dispBC.or.BCType.eq.forceBC) then
    do iNode=1,meshdatactx%nvglobal
        do iDof =1,3
            if(modeldatactx%IBONDglobal(iDof,iNode).eq.BCType) then
                BCarray(iDof,iNode) = modeldatactx%BOND(iDof,iNode)
            endif
        enddo
    enddo
else if (BCType.eq.slideInPlaneBC) then
    do iNode=1,meshdatactx%nvglobal
        if      (modeldatactx%IBONDglobal(1,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(1,iNode).eq.0             .and.&
                 modeldatactx%IBONDglobal(2,iNode).eq.0     .and.&
                 modeldatactx%IBONDglobal(3,iNode).eq.0)     then
                 ! normal to slide plane in x direction
                 BCarray(1,inode) = scaleFactor
        else if (modeldatactx%IBONDglobal(1,iNode).eq.0     .and.&
                 modeldatactx%IBONDglobal(2,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(2,iNode).eq.0             .and.&
                 modeldatactx%IBONDglobal(3,iNode).eq.0)     then
                 ! normal to slide plane in y direction
                 BCarray(2,inode) = scaleFactor
        else if (modeldatactx%IBONDglobal(1,iNode).eq.0     .and.&
                 modeldatactx%IBONDglobal(2,iNode).eq.0     .and.&
                 modeldatactx%IBONDglobal(3,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(3,iNode).eq.0)            then
                 ! normal to slide plane in z direction
                 BCarray(3,inode) = scaleFactor
        else
                 ! this node is not a plane slider. Do nothing
        endif
    enddo
else if (BCType.eq.slideOnLineBC) then
    do iNode=1,meshdatactx%nvglobal
        if      (modeldatactx%IBONDglobal(1,iNode).eq.0     .and.&
                 modeldatactx%IBONDglobal(2,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(2,iNode).eq.0             .and.&
                 modeldatactx%IBONDglobal(3,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(3,iNode).eq.0)            then
                 ! slide on line in x direction
                 BCarray(1,inode) = scaleFactor
        else if (modeldatactx%IBONDglobal(1,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(1,iNode).eq.0             .and.&
                 modeldatactx%IBONDglobal(2,iNode).eq.0      .and.&
                 modeldatactx%IBONDglobal(3,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(3,iNode).eq.0)            then
                 ! slide on line in y direction
                 BCarray(2,inode) = scaleFactor
        else if (modeldatactx%IBONDglobal(1,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(1,iNode).eq.0             .and.&
                 modeldatactx%IBONDglobal(2,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(2,iNode).eq.0            .and.&
                 modeldatactx%IBONDglobal(3,iNode).eq.0)     then
                 ! slide on line in z direction
                 BCarray(3,inode) = scaleFactor
        else 
                 ! this node is not a 1D roller. Do nothing
        endif
    enddo
else if (BCType.eq.fixedPointBC) then
    do iNode=1,meshdatactx%nvglobal
        if      (modeldatactx%IBONDglobal(1,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(1,iNode).eq.0             .and.&
                 modeldatactx%IBONDglobal(2,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(2,iNode).eq.0            .and.&
                 modeldatactx%IBONDglobal(3,iNode).eq.dispBC .and.&
                 modeldatactx%BOND(3,iNode).eq.0)            then
                 ! Where are you going? Nowhere! That's Right... you're going nowhere.
                 BCarray(1,inode) = scaleFactor
                 BCarray(2,inode) = scaleFactor
                 BCarray(3,inode) = scaleFactor
        endif
    enddo
else if (BCType.eq.EulerAnglesBC) then
    do iNode=1,meshdatactx%nvglobal
        if      (modeldatactx%SKEW(1,iNode).ne.0d0 .or. &
                 modeldatactx%SKEW(2,iNode).ne.0d0) then
                  BCarray(1,inode) = scaleFactor
        else
                ! no rotation; vector remains 0
        endif
    enddo
else
    ! no BC, do nothing
endif

! Rotate all the arrows in BCarray if they are subjected to Euler angles        
do iNode =1,meshdatactx%nvglobal
               if(modeldatactx%SKEW(1,iNode).ne.0d0 .or. &
                  modeldatactx%SKEW(2,iNode).ne.0d0) then
                call rotatePoint(modeldatactx%SKEW(:,iNode), 3, BCarray(:,inode), rotated)
!               write(*,*) "rotating veloBC", iNode, "from", BCarray(:,inode), "to", rotated 
                BCarray(:,inode) = rotated
               endif
enddo


if (BCType.eq.EulerAnglesBC) then
    call  writeVTUFloatBlock(fp, "Euler Angles", 12, 3, BCarray, nAllNodes, .true.)
endif


if (NUMSLPglobal.gt.0 .or. NUMFNglobal.gt.0) then
    if (BCType.eq.veloBC) then
        call  writeVTUFloatBlock(fp, "Velocity BC", 11, 3, BCarray, nAllNodes, .true.)
    else if (BCType.eq.dispBC) then
        call  writeVTUFloatBlock(fp, "Displacement BC", 15, 3, BCarray, nAllNodes, .true.)
    else if (BCType.eq.forceBC) then
        call  writeVTUFloatBlock(fp, "Force BC", 8, 3, BCarray, nAllNodes, .true.)
    else if (BCType.eq.fixedPointBC) then
        call  writeVTUFloatBlock(fp, "Fixed points", 12, 3, BCarray, nAllNodes, .true.)
    else if (BCType.eq.slideOnLineBC) then
        call  writeVTUFloatBlock(fp, "Slide on line", 13, 3, BCarray, nAllNodes, .true.)
    else if (BCType.eq.slideInPlaneBC) then
        call  writeVTUFloatBlock(fp, "Slide in plane", 14, 3, BCarray, nAllNodes, .true.)
    else if (BCType.eq.EulerAnglesBC) then
        ! do nothing, this has been handled above
    else
         write(*,*) "did not recognize BC type", BCType
         stop "Contact model support. Leaving Plt3d"
    endif
else
    if (BCType.eq.veloBC) then
        call  writeVTUFloatBlock(fp, "Velocity BC", 11, 3, BCarray, meshdatactx%nvglobal, .true.)
    else if (BCType.eq.dispBC) then
        call  writeVTUFloatBlock(fp, "Displacement BC", 15, 3, BCarray, meshdatactx%nvglobal, .true.)
    else if (BCType.eq.forceBC) then
        call  writeVTUFloatBlock(fp, "Force BC", 8, 3, BCarray, meshdatactx%nvglobal, .true.)
    else if (BCType.eq.fixedPointBC) then
        call  writeVTUFloatBlock(fp, "Fixed points", 12, 3, BCarray, meshdatactx%nvglobal, .true.)
    else if (BCType.eq.slideOnLineBC) then
        call  writeVTUFloatBlock(fp, "Slide on line", 13, 3, BCarray, meshdatactx%nvglobal, .true.)
    else if (BCType.eq.slideInPlaneBC) then
        call  writeVTUFloatBlock(fp, "Slide in plane", 14, 3, BCarray, meshdatactx%nvglobal, .true.)
    else if (BCType.eq.EulerAnglesBC) then
        ! do nothing, this has been handled above
    else
        write(*,*) "did not recognize BC type", BCType
        stop "Contact model support. Leaving Plt3d"
    endif

endif

deallocate(BCarray)

end subroutine




!--------------------------------------------------------------------------------

subroutine writeConnectivityToVTU()

use meshdatamodule,  only: meshdatactx
use modeldefinition, only: NUMFNglobal, &
                           NUMSLPglobal
use base64Encode,    only: EncodeArrayToFile

implicit none

integer :: i
integer :: nNodes

if (meshAdjusted) then
    ! the connectivity has been adapted into ienDuplicate
    call writeVTUIntegerBlock(fp,"connectivity",12,4, ienDuplicate-1, meshdatactx%neglobal,.false.)
else
    ! the connectivity is in the original ien array
    call writeVTUIntegerBlock(fp,"connectivity",12,4, meshdatactx%ien-1, meshdatactx%neglobal,.false.)
endif

! offsets are relatively small, no need for double precision encoding
if (plotthis%writeASCII) then
    write(fp,*) '            <DataArray type="UInt32" Name="offsets" Format="ascii">'

    do i=1,meshdatactx%neglobal
        write(fp,*) 4*i
    enddo
else
    write(fp,*) '            <DataArray type="UInt32" Name="offsets" Format="binary">'

    allocate(offsets(meshdatactx%neglobal))
    do i=1,meshdatactx%neglobal
        offsets(i) = 4*i
    enddo
    call EncodeArrayToFile(fp, meshdatactx%mat, plotthis%precision)
endif

write(fp,*) '            </DataArray>'

! write the types of the elements. In this case they are all type 10, for tetrahedrons
! because the numbers stay small, it is hardly beneficial to write this in binary
! Hence, we will keep writing this always in ascii 
write(fp,*) '            <DataArray type="Int32" Name="types" Format="ascii">'
do i=1,meshdatactx%neglobal-1
! 10 is paraviewinese for 'tetrahedron'.
    write(fp,'(a3)',advance="no") '10 ' ! write them all on one line to keep it compact
enddo
write(fp,'(a2)') '10'

write(fp,*) '            </DataArray>'

end subroutine

!--------------------------------------------------------------------------------

subroutine writePointCoordsToVTU(deformed)

use modeldefinition, only: NUMFNglobal, &
                           NUMSLPglobal
use  meshdatamodule, only:  meshdatactx
use modeldatamodule, only: modeldatactx
use postProcessData, only: deffac ! deformation factor.
use iomodule,        only: stderr

implicit none

logical  :: deformed

integer  :: iPoint, nodeID, iElem
integer  :: originalElem
integer  :: faultEntry

double precision, allocatable :: newDisplacedCoords(:,:)
double precision              :: rotatedVector(3)

integer :: originalID

if (deformed) then
    ! a plot of the deformed mesh is desired.
    if (NUMSLPglobal.eq.0 .and. NUMFNglobal.eq.0) then
        ! If the mesh does not have any slip of fault, this is fairly simple:  
        call  writeVTUFloatBlock(fp, "Vertex coordinates", 18, 3, meshdatactx%X + modeldatactx%D, nAllNodes, .true.)
    else
        ! if the there are slippery of faulted interfaces, the new ienDuplicate
        ! does not correspond to point numbers in modeldatactx%D,
        ! and the lookup table must be used to ensure that the proper value goes to the proper place.
        allocate(newDisplacedCoords(3, nAllNodes))

        ! first fill with the original coordinates (including the duplicates)
        newDisplacedCoords = allNodeCoords 

        ! add the regular displacements from the original node in modeldatactx%D
        do iPoint = 1, nAllNodes

             call originalNode(iPoint, originalID)

             newDisplacedCoords(:, iPoint) = newDisplacedCoords(:, iPoint) + &
               modeldatactx%D(:, originalID)
        enddo

        if (NUMSLPglobal .gt. 0) then
            ! and then add the slippery displacements. Get this from DXE

            do iElem = 1, meshdatactx%neglobal
                do iPoint = 1, 4
                    nodeID = ienDuplicate(iPoint, iElem)
                    newDisplacedCoords(:, nodeID) = newDisplacedCoords(:, nodeID) + &
                     modeldatactx%DXE(:, iPoint, iElem)
                enddo
            enddo
        endif

        if (NUMFNglobal .gt. 0) then
            ! and then add the regular displacements
            do iPoint = 1, nAllNodes
                faultEntry = faultEntryLookup(iPoint)
                if (faultEntry .gt. 0) then
                    ! this node is faulted.

                        call originalNode(iPoint, originalID)

                    ! take possible Euler angles into account
                    call rotateVector(modeldatactx%tfault(:, faultEntry), &
                                      originalID, &
                                      rotatedVector)


                    newDisplacedCoords(:, iPoint) = newDisplacedCoords(:, iPoint) + &
                        rotatedVector
                endif
            enddo
        endif

        call  writeVTUFloatBlock(fp, "Vertex coordinates", 18, 3, newDisplacedCoords, nAllNodes, .true.)
        
        deallocate(newDisplacedCoords)
    endif

else
    if (NUMSLPglobal.gt.0 .or. NUMFNglobal.gt.0) then
        call  writeVTUFloatBlock(fp, "Vertex coordinates", 18, 3, allNodeCoords, nAllNodes, .true.)
    else
        call  writeVTUFloatBlock(fp, "Vertex coordinates", 18, 3, meshdatactx%X, meshdatactx%nvglobal, .true.)
    endif

endif

end subroutine

!--------------------------------------------------------------------------------

subroutine writeVTUIntegerBlock(fp, name, nNameChars, nComponents, dataArray, nEntries, mentionNrComponents)

use constants,    only: int, float, single, double
use base64Encode, only: EncodeArrayToFile

implicit none

!-pass
integer :: fp        ! pointer to the VTU file
character(len=*) :: name ! name of the variable. Will appear in paraview list
integer :: nNameChars  ! length of name
integer :: nComponents ! 1 for scalar; 3 for vectors; 4 for connectivity; 6 for tensors
integer :: nEntries    ! can be equal to:
                       ! a) numnp, for nodal quantities with no slippery of faulted nodes
                       ! b) nElems * 4, for nodal quantities in models slippery of faulted nodes, or for deformed meshes
                       !      (should be changed into only separate for slippery nodes)
                       ! c) nElems, for element based quantities
integer :: dataArray (nComponents, nEntries)
logical :: mentionNrComponents ! some entries in the VTU file, such as connectivity, require the number of components to not be mentioned.
                               ! however, the subroutine still needs it in order to print the data correctly. Hence an extra flag

!-local
integer :: iComponent, iEntry

!call writeVTUBlockHeader(fp, name, nNameChars, nComponents, mentionNrComponents)

character(len=1)  :: nComponentsString
character(len=10) :: varTypeString
character(len=10) :: formatString

write (nComponentsString, "(i1)") nComponents

if   (plotthis%precision .eq. single) then
    varTypeString="UInt32"
else
    varTypeString="Uint64"
endif

if (plotthis%writeASCII) then
    formatString = "ascii"
else
    formatString = "binary"
endif

! ***** now write the actual block ******

if (mentionNrComponents) then
    write(fp,*) '            <DataArray type="'//trim(varTypeString)// \
                            '" NumberOfComponents="'//nComponentsString// \
                            '" Name="'// trim(name(1:nNameChars))// \
                            '" Format="' // trim(formatString) // '">'
else
    write(fp,*) '            <DataArray type="'//trim(varTypeString)// \
                            '" Name="'// trim(name(1:nNameChars))// \
                            '" Format="' // trim(formatString) // '">'
endif

! content
if (plotthis%writeASCII) then
    do iEntry=1,nEntries
        do iComponent=1,nComponents-1
            write(fp,'(i20,a1)', advance="no") dataArray(iComponent, iEntry), " "
        enddo
        ! write the last component with a line break
        write(fp,*) dataArray(nComponents, iEntry)
!        write(fp,*) dataArray(iEntry)
    enddo
else
    call EncodeArrayToFile(fp, dataArray, plotthis%precision)
endif

! footer
write(fp,*) '            </DataArray>'


end subroutine

!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------

subroutine writeVTUFloatBlock(fp, name, nNameChars, nComponents, dataArray, nEntries, mentionNrComponents)

use constants,       only: int, float, single, double
use base64Encode,    only: EncodeArrayToFile
use meshdatamodule,  only: meshdatactx

implicit none

!-pass
integer :: fp        ! pointer to the VTU file
character(len=*) :: name ! name of the variable. Will appear in paraview list
integer :: nNameChars  ! length of name
integer :: nComponents ! 1 for scalar; 3 for vectors; 6 for tensors
integer :: nEntries    ! can be equal to:
                       ! a) numnp, for nodal quantities with no slippery of faulted nodes
                       ! b) nElems * 4, for nodal quantities in models slippery of faulted nodes, or for deformed meshes
                       !      (should be changed into only separate for slippery nodes)
                       ! c) nElems, for element based quantities
double precision :: dataArray (nComponents, nEntries)
logical :: mentionNrComponents ! some entries in the VTU file, such as connectivity, require the number of components to not be mentioned.
                               ! however, the subroutine still needs it in order to print the data correctly. Hence an extra flag
logical :: hasDuplicateNodes   ! has slippery nodes and must handle this by selecting the proper ones

!-local
integer :: iComponent, iEntry

character(len=1)  :: nComponentsString
character(len=10) :: varTypeString
character(len=10) :: formatString

write (nComponentsString, "(i1)") nComponents

if   (plotthis%precision .eq. single) then
    varTypeString="Float32"
else
    varTypeString="Float64"
endif

if (plotthis%writeASCII) then
    formatString = "ascii"
else
    formatString = "binary"
endif

! ***** now write the actual block ******

if (mentionNrComponents) then
    write(fp,*) '            <DataArray type="'//trim(varTypeString)// \
                            '" NumberOfComponents="'//nComponentsString// \
                            '" Name="'// trim(name(1:nNameChars))// \
                            '" Format="' // trim(formatString) // '">'
else
    write(fp,*) '            <DataArray type="'//trim(varTypeString)// \
                            '" Name="'// trim(name(1:nNameChars))// \
                            '" Format="' // trim(formatString) // '">'
endif


! content
if (plotthis%writeASCII) then
    do iEntry=1,nEntries  ! usually the number of nodes/elements/thingies

        do iComponent=1,nComponents-1  ! usually iDOF, or number of tensor components or so
           write(fp,'(e20.8,a1)', advance="no") dataArray(iComponent, iEntry), " "
        enddo
        ! write the last component with a line break
        write(fp,*) dataArray(nComponents, iEntry)
    enddo
else
    call EncodeArrayToFile(fp, dataArray, plotthis%precision)
endif

! footer
write(fp,*) '            </DataArray>'

end subroutine

!--------------------------------------------------------------------

subroutine nodalDoubleQtyToAllNodesQty(QtyOnNodes, QtyOnAllNodes)

use  meshdatamodule, only:  meshdatactx

implicit none

double precision :: QtyOnNodes(3, meshdatactx%nvglobal)
double precision :: QtyOnAllNodes(3, nAllNodes)

integer :: iPoint


! copy all values from original nodes directly
do iPoint = 1, nAllNodes           
    if (iPoint .le. meshdatactx%nvglobal) then
        QtyOnAllNodes(1:3,iPoint) = QtyOnNodes(1:3,iPoint)
    else
        QtyOnAllNodes(1:3,iPoint) = &
          QtyOnNodes(1:3,referenceNode(iPoint-meshdatactx%nvglobal))
    endif
enddo

end subroutine

!--------------------------------------------------------------------


subroutine nodalIntegerQtyToAllNodesQty(QtyOnNodes, QtyOnAllNodes)

use  meshdatamodule, only:  meshdatactx

implicit none

integer :: QtyOnNodes(meshdatactx%nvglobal)
integer :: QtyOnAllNodes(nAllNodes)

integer :: iPoint

! copy all values from original nodes directly
do iPoint = 1, nAllNodes
    if (iPoint .le. meshdatactx%nvglobal) then
        QtyOnAllNodes(iPoint) = QtyOnNodes(iPoint)
    else
        QtyOnAllNodes(iPoint) = QtyOnNodes(referenceNode(iPoint-meshdatactx%nvglobal))
    endif
enddo

end subroutine

!--------------------------------------------------------------------


subroutine adjustMeshForSlipAndFaults()

    ! mesh needs adaptation
    ! 1: check how much extra room we need, to allocate arrays.
    !    How many extra nodes will be created.
    !    a) allocate an array to count the the slippery nodes
    !    b) walk to the slippery and through the faulted arrays add all slippery combinations
    !    c) reallocate the nodal point array

    ! 2: add coordinates and modify connectivity, by walking over the elements with slippery nodes.
    !    At the same time, add new coordinates for the new nodal points to the nodal point array
    !     Because the number of elements does not change, the connectivity array does not need to be resized.
    !     (MATT: BUT IT DOES NEED TO BE REDEFINED TO HAVE NODES MAKING UP ELEMENTS CONSIST OF DUPLICATED NODES)
    !     (Lukas: Yes, and that ius what this code below is all bout)

    ! 3: The array for new nodal point quantities (displacement and velocity) will
    !    be created when they are needed, because they are not necessarily needed, and they are
    !    not needed at the same time. We do not want both to take up memory at the same time, in case they are very large
    !    We will use the same 'scratch' array for this, which is passed along to the VTU blockwriter

use  meshdatamodule, only:  meshdatactx
use modeldatamodule, only: modeldatactx
use     debugmodule, only: iecho, debug
use modeldefinition, only: numfnglobal, NUMSLPglobal

implicit none

integer, allocatable :: nMultiNodeLookup(:)
integer              :: nodeOfEntry, thisNode
integer              :: offset
integer, allocatable :: firstDuplicateOffset(:)

integer              :: iNode, iElem, iFaultEntry, iEntry, iDOF, elemNR, newNodeID, nodeNR
integer              :: iTime
integer              :: sequenceID


    ! 1a

    ! nMultiNodeLookup: the number of extra nodes that are added at each node
    allocate(nMultiNodeLookup(meshdatactx%nvglobal))
    ! firstReferenceElem: the first element connected to each node
    allocate(firstReferenceElem(meshdatactx%nvglobal))
    nMultiNodeLookup = 0
    firstReferenceElem = 0

    ! 1b

    ! walk through the slippery array
    if (NUMSLPglobal.gt.0) then
        do iEntry=1,NUMSLPglobal
            ! nodeOfEntry: the node number associated with this node-element pair
            nodeOfEntry = modeldatactx%NSLIPglobal(2,iEntry)
            if (nMultiNodeLookup(nodeOfEntry).eq.0) then
!                write(*,*) "first slip entry for node", nodeOfEntry, &
!                           "with elem", modeldatactx%NSLIPglobal(1,iEntry)
                firstReferenceElem(nodeOfEntry) = modeldatactx%NSLIPglobal(1,iEntry)
           endif
            nMultiNodeLookup(nodeOfEntry) = \
            nMultiNodeLookup(nodeOfEntry) + 1
        enddo
    endif

    ! and through the faulted array
    if (NUMFNglobal.gt.0) then

        if (debug) write(*,*) "making space for faulted nodes", NUMFNglobal

        do iEntry=1,NUMFNglobal

            ! nodeOfEntry: the node number associated with this node-element pair
            nodeOfEntry = modeldatactx%NFAULTglobal(2,iEntry)

            if (debug) write(*,*) "faulted node entry ", iEntry, "is node", nodeOfEntry

            if (nMultiNodeLookup(nodeOfEntry).eq.0) then
!                write(*,*) "first faulted entry for node", nodeOfEntry, &
!                           "with elem", modeldatactx%NSLIPglobal(1,iEntry)
                firstReferenceElem(nodeOfEntry) = modeldatactx%NFAULTglobal(1,iEntry)
            endif

            nMultiNodeLookup(nodeOfEntry) = \
            nMultiNodeLookup(nodeOfEntry) + 1
        enddo
    endif

    ! of all the multiplied nodes, remove one to correct for the original node
    do iNode=1,meshdatactx%nvglobal
        if (nMultiNodeLookup(iNode).gt.0) then
            nMultiNodeLookup(iNode) = \
            nMultiNodeLookup(iNode) - 1
        endif
    enddo

    ! sum to get the extra nodes
    ! nExtraNodes: total number of extra nodes
    nExtraNodes = sum(nMultiNodeLookup)

    ! nAllNodes: total number of original and extra nodes
    nAllNodes = nExtraNodes + meshdatactx%nvglobal

    if (iecho.eq.2) then
        write(*,*) "nNodes:", meshdatactx%nvglobal
        write(*,*) "added:", nExtraNodes
        write(*,*) "new total:", nAllNodes
    endif


    ! 1c, allocate data for the new amout of nodes, nAllNodes

    if (.not. allocated(allNodeCoords)) then
        allocate(allNodeCoords(3,nAllNodes))
        allNodeCoords = 0d0
    endif
    ! referenceNode: the original node at each of the extra nodes
    if (.not. allocated(referenceNode)) then
        allocate(referenceNode(nExtraNodes))
        referenceNode = 0
    endif
    ! referenceElem: the element for each (duplicated) node
    if (.not. allocated(referenceElem)) then
        allocate(referenceElem(nAllNodes))
        referenceElem = 0
    endif

    ! copy all the ordinary nodes
    do iNode=1,meshdatactx%nvglobal
        do iDOF = 1,3
            allNodeCoords(iDOF,iNode) = meshdatactx%X(iDOF,iNode)
        enddo
    enddo

    ! 2

    ! and add the coordinates of the duplicate nodes
    ! firstDuplicateOffset: location of the first duplicate in the total node array
    allocate(firstDuplicateOffset(meshdatactx%nvglobal))
    firstDuplicateOffset = 0
    offset = meshdatactx%nvglobal
    do iNode=1,meshdatactx%nvglobal
        if (nMultiNodeLookup(iNode).gt.0) then
            ! this node occurd multiple times
            firstDuplicateOffset(iNode) = offset+1
            do iTime=1,nMultiNodeLookup(iNode)
                offset = offset + 1
                referenceNode(offset- meshdatactx%nvglobal) = iNode ! lookup table to be used elsewhere
                do iDOF = 1,3
                    allNodeCoords(iDOF,offset) = meshdatactx%X(iDOF,iNode)
                enddo
            enddo
        endif
    enddo

    ! Define a new IEN array, with points corresponding to modified mesh instead of original mesh
    ! Initialize it with values from the original IEN array and modify when a duplicated node is found.
    ! Note that while nodes are being added, the number of elements remains the same,
    ! so the ienDuplicate array can also be the same size as the original IEN array.
    if (.not.allocated(ienDuplicate)) then
        write(*,*) "allocating ienDuplicate", 4,meshdatactx%neglobal
        allocate(ienDuplicate(4,meshdatactx%neglobal))
    endif
    ienDuplicate = meshdatactx%ien

    ! We need a mechanism to avoid changing the original node as well.
    ! To avoid needing extra memory, we add a minus sign.
    ! After the first occurence, we remove the minus sign and continue replacing the node by its duplicates,
    ! leaving the original intact
    firstDuplicateOffset = -firstDuplicateOffset

    ! Modify the connectivity by splitting the slippery and faulted nodes into separate nodes.
    ! Walk through the array, and see if we get a slippery/split node.
    do iElem = 1,meshdatactx%neglobal
        do iNode = 1,4
            thisNode = meshdatactx%ien(iNode, iElem)
            if (firstDuplicateOffset(thisNode).ne.0 ) then
                ! Yes, this is a slippery/split node

                if (firstDuplicateOffset(thisNode).lt.0) then
                    ! First occurence. Remove minus sign and do nothing else.
                    ! Already been considered above when defining firstReferenceElem
                    ! Added to referenceElem below
                    firstDuplicateOffset(thisNode) = -firstDuplicateOffset(thisNode)
                    referenceElem(thisNode) = iElem
                else
                    ! It is a slippery or faulted node

                    ! Enter the duplicates in the reference element list, 
                    ! to determine which differential disp/velo
                    ! belongs to this particular duplicate point.
                    referenceElem(firstDuplicateOffset(thisNode)) = iElem

                    ! Update the connectivity
                    ienDuplicate(iNode,iElem) = firstDuplicateOffset(thisNode)
                    ! Matt says: I don't think we should mess with this array;
                    ! seems pretty fundamental
!                    meshdatactx%ien(iNode, iElem) = firstDuplicateOffset(thisNode)
                    ! Lukas replies: You are probably right. It is safe when the writing 
                    ! of the VTU file is the last thing that happens during post processing, 
                    ! which is typically the case, but not necessarily so!
                    ! I wrote the new connectivity to the same IEN array
                    ! to save memory, as the IEN array can become large for large meshes,
                    ! but we can always change it when we run into the limit,
                    ! and when 64 GB is becoming standard, that limit is far away.
                    ! (I imagine the programmers of the year 2050 reading this comment,
                    ! and thinking: 'whoa, those poor bastards had only 64 GB? How could
                    ! they do *anything*?' )

                    firstDuplicateOffset(thisNode) = \
                    firstDuplicateOffset(thisNode) + 1

                endif
            endif
        enddo
   enddo

    ! using this update in the connectivity we can create a new nfault and tfault array,
    ! that contains the new node IDs, so there should not be any duplicate nodes in there.
    ! We then have a unique node number linked to each differential faulted displacement.

    allocate(faultEntryLookup(nAllNodes))
    faultEntryLookup = 0

    do iFaultEntry = 1, NUMFNglobal

        elemNR = modeldatactx%NFAULTglobal(1, iFaultEntry)
        nodeNR = modeldatactx%NFAULTglobal(2, iFaultEntry)

        ! check in IEN the position of this node in the element
        sequenceID = 0
        do iNode=1,4
            if (meshdatactx%ien(iNode, elemNr) .eq. nodeNR) then
                sequenceID = iNode
            endif
        enddo

        newNodeID = ienDuplicate(sequenceID,elemNr)

        faultEntryLookup(newNodeID) = iFaultEntry

    enddo


    ! Now the connectivity and the coordinates have been adapted.
    ! if we need to write a quantity to file, we wish the VTU block writer to do this whenever needed.

    ! add the first reference Elements, so that we may remove that array
    !do iNode=1,meshdatactx%nvglobal
    !    referenceElem(iNode) = referenceElem(iNode) + firstReferenceElem(iNode)
    !enddo
    deallocate(firstReferenceElem)

end subroutine

!-----------------------------------------------------------------------------

subroutine originalNode(newNode, originalID)

use  meshdatamodule, only:  meshdatactx

implicit none

integer :: newNode, originalID

! this gets the old node ID back from the list of new duplicated nodes

originalID = newNode

if (newNode .gt. meshdatactx%nvglobal) then
    originalID = referenceNode(newNode - meshdatactx%nvglobal)
endif

end subroutine

!-----------------------------------------------------------------------------

subroutine getLocalElementCoordinates(iElem, deformed, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

use meshdatamodule,  only: meshdatactx
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: modeldatactx
#endif
use modeldefinition, only: NUMSLPglobal

implicit none

!-pass
integer,          intent(in)  :: iElem
logical,          intent(in)  :: deformed
double precision, intent(out) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4

!-local
integer :: p1, p2, p3, p4
integer :: slipperySequenceNr

p1 = meshdatactx%ien(1,iElem)
p2 = meshdatactx%ien(2,iElem)
p3 = meshdatactx%ien(3,iElem)
p4 = meshdatactx%ien(4,iElem)

x1 = meshdatactx%x(1,p1)
y1 = meshdatactx%x(2,p1)
z1 = meshdatactx%x(3,p1)

x2 = meshdatactx%x(1,p2)
y2 = meshdatactx%x(2,p2)
z2 = meshdatactx%x(3,p2)

x3 = meshdatactx%x(1,p3)
y3 = meshdatactx%x(2,p3)
z3 = meshdatactx%x(3,p3)

x4 = meshdatactx%x(1,p4)
y4 = meshdatactx%x(2,p4)
z4 = meshdatactx%x(3,p4)


if (deformed) then

    ! add regular displacement
    x1 = meshdatactx%X(1,p1) + modeldatactx%D(1,p1)
    y1 = meshdatactx%X(2,p1) + modeldatactx%D(2,p1)
    z1 = meshdatactx%X(3,p1) + modeldatactx%D(3,p1)

    x2 = meshdatactx%X(1,p2) + modeldatactx%D(1,p2)
    y2 = meshdatactx%X(2,p2) + modeldatactx%D(2,p2)
    z2 = meshdatactx%X(3,p2) + modeldatactx%D(3,p2)

    x3 = meshdatactx%X(1,p3) + modeldatactx%D(1,p3)
    y3 = meshdatactx%X(2,p3) + modeldatactx%D(2,p3)
    z3 = meshdatactx%X(3,p3) + modeldatactx%D(3,p3)

    x4 = meshdatactx%X(1,p4) + modeldatactx%D(1,p4)
    y4 = meshdatactx%X(2,p4) + modeldatactx%D(2,p4)
    z4 = meshdatactx%X(3,p4) + modeldatactx%D(3,p4)

    ! and add slippery displacement

    if (NUMSLPglobal.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc

        slipperySequenceNr = elementSlipperyID(iElem,142)

        if(slipperySequenceNr.gt.0) then
            ! yes, this element has slippery nodes
            x1 = x1 + modeldatactx%DXE(1, 1, slipperySequenceNr)
            y1 = y1 + modeldatactx%DXE(2, 1, slipperySequenceNr)
            z1 = z1 + modeldatactx%DXE(3, 1, slipperySequenceNr)
            x2 = x2 + modeldatactx%DXE(1, 2, slipperySequenceNr)
            y2 = y2 + modeldatactx%DXE(2, 2, slipperySequenceNr)
            z2 = z2 + modeldatactx%DXE(3, 2, slipperySequenceNr)

            x3 = x3 + modeldatactx%DXE(1, 3, slipperySequenceNr)
            y3 = y3 + modeldatactx%DXE(2, 3, slipperySequenceNr)
            z3 = z3 + modeldatactx%DXE(3, 3, slipperySequenceNr)

            x4 = x4 + modeldatactx%DXE(1, 4, slipperySequenceNr)
            y4 = y4 + modeldatactx%DXE(2, 4, slipperySequenceNr)
            z4 = z4 + modeldatactx%DXE(3, 4, slipperySequenceNr)
       endif

#else
       x1 = x1 + modeldatactx%DXE(1, 1, iElem)
       y1 = y1 + modeldatactx%DXE(2, 1, iElem)
       z1 = z1 + modeldatactx%DXE(3, 1, iElem)

       x2 = x2 + modeldatactx%DXE(1, 2, iElem)
       y2 = y2 + modeldatactx%DXE(2, 2, iElem)
       z2 = z2 + modeldatactx%DXE(3, 2, iElem)

       x3 = x3 + modeldatactx%DXE(1, 3, iElem)
       y3 = y3 + modeldatactx%DXE(2, 3, iElem)
       z3 = z3 + modeldatactx%DXE(3, 3, iElem)

       x4 = x4 + modeldatactx%DXE(1, 4, iElem)
       y4 = y4 + modeldatactx%DXE(2, 4, iElem)
       z4 = z4 + modeldatactx%DXE(3, 4, iElem)

#endif
    endif
endif

end subroutine

!-----------------------------------------------------------------------------

subroutine VTKInititialisePlotoptions

use timestepmodule, only: timestepdatactx, &
                          nintg ! number of time step groups

implicit none
!    type(plotoptions plotthis

plotthis%prestress       = .false.
plotthis%materials       = .false.
plotthis%slipperyMarkers = .false.
plotthis%slipperyVectors = .false.
plotthis%slipperyWeights = .false.
plotthis%faultedDisplacement = .false.
plotthis%velocity        = .true.
plotthis%velocityBC      = .false.
plotthis%displacement    = .true.
plotthis%displacementBC  = .false.
plotthis%slipperyDisplacement = .false.
plotthis%forceBC         = .false.
plotthis%EulerAnglesBC   = .false.
plotthis%stress          = .false.
plotthis%strain          = .false.
plotthis%strainrate      = .false.
plotthis%viscosity       = .false.
#ifdef EXPERIMENTALTRANSIENTDISCREEP
plotthis%taylorstress    = .false.
#endif
plotthis%volume          = .false.
plotthis%elempartition   = .false.
plotthis%vertexpartition = .false.
plotthis%writeASCII      = .true.
plotthis%precision       = 1
plotthis%elemQualitySide     = .false.
plotthis%elemQualitySurface  = .false.
plotthis%elemQualityVolume   = .false.
plotthis%error           = .false.

end subroutine





end module
