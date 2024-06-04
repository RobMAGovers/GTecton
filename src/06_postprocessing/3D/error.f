subroutine sumFNerror3D ()

! for each element,
! COmpute Sum over F . n
! With F the displacement gradient tensor,
! and the n the outgoing normal vector of every side,
! 
! check the difference with the incoming, summer over the neighbors.

use meshdatamodule,  only: meshdatactx, &
                           elementNeighborTableComplete, &
                           BuildElementNeighbourTable3D
use modeldatamodule, only: modeldatactx
use modeldefinition, only: lgdef
use algebra,         only: vmprd, vectorLength, angleBetweenTwoVectors3D
use modeltopology,   only: ndof, nsd, nen
use constants,       only: useDXE, one

implicit none

integer :: iElement, iSide, iPoint, neighborID
integer :: fileID
double precision, allocatable :: SumFnDiff(:,:)
logical :: elemIsQuadrilateral

double precision :: BP(ndof, nsd) ! ndof and nsd are not the same in axi, opn(?)

double precision :: sideVec(nsd), normalVec(nsd)
double precision :: xl(nsd,nen), dl(nsd,nen)

double precision :: contributionVec(nsd)

double precision :: elemCenter(3)

integer          :: ierr


!call testNormals()

if (.not. elementNeighborTableComplete) then
    call BuildElementNeighbourTable3D()
endif

allocate(SumFnDiff(3,meshdatactx%neglobal))

do iElement = 1, meshdatactx%neglobal

    call LCOORD (meshdatactx%X,xl,meshdatactx%IEN(:,iElement),iElement)
    call LDISP  (dl,modeldatactx%D,meshdatactx%IEN(:,iElement),NDOF,NEN)

    call AddFaultDisplacement (dl,iElement,modeldatactx%TFAULT,NDOF,NEN,35)
    call ADDSNE (dl,iElement,NDOF,NEN,useDXE)
    call REZONE(xl,dl,ONE)

    elemCenter = 0.25d0 * (xl(:,1) + xl(:,2) + xl(:,3) + xl(:,4))

    ! elemIsQuadrilateral is not used in 3D,
    ! but the subroutine requires the argument.
    call BPMATRIX (XL,DL,elemIsQuadrilateral,BP,ierr)

    do iSide = 1,4
        if (meshdatactx%ElementNeighbours(iElement, iSide) .eq. 0) then
            ! the element is on edge of the domain, it has no
            ! neighbor on this side. Skip it.
        else
            ! there is neighbor here. 
            ! add the contribution of the normal vector, and subtract is from the neighbor.

            if (iSide.eq.1) then ! side from point 1, 2, 4
                call getNormal3D (xl(:,1), xl(:,2), xl(:,4), elemCenter, normalVec)

            else if (iSide.eq.2) then ! side from point 1, 3, 2
                call getNormal3D (xl(:,1), xl(:,3), xl(:,2), elemCenter, normalVec)

            else if (iSide.eq.3) then ! side from point 1, 4, 3
                call getNormal3D (xl(:,1), xl(:,4), xl(:,3), elemCenter, normalVec)

            else ! (iSide.eq.4) ! side from point 2, 3, 4
                call getNormal3D (xl(:,2), xl(:,3), xl(:,4), elemCenter, normalVec)
                
            endif

            call VMPRD(BP, normalVec, contributionVec, 3, 3)

            modeldatactx%elemError(iElement,:) = &
            modeldatactx%elemError(iElement,:) + contributionVec

            neighborID = meshdatactx%ElementNeighbours(iElement, iSide)

            modeldatactx%elemError(neighborID,:) = &
            modeldatactx%elemError(neighborID,:) + contributionVec

        endif
    enddo

enddo

end subroutine


subroutine testNormals()

use algebra, only: angleBetweenTwoVectors3D

implicit none

! take unit tetrahedron

double precision :: xl(3,4), elemCenter(3), normalVec(3)
integer          :: iSide

double precision :: v1(3), v2(3)

xl(:,1) = (/0,0,0/)
xl(:,2) = (/1,0,0/)
xl(:,3) = (/0,1,0/)
xl(:,4) = (/0,0,1/)

elemCenter = 0.25d0 * (xl(:,1) + xl(:,2) + xl(:,3) + xl(:,4))


v1 = (/0d0,2d0,2d0/)
v2 = (/0d0,0d0,3d0/)

write(*,*) "angle1 ", angleBetweenTwoVectors3D(v1, v2)

v1 = (/-1d0,-1d0,-1d0/)
v2 = (/-3d0,-3d0,-3d0/)

write(*,*) "angle2 ", angleBetweenTwoVectors3D(v1, v2)



write(*,*) "normal of side 1 shouild be:  0, -1,  0"
write(*,*) "normal of side 2 shouild be:  0,  0, -1"
write(*,*) "normal of side 3 shouild be: -1,  0,  0"
write(*,*) "normal of side 4 shouild be: (1,  1,  1) / sqrt3)"


do iSide = 1,4

            if (iSide.eq.1) then ! side from point 1, 2, 4
                call getNormal3D (xl(:,1), xl(:,2), xl(:,4), elemCenter, normalVec)

            else if (iSide.eq.2) then ! side from point 1, 3, 2
                call getNormal3D (xl(:,1), xl(:,3), xl(:,2), elemCenter, normalVec)

            else if (iSide.eq.3) then ! side from point 1, 4, 3
                call getNormal3D (xl(:,1), xl(:,4), xl(:,3), elemCenter, normalVec)

            else ! (iSide.eq.4) ! side from point 2, 3, 4
                call getNormal3D (xl(:,2), xl(:,3), xl(:,4), elemCenter, normalVec)

            endif

    write(*,*) "Side", iSide, "has 3D normal:", normalVec

enddo


stop "Normals test completed"

end subroutine


subroutine getNormal3D (p1, p2, p3, elemCenter, normal)
! normal of a tetrahedron is not entirely triavial.
! tested and OK.

use algebra, only: vectorLength, angleBetweenTwoVectors3D, crossp, dot
use constants, only: halfpi, pi

implicit none
double precision :: p1(3), p2(3), p3(3), elemCenter(3), normal(3)

double precision :: vec1(3), vec2(3), vecToCenter(3), centerOfPlane(3)

double precision :: angle, toacos

! Not entirely sure whether the order of the sides is also fixed.
! So we compute a vector that is orhtogongal to the plane spanned
! by the two vectors vec1 and vec2, compare its angle to the center,
! and if it is smaller than 90 degrees, we have accidentally
! computed the inside normal and we switch sides.

centerOfPlane = (p1 + p2 + p3) / 3d0
vecToCenter = elemCenter - centerOfPlane 

vec1 = p3 - p2
vec2 = p3 - p1

call crossp(vec1, vec2, normal)

toacos = dot(vecToCenter, normal, 3) / &
                 (vectorLength(normal, 3) * vectorLength(vecToCenter, 3))

! due to numerical rounding in the division above, 
! this term can reach outside the [-1, 1] interval,
! which would crash the acos operation.
! Correct these values:
if (toacos.lt.-1d0) then
    toacos=-1d0
else if (toacos.gt.1d0) then
    toacos=1d0
endif

angle = acos(toacos)

if (angle .lt. halfpi .or. angle .gt. pi + halfpi) then
    normal = -normal
endif

! normalize the normal
normal = normal / vectorLength(normal,3)

end subroutine

