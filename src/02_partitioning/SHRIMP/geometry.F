module geometry

double precision, parameter :: deg2rad = 4d0 * atan(1d0) / 180d0


! contains subroutines separately, so that they can be properly tested

contains

subroutine angle2abc(angle, x0, y0, a, b, c)

! converts a line knows as going through x0, y0 and with angle alpha
! into ax + by + c = 0

implicit none

double precision :: angle, x0, y0
double precision :: a, b, c

double precision :: slope

! determine slope 

if (angle .lt. -89.9 .or. angle .gt. 89.9) then
    stop "angle2abc gets too steep angle. This is not supposed to happen. Leaving SHRIMP."
endif

!a = tan(angle * deg2rad)
a = tan(angle)
b = -1d0

! to cmpute c, determine where it goes through x = 0

c = y0 - x0 * a

end subroutine

!-------------------------------------------------------------------------------

subroutine angles2abcd(phi, theta, x0, y0, z0, a, b, c, d)

! converts a line knows as going through x0, y0 and with angle alpha
! into ax + by + c = 0

implicit none

double precision :: phi, theta, x0, y0, z0
double precision :: a, b, c, d

double precision :: p(3), q(3), u(3), v(3)
double precision :: vOrthogonal(3)

double precision :: slope



! determine slope 

!if (phi .lt. -89.9 .or. phi .gt. 89.9) then
!    stop "angle2abc gets too steep angle. This is not supposed to happen. Leaving PRIMP."
!endif


! we have two points, p at (1,0,0) and q at (0,1,0).
! These are rotated over the angles given.
! They are seen as base vectors spanning a plane.
! And an orthogonal vector is computed to define the plane 
! in the form of ax + by + cy + z = 0

p = (/1,0,0/)
q = (/0,1,0/)

call rotatePointIn3D(phi, theta, p, u)
call rotatePointIn3D(phi, theta, q, v)


call findVectorOrthoganalToTwo(vOrthogonal, u, v)

a = vOrthogonal(1)
b = vOrthogonal(2)
c = vOrthogonal(3)

! to compute d, determine where it goes through the midpoint.

d = - x0 * a - y0 * b - z0 * c

end subroutine




!-------------------------------------------------------------------------------

subroutine rotatePointIn3D(phi, theta, oldCoords, newCoords)

implicit none

double precision :: phi, theta
double precision :: oldCoords(3), newCoords(3)

! use a combination of rotation on x (phi),
! and a rotation on y (theta).

!      [ 1    0          0      ]
! Rx = [ 0  cos(phi)  -sin(phi) ]
!      [ 0  sin(phi)   cos(phi) ]

!      [  cos(theta)  0  sin(theta) ]
! Ry = [     0        1      0      ]
!      [ -sin(theta)  0  cos(theta) ]

!         [ cos(theta)                   0             sin(theta)      ]
! Rx Ry = [ sin(phi) * sin(theta)      cos(phi)   -sin(phi)*cos(theta) ]
!         [ -cos(phi) * sin(theta)     sin(phi)   cos(phi)*cos(theta)  ]


newCoords(1) =  cos(theta) *            oldCoords(1) + &
                sin(theta) *            oldCoords(3)

newCoords(2) =  sin(phi) * sin(theta) * oldCoords(1) + &
                cos(phi) *              oldCoords(2) - &
                sin(phi) * cos(theta) * oldCoords(3)

newCoords(3) = -cos(phi) * sin(theta) * oldCoords(1) + &
                sin(phi) *              oldCoords(2) + &
                cos(phi) * cos(theta) * oldCoords(3)

end subroutine

!-------------------------------------------------------------------------------

subroutine findVectorOrthoganalToTwo(vOrthogonal, u, v)

! use cross product to find an orthoginal vector

implicit none

double precision :: u(3), v(3)
double precision :: vOrthogonal(3)

vOrthogonal(1) = u(2)*v(3) - v(2)*u(3)
vOrthogonal(2) = u(1)*v(3) - v(1)*u(3)  
vOrthogonal(3) = u(1)*v(2) - v(1)*u(2)  

end subroutine



!-------------------------------------------------------------------------------

double precision function DistanceBetweenPointAndLine(xp, yp, a, b, c)

implicit none

! computes the distance between point coordsPoint and line with a 

! the line is described as ax + by + c = 0

double precision :: xp, yp
double precision :: a, b, c

! shamelessly following wikipedia, we get, the distance between a line
! al x + bl y + cl = 0
! and a point xp, yp is
!  | al xp + bl yp + cl |
! ---------------------
!    sqrt(a^2 + b^2)

DistanceBetweenPointAndLine = abs(a * xp + b * yp + c)/sqrt(a**2 + b**2)

end function

!-------------------------------------------------------------------------------

double precision function DistanceBetweenPointAndPlane(xp, yp, zp, a, b, c, d)

implicit none

! computes the distance between point coordsPoint and plane

! a plane is described as ax + by + cz + d = 0

integer :: nDimensions
double precision :: xp, yp, zp
double precision :: a,b,c,d

DistanceBetweenPointAndPlane = abs(a * xp + b * yp + c * zp + d)/sqrt(a**2 + b**2 + c**2)

end function

!-------------------------------------------------------------------------------

integer function nElementsWithinRangeOfLine(nLocalPoints,localCoords,a,b,c,maxDistance)

implicit none

integer          :: nLocalPoints
double precision :: localCoords(2,nLocalPoints)
double precision :: a,b,c
double precision :: maxDistance, thisDistance

integer          :: iPoint

nElementsWithinRangeOfLine = 0

do iPoint = 1, nLocalPoints
	thisDistance = DistanceBetweenPointAndLine(localCoords(1,iPoint), localCoords(2,iPoint), a, b, c)

!	write(*,*) iPoint, localCoords(:,iPoint), thisDistance, maxDistance

	if (thisDistance .le. maxDistance) then
!		write(*,*) "yay!"
		nElementsWithinRangeOfLine = nElementsWithinRangeOfLine + 1
	endif

enddo

end function

!-------------------------------------------------------------------------------

double precision function findBestRange(nLocalPoints, localCoords, nDimensions)

! determine extrema of the domain,
! pick the distance one 20-th of it.
! Works great for squares.

use moduleMPI

implicit none

integer          :: nLocalPoints, nDimensions
double precision :: localCoords(nDimensions, nLocalPoints)


double precision :: coordsMinLocal(nDimensions)
double precision :: coordsMaxLocal(nDimensions)
double precision :: coordsMinGlobal(nDimensions)
double precision :: coordsMaxGlobal(nDimensions)

integer          :: iPoint, iDimension
integer	       	 :: iStart, iEnd, stepsize
integer	       	 :: iError

if (nLocalPoints .lt. 10000) then
	iStart   = 1
	iEnd     = nLocalPoints
	stepsize = 1
else
	call setRange(10000, nLocalPoints, iStart, iEnd, stepsize)
endif

coordsMinLocal =  9999999999d0
coordsMaxLocal = -9999999999d0


do iPoint = iStart, iEnd, stepsize
	do iDimension = 1,nDimensions
		if (localCoords(iDimension,iPoint) .lt. coordsMinLocal(iDimension)) then
			coordsMinLocal(iDimension) = localCoords(iDimension,iPoint)
		endif
        if (localCoords(iDimension,iPoint) .gt. coordsMaxLocal(iDimension)) then
    	    coordsMaxLocal(iDimension) = localCoords(iDimension,iPoint)
    	endif
	enddo
enddo

! allreduce to a global mimimum and maximum


call MPI_Allreduce( &
    coordsMinLocal, &
    coordsMinGlobal, &
    nDimensions, &
    MPI_Double, &
    MPI_MIN, &
    MPI_COMM_WORLD, &
    iError)

call MPI_Allreduce( &
    coordsMaxLocal, &
    coordsMaxGlobal, &
    nDimensions, &
    MPI_Double,	&
    MPI_MAX, &
    MPI_COMM_WORLD, & 
    iError) 

findBestRange = 0.05 * sqrt(sum((coordsMaxGlobal - coordsMinGlobal)**2))



end function

!-------------------------------------------------------------------------------

subroutine setRange(selectionSize, totalSize, iStart, iEnd, stepsize)

implicit none

! let our local node have a few million nodes. (totalSize)
! We wish to find the extrema, but we only wish to sample a selection
! of these nodes (selectionSize,), for speed.

integer :: selectionSize, totalSize

integer :: iStart, iEnd, stepsize

stepsize = (totalSize - modulo(totalSize, selectionSize)) / selectionSize

iStart = stepsize
iEnd   = stepsize * selectionSize

end subroutine

!-------------------------------------------------------------------------------

subroutine neighborKnown(neighborListLength, &
                         neighborList, &
                         blockSize, &
                         neighborOf, &   ! a neighbor of which point? (Local nr!)
                         neighborID, &
                         neighborIsKnown, &
                         highestCheckedIndex)

implicit none

integer :: neighborListLength
integer :: neighborList(neighborListLength)
integer	:: blockSize
integer :: neighborOf
integer :: neighborID
logical :: neighborIsKnown
integer :: highestCheckedIndex

integer :: checkIndex



neighborIsKnown = .false.

! walk through the list, and 

checkIndex = blockSize * (neighborOf-1) + 1

do while (neighborList(checkIndex) .ne. 0)
	if (neighborList(checkIndex) .eq. neighborID) then
		neighborIsKnown	= .true. 
		return
	else
		! maybe the next
		checkIndex = checkIndex + 1
		
		! handle the end of a block
		if (modulo(checkIndex, blockSize) .eq. 0) then
			! shift to the next block
			checkIndex = neighborList(checkIndex)
		endif
	endif
enddo

! we did not find the neighbor. New one. Yay :-)


end subroutine

!-------------------------------------------------------------------------------


end module
