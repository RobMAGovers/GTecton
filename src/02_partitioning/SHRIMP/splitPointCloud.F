recursive subroutine splitPointCloud(nSteps,&           ! calling depth, decreasing for deeper recursion 
                                     nPointsLocal, &    ! number of points on this call on this thread
                                     nPointsGlobal, &   ! number of points on this call on all threads
                                     pointCoords, &     ! coordinates of local points
                                     partitionOfPoints, &
                                     partitionNumber)   ! a number that increases for deeper recursion,
                                                        ! and is applied to partitionOfPoints at the end
                                                        ! of the recursion tree.


use meshdata,  only: nDimensions, &
                     PointCountPerPartitionThread
                     
use moduleMPI, only: nProcessors, &
                     myProcessorID, &
                     parallelLog

use geometry,  only: deg2rad, &
                     distanceBetweenPointAndLine, &
                     distanceBetweenPointAndPlane, &
                     angle2abc, &
                     angles2abcd, &
                     setRange, &
                     findBestRange

implicit none

#include "mpif.h"


!-arguments

integer :: nSteps
integer :: nPointsLocal ! nr of points in this cloud in this particular thread
integer :: nPointsGlobal ! nr points in this entire point cloud
                         ! sum of all nPointsLocal equals nPointsGlobal.
double precision :: pointCoords(nDimensions, nPointsLocal)
integer :: partitionOfPoints(nPointsLocal)

!-internal variables

double precision, parameter   :: pi = 4d0 * atan(1d0)

double precision, allocatable :: CoordMidLocal(:), CoordMidGlobal(:)

double precision              :: weightOfMyAverage

integer, allocatable          :: nPointsLocalAllThreads(:)

! loop iterators
integer :: iThread, iPoint, jPoint

! error for the MPI calls
integer :: iError

integer                       :: myReceiveOffset
integer, allocatable          :: receiveOffset(:)
double precision, allocatable :: sampleCoordinates(:,:)
integer                       :: sampleInterval ! in case the mesh is too large to take all
integer                       :: nPointsPerThread
integer                       :: firstPointFromThread
integer                       :: lastPointFromThread
integer                       :: nsamplepointsfromthread
integer, allocatable          :: nPointsInLowerThreads(:)
integer, allocatable          :: nSamplePointsFromEachThread(:)

integer                       :: lengthWithoutHead
integer                                :: nMainBody
integer                                :: tailLength

double precision              :: majorAxis(nDimensions)

!- variables to pass on the next iteration
integer                       :: nPointsAbove, nPointsAboveGlobal
integer                       :: nPointsBelow, nPointsBelowGlobal
integer, allocatable          :: pointIDsAbove(:)
integer, allocatable          :: pointIDsBelow(:)
double precision, allocatable :: coordsAbove(:,:)
double precision, allocatable :: coordsBelow(:,:)
integer, allocatable          :: partitionsAbove(:)
integer, allocatable          :: partitionsBelow(:)

integer, allocatable          :: pointLookUp(:)

! variables to define the dividing plane
double precision              :: a,b,c,d
double precision              :: maxLength, thisLength
integer                       :: nSamplePoints
! support to make sure the lines do not become too long.
double precision              :: thisX, thisY, thisZ

integer                       :: partitionNumber
integer,          allocatable :: nCLosePointsatAngle2DLocal(:)
integer,   	  allocatable :: nCLosePointsatAngle2DGlobal(:)
integer,          allocatable :: nCLosePointsatAngle2DGlobalSmoothed(:)
integer,          allocatable :: nCLosePointsatAngle3DLocal(:,:)
integer,          allocatable :: nCLosePointsatAngle3DGlobal(:,:)
integer,          allocatable :: nCLosePointsatAngle3DGlobalSmoothed(:,:)

double precision, allocatable :: distanceAtAngle3D(:,:)

integer,          parameter   :: nAnglesteps = 178  ! 1 degree steps for now
double precision              :: separationSteps(nAnglesteps)

double precision, parameter   :: angleStart =  -89.0 * deg2rad
double precision, parameter   :: angleEnd   =   89.0 * deg2rad
double precision              :: angleStepSize
double precision              :: angle, phi, theta
integer                       :: iAngle, jAngle
double precision              :: smallestDistance
double precision              :: smallestAngle
integer                       :: iStart, iEnd, stepsize ! sample interval parameters
double precision              :: maxDistance

double precision              :: bmin, cmin
double precision              :: angleWithLowestProximityCount
double precision              :: phiWithLowestProximityCount
double precision              :: thetaWithLowestProximityCount


integer                       :: lowestProximityCount

integer                       :: sampleUpToPoint



! define the angle interval


do iAngle = 1,nAnglesteps
    separationSteps(iAngle) = -0.5 * pi + dble(iAngle) * pi / dble(nAnglesteps+1)
    separationSteps(iAngle) = tan(separationSteps(iAngle))
enddo

!write(*,*) "angle steps", separationSteps

!call parallelLog("Entering a new round of splitcloud")

allocate(CoordMidLocal(nDimensions))
allocate(CoordMidGlobal(nDimensions))

!write(*,*) "thread", myProcessorID, "has total cloud size", nPointsGlobal

weightOfMyAverage = dble(nPointsLocal) / dble(nPointsGlobal)



!write(*,*) "thread", myProcessorID, "has called split with nSteps", nSteps, nDimensions


!write(*,*) "Thread", myProcessorID, "has nPointsLocal", nPointsLocal

!write(*,*) "Thread", myProcessorID, "has pointCoords", pointCoords


!write(*,*) "Thread", myProcessorID, "has partitionOfPoints:",partitionOfPoints
!write(*,*) "Thread", myProcessorID, "has partitionNumber:", partitionNumber


if (nPointsLocal.gt.0) then
    ! determine absolute midpoint of the point cloud
    CoordMidLocal(1) = sum(pointCoords(1,:)) / nPointsLocal
    CoordMidLocal(2) = sum(pointCoords(2,:)) / nPointsLocal
    if (nDimensions.eq.3) then
        CoordMidLocal(3) = sum(pointCoords(3,:)) / nPointsLocal
    endif
else
    CoordMidLocal(1) = 0d0
    CoordMidLocal(2) = 0d0
    if (nDimensions.eq.3) then
        CoordMidLocal(3) = 0d0
    endif
endif

!write(*,*) "Thread", myProcessorID, "has CoordMidLocal", CoordMidLocal
!write(*,*) "Thread", myProcessorID, "has weightOfMyAverage", weightOfMyAverage

!write(*,*) "Thread", myProcessorID, "assembles local midpoints ", CoordMidLocal * weightOfMyAverage

call MPI_Allreduce( &
    CoordMidLocal * weightOfMyAverage, &  ! sending from
    CoordMidGlobal, & ! sending to
    nDimensions, &    ! this many numbers in the array
    MPI_Double, &     ! data type
    MPI_SUM, &        ! Summation of weighted gives the average
    MPI_COMM_WORLD, iError)

!write(*,*) "Thread", myProcessorID, "has global midpoint", CoordMidGlobal

! CoordMidGlobal now has the midpoint for this point cloud.
! for this we must assemble a thousand random points
! (or all of the partition)
! However, not all partitions necessarily have nodes in this cloud.
! So first determine which thread is going to send how many random representatives. 

! There are three possibilities, written from easy to hard:
! 1: there are less that a thousand points in the partition.
!    In that case, take them all
! 2: There are more than thousand points in the point cloud,
!    And all partition have enough to contribute their part
! 3: There are more than thousand points in the point cloud,
!    But some partition do not have enough.

allocate(nPointsLocalAllThreads(nProcessors))

!making sure that the 0 processor knows how many 
call MPI_Allgather(nPointsLocal, &  ! data being sent
                1 ,           &  ! number of variables being sent. It is only a single integer
                MPI_Integer,  &  ! send type
                nPointsLocalAllThreads, & ! receive data here, +1 for counting starting at 0 (MPI threads) vs. 1 (Fortran)
                1 ,           &  ! receive this many per thread
                MPI_Integer,  &  ! yes yes, type does not change halfway on the road
                MPI_COMM_WORLD, &
                iError)  !

if (iError .ne. 0) then
    write(*,*) "Failed to gather local number of points. Error: ", iError
    stop "Leaving SHRIMP..."
endif


!write(*,*) "thread", myProcessorID, "has localAll", nPointsLocalAllThreads
!write(*,*) "thread",  myProcessorID, "has nPointsGlobal: ", nPointsGlobal



maxDistance = findBestRange(nPointsLocal, pointCoords, nDimensions)



if (nDimensions.eq.2) then

    allocate(nCLosePointsAtAngle2DLocal(nAnglesteps))
    allocate(nCLosePointsAtAngle2DGlobal(nAnglesteps))
    allocate(nCLosePointsatAngle2DGlobalSmoothed(nAnglesteps))

    nCLosePointsatAngle2DLocal = 0
    nCLosePointsatAngle2DGlobal = 0
	nCLosePointsatAngle2DGlobalSmoothed = 0

    if (nPointsLocal .ge. 1) then

        ! determine the range of points we will sample over.
        if (nPointsLocal .lt. 10000) then
            ! not very many points. Sample all of them
            iStart   = 1
            iEnd     = nPointsLocal
            stepsize = 1
        else
            ! many points present. Sample over all of them.
            call setRange(10000, nPointsLocal, iStart, iEnd, stepsize)
        endif


!        write(*,*) "thread doing smpling: ", iStart, iEnd, stepsize

        ! compute the angle at which a line through the points crosses the fewest points.




!*****************

!        angleStepSize = (angleEnd - angleStart) / dble(nAnglesteps-1)

!        do iAngle = 1, nAnglesteps
!            angle = angleStart + (iAngle-1) * angleStepSize

!            call angle2abc(angle, CoordMidGlobal(1), CoordMidGlobal(2), a, b, c)

        a = 1d0
        do iAngle = 1, nAnglesteps


!            call angle2abc(angle, CoordMidGlobal(1), CoordMidGlobal(2), a, b, c)
!            write(*,*) "old: ", a,b,c

            a = 1d0
            b = separationSteps(iAngle)
            ! use midpoint and ax+by+c=0 to determine c
            c = -a*CoordMidGlobal(1) - b*CoordMidGlobal(2)

!            write(*,*) "new: ",    a,b,c
            do iPoint = iStart, iEnd, stepsize

                if(DistanceBetweenPointAndLine(pointCoords(1,iPoint), &
                                               pointCoords(2,iPoint), &
                                               a, b, c) &
                      .lt. maxDistance) then
                    nCLosePointsatAngle2DLocal(iAngle) = &
                    nCLosePointsatAngle2DLocal(iAngle) + 1
                endif

            enddo

        enddo

!*********************

!        write(*,*) "finished angle loop", nCLosePointsatAngle2DLocal

!        stop "did point loop"

    else
        ! This thread does not have a point in this partition. 
        ! No distance can be determined.
        ! Return a zero-vector
    endif

    ! add the local proximity count into its global equivalent

    call MPI_Allreduce( &
        nCLosePointsAtAngle2DLocal, & ! sending from
        nCLosePointsAtAngle2DGlobal, &  ! sending to
        nAngleSteps, &                  ! this many numbers in the array
        MPI_Integer, &                  ! data type
        MPI_SUM, &                      ! Summation of weighted gives the average
        MPI_COMM_WORLD, iError)

    ! take moving window average of the global count
    ! Eventually. Take just the smallest for now, though.

    nCLosePointsatAngle2DGlobalSmoothed = nCLosePointsAtAngle2DGlobal



else if (nDimensions .eq. 3) then

    allocate(nCLosePointsAtAngle3DLocal(nAnglesteps,nAnglesteps))
    allocate(nCLosePointsAtAngle3DGlobal(nAnglesteps,nAnglesteps))
    allocate(nCLosePointsatAngle3DGlobalSmoothed(nAnglesteps,nAnglesteps))

    nCLosePointsatAngle3DLocal = 0
    nCLosePointsatAngle3DGlobal = 0
	nCLosePointsatAngle3DGlobalSmoothed = 0

    if (nPointsLocal .ge. 1) then

        ! determine the range of points we will sample over.
        if (nPointsLocal .lt. 10000) then
            ! not very many points. Sample all of them
            iStart   = 1
            iEnd     = nPointsLocal
            stepsize = 1
        else
            ! many points present. Sample over all of them.
            call setRange(10000, nPointsLocal, iStart, iEnd, stepsize)
        endif


!        write(*,*) "thread doing smpling: ", iStart, iEnd, stepsize

        ! compute the angle at which a line through the points crosses the fewest points.
        angleStepSize = (angleEnd - angleStart) / dble(nAnglesteps-1)

!        do iAngle = 1, nAnglesteps
!            phi = angleStart + (iAngle-1) * angleStepSize
!            do jAngle = 1, nAnglesteps
!                theta = angleStart + (jAngle-1) * angleStepSize

!                call angles2abcd(phi, theta, CoordMidGlobal(1), &
!                                             CoordMidGlobal(2), &
!                                             CoordMidGlobal(3), a, b, c, d)
        a = 1d0
        do iAngle = 1, nAnglesteps
            b = separationSteps(iAngle)
            do jAngle = 1, nAnglesteps
                c = separationSteps(jAngle)
                d = - a * CoordMidGlobal(1) - &
                      b * CoordMidGlobal(2) - &
                      c * CoordMidGlobal(3)


                do iPoint = iStart, iEnd, stepsize

                    if(DistanceBetweenPointAndPlane(pointCoords(1,iPoint), &
                                                    pointCoords(2,iPoint), &
                                                    pointCoords(3,iPoint), &
                                                    a, b, c, d) &
                          .lt. maxDistance) then
                        nCLosePointsatAngle3DLocal(iAngle, jAngle) = &
                        nCLosePointsatAngle3DLocal(iAngle, jAngle) + 1
                    endif

                enddo

!               write(*,*) "plane has max Distance", maxDistance

            enddo
        enddo

    else
        ! This thread does not have a point in this partition. 
        ! No distance can be determined.
        ! Return a zero-vector
    endif

    call MPI_Allreduce( &
        nCLosePointsAtAngle3DLocal, &   ! sending from
        nCLosePointsAtAngle3DGlobal, &  ! sending to
        nAngleSteps**2, &               ! this many numbers in the array
        MPI_Integer, &           	! data type
        MPI_SUM, &                      ! Summation of weighted gives the average
        MPI_COMM_WORLD, iError)



    nCLosePointsatAngle3DGlobalSmoothed = nCLosePointsAtAngle3DGlobal


else
	write(*,*) "Found n Dimensions to be", nDimensions
    stop "n Dimensions not 2 or 3. This is not supposed to happen. Leaving SHRIMP..."
endif



! all thread have constrcted their local angle table.
! add results

!write(*,*) "sharing angle data", nCLosePointsAtAngle3DGlobal

!stop "done"


! Determine plane that separates the two halves.
! This is the plane that is orthogonal to the connecting line determined above,
! and passes through the point

! It is possible that the line would be exactly horizontal, in which case 
! the diving plane would be exactly vertical, which could not be described by an
! equation of the form y = ax + b or z = ax + by + c

! To prevent this, we will give the line a slight angle off the horizontal.


if (nDimensions .eq. 2) then

    angleWithLowestProximityCount = 0d0
    bmin = 0d0
    lowestProximityCount = 99999999


    a = 1d0
    do iAngle = 1, nAnglesteps
        b = separationSteps(iAngle)
        ! use midpoint and ax+by+c=0 to determine c
        c = -a*CoordMidGlobal(1) - b*CoordMidGlobal(2)


!    do iAngle = 1, nAnglesteps

!        angle = angleStart + (iAngle-1) * angleStepSize
!        write(*,*) "angle check", iAngle, angle, nCLosePointsatAngle2DGlobalSmoothed(iAngle), angleWithLowestProximityCount


        if (nCLosePointsatAngle2DGlobalSmoothed(iAngle) .lt. lowestProximityCount) then
!            angleWithLowestProximityCount = angle
            bmin                          = b
            lowestProximityCount = nCLosePointsatAngle2DGlobalSmoothed(iAngle)
!             write(*,*) "angle check YES" , lowestProximityCount, angleWithLowestProximityCount
        endif

    enddo
!    write(*,*) "best line has angle: angleWithLowestProximityCount", angleWithLowestProximityCount, CoordMidGlobal


!    call angle2abc(angleWithLowestProximityCount, CoordMidGlobal(1), CoordMidGlobal(2), a, b, c)
    b = bmin
    c = -a*CoordMidGlobal(1) - b*CoordMidGlobal(2)


!    write(*,*) "thread", myProcessorID, "has a,b,c of partition divider", a,b,c

else
      phiWithLowestProximityCount = 0d0
    thetaWithLowestProximityCount = 0d0
    bmin = 0d0
    cmin = 0d0


    lowestProximityCount = 99999999


    do iAngle = 1, nAnglesteps
        b = separationSteps(iAngle)
        do jAngle = 1, nAnglesteps
            c = separationSteps(jAngle)
            d = - a * CoordMidGlobal(1) - &
                  b * CoordMidGlobal(2) - &
                  c * CoordMidGlobal(3)



!    do iAngle = 1, nAnglesteps
!        phi = angleStart + (iAngle-1) * angleStepSize
!        do jAngle = 1, nAnglesteps
!            theta = angleStart + (jAngle-1) * angleStepSize


!        write(*,*) "angle check", iAngle, angle, nCLosePointsatAngle2DGlobalSmoothed(iAngle), angleWithLowestProximityCount

            if (nCLosePointsatAngle3DGlobalSmoothed(iAngle, jAngle) .lt. lowestProximityCount) then
!                  phiWithLowestProximityCount = phi
!                thetaWithLowestProximityCount = theta
                bmin = b
                cmin = c


                lowestProximityCount = nCLosePointsatAngle3DGlobalSmoothed(iAngle, jAngle)
!             write(*,*) "angle check YES" , lowestProximityCount, angleWithLowestProximityCount
            endif

        enddo
    enddo

    b = bmin
    c = cmin
    d = - a * CoordMidGlobal(1) - &
          b * CoordMidGlobal(2) - &
          c * CoordMidGlobal(3)

!    write(*,*) "best plane has angles: ", phiWithLowestProximityCount, &
!                                          thetaWithLowestProximityCount, &
!                                          lowestProximityCount, &
!                                          CoordMidGlobal


!    call angles2abcd(  phiWithLowestProximityCount,& 
!                     thetaWithLowestProximityCount,&
!                      CoordMidGlobal(1), &
!                      CoordMidGlobal(2), &
!                      CoordMidGlobal(3), &
!                      a, b, c, d)

!    write(*,*) "has a, b, c, d", a, b, c, d


endif




! split the point cloud into two, split by the plane.
! And when necessary, 

! allocate the Above and Below spaces to accomodate all the points.
! This is too much, but it saves on time for reallocating the
! array when more space is needed.
! Shaving off space will eventually lead to shaving off further,
! until this one model where it runs out fo space.

! Only the required amount of data will be passed on to the next iteration,
! so not much memory is lost.

!write(*,*) "thread", myProcessorID, "allocates splitting room"
nPointsAbove = 0
nPointsBelow = 0
allocate(pointIDsAbove(nPointsLocal))
allocate(pointIDsBelow(nPointsLocal))
allocate(coordsAbove(nDimensions,nPointsLocal))
allocate(coordsBelow(nDimensions,nPointsLocal))
allocate(partitionsAbove(nPointsLocal))
allocate(partitionsBelow(nPointsLocal))
pointIDsAbove = 0
pointIDsBelow = 0
coordsAbove = 0d0
coordsBelow = 0d0
partitionsAbove = 0
partitionsBelow = 0


! lookup table to recombine the partitions from the lower and upper halves.
! Nodes from the upper half will have their sequence there with a + sign,
! Nodes from the lower half will have their sequence there with a - sign
allocate(pointLookUp(nPointsLocal))
pointLookUp = 0

if (nDimensions .eq. 2) then

    do iPoint = 1, nPointsLocal
        thisX = pointCoords(1, iPoint)
        thisY = pointCoords(2, iPoint)
        ! compare with line ax + by + c = 0
!        write(*,*) "comparing point", iPoint, thisX, thisY, a,b,c,thisX * a + thisY * b + c
        if (thisX * a + thisY * b + c .gt. 0d0) then
!            write(*,*) "thread", myProcessorID, "has point above", thisX, thisY, a, b
            nPointsAbove                = nPointsAbove + 1
            pointIDsAbove(nPointsAbove) = iPoint
            coordsAbove(:,nPointsAbove) = pointCoords(:,iPoint)
            pointLookUp(iPoint)         = nPointsAbove
        else
!            write(*,*) "thread", myProcessorID, "has point below", thisX, thisY, a, b

            nPointsBelow                = nPointsBelow + 1
            pointIDsBelow(nPointsBelow) = iPoint
            coordsBelow(:,nPointsBelow) = pointCoords(:,iPoint)
            pointLookUp(iPoint)         = -nPointsBelow
        endif

    enddo

!    write(*,*) "thread", myProcessorID, "has coordsAbove", coordsAbove
!    write(*,*) "thread", myProcessorID, "has coordsBelow", coordsBelow

else

    do iPoint = 1, nPointsLocal
        thisX = pointCoords(1, iPoint)
        thisY = pointCoords(2, iPoint)
        thisZ = pointCoords(3, iPoint)

        if (thisX * a + thisY * b + thisZ * c + d .gt. 0) then
            nPointsAbove = nPointsAbove + 1
            pointIDsAbove(nPointsAbove) = iPoint
            coordsAbove(:,nPointsAbove) = pointCoords(:,iPoint)
            pointLookUp(iPoint)    = nPointsAbove
        else
            nPointsBelow = nPointsBelow + 1
            pointIDsBelow(nPointsBelow) = iPoint
            coordsBelow(:,nPointsBelow) = pointCoords(:,iPoint)
            pointLookUp(iPoint) = -nPointsBelow
        endif

    enddo



endif

!write(*,*) "thread", myProcessorID, "has split the points"


!write(*,*) "thread", myProcessorID, "reduce local count" , nPointsAbove, nPointsBelow

! gather the points split up above and below, to get the global number, for the recursive calls.

call MPI_Allreduce( &
    nPointsAbove, &  ! sending from
    nPointsAboveGlobal, & ! sending to
    1, &    ! this many numbers in the array
    MPI_Integer, &     ! data type
    MPI_SUM, &        ! Summation of weighted gives the average
    MPI_COMM_WORLD, iError)

call MPI_Allreduce( &
    nPointsBelow, &  ! sending from
    nPointsBelowGlobal, & ! sending to
    1, &    ! this many numbers in the array
    MPI_Integer, &     ! data type
    MPI_SUM, &        ! Summation of weighted gives the average
    MPI_COMM_WORLD, iError)


!write(*,*) "thread", myProcessorID, "has reduced the global point count: ", nPointsAboveGlobal, nPointsBelowGlobal

if (nSteps.gt.1) then
    ! Not yet reached bottom of recursion.
    ! Split each half that we just created into two smaller clouds.

!    write(*,*) "thread", myProcessorID, "calls next split with nSteps", nSteps, nDimensions
!    write(*,*) "thread", myProcessorID, "sends coordsAbove", coordsAbove

!write(*,*) "*******************************************************************************"
!write(*,*) "*******************************************************************************"
!write(*,*) "***  Calling Above recursion **************************************************"
!write(*,*) "*******************************************************************************"

!write(*,*) "*******************************************************************************"


!    call parallelLog("calling above recursion")


    call MPI_Barrier(MPI_COMM_WORLD, iError)

 !   write(*,*) "thread", myProcessorID,	"calls above with part", partitionNumber, partitionNumber + 2**(nSteps-1)


    ! points above the dividing plane
    call splitPointCloud(nSteps-1, &
                         nPointsAbove, &
                         nPointsAboveGlobal, &
                         coordsAbove(:,1:nPointsAbove), &
                         partitionsAbove(1:nPointsAbove), &
                         partitionNumber + 2**(nSteps-1))


    !call parallelLog("finished above recursion")



!write(*,*) "*******************************************************************************"
!write(*,*) "*******************************************************************************"
!write(*,*) "***  Calling Below recursion **************************************************"
!write(*,*) "*******************************************************************************"
!write(*,*) "*******************************************************************************"

!    write(*,*) "thread", myProcessorID, "sends coordsBelow", coordsBelow

!    call parallelLog("calling below recursion")

!	write(*,*) "thread", myProcessorID, "calls below with part", partitionNumber, partitionNumber


    ! point underneath the dividing plane
    call splitPointCloud(nSteps-1, &
                         nPointsBelow, &
                         nPointsBelowGlobal, &
                         coordsBelow(:,1:nPointsBelow), &
                         partitionsBelow(1:nPointsAbove), &
                         partitionNumber)

    !call parallelLog("finished a below recursion")

!write(*,*) "*******************************************************************************"
!write(*,*) "*******************************************************************************"
!write(*,*) "******* done recursive calls **************************************************"
!write(*,*) "*******************************************************************************"
!write(*,*) "*******************************************************************************"

    ! recombine partitionsAbove and partitionsBelow in the proper sequence.
    do iPoint = 1,nPointsLocal
        if (pointLookUp(iPoint) .lt. 0) then
            ! compensate for the negative numbers.
            partitionOfPoints(iPoint) = partitionsBelow(-pointLookUp(iPoint))
        else
            partitionOfPoints(iPoint) = partitionsAbove(pointLookUp(iPoint))
        endif
    enddo


else

	! we are the bottom of the recursion.
	! Assign the partition numbers.

	do iPoint = 1,nPointsLocal
		if (pointLookUp(iPoint) .lt. 0) then
            ! compensate for the negative numbers.
            partitionOfPoints(iPoint) = partitionNumber 

			PointCountPerPartitionThread(partitionNumber) = &
            PointCountPerPartitionThread(partitionNumber) + 1

        else
            partitionOfPoints(iPoint) = partitionNumber + 1

    	    PointCountPerPartitionThread(partitionNumber + 1) = &
      	    PointCountPerPartitionThread(partitionNumber + 1) + 1


        endif

!		write(*,*) "thread", myProcessorID, "assigned point partition ", partitionOfPoints(iPoint)

    enddo

endif


!call parallelLog("Reached very end of splitcloud")

end subroutine


