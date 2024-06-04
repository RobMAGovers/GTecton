subroutine algo3()

use globals

implicit none

integer :: iElem, iPoint, iInterval
integer :: matchID(3)
logical :: match(3)
integer :: sideID
double precision  :: centerCoords(2)
double precision, allocatable :: tmp(:,:)
integer, external :: numInArray
integer, external :: setWeight
integer :: p1match, p2match
integer          :: weight
integer          :: p1, p2
integer, allocatable :: weights(:)

allocate(weights(nelems))
weights = 0
allocate(tmp(2,numloop))
tmp = 0d0

if (ndimensions .eq. 3) then
    write(*,*) "algorithm 3 only works in two dimensions."
else
    if (iecho.eq.1) then
        write(*,*) "algo3 says: hello world"
    endif
endif

! for this algorithm,
! the nodelist is ordered spatially using the polysort tool.
! we step through the nodes

if (iecho.eq.1) then
    write(*,*) "# of nodes = ",nvertices," interpreteds as # of interface intervals = ",nvertices-1
endif

do iInterval = 1, nvertices-1

    p1 = vertices(iInterval)
    p2 = vertices(iInterval+1)

    if (iecho.eq.1) then
        write(*,*) "----------------------------------------------------------------------------"
        write(*,*) "Next interface section, from node ", p1, "to node", p2
        write(*,*) "----------------------------------------------------------------------------"
    endif

    do iElem = 1, nelems

        p1match = numInArray(p1, v(1:3,iElem),3)
        p2match = numInArray(p2, v(1:3,iElem),3)

        if (iecho.eq.1) then
            if (p1match.gt.0 .and. p2match.gt.0) then
                write(*,*) "element",iElem," has two nodes on the interface: ",p1, p2
            else if (p1match.gt.0 .and. p2match.eq.0) then
                write(*,*) "element",iElem," has one node on the interface: ",p1
            else if (p1match.eq.0 .and. p2match.gt.0) then
                write(*,*) "element",iElem," has one node on the interface: ",p2
            endif
        endif

        ! find elements that have a side on the edge. Determine the edge and sign
        if (p1match .gt. 0 .and. p2match .gt. 0) then
            if (p1match + p2match .eq. 3) then
                ! 1 and 2
                sideID = 1
            else if (p1match + p2match .eq. 5) then
                ! 2 and 3
                sideID = 2
            else if (p1match + p2match .eq. 4) then
                ! 1 and 3
                sideID = 3
            else
                write(0,*) "unknown side ID", p1match, p2match
                write(0,*) "This is not supposed to happen. Please contact model support."
                stop "Leaving elmside..."
            endif
            if (iecho.eq.1) then
                write(*,*) "Element",iElem,"has side ",sideID," on the interface"
            endif
            call elemCenter(iElem, centerCoords)
            call pnpoly(centerCoords,loopcoords,numloop,tmp,weight)
            if (iecho.eq.1) then
                write(*,*) "element center coordinates: ",centerCoords
            endif
            if (weight.eq.0) then
                write(0,*) "Element center of element ",iElem," lies on the loop"
                stop "leaving elmside"
            endif
            if (weights(iElem) .eq. 0) then
                weights(iElem) = weight
            else if (weights(iElem).eq.weight) then
                continue
            else
                write(0,*) "CONFLICT: Element",iElem,"has the side between", p1, p2, &
                  "on the interface, but was already assigned an opposite weight"
            endif
            call output2D(iElem, p1, p2, 1d0, 1d0, sideID, weight)
        endif

        ! if necessary, also find elements that have just one single point on the interface,
        ! and not a side in common with it.
        if (.not. sideonly) then
            if (p1match .gt. 0 .and.  p2match .eq. 0) then
                ! element only borders on the first point of the interval
                call elemCenter(iElem, centerCoords)
                call pnpoly(centerCoords,loopcoords,numloop,tmp,weight)
                if (iecho.eq.1) then
                    write(*,*) "element center coordinates: ",centerCoords
                endif
                if (weight.eq.0) then
                    write(0,*) "Element center of element ",iElem," lies on the loop"
                    stop "leaving elmside"
                endif

                if (weights(iElem) .eq. 0) then
                    weights(iElem) = weight
                else if (weights(iElem) .eq. weight) then
                    continue
                else
                    write(0,*) "CONFLICT: Element",iElem,"has node ", p1, &
                      "on the interface, but was already assigned an opposite weight"
                endif
                call output2D(iElem, p1, 0, 1d0, 1d0, 0, weight)
            else if (iInterval .eq. nvertices-1  .and. p1match .eq. 0  .and. p2match .gt. 0 ) then
                ! on the very last interval, we also need points that connect to the last point of the line
                call elemCenter(iElem, centerCoords)
                call pnpoly(centerCoords,loopcoords,numloop,tmp,weight)
                if (iecho.eq.1) then
                    write(*,*) "element center coordinates: ",centerCoords
                endif
                if (weight.eq.0) then
                    write(0,*) "Element center of element ",iElem," lies on the loop"
                    stop "leaving elmside"
                endif
                if (weights(iElem) .eq. 0) then
                    weights(iElem) = weight
                else if (weights(iElem) .eq. weight) then
                    continue
                else
                    write(0,*) "CONFLICT: Element",iElem,"has node ", p1, &
                      "on the interface, but was already assigned an opposite weight"
                endif
                call output2D(iElem, p2, 0, 1d0, 1d0, 0, weight)
            else
                ! the element does not have a point in common with this interface.
                continue
            endif
        endif
    enddo
enddo

deallocate(weights)
deallocate(tmp)

end subroutine
!-------------------------------------------------------------------------------
subroutine PNPOLY(P,LOOP,N,TMP,INOUT)
 
! Determines whether a point is inside a polygon 
! The vertices may be listed clockwise or anticlockwise.
! The first may optionally be repeated; if so N may 
! optionally be increased by 1. 
! The input polygon may be a compound polygon consisting 
! of several separate subpolygons. If so, the first vertex 
! of each subpolygon must be repeated, and when calculating 
! N, these first vertices must be counted twice. 
! INOUT is the only parameter whose value is changed. 
! written by Randolph Franklin, university of ottawa, 7/70. 
! Method: a vertical line is drawn thru the point in question. If it 
! crosses the polygon an odd number of times, then the 
! point is inside of the polygon.
 
IMPLICIT NONE
!-pass
INTEGER :: N       ! Number of vertices of polygon
INTEGER :: INOUT   ! Return value: -1 outside, 0 = on edge, 1 = inside
DOUBLE PRECISION :: LOOP(2,N) ! polygon coordinates
DOUBLE PRECISION :: P(2)      ! point coordinate
!-local
LOGICAL mx,my,nx,ny
INTEGER i,j
DOUBLE PRECISION :: TMP(2,N), testMe
 
do i=1,N
    tmp(1,i) = LOOP(1,i) - P(1)
    tmp(2,i) = LOOP(2,i) - P(2)
enddo

INOUT = -1 
do i=1,N
    j = 1 + MOD(i,N)
    MX = TMP(1,i).GE.0d0
    NX = TMP(1,j).GE.0d0 
    MY = TMP(2,i).GE.0d0
    NY = TMP(2,j).GE.0d0 

    if (.not.((my.or.ny).and.(mx.or.nx)).or.(mx.and.nx)) then
        cycle
    endif

    if (.not.(my.and.ny.and.(mx.or.nx).and..not.(mx.and.nx))) then
        goto 3 
    endif

    INOUT=-INOUT 
    cycle

3   testMe = (TMP(2,i)*TMP(1,j)-TMP(1,i)*TMP(2,j))/(TMP(1,j)-TMP(1,i))

    if (testme.lt.0d0) then
        cycle
    else if (testme.eq.0d0) then
        goto 4
    else
        goto 5
    endif

4   INOUT = 0
    return
5   INOUT = -INOUT
enddo

return
end subroutine
