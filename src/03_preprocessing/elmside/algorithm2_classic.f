subroutine algo2()

use globals

implicit none

integer :: iElem, iPoint, iInterval

integer :: matchID(3)
logical :: match(3)

integer :: sideID

double precision :: centerCoords(2)

double precision, external :: crossProduct
integer, external :: numInArray
integer, external :: setWeight

integer :: p1match, p2match

double precision :: vecToCenter(2), vecAlongInterface(2)
integer          :: weight
integer          :: p1, p2

integer, allocatable :: weights(:)

allocate(weights(nelems))
weights = 0

if (ndimensions .eq. 3) then
    write(*,*) "algorithm 2 only works in two dimensions."
else
    if (iecho.eq.1) then
        write(*,*) "algo2 says: hello world"
    endif
endif


! for this algorithm,
! the nodelist is ordered spatially using the polysort tool.
! we step through the nodes


    ! again let us rephrase the convention for sides
    ! 2D:
    !   face 1: vertices 1 2
    !        2:          2 3
    !        3:          3 1 

if (iecho.eq.1) then
    write(*,*) "n intervals", nvertices-1
endif

do iInterval = 1, nvertices-1


    p1 = vertices(iInterval)
    p2 = vertices(iInterval+1)


    if (iecho.eq.1) then
        write(*,*) "----------------------------------------------------------------------------"
        write(*,*) "----------  Next interface section, from ", p1, "to", p2,"-----------"
        write(*,*) "----------------------------------------------------------------------------"
    endif

    do iElem = 1, nelems

        p1match = numInArray(p1, v(1:3,iElem),3)
        p2match = numInArray(p2, v(1:3,iElem),3)

!        write(*,*) "check",iElem, &
!                   "with verts", v(1:3,iElem),&
!                   " on interface, with: ", p1, p2


        if (iecho.eq.1 .and. (p1match .gt. 0 .or. p2match .gt. 0)) then
            write(*,*) "got element",iElem," on interface, with: ", p1match, p2match, &
                       "which are: ", p1, p2

        endif

        ! find elements that have a side on the edge. Determine the edge and sign
        if (p1match .gt. 0 .and. p2match .gt. 0) then

            if (weights(iElem) .eq. 0) then

                if (iecho.eq.1) then
                    write(*,*) "Element",iElem, &
                               "has the side between", p1, p2, &
                               "on the interface"
                endif


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
                    write(*,*) "unknown side ID", p1match, p2match
                    write(*,*) "This is not supposed to happen. Please contact model support."
                    stop "Leaving elmside..."
                endif
    
                call elemCenter(iElem, centerCoords)
    
                vecToCenter = centerCoords - coords(:,p1)
                vecAlongInterface = coords(:,p2) - coords(:,p1)
    
                weight = setWeight(crossProduct(vecToCenter, vecAlongInterface))
                weights(iElem) = weight

                call output2D(iElem, p1, p2, 1d0, 1d0, sideID, weight)

            else

                if (iecho.eq.1) then
                    write(*,*) "Element",iElem, &
                               "has the side between", p1, p2, &
                               "on the interface, but was already assigned weight", weights(iElem)
                endif


            endif

        endif

        ! if necessary, also find elements that have just one single point on the interface,
        ! and not a side in common with it.
        if (.not. sideonly) then

            if (p1match .gt. 0 .and.  p2match .eq. 0) then

                if (weights(iElem) .eq. 0) then

                    if (iecho.eq.1) then 
                        write(*,*) "Element",iElem, &
                                        "has    just a point", p1, &
                                        "on the interface"
                       endif
    
                    ! element only borders on the first point of the interval
                    call elemCenter(iElem, centerCoords)

                    vecToCenter    = centerCoords - coords(:,p1)
                    vecAlongInterface =    coords(:,p2) - coords(:,p1)

                    weight = setWeight(crossProduct(vecToCenter, vecAlongInterface))
                    weights(iElem) = weight

                    call output2D(iElem, p1, 0, 1d0, 1d0, 0, weight)
                else
    
                    if (iecho.eq.1) then
                        write(*,*) "Element",iElem, &
                               "has a single point ", p1, &
                               "on the interface, but was already assigned weight", weights(iElem)
                    endif


                endif



            else if (iInterval .eq. nvertices-1  .and. p1match .eq. 0  .and. p2match .gt. 0 ) then

                if (weights(iElem) .eq. 0) then


                    if (iecho.eq.1) then
                        write(*,*) "Element",iElem, &
                               "has just a point", p2, &
                               "on the last point of the interface"
                    endif

                    ! on the very last interval, we also need points that connect to the last point of the line
                    call elemCenter(iElem, centerCoords)

                    vecToCenter    = centerCoords - coords(:,p1)
                    vecAlongInterface =    coords(:,p2) - coords(:,p1)

                    weight = setWeight(crossProduct(vecToCenter, vecAlongInterface))
                    weights(iElem) = weight

                    call output2D(iElem, p2, 0, 1d0, 1d0, 0, weight)

                else
    
                    if (iecho.eq.1) then
                        write(*,*) "Element",iElem, &
                               "has a single point ", p1, &
                               "on the last point of the interface but was already assigned weight", weights(iElem)
                    endif

                endif


            else
                ! the element does not have a point in common with this interface.
            endif

        endif



    enddo

enddo


end subroutine


subroutine elemCenter(elemID, center)

use globals

implicit none

integer          :: elemID
double precision :: center(2)

! take the average of three points of triangle.
center(:) = (coords(1:2,v(1,elemID)) + &
             coords(1:2,v(2,elemID)) + &
             coords(1:2,v(3,elemID))) / 3d0

end subroutine


function crossProduct(vec1, vec2)

implicit none

double precision :: vec1(2), vec2(2)
double precision :: crossProduct

crossProduct = vec2(2) * vec1(1) - vec2(1) * vec1(2)

end function

function setWeight(double)
! native function 'sign' needs same variable type to set sign,
! and here we set sign from a real to an integer, hence separate function
implicit none

integer :: setWeight
double precision :: double

if (double.gt.0d0) then
    setWeight = 1
else
    setWeight =    -1
endif

end function



