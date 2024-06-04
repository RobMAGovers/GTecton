! first fill the ElementNeighbours
! then this function can be called.

subroutine findBoundaries(nodeSpace, fromNode, toNode)

use meshdatamodule, only: meshdatactx                  

implicit none

integer :: nodeSpace
integer :: iNode, iBoundary, iSide, iElem
integer :: checkElem
integer :: v1, v2, v3
integer :: fromNode(nodeSpace), toNode(nodeSpace)



iBoundary = 1

do iElem=1,meshdatactx%neglobal

    v1 = meshdatactx%ien(1,iElem)
    v2 = meshdatactx%ien(2,iElem)
    v3 = meshdatactx%ien(3,iElem)


    do iSide = 1,3
        checkElem = meshdatactx%ElementNeighbours(iElem, iSide)
        if (checkElem.eq.0) then
            ! we have an element that borders on the edge of the domain.
            ! include it anyway, because it also the edge of a partition.
            if (iSide.eq.1) then
                ! boundary on side 1, node 1 and 2
                fromNode(iBoundary) = v1
                  toNode(iBoundary) = v2
            else if (iSide.eq.2) then
                 ! boundary on side 2, node 2 and 3
                 fromNode(iBoundary) = v2
                   toNode(iBoundary) = v3
            else if (iSide.eq.3) then
                 ! boundary on side 3, node 1 and 3
                 fromNode(iBoundary) = v1
                   toNode(iBoundary) = v3
            endif
            iBoundary = iBoundary + 1
        else
            ! not on the domain boundary, we have to check partitions
            if (meshdatactx%e2p(iElem) .ne. meshdatactx%e2p(checkElem)) then
                ! we have two adjacent elements of different partition
                if (iSide.eq.1) then
                    ! boundary on side 1, node 1 and 2
                    fromNode(iBoundary) = v1
                      toNode(iBoundary) = v2
                else if (iSide.eq.2) then
                     ! boundary on side 2, node 2 and 3
                     fromNode(iBoundary) = v2
                       toNode(iBoundary) = v3
                else if (iSide.eq.3) then
                     ! boundary on side 3, node 1 and 3
                     fromNode(iBoundary) = v1
                       toNode(iBoundary) = v3
                endif
                iBoundary = iBoundary + 1
            endif
        endif

        if (iBoundary.gt.nodeSpace) then
            write(*,*) "Error: More boundaries than allocated space"
            return
        endif

    enddo
enddo


end subroutine findBoundaries
