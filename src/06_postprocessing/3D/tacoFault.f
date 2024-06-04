subroutine LOADFT (IEN,NFAULT,FAULT,TFAULT,NDOF,NUMFN,NEN,NUMEL)

USE TIMESTEPMODULE,  only: nstep
USE ALGEBRA
use modeldatamodule, only: modeldatactx
!
!    set total faulted node displacements TFAULT 
!
!    logical VELBC (also FORMF and RFAULT) is used to change interpretation
!    of  NFAULT(3,.) as follows:
!    VELBC = .true.: NFAULT(3,.) < 0 means constant velocity
!    VELBC = .false.: NFAULT(3,.) < 0 displacement @ time step ABS(NFAULT(3,.)
!
! taco broerse 2016
implicit none
 
! following has been copied from LOADF, but its meanings are unknown, therefore
! commented out
! logical VELBC
!#ifdef FAULTVELO
!parameter (VELBC=.true.) ! consistent in FORMF and RFAULT
!#else
!parameter (VELBC=.false.) ! consistent in FORMF and RFAULT
!#endif

!-pass
integer NDOF,NUMFN,NUMEL,NEN
integer NFAULT,IEN
double precision FAULT,TFAULT
dimension NFAULT(3,*),FAULT(NDOF,*),TFAULT(NDOF,*),IEN(NEN,NUMEL)
!-locl
integer j,k,i
integer nodeID,elemID,vertexID
logical NoMatch
#include "io.i"
! check if there are faulted nodes
if (NUMFN.le.0) then
    return
endif
!
call CLEAR(TFAULT,NDOF*NUMFN)
!
! write(stderr,*) 'numfn in loadft',numfn
do k=1,NUMFN
  
    if (NFAULT(3,k).lt.0) then
        write(stderr,*) 'application mode not catered for currently'
        write(stderr,*) 'does not handle cyclic fault displacements'
        write(stderr,*) k,NFAULT(:,k)
        return
        ! if this subroutine proves valuable, cyclic fault displacements may
        ! be introduced as well
    endif

    ! if time of application is equal or less than current time, add to TFAULT
    if (NFAULT(3,k).le.NSTEP) then

		! check on whether node is in element
	    NoMatch = .true.
	    ! node number
	    nodeID = NFAULT(2,k)
	    ! element number
	    elemID = NFAULT(1,k)
	    ! now find which vertex (1,2,..,NEN) this is of the current element
	    do j=1,NEN
	        if (IEN(j,elemID).eq.nodeID) then
	            vertexID = j
	            NoMatch = .false.
	        endif
	    enddo

	    if (NoMatch) then
	        write(*,*) 'LOADFT: node ',nodeID,' doesn''t occur in', &
	                   'element ',elemID
	        call exitp(1)
	    endif
    
	    ! this index is identical to the last entry of NFAULT with this element-node
	    ! combination, and is consistent with AddFaultDisplacement subroutine
	    i = modeldatactx%LMF(1,vertexID,elemID)
	    do j=1,NDOF
	        TFAULT(j,i) = TFAULT(j,i) + FAULT(j,k)
	    enddo
    endif

enddo

return
end subroutine
!-------------------------------------------------------------------------------

