program testDistance()

implicit none




end program


double precision function DistancePointAndLine(coordsPoint, a,b,c, nDimensions)

imlicit none

! computes the distance between point coordsPoint and line with a 
! basepoint of baseLine and a direction of vectorLine.

! the line is described as y = ax + b
! a plane is described as z = ax + by + c
! nDimensions is 2 or 3


integer :: nDimensions
double precision :: coordsPoint(nDimensions)
double precision :: a,b,c

if (nDimensions.eq.2) then
	
else if (nDimensions.eq.3) then
	! we follow https://mathinsight.org/distance_point_plane
    ! accessed on 16 jun 2020

    DistancePointAndLine = 

else
    stop "nDimensions can only be 2 or 3, but is"// str(nDimensions)
endif

end function

