program trialMesh2D

implicit none


! This subroutine creates the tecin.dat.partf.elm en tecin.dat.partf.nps
! to test the parallel numbering scheme, and the neighbor search.

! We force our own partition, because a METIS partitioning
! is freuently irregular in shape, and needs to be drawn out by hand, 
! which is a lot of work.

! The mesh is square, around the origin, with sides
! [-squareSide/2, squareSide/2].
! There is a circle section starting at 12 o'clock,
! and running in clockwise direction, of radius circleRadius
! and it spans an angle of circleSegmentAngle

! The points on the circle have marker 5
! The points on the surrounding square have marker 6


! parameters
double precision,  parameter :: squareSide = 2d0

double precision,  parameter :: circleRadius       = 0.5d0
double precision,  parameter :: circleSegmentAngle = 135d0  ! in degrees

integer,           parameter :: nPointsOnSegment    = 10

character(len=19), parameter :: domainFile   = "domain.poly"


double precision,  parameter :: pi = 3.1415926535

! derived parameters

double precision,  parameter :: halfSquareSide = squareSide * 0.5

! other values

integer :: iPoint
double precision             ::  angleProgress



open(unit=42, file="domain.poly")

write(42,*) "# points #"
write(42,*) nPointsOnSegment + 4," 2 0 1"
write(42,*) " 0: ", -halfSquareSide,  -halfSquareSide, 4
write(42,*) " 1: ",  halfSquareSide,  -halfSquareSide, 4
write(42,*) " 2: ",  halfSquareSide,   halfSquareSide, 4
write(42,*) " 3: ", -halfSquareSide,   halfSquareSide, 4

do iPoint=1,nPointsOnSegment
    angleProgress = circleSegmentAngle * pi / 180.0 * (iPoint-1) / (nPointsOnSegment-1)
    write(42,*)  3+iPoint,":",&
                 circleRadius * sin(angleProgress), &
                 circleRadius * cos(angleProgress), &
                 5
enddo

write(42,*) "# domain boundary #"
write(42,*) nPointsOnSegment + 3," 1"
! outer boundary

write(42,*) 0," :", 0, 1, 4
write(42,*) 1," :", 1, 2, 4
write(42,*) 2," :", 2, 3, 4
write(42,*) 3," :", 3, 0, 4

do iPoint=1,nPointsOnSegment-1
    write(42,*) 3+iPoint," :", iPoint+3, iPoint+4, 5
enddo

write(42,'(a)') "#"
write(42,'(a)') "# no gaps #"
write(42,'(a)') "0"
write(42,'(a)') "#"
write(42,'(a)') "# material markers #"
write(42,'(a)') "1"
write(42,'(a,2g14.6,a)') "0: ", 0, 0, "  1"

close(42)


end program
