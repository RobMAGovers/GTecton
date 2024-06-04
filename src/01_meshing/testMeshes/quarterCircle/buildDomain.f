program buildDomain

! build quarter circle domain.poly

implicit none

double precision,  parameter :: pi = 3.1415926535
double precision,  parameter :: radius = 1d0
integer,           parameter :: nCirclePoints = 5


double precision             :: angle, x, y
integer                      :: iPoint

open(unit=42, file="domain.poly")

write(42,*) "# points #"
write(42,*) nCirclePoints + 1, " 2 0 1"
write(42,*) 0,0,0

do iPoint = 1, nCirclePoints

    angle = 0.5 * pi * dble(iPoint-1) / (dble(nCirclePoints-1))

    x = radius * sin(angle)
    y = radius * cos(angle)

    write(42,*) iPoint, x, y, 1000
enddo

write(42,*) "# domain boundary #"
write(42,*) nCirclePoints+1, 1

write(42,*) "0: 0", nCirclePoints, "1" 
do iPoint = 1, nCirclePoints          
	write(42,*) iPoint, ":",iPoint-1, iPoint, "1"
enddo

write(42,'(a)') "#"
write(42,'(a)') "# no gaps #"
write(42,'(a)') "0"
write(42,'(a)') "#"
write(42,'(a)') "# material markers #"
write(42,*) 1
write(42,*) "0:", 0.5 * radius, 0.5 * radius, 1

close(42)

end program
