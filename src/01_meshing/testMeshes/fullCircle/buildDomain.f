program buildDomain

! build full circle domain.poly.
! The central point (0,0) is generated
! and connected with a line to the North Pole (0,1000)

implicit none

double precision,  parameter :: pi = 3.1415926535
double precision,  parameter :: radius = 1000d0
integer,           parameter :: nCirclePoints = 400


double precision             :: angle, x, y
integer                      :: iPoint

open(unit=42, file="domain.poly")


write(42,*) "# points #"
write(42,*) nCirclePoints+1, " 2 0 1"
write(42,*) 0,0,0

do iPoint = 1, nCirclePoints

    angle = 2.0 * pi * dble(iPoint) / (dble(nCirclePoints))

    x = radius * sin(angle)
    y = radius * cos(angle)

    write(42,*) iPoint, x, y, 99

enddo

write(42,*) "# domain boundary #"
write(42,*) nCirclePoints+1, 1


write(42,*) "0: 0", nCirclePoints, "1" 
write(42,*) "1: 1", nCirclePoints, 99


do iPoint = 1, nCirclePoints-1
	write(42,*) iPoint+1, ":",iPoint, iPoint+1,99
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
