program setmaterials

implicit none
character (len=1) :: a
integer :: NUMNP,i,node,label
double precision, allocatable :: X(:,:),XELM(:,:)
double precision :: E,v
integer :: NUMEL,element
integer, allocatable :: MAT(:),IEN(:,:)

interface 
    subroutine PREM(MAT,XELM,E,v)
        integer, intent(in) :: MAT
        double precision, intent(in) :: XELM(3)
        double precision, intent(out) :: E,v
    end subroutine PREM

    subroutine ELCTR(X,IEN,XELM)
        integer, intent(in) :: IEN(4)
        double precision, intent(in) :: X(3,*)
        double precision, intent(out) :: XELM(3)
    end subroutine
end interface

open (unit=10,file='tecin.dat.nps',status='old')
rewind (10)

i = 0
do
   read(10,'(a1)',end=100,err=1000) a
   i = i + 1
enddo
100 NUMNP = i-1
allocate (X(3,NUMNP))
rewind(10)
do i=1,NUMNP
    read(10,1,err=1000,end=1100) node,label,X(1,i),X(2,i),X(3,i)
1   format(2I12,3E25.0)
enddo
close(10)

open (unit=10,file='tecin.dat.elm',status='old')
rewind (10)
i = 0
do
   read(10,'(a1)',end=200,err=2000) a
   i = i + 1
enddo
200 NUMEL = i-1
allocate (IEN(4,NUMEL))
allocate (MAT(NUMEL))
rewind(10)
do i=1,NUMEL
    read(10,2,err=2000,end=2100) element,MAT(i),IEN(1,i),IEN(2,i),IEN(3,i),IEN(4,i)
2   format(6I12)
enddo
close(10)

allocate (XELM(3,NUMEL))
do i=1,NUMEL
    call ELCTR(X,IEN(1,i),XELM(1,i))
enddo

deallocate (X)
deallocate (IEN)

open (unit=10,file='tecin.dat.elasticity')
rewind(10)

! assign element property in subroutine (potentially based on material number)
do i=1,NUMEL
    call PREM(MAT(i),XELM(1,i),E,v)
    write(10,3) MAT(i),E,v
3   format(I12,2E25.12)
enddo

close(10)
deallocate (MAT)
deallocate (XELM)
stop

1000 write(0,*) 'setmaterials: read error tecin.dat.nps record=',i
stop
1100 write(0,*) 'setmaterials: EOF error tecin.dat.nps'
stop
2000 write(0,*) 'setmaterials: read error tecin.dat.elm record=',i
stop
2100 write(0,*) 'setmaterials: EOF error tecin.dat.elm'
stop
end program
!-------------------------------------------------------------------------------
subroutine PREM(MAT,XELM,E,v)

! Compute elastic Young's modulus and Poisson ratio from PREM. Material number is irrelevant.
! Assumption is that the third argument of XELM represents distance (height or depth)
! in meters from the free surface.

implicit none
integer, intent(in) :: MAT
double precision, intent(in) :: XELM(3)
double precision, intent(out) :: E,v    ! Young's modulus (Pa) and Poisson ratio

integer :: i
double precision :: Depth(13),BulkModulus(13),ShearModulus(13)
double precision :: distance,dKdz,Ki,dGdz,Gi

! initialize PREM depth (km), K (GPa), G (GPa)
Depth = (/ 0d0, 12d0, 40d0, 68d0, 170d0, 200d0, 220d0, 271d0, 371d0, 400d0, &
  471d0, 571d0, 670d0 /) 
BulkModulus = (/ 52d0, 75.3d0, 131.5d0, 130.4d0, 128.1d0, 128.1d0, 152.9d0, &
  158.6d0, 170.1d0, 189.9d0, 209.7d0, 239.7d0, 255.6d0 /)
ShearModulus = (/ 26.6d0, 44.1d0, 68.2d0, 67.6d0, 66.2d0, 66.2d0, 74.1d0, &
  75.9d0, 79.5d0, 90.6d0, 100.7d0, 116.2d0, 123.9d0 /)

distance = ABS(XELM(3)*1D-3)

i = 1
do while (i.lt.13)
    if ( (Depth(i)-distance)*(Depth(i+1)-distance).le.0d0 ) then
        dKdz = (BulkModulus(i+1)-BulkModulus(i))/(Depth(i+1)-Depth(i))
        Ki = BulkModulus(i) + dKdz*(distance-Depth(i))
        dGdz = (ShearModulus(i+1)-ShearModulus(i))/(Depth(i+1)-Depth(i))
        Gi = ShearModulus(i) + dGdz*(distance-Depth(i))
        E = 9d0*Ki*Gi/(3d0*Ki+Gi) * 1d9
        v = (3d0*Ki-2d0*Gi)/(6d0*Ki+2d0*Gi)
        return
    endif
    i = i + 1
enddo

write(0,*) 'PREM: ',distance,' km not in table'
stop

end subroutine
!-------------------------------------------------------------------------------
subroutine ELCTR(X,IEN,XELM)

! compute element center coordinate

implicit none
integer, intent(in) :: IEN(4)
double precision, intent(in) :: X(3,*)
double precision, intent(out) :: XELM(3)
double precision ::XL(3,4)
integer :: k

k = IEN(1)
XL(1,1) = X(1,k)
XL(2,1) = X(2,k)
XL(3,1) = X(3,k)
k = IEN(2)
XL(1,2) = X(1,k)
XL(2,2) = X(2,k)
XL(3,2) = X(3,k)
k = IEN(3)
XL(1,3) = X(1,k)
XL(2,3) = X(2,k)
XL(3,3) = X(3,k)
k = IEN(4)
XL(1,4) = X(1,k)
XL(2,4) = X(2,k)
XL(3,4) = X(3,k)
XELM(1) = (xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))*0.25d0
XELM(2) = (xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))*0.25d0
XELM(3) = (xl(3,1)+xl(3,2)+xl(3,3)+xl(3,4))*0.25d0

return
end subroutine
