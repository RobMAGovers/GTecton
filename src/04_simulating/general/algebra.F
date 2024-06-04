MODULE ALGEBRA

!******************************************************************
! This model contains basic algebraic operations
! that do not depend on geophysical context.
! Here you will find vector and matrix operations and the like.
!
! This module does only depend on the constants module, but many depend
! on it, hence it is best compiled first in the module chain.
!
! Note that it does not include the parallel PETSC algebra
! Here are only mere mortal serial routines
!******************************************************************

implicit none

! invert a 3x3 matrix
PUBLIC :: INVMT3  
! compute a cross product
PUBLIC :: CROSSP

! construct a rotation matrix, both for 2D and 3D
PUBLIC :: FORMRT
public :: rotatePoint

PUBLIC :: SWAP
PUBLIC :: LINLIN
PUBLIC :: TENLOG

! multiply two matrices
PUBLIC :: MATPRD
PUBLIC :: PUTMAT
PUBLIC :: PUTVEC
! vector matrix product
PUBLIC :: VMPRD
! copy lower triangle into upper triangle,
! making matrix symmetrical
PUBLIC :: SYMMET
! transpose a matrix
PUBLIC :: TRANSP
! compute dot product
PUBLIC :: DOT
! invert any matrix
!PUBLIC :: INVMAT   ! TODO This and the solve routines should move from elastc to here.
! find position of diagonal elements when matrix is
! represented as a vector
PUBLIC :: IDFORM
! rotate a basis around the Z/X-axis
PUBLIC :: rotz
PUBLIC :: rotx
! rotate a vector by multiplication with a rotation matrix
PUBLIC :: rotvec


! fill a vector with zeroes
PUBLIC :: CLEAR

! TODO: encapsulate those two in an interface
PRIVATE :: CLEAR_double3D
PRIVATE :: CLEAR_double2D
PRIVATE :: CLEAR_double1D
PRIVATE :: CLEAR_real2D  ! specifically for the XPT array which holds single precision
PRIVATE :: CLEAR_real1D
PRIVATE :: CLEAR_int3D
PRIVATE :: CLEAR_int2D
PRIVATE :: CLEAR_int1D


! clear a matrix       
public :: clear_matrix
public :: copy_matrix

! add two vectors
PUBLIC :: UPDATE
PUBLIC :: UPDATE_real
! TODO: encapsulate those two in an interface
!        PRIVATE :: UPDATE_real
!        PRIVATE :: UPDATE_double

! Copy entries of one vector into another
PUBLIC :: EQUATE
! TODO: encapsulate those two in an interface
PRIVATE :: EQUATE_int3D
PRIVATE :: EQUATE_int2D
PRIVATE :: EQUATE_int1D
PRIVATE :: EQUATE_double3D
PRIVATE :: EQUATE_double2D
PRIVATE :: EQUATE_double1D

! multiply... thingy...
PUBLIC :: BMULT
PUBLIC :: BMULT_real
! TODO: encapsulate those two in an interface
!        PRIVATE :: BMULT_real
!        PRIVATE :: BMULT_double

public :: sort4vec

public :: determinant3x3
public :: determinant4x4


public :: arraycontains


public :: twoPointsOnSameSideOfPlane

public :: vectorLength
public :: areaOfTriangleIn3D

public :: distanceBetweenTwoPoints
public :: distanceBetweenPointAndPlane
public :: distanceBetweenTwoPointsOnASphere

public :: radiusOfSphereInsideTetrahedron
public :: radiusOfSphereAroundTetrahedron

#ifdef SPARSE
public :: globalSumSquare
#endif

public :: angleBetweenTwoVectors3D

!interface formrt
!    module procedure formrt2d, formrt3d
!end interface


! The module procedure CLEAR_double can only take
! and argument of one dimensionality.
! So to handle all the various array shapes, multiple
! procedures are needed.

! This seems like an unnecessary hassle.
! The great advantage is that it become possible to tag
! a certain call with an extra parameter, and only needing
! to make an extra procedure, without having to change
! every single to any clear function.

! This has already been used. :-)

interface CLEAR
    module procedure CLEAR_double3D, &
                     CLEAR_double2D, &
                     CLEAR_double1D, &
                     CLEAR_real2D, &
                     CLEAR_real1D, &
                     CLEAR_int3D, &
                     CLEAR_int2D, &
                     CLEAR_int1D

end interface

interface EQUATE
    module procedure EQUATE_double3D, &
                     EQUATE_double2D, &
                     EQUATE_double1D, &
                     EQUATE_int3D, &
                     EQUATE_int2D, &
                     EQUATE_int1D

end interface 

! interface UPDATE
!    module procedure UPDATE_double, UPDATE_real
! end interface 

! interface BMULT
!    module procedure BMULT_double, BMULT_real
! end interface 


!******************************************************************
contains

double precision function angleBetweenTwoVectors3D(vec1, vec2)

implicit none

double precision :: vec1(3), vec2(3)

!write(*,*) "a ", vec1
!write(*,*) "b ", vec2
!write(*,*) "c ", dot(vec1, vec2, 3)
!write(*,*) "d ", vectorLength(vec1, 3)
!write(*,*) "e ", vectorLength(vec2, 3)
!write(*,*) "f ", vectorLength(vec1, 3) * vectorLength(vec2, 3)
!write(*,*) "g ",  1d0 / vectorLength(vec1, 3) * vectorLength(vec2, 3)
!write(*,*) "h " 

angleBetweenTwoVectors3D = acos((dot(vec1, vec2, 3)) / & 
        (vectorLength(vec1, 3) * vectorLength(vec2, 3)))

end function

!-----------------------------------------------------------------------


#ifdef SPARSE
! sparse because the MPI call
double precision function globalSumSquare(vector, nLocal)

use petscksp  ! for the MPI_Allreduce

implicit none

integer          :: nLocal
double precision :: vector(nLocal)
double precision :: localSumSquare

integer          :: iEntry, iError

localSumSquare = 0d0
do iEntry=1,nLocal
    localSumSquare = \
    localSumSquare + vector(iEntry) * vector(iEntry)
enddo

call mpi_Allreduce(localSumSquare, &   ! send
                   globalSumSquare, &  ! receive
                   1,&                 ! number of variables in array
                   MPI_double, &       ! type
                   MPI_sum, &          ! operator
                   MPI_COMM_WORLD, &
                   iError)
end function
# endif

!-----------------------------------------------------------------------


subroutine INVMT3 (A,IERR)

use constants, only: eps

implicit none
!
! inverts 3x3 matrix A (Cramers rule)
!
double precision, parameter :: ONE=1D0
integer :: IERR
!-pass
 double precision :: A(3,3)
!-locl
 double precision Xi11,Xi21,Xi31,Xi12,Xi22,Xi32,Xi13,Xi23,Xi33,det
!
det = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2) &
     -A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2) &
     +A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2)

if (ABS(det).le.EPS) then
    IERR = 1
    return
endif
!
det=ONE/det

Xi11=(A(2,2)*A(3,3)-A(3,2)*A(2,3))*det
Xi21=(A(2,3)*A(3,1)-A(2,1)*A(3,3))*det
Xi31=(A(2,1)*A(3,2)-A(2,2)*A(3,1))*det
Xi12=(A(1,3)*A(3,2)-A(1,2)*A(3,3))*det
Xi22=(A(1,1)*A(3,3)-A(1,3)*A(3,1))*det
Xi32=(A(1,2)*A(3,1)-A(1,1)*A(3,2))*det
Xi13=(A(1,2)*A(2,3)-A(1,3)*A(2,2))*det
Xi23=(A(1,3)*A(2,1)-A(1,1)*A(2,3))*det
Xi33=(A(1,1)*A(2,2)-A(1,2)*A(2,1))*det

A(1,1)=Xi11
A(2,1)=Xi21
A(3,1)=Xi31
A(1,2)=Xi12
A(2,2)=Xi22
A(3,2)=Xi32
A(1,3)=Xi13
A(2,3)=Xi23
A(3,3)=Xi33
!
IERR = 0
return
end subroutine

!-----------------------------------------------------------------------

subroutine rotatePoint(angles, ndof, original, rotated)

implicit none

integer          :: ndof
double precision :: angles(2)
double precision :: matrix(ndof, ndof)

double precision :: original(ndof)
double precision :: rotated(ndof)

call FORMRT (angles, matrix,ndof)
call VMPRD(matrix,original,rotated,ndof,ndof)

end subroutine

!-----------------------------------------------------------------------

subroutine FORMRT (SKEW,ROT,ndof)
!
!    Constructs a two or three dimensional rotation matrix ROT
!    from the Euler angles contained in SKEW
!
implicit none
!-pass
integer          :: ndof
double precision :: SKEW(2), ROT(3,3)
!-local
double precision :: cosb, sinb

ROT(1,1) = COS(SKEW(1))
ROT(1,2) = SIN(SKEW(1))
ROT(2,1) = -ROT(1,2)
ROT(2,2) =  ROT(1,1)

! if 2d, then skew(2) = 0, so COSB = 1 and SINB=0
! making the 3D complement just as if it were in 2D.
SINB = SIN(SKEW(2))
COSB = COS(SKEW(2))
ROT(3,1) = -ROT(1,1)*SINB
ROT(1,1) =  ROT(1,1)*COSB
ROT(1,2) =  ROT(1,2)*COSB
ROT(1,3) =  SINB
ROT(3,2) =  ROT(2,1)*SINB
ROT(2,3) =  0d0
ROT(3,3) =  COSB

return

end subroutine


!-------------------------------------------------------------------------------

subroutine CROSSP (A,B,C)
!
! computes the cross product C = A x B
!
implicit none

!-pass
double precision :: A(3),B(3),C(3)

C(1) = A(2)*B(3)-A(3)*B(2)
C(2) = A(3)*B(1)-A(1)*B(3)
C(3) = A(1)*B(2)-A(2)*B(1)

return
end subroutine
!-------------------------------------------------------------------------------
subroutine SWAP (A,B,N)
!
implicit none
!-pass
double precision A(*),B(*)
!-locl
integer i,N
double precision tmp
!
do i=1,N
    tmp = A(i)
    A(i) = B(i)
    B(i) = tmp
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine LINLIN (x1,y1,x2,y2,a,b,mode)

use constants, only: dfmin

! Calculates line through x1,y1 and x2,y2. Mode=0: y=ax+b, Mode<>0: x=b

implicit none

integer          :: mode
double precision :: x1, y1
double precision :: x2, y2
double precision :: a , b

double precision :: dx

dx = x2 - x1
if (ABS(dx).le.dfmin) then
    mode = 1
    b = x1
else
    mode = 0
    a = (y2-y1)/dx
    b = (y1*x2-y2*x1)/dx
endif
return

end subroutine 
!-------------------------------------------------------------------------------
double precision function TENLOG(x)

use constants, only: dfmin,dfmax,dlogmax
use iomodule,  only: stderr

implicit none

! takes the 10Log of x

double precision :: x, zero

!-init
zero = 0d0

if (x.lt.zero) then
    write(stderr,1)
 1       format(1x,'LOG10: negative argument error')
    TENLOG = dfmax
elseif (x.le.dfmin) then
    TENLOG = -dlogmax
else
    TENLOG = LOG10(x)
endif

return
end function
!-----------------------------------------------------------------------

subroutine BMULT (B,N,F)
!
use iomodule, only: stderr

implicit none

integer N,i
double precision :: B(N)
double precision :: F
!
if (N.eq.0) then
    write(stderr,*) "Vulnerability: BMULT wants to scale a vector of length 0."
    write(stderr,*) "Please contact model support."
endif

do i=1,N
    B(i) = B(i) * F
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine BMULT_real (B,N,F)
!
use iomodule, only: stderr

implicit none
integer :: N,i
real :: F, B(N)

if (N.eq.0) then
    write(stderr,*) "Vulnerability: BMULT_real wants to scale a vector of length 0."
    write(stderr,*) "Please contact model support."
endif

do i=1,N
    B(i) = B(i) * F
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine CLEAR_int1D (A,N, arrayName)

use iomodule, only: stderr

implicit none

integer, parameter :: ZERO=0
!-pass
integer N
integer :: A(N)
character(len=*) arrayName
!dimension A(N)
!-locl
integer i

if (N.eq.0) then
    write(stderr,*) "Vulnerability: CLEAR_int wants to clear array ", trim(arrayName), "."
    write(stderr,*) "This array has length 0. This is not supposed to happen"
    write(stderr,*) "Please contact model support."
endif

do i=1,N
    A(i) = ZERO
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine CLEAR_int2D (A,N, arrayName)
!
implicit none

!-pass
integer N
integer :: A(:,:)
character(len=*) arrayName


call CLEAR_int1D(A,N, arrayName)

return
end subroutine
!-----------------------------------------------------------------------
subroutine CLEAR_int3D (A,N, arrayName)
!
implicit none
!-pass
integer N
integer :: A(:,:,:)
character(len=*) arrayName

call CLEAR_int1D(A,N, arrayName)

return
end subroutine
!-----------------------------------------------------------------------
subroutine CLEAR_double1D (A,N, arrayName)
!
use iomodule, only: stderr

implicit none
double precision, parameter :: ZERO=0d0
!-pass
integer N
double precision :: A(N)
character(len=*) arrayName
! dimension A(*)
!-locl
 integer i
!
if (N.eq.0) then
    write(stderr,*) "Vulnerability: CLEAR_double wants to clear array ", trim(arrayName), "."
    write(stderr,*) "This array has length 0. This is not supposed to happen"
    write(stderr,*) "Please contact model support."
endif

do i=1,N
    A(i) = ZERO
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine CLEAR_double2D (A,N, arrayName)
!
implicit none
!-pass
integer :: N
double precision :: A(:,:)
character(len=*) arrayName

call CLEAR_double1D(A,N, arrayName)

end subroutine
!-----------------------------------------------------------------------
subroutine CLEAR_double3D (A,N, arrayName)
!
implicit none
!-pass
integer :: N
double precision :: A(:,:,:)
character(len=*) arrayName

call CLEAR_double1D(A,N, arrayName)

end subroutine
!-----------------------------------------------------------------------
subroutine CLEAR_real1D (A,N, arrayName)

use iomodule, only: stderr

implicit none
real, parameter :: ZERO=0.0
!-pass
integer N
real :: A(N)
character(len=*) arrayName
! dimension A(*)
!-locl
 integer i
!
if (N.eq.0) then
    write(stderr,*) "Vulnerability: CLEAR_real wants to clear array ", trim(arrayName), "."
    write(stderr,*) "This array has length 0. This is not supposed to happen"
    write(stderr,*) "Please contact model support."
endif

do i=1,N
    A(i) = ZERO
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine CLEAR_real2D (A,N, arrayName)
!
implicit none
!-pass
integer :: N
real :: A(:,:)
character(len=*) arrayName

call CLEAR_real1D(A,N, arrayName)

end subroutine

!-------------------------------------------------------------------------------
subroutine MATPRD (a,b,c,l,m,n)
!
! performs matrix product C(L,N) = sum over M { A(L,M)*B(M,N) }
!
use iomodule, only: stderr

implicit none

integer          :: l,m,n, i,j,k
double precision :: a(l,m), b(m,n), c(l,n)

if (l.eq.0 .or. m.eq.0 .or. n.eq.0) then
    write(stderr,*) "Vulnerability: MATPRD wants to multiply a matrices of size 0."
    write(stderr,*) "Please contact model support."
endif

do i=1,l
    do j=1,n
        c(i,j) = 0d0
        do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
        enddo
    enddo
enddo

return
end subroutine
!-------------------------------------------------------------------------------
 subroutine PUTMAT (A,N,M,lu)

 implicit none
!-pass
 integer N,M, lu, i, j
double precision :: A(N,M)
!-locl
 character(len=20) f

 write(f,1) M
 1   format('(1x,',I3,'(1PG20.8))')
 write(lu,fmt=f) ((A(i,j),j=1,M),i=1,N)

 return

 end subroutine
!-------------------------------------------------------------------------------
 subroutine PUTVEC (A,N,lu)

 implicit none
!-pass
 integer N, lu, i
double precision :: A(N)
!-locl
 character(len=20) f

 write(f,1) N
 1   format('(1x,',I3,'(1PG20.8))')
 write(lu,fmt=f) (A(i),i=1,N)
 return
 end subroutine
!-----------------------------------------------------------------------
subroutine VMPRD(A,B,C,NDIM,NA)
!
!    computes the product of matrix A and column vector B of height NDIM:
!    C(I) = A(I,J)*B(J)
!
use iomodule, only: stderr

implicit none

integer                      :: NDIM, NA, i,j
double precision,  parameter :: ZERO=0d0
double precision             :: A(NA,NA),B(NDIM),C(NDIM)

double precision             :: sum

if (NDIM.eq.0) then
    write(stderr,*) "Vulnerability: MATPRD wants to multiply a matrix or vector of length 0."
    write(stderr,*) "Please contact model support."
endif

!write(*,*) "vmprd: A", NDIM, A
!write(*,*) "vmprd: A bits:", A(1,1), A(1,2), A(2,1), A(2,2)


do i=1,NDIM
    sum = ZERO
    do j=1,NDIM
        sum = sum + A(i,j)*B(j)
    enddo
    C(i) = sum
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine SYMMET(S,NEE)
!
!    subroutine to symmetrize a matrix S
!
use iomodule, only: stderr

implicit none

integer          :: NEE, i,j
double precision :: S(NEE,NEE)

if (NEE.eq.0) then
    write(stderr,*) "Vulnerability: SYMMET wants to edit a matrix of size 0."
    write(stderr,*) "Please contact model support."
endif


do i=2,NEE
    do  j=1,i-1
        S(i,j) = S(j,i)
    enddo
enddo
 
return
end subroutine
!-----------------------------------------------------------------------
subroutine TRANSP(AMAT,NDIM)
!
!    replaces the square matrix AMAT by its transpose
!
use iomodule, only: stderr

implicit none
integer NDIM, i, j
double precision :: AMAT(NDIM,NDIM)
double precision :: temp

if (NDIM.eq.0) then
    write(stderr,*) "Vulnerability: TRANSP wants to transpose a matrix of size 0."
    write(stderr,*) "Please contact model support."
endif


do  i=1,NDIM-1
 do  j=i+1,NDIM
temp = AMAT(i,j)
AMAT(i,j) = AMAT(j,i)
AMAT(j,i) = TEMP
 enddo
enddo
!
 return
 end subroutine
!-------------------------------------------------------------------------------
double precision function DOT (A,B,N)

!    Routine to perform the dot product of two vectors

implicit none

!-pass
integer N, i
double precision :: A(N),B(N)

!-init
DOT = 0d0
do i=1,N
    DOT = DOT + A(i)*B(i)
enddo
return
end function

!-------------------------------------------------------------------------------
subroutine clear_matrix(matrix, a, b)

use iomodule, only: stderr

implicit none

integer :: a, b, i, j
double precision :: matrix(a,b)

if (a.eq.0 .or. b.eq.0) then
    write(stderr,*) "Vulnerability: clear_matrix wants to clear a matrix of size 0."
    write(stderr,*) "Please contact model support."
endif

call CLEAR_double1D(matrix, a*b, "unknown non-square matrix")

end subroutine
!-------------------------------------------------------------------------------
subroutine copy_matrix(from, to, a, b)

use iomodule, only: stderr

implicit none

integer :: a, b, i, j
double precision :: from(a,b), to(a,b)

if (a.eq.0 .or. b.eq.0) then
    write(stderr,*) "Vulnerability: copy_matrix wants to copy a matrix of size 0."
    write(stderr,*) "Please contact model support."
endif

do i=1,a
    do j=1,b
        to(i,j) = from(i,j)
    enddo
enddo


 end subroutine

!-------------------------------------------------------------------------------
subroutine IDFORM (IDMAT,NDIM)
!
!  Define values of IDMAT, which contains the addresses of the diagonals of a
!  square matrix. Used for INVMAT only
!
use iomodule, only: stderr

implicit none
!-pass
integer NDIM, i
integer :: IDMAT(NDIM)
!
if (NDIM.eq.0) then
    write(stderr,*) "Vulnerability: IDFORM want to fidget with a matrix of size 0."
    write(stderr,*) "Please contact model support."
endif

do i=1,NDIM
    IDMAT(i) = (i*(i+1))/2
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine EQUATE_int1D (TEMP,AMAT,NDIM)

! Copies the vector (matrix) AMAT into TEMP

use iomodule, only: stderr

implicit none
!-pass
integer i, NDIM
integer :: TEMP(NDIM),AMAT(NDIM)

if (NDIM.eq.0) then
    write(stderr,*) "Vulnerability: EQUATE_int wants to equate data of size 0."
    write(stderr,*) "Possible cause: a partition has no degrees of freedom."
    write(stderr,*) "If this is not the cause, please contact model support."
endif

do i=1,NDIM
    TEMP(i) = AMAT(i)
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine EQUATE_int2D (TEMP,AMAT,NDIM)

! Copies the vector (matrix) AMAT into TEMP

implicit none

!-pass
integer NDIM
integer :: TEMP(:,:),AMAT(:,:)

call EQUATE_int1D(TEMP,AMAT,NDIM)

end subroutine
!-------------------------------------------------------------------------------
subroutine EQUATE_int3D (TEMP,AMAT,NDIM)

! Copies the vector (matrix) AMAT into TEMP

implicit none

!-pass
integer NDIM
integer :: TEMP(:,:,:),AMAT(:,:,:)

call EQUATE_int1D(TEMP,AMAT,NDIM)

end subroutine

!-------------------------------------------------------------------------------
subroutine EQUATE_double1D (TEMP,AMAT,NDIM)

! Copies the vector (matrix) AMAT into TEMP

use iomodule, only: stderr

implicit none
!-pass
integer i, NDIM
double precision :: TEMP(NDIM),AMAT(NDIM)


if (NDIM.eq.0) then
    write(stderr,*) "Vulnerability: EQUATE_double wants to equate data of size 0."
    write(stderr,*) "Possible cause: a partition has no degrees of freedom."
    write(stderr,*) "If this is not the cause, please contact model support."
endif

do i=1,NDIM
    TEMP(i) = AMAT(i)
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine EQUATE_double2D (TEMP,AMAT,NDIM)

! Copies the vector (matrix) AMAT into TEMP

implicit none

!-pass
integer NDIM
double precision :: TEMP(:,:),AMAT(:,:)

call EQUATE_double1D(TEMP,AMAT,NDIM)

end subroutine
!-------------------------------------------------------------------------------
subroutine EQUATE_double3D (TEMP,AMAT,NDIM)

! Copies the vector (matrix) AMAT into TEMP

implicit none

!-pass
integer NDIM
double precision :: TEMP(:,:,:),AMAT(:,:,:)

call EQUATE_double1D(TEMP,AMAT,NDIM)

end subroutine
!-------------------------------------------------------------------------------
subroutine UPDATE(A,DA,NA)
!
! Update array A with increment DA
!
use iomodule, only: stderr

implicit none

!-pass
integer NA, i
double precision :: A(NA),DA(NA)
double precision :: temp
!
if (NA.eq.0) then
    write(stderr,*) "Vulnerability: UPDATE wants to modify vectors of length 0."
    write(stderr,*) "Possible cause: a partition has no degrees of freedom."
    write(stderr,*) "If this is not the cause, please contact model support."
endif

do i=1,NA
    temp = A(i)
    A(i) = A(i) + DA(i)
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine UPDATE_real(A,DA,NA)

! Update array A with increment DA

use iomodule, only: stderr

implicit none
!-pass
integer NA, i
real A(NA),DA(NA)

if (NA.eq.0) then
    write(stderr,*) "Vulnerability: UPDATE_real wants to modify vectors of length 0."
    write(stderr,*) "Possible cause: a partition has no degrees of freedom."
    write(stderr,*) "If this is not the cause, please contact model support."
endif


do i=1,NA
 A(i) = A(i) + DA(i)
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine rotz(R,rangle)

USE CONSTANTS, only: deg2rad

implicit none

double precision :: R(3,3)
double precision :: sinA, cosA, rangle, angle

angle = rangle * deg2rad
SINA = SIN(angle)
COSA = COS(angle)

R(1,1) = COSA
R(2,2) = COSA
R(1,2) = -SINA
R(2,1) = SINA
R(3,3) = 1d0
R(1,3) = 0d0
R(2,3) = 0d0
R(3,1) = 0d0
R(3,2) = 0d0

return
end subroutine
!-------------------------------------------------------------------------------
subroutine rotx(R,rangle)

USE CONSTANTS, only: deg2rad

implicit none

double precision :: R(3,3)
double precision :: sinA, cosA, rangle, angle

angle = rangle * deg2rad
SINA = SIN(angle)
COSA = COS(angle)


R(1,1) = 1d0
R(1,2) = 0d0
R(1,3) = 0d0
R(2,1) = 0d0
R(2,2) = COSA
R(2,3) = -SINA
R(3,1) = 0d0
R(3,2) = SINA
R(3,3) = COSA

return
end subroutine

!-------------------------------------------------------------------------------

subroutine rotvec (a,R)

implicit none

double precision :: a(3),R(3,3)
double precision :: c(3)
integer          :: i,j

do i=1,3
    c(i) = 0d0
    do j=1,3
        c(i) = c(i) + R(i,j)*a(j)
    enddo
enddo

do i=1,3
    a(i) = c(i)
enddo

return
end subroutine

!-------------------------------------------------------------------------------


function determinant3x3(A)

implicit none

double precision :: A(3,3), determinant3x3

determinant3x3 = A(1,1) * (A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
               - A(2,1) * (A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
               + A(3,1) * (A(1,2)*A(2,3)-A(1,3)*A(2,2))

end function

!-------------------------------------------------------------------------------

function determinant4x4(A)

implicit none

double precision :: A(4,4), determinant4x4

double precision :: a3x3(3,3)
double precision :: D33,sum


a3x3(1,1) = A(2,2)
a3x3(1,2) = A(2,3)
a3x3(1,3) = A(2,4)
a3x3(2,1) = A(3,2)
a3x3(2,2) = A(3,3)
a3x3(2,3) = A(3,4)
a3x3(3,1) = A(4,2)
a3x3(3,2) = A(4,3)
a3x3(3,3) = A(4,4)

D33 = determinant3x3(a3x3)

sum = A(1,1)*D33
a3x3(1,1) = A(1,2)
a3x3(1,2) = A(1,3)
a3x3(1,3) = A(1,4) 
!    a21 = A(3,2)
!    a22 = A(3,3)
!    a23 = A(3,4)
!    a31 = A(4,2)
!    a32 = A(4,3)
!    a33 = A(4,4)

D33 = determinant3x3(a3x3)

sum = sum - A(2,1)*D33
!    a11 = A(1,2)
!    a12 = A(1,3)
!    a13 = A(1,4)
a3x3(2,1) = A(2,2) 
a3x3(2,2) = A(2,3)
a3x3(2,3) = A(2,4)
!    a31 = A(4,2)
!    a32 = A(4,3)
!    a33 = A(4,4)

D33 = determinant3x3(a3x3)

sum = sum + A(3,1)*D33
!    a11 = A(1,2)
!    a12 = A(1,3)
!    a13 = A(1,4)
!    a21 = A(2,2)
!    a22 = A(2,3)
!    a23 = A(2,4)
a3x3(3,1) = A(3,2)
a3x3(3,2) = A(3,3)
a3x3(3,3) = A(3,4)

D33 = determinant3x3(a3x3)

determinant4x4 = sum - A(4,1)*D33

end function

!-------------------------------------------------------------------------------

double precision function VolumeOfTetrahedron(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)

implicit none

! let T be tetrahedron with 4 vertices with coordinates

! x1, y1, z1
! x2, y2, z2
! x3, y3, z3
! x4, y4, z4

! which is moved with the fourth points to be on the origin,
! so that we get coordinates:

! x1-x4, y1-y4, z1-z4
! x2-x4, y2-y4, z2-z4
! x3-x4, y3-y4, z3-z4
!     0,     0,     0

! six times the volume of the element is the determinant of matrix

! x1-x4, y1-y4, z1-z4
! x2-x4, y2-y4, z2-z4
! x3-x4, y3-y4, z3-z4

double precision :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
double precision :: spanningVectors(3,3)

spanningVectors(1,1) = x1-x4
spanningVectors(2,1) = x2-x4
spanningVectors(3,1) = x3-x4

spanningVectors(1,2) = y1-y4
spanningVectors(2,2) = y2-y4
spanningVectors(3,2) = y3-y4

spanningVectors(1,3) = z1-z4
spanningVectors(2,3) = z2-z4
spanningVectors(3,3) = z3-z4


VolumeOfTetrahedron = -determinant3x3(spanningVectors) / 6d0

end function

!-------------------------------------------------------------------------------


integer function IDIG (number,length,i)
!
! returns digit number I of integer NUMBER with LENGTH positions
!

implicit none

integer number, length, i
integer j, jdig, n


n = IABS(number)
if (length.gt.i) then
    do j=1,length-i
        jdig = MOD(n,10)
        n = (n-jdig)/10
    enddo
    IDIG = MOD(n,10)
else if (length.eq.i) then
    IDIG = MOD(n,10)
else
    IDIG = 0
endif

return
end function

!-------------------------------------------------------------------------------


! below here are helper functions of the sampling in 3D
! We have to determine which element is closer to the sampling point,
! from a given element.


function twoPointsOnSameSideOfPlane(vecSpan1x, vecSpan1y, vecSpan1z, &
                                    vecSpan2x, vecSpan2y, vecSpan2z, &
                                      vecPt1x,   VecPt1y,   VecPt1z, &
                                      vecPt2x,   VecPt2y,   VecPt2z)

! given is a plane in 3D spanned up by two vectors.
! two points are given my means of vectors from the same
! origin as the spanning vectros. 
! Returns true when they are on the same side of the plane.
! Returns false otherwise.

! This is checked by taking the determinants of the spanning vectors with
! each of the pointing vectors and check whether the signs are equal (same side)
! or not (opposite sides)

use constants, only: eps

implicit none

logical          :: twoPointsOnSameSideOfPlane
double precision :: vecSpan1x, vecSpan1y, vecSpan1z
double precision :: vecSpan2x, vecSpan2y, vecSpan2z
double precision ::   vecPt1x,   VecPt1y,   VecPt1z
double precision ::   vecPt2x,   VecPt2y,   VecPt2z

double precision :: DeterminantMatrix(3,3)

double precision :: detPoint1, detPoint2

DeterminantMatrix(1,1) = vecSpan1x
DeterminantMatrix(2,1) = vecSpan1y
DeterminantMatrix(3,1) = vecSpan1z

DeterminantMatrix(1,2) = vecSpan2x
DeterminantMatrix(2,2) = vecSpan2y
DeterminantMatrix(3,2) = vecSpan2z

DeterminantMatrix(1,3) = vecPt1x
DeterminantMatrix(2,3) = vecPt1y
DeterminantMatrix(3,3) = vecPt1z

detPoint1 = determinant3x3(DeterminantMatrix)


DeterminantMatrix(1,3) = vecPt2x
DeterminantMatrix(2,3) = vecPt2y
DeterminantMatrix(3,3) = vecPt2z

detPoint2 = determinant3x3(DeterminantMatrix)


! -eps, we give it the benefit of the doubt.
! in the worst case scenario it will take one extra 
! step to approach the sample point.
if ((detPoint1 * detPoint2) .gt. -eps) then
    twoPointsOnSameSideOfPlane = .true.
else
    twoPointsOnSameSideOfPlane = .false.
endif


return


end function

!-------------------------------------------------------------------------------

double precision function areaOfTriangleIn3D(xa, ya, za, xb, yb, zb, xc, yc, zc)

! area of a 3D-triangle with points A, B an C
! define sides AB and AC as vectors from point A to B and C
! find angle between AB and BC: AB . AC = |AB||AC| cos theta
! now area = 1/2 |AB||AC| sin theta



implicit none

double precision :: xa, ya, za
double precision :: xb, yb, zb
double precision :: xc, yc, zc

double precision :: ab(3)
double precision :: ac(3)

double precision :: abDOTac
double precision :: abLength
double precision :: acLength
double precision :: theta

!double precision, external :: vectorLength
ab(1) = xb - xa
ab(2) = yb - ya
ab(3) = zb - za

ac(1) = xc - xa
ac(2) = yc - ya
ac(3) = zc - za

abDOTac  = dot(ab, ac, 3)

abLength = vectorLength(ab,3)
acLength = vectorLength(ac,3)

theta = acos(abDOTac / (abLength * acLength))

areaOfTriangleIn3D = 0.5d0 * sin(theta) * abLength * acLength

end function

!-------------------------------------------------------------------------------

double precision function vectorLength(vector, length)

implicit none

integer          :: length
double precision :: vector(length)

integer          :: iEntry

vectorLength = 0d0
do iEntry=1,length
    vectorLength = vectorLength + vector(iEntry)**2
enddo

vectorLength = sqrt(vectorLength)

end function


!-------------------------------------------------------------------------------

double precision function distanceBetweenTwoPoints(x1, y1, z1, x2, y2, z2)

implicit none

double precision :: x1, y1, z1, x2, y2, z2

distanceBetweenTwoPoints = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )

end function
!-------------------------------------------------------------------------------

double precision function distanceBetweenTwoPointsOnASphere(x1, y1, x2, y2, r)

use constants, only: one


implicit none

double precision :: x1, y1, x2, y2, r
double precision :: distance

double precision :: A(3), B(3), C(3)

! first set the vector of the first nodal location on a shell
call cartesian(one, x1, y1, A)

! then set the vector of the second nodal location on a shell
call cartesian(one, x2, y2, B)

! next calculate the area
call CROSSP(A,B,C)

distance = DOT(C,C,3)              ! compute horizontal vector length
distance = SQRT(distance)
distance = ASIN(distance)             ! |a x b| = |a|.|b| sin(alpha)

distanceBetweenTwoPointsOnASphere = distance*r

end function

!******************************
subroutine cartesian (R,lon,lat,vec)

USE CONSTANTS, only: deg2rad
! subroutine to set normal vectors on a unit sphere
! useful for calculating distances on a sphere
implicit none

double precision :: R,lon,lat,vec(3)

lon = lon * deg2rad
lat = lat * deg2rad

vec(1) = R*COS(lat)*COS(lon)
vec(2) = R*COS(lat)*SIN(lon)
vec(3) = R*SIN(lat)

return

end subroutine

!-------------------------------------------------------------------------------




function distanceBetweenPointAndPlane(x1, y1, z1, x2, y2, z2, x3, y3, z3, xp, yp, zp)

! given is a plane in 3D spanned up by three points

use constants, only: eps

implicit none

!-pass
double precision :: distanceBetweenPointAndPlane
double precision :: x1, y1, z1, x2, y2, z2, x3, y3, z3, xp, yp, zp

!-local
double precision :: DeterminantMatrix(3,3)

double precision :: vecSpan1x, vecSpan1y, vecSpan1z
double precision :: vecSpan2x, vecSpan2y, vecSpan2z

vecSpan1x = x1 - x2
vecSpan1y = y1 - y2
vecSpan1z = z1 - z2

vecSpan2x = x1 - x3
vecSpan2y = y1 - y3
vecSpan2z = z1 - z3

DeterminantMatrix(1,1) = vecSpan1x
DeterminantMatrix(2,1) = vecSpan1y
DeterminantMatrix(3,1) = vecSpan1z

DeterminantMatrix(1,2) = vecSpan2x
DeterminantMatrix(2,2) = vecSpan2y
DeterminantMatrix(3,2) = vecSpan2z

DeterminantMatrix(1,3) = xp
DeterminantMatrix(2,3) = yp
DeterminantMatrix(3,3) = zp

distanceBetweenPointAndPlane = determinant3x3(DeterminantMatrix)

return

end function

!-------------------------------------------------------------------------------

function radiusOfSphereInsideTetrahedron(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4) 

implicit none

!-pass
double precision :: radiusOfSphereInsideTetrahedron
double precision :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4

!-local
double precision :: side1, side2, side3, side4
double precision :: totalVolume

! We divide the tetrahedron into 4 smaller ones, one with each of the sides as its base
! and a point in the middle (the center of the inscribed sphere) as the tip for each of the tetrahedrons.

totalVolume = VolumeOfTetrahedron(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

side1 = areaOfTriangleIn3D(x1, y1, z1, x2, y2, z2, x3, y3, z3)
side2 = areaOfTriangleIn3D(x1, y1, z1, x2, y2, z2, x4, y4, z4)
side3 = areaOfTriangleIn3D(x1, y1, z1, x3, y3, z3, x4, y4, z4)
side4 = areaOfTriangleIn3D(x2, y2, z2, x3, y3, z3, x4, y4, z4)

! sum of the volumes of the individual pyramids equals that of the whole tetrahedron

! r * side1   r * side2   r * side3   r * side4
! --------- + --------- + --------- + ---------  = totalVolume 
!     3           3           3           3

! r * side1   r * side2   r * side3   r * side4
! --------- + --------- + --------- + ---------  = totalVolume 
!     3           3           3           3

! r * side1 + r * side2 + r * side3 + r * side4  = 3 totalVolume

!              3 totalVolume
! r = -------------------------------
!      side1 + side2 + side3 + side4

radiusOfSphereInsideTetrahedron = 3d0 * totalVolume / (side1 + side2 + side3 + side4)

end function

!-------------------------------------------------------------------------------

function radiusOfSphereAroundTetrahedron(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

! we follow Michael Woltermann
! http://www2.washjeff.edu/users/mwoltermann/Dorrie/70.pdf

implicit none

!-pass
double precision :: radiusOfSphereAroundTetrahedron
double precision :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4

!-local
double precision :: e,f,g
double precision :: a,p,b,q,c,r
double precision :: j
double precision :: totalVolume

double precision :: peri !meter of triangle

a = distanceBetweenTwoPoints(x3, y3, z3, x4, y4, z4)
b = distanceBetweenTwoPoints(x4, y4, z4, x2, y2, z2)
c = distanceBetweenTwoPoints(x3, y3, z3, x2, y2, z2)
p = distanceBetweenTwoPoints(x1, y1, z1, x2, y2, z2)
q = distanceBetweenTwoPoints(x1, y1, z1, x3, y3, z3)
r = distanceBetweenTwoPoints(x1, y1, z1, x4, y4, z4)

e = a * p
f = b * q
g = c * r

! use Heron's formula to compute area of triangle
peri = 0.5d0 * (e + f + g)
j = sqrt(peri * (peri - e) * (peri - f) * (peri - g))

totalVolume = VolumeOfTetrahedron(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

radiusOfSphereAroundTetrahedron = j / (6d0 * totalVolume)


end function

!-------------------------------------------------------------------------------

subroutine arraycontains(array, length, value, pos)

implicit none

integer array(*)
integer length, value, pos,i


do i=1,length
    if (array(i).eq.value) then
        pos = i
        return
    endif
enddo
pos = length + 1

end subroutine arraycontains

!-------------------------------------------------------------------------------

function threeInFour(threeArray, fourArray)

implicit none
logical :: threeInFour

integer :: threeArray(3)
integer :: fourArray(4)

integer :: i, pos

threeInFour = .true.
do i=1,3
    call arraycontains(fourArray, 4, threeArray(i), pos)
    if (pos.eq.5) then
        threeInFour =  .false.
    endif
enddo

end function

!-------------------------------------------------------------------------------

subroutine sort4vec(vec)
! sort result from small to big, using the Shell method

implicit none
integer vec(4)

if (vec(1).gt.vec(3)) then
    call swapInts(vec(1), vec(3))
endif

if (vec(2).gt.vec(4)) then
    call swapInts(vec(2), vec(4))
endif

if (vec(1).gt.vec(2)) then
    call swapInts(vec(1), vec(2))
endif

if (vec(3).gt.vec(4)) then
    call swapInts(vec(3), vec(4))
endif

if (vec(2).gt.vec(3)) then
    call swapInts(vec(2), vec(3))
endif


end subroutine

!-------------------------------------------------------------------------------

subroutine swapInts(a,b)
implicit none
integer a,b,tmp
tmp = a
a = b
b = tmp
end subroutine

!-------------------------------------------------------------------------------

END MODULE ALGEBRA
!-------------------------------------------------------------------------------



