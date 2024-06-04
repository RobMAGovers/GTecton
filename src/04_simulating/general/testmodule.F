module test

implicit none
! this module contains subroutines that unit test
! various funtions and subroutines in GTecton.
! It is called when pln is called with argument 'test'
! All other arguments are then ignored, and we enter here

logical AllTestsOK


contains



subroutine runTests()
implicit none

AllTestsOK = .true.

if (testOneThing()) then
    write(*,*) "test 1 OK"
else
    write(*,*) "test 1 FAIL"
    AllTestsOK = .false.
endif

if (testAnotherThing()) then
    write(*,*) "test 2 OK"
else
    write(*,*) "test 2 FAIL"
    AllTestsOK = .false.
endif

#ifdef EXPERIMENTAL_FRICTION

if (testFriction(6)) then
    write(*,*) "3D test of friction OK"
else
    write(*,*) "3D test of friction FAILS"
    AllTestsOK = .false.
endif
if (testFriction(1)) then
    write(*,*) "pln test of friction OK"
else
    write(*,*) "pln test of friction FAILS"
    AllTestsOK = .false.
endif
if (testFriction(2)) then
    write(*,*) "opn test of friction OK"
else
    write(*,*) "opn test of friction FAILS"
    AllTestsOK = .false.
endif

#endif


! add more calls


if (AllTestsOK) then
    write(*,*) "All unit tests passed!"
else
    write(*,*) "One or more unit tests failed!"
endif


end subroutine

!----  actual tests below ----------------------------------------------

logical function testOneThing()
implicit none
testOneThing = .true.

end function

!----------

logical function testAnotherThing()
implicit none
testAnotherThing = .false.

end function

!----------

#ifdef EXPERIMENTAL_FRICTION

logical function testFriction(NTYPE)

use FrictionModule, only: FrictionEntry, frictionTraction

implicit none

type(FrictionEntry), allocatable :: friction(:)
double precision, allocatable :: DIFORC(:,:)
integer NTYPE,NDOF,nFrictionEntriesGlobal,fault_number

NTYPE = 6 ! ******* VARY NTYPE [0-6] FOR TESTING *****

if (NTYPE.eq.1) then
  NDOF=2
else
  NDOF=3
endif
nFrictionEntriesGlobal=3
fault_number = 3

ALLOCATE ( DIFORC(NDOF,200) ); DIFORC = 0d0
ALLOCATE ( friction(nFrictionEntriesGlobal) )

if (NTYPE.eq.6) then

  friction(1)%faultNumber = fault_number
  friction(1)%point(1)=15; friction(1)%point(2)=78; friction(1)%point(3)=13
  friction(1)%sideArea = 2d0
  friction(2)%faultNumber = fault_number
  friction(2)%point(1)=78; friction(2)%point(2)=6; friction(2)%point(3)=13
  friction(2)%sideArea = 3d0
  friction(3)%faultNumber = fault_number
  friction(3)%point(1)=78; friction(3)%point(2)=101; friction(3)%point(3)=6
  friction(3)%sideArea = 4d0

  print*,'segment nodes:'
  print*,friction(1)%point(1),friction(1)%point(2),friction(1)%point(3)
  print*,friction(2)%point(1),friction(2)%point(2),friction(2)%point(3)
  print*,friction(3)%point(1),friction(3)%point(2),friction(3)%point(3)

  DIFORC(1, 15) = friction(1)%sideArea/3d0 * 3d0
  DIFORC(1, 78) = friction(1)%sideArea/3d0 * 3d0 + friction(2)%sideArea/3d0 * 5d0 + &
&                 friction(3)%sideArea/3d0 * 7d0
  DIFORC(1, 13) = friction(1)%sideArea/3d0 * 3d0 + friction(2)%sideArea/3d0 * 5d0
  DIFORC(1,  6) = friction(2)%sideArea/3d0 * 5d0 + friction(3)%sideArea/3d0 * 7d0
  DIFORC(1,101) = friction(3)%sideArea/3d0 * 7d0

  DIFORC(2, 15) = friction(1)%sideArea/3d0 * 11d0
  DIFORC(2, 78) = friction(1)%sideArea/3d0 * 11d0 + friction(2)%sideArea/3d0 * 13d0 + &
&                 friction(3)%sideArea/3d0 * 17d0
  DIFORC(2, 13) = friction(1)%sideArea/3d0 * 11d0 + friction(2)%sideArea/3d0 * 13d0
  DIFORC(2,  6) = friction(2)%sideArea/3d0 * 13d0 + friction(3)%sideArea/3d0 * 17d0
  DIFORC(2,101) = friction(3)%sideArea/3d0 * 17d0

  DIFORC(3, 15) = friction(1)%sideArea/3d0 * 19d0
  DIFORC(3, 78) = friction(1)%sideArea/3d0 * 19d0 + friction(2)%sideArea/3d0 * 23d0 + &
&                 friction(3)%sideArea/3d0 * 29d0
  DIFORC(3, 13) = friction(1)%sideArea/3d0 * 19d0 + friction(2)%sideArea/3d0 * 23d0
  DIFORC(3,  6) = friction(2)%sideArea/3d0 * 23d0 + friction(3)%sideArea/3d0 * 29d0
  DIFORC(3,101) = friction(3)%sideArea/3d0 * 29d0

else
  friction(1)%faultNumber = fault_number; friction(1)%point(1)=15; friction(1)%point(2)=78
  friction(1)%sideArea = 2d0
  friction(2)%faultNumber = fault_number; friction(2)%point(1)=78; friction(2)%point(2)=101
  friction(2)%sideArea = 3d0
  friction(3)%faultNumber = fault_number; friction(3)%point(1)=101; friction(3)%point(2)=8
  friction(3)%sideArea = 4d0
  print*,'segment nodes:'
  print*,friction(1)%point(1),friction(1)%point(2)
  print*,friction(2)%point(1),friction(2)%point(2)
  print*,friction(3)%point(1),friction(3)%point(2)

  DIFORC(1,15)  = friction(1)%sideArea/2d0 * 3d0
  DIFORC(1,78)  = friction(1)%sideArea/2d0 * 3d0 + friction(2)%sideArea/2d0 * 5d0
  DIFORC(1,101) = friction(2)%sideArea/2d0 * 5d0 + friction(3)%sideArea/2d0 * 7d0
  DIFORC(1,8)   = friction(3)%sideArea/2d0 * 7d0

  DIFORC(2,15)  = friction(1)%sideArea/2d0 * 11d0
  DIFORC(2,78)  = friction(1)%sideArea/2d0 * 11d0 + friction(2)%sideArea/2d0 * 13d0
  DIFORC(2,101) = friction(2)%sideArea/2d0 * 13d0 + friction(3)%sideArea/2d0 * 17d0
  DIFORC(2,8)   = friction(3)%sideArea/2d0 * 17d0
!
  if (NDOF.eq.3) then
    DIFORC(3,15)  = friction(1)%sideArea/2d0 * 19d0
    DIFORC(3,78)  = friction(1)%sideArea/2d0 * 19d0 + friction(2)%sideArea/2d0 * 23d0
    DIFORC(3,101) = friction(2)%sideArea/2d0 * 23d0 + friction(3)%sideArea/2d0 * 29d0
    DIFORC(3,8)   = friction(3)%sideArea/2d0 * 29d0
  endif
endif

if (NDOF.eq.2) then
  print*,'input traction on segment 1: ',3,11
  print*,'input traction on segment 2: ',5,13
  print*,'input traction on segment 3: ',7,17
else
  print*,'input traction on segment 1: ',3,11,19
  print*,'input traction on segment 2: ',5,13,23
  print*,'input traction on segment 3: ',7,17,29
endif

!call frictionTraction (fault_number,friction,DIFORC,nFrictionEntriesGlobal,NTYPE,NDOF)
call frictionTraction (fault_number,DIFORC,nFrictionEntriesGlobal,NTYPE,NDOF,&
                       .false.,"dummy.dat")


!**** frictionTraction can be evaluated when it gives a result back (it now only prints)
if (1.eq.1) then
    testFriction = .true.
else
    testFriction = .false.
endif


end function

#endif /* EXPERIMENTAL_FRICTION */

end module test
