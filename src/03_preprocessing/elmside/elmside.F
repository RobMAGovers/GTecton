! the goal of this program is to replace the elmside in C.
! The reason for doing this in Fortran was that the routine
! for reading the tecin.dat.partf.nps already existed in Fortran.

! moreover, fortran has logicals, it has allocate which is easier than malloc.
! the only drawback is that is a bit more iffy with reading command line args,
! but alas...

! for a quick overview of the command line arguments:
!            call stringfromarg(arg, 'n', tecindatpartfnps)
!            call stringfromarg(arg, 'e', tecindatpartfelm)
!            call stringfromarg(arg, 'l', vertexlist)

!            call intfromarg(arg, 'd', ndimensions)

!            call logicalfromarg(arg, 'w', printweights)
!            call logicalfromarg(arg, 'c', printcoords)
!            call logicalfromarg(arg, 'o', opn)
!            call logicalfromarg(arg, 's', sideonly)
!            call logicalfromarg(arg, 'i', invertsign)
!            call logicalfromarg(arg, 'p', endpoints)



! Just for rememberance, the side convention:

    ! 2D:
    !   face 1: vertices 1 2 
    !        2:          2 3 
    !        3:          3 1 
    ! 3D:
    !   face 1: vertices 1 2 4
    !        2:          1 3 2
    !        3:          1 4 3
    !        4:          2 3 4 


!***********************************************************
!           Main program  
!***********************************************************

program elmside

    use globals  ! literally everything, so no only statement needed

    implicit none

! internal stuff
    integer, dimension(4)              :: selected    ! check, vertex of element is in selection list
                                                      ! in 2D, only the first 3 are used. in 3D, all are used.
    integer, dimension(4)              :: isselected  ! same as selected, but holds 0/1 in stead of position

    integer :: i,j, ierr
    integer :: nmatch, numinarray
    logical                            :: publish

! sign/weight determination
    integer                            :: b1, b2, b3, thissign
    integer                            :: elm_sign    

! coordinates for output
    double precision                   :: xmidpoint, ymidpoint

! functions
    double precision                   :: triangle_area

! debuging:
    integer                            :: foo

! for Nicolai
    integer, external                  :: signOfElemWithOneNodeAdjacent2D

    eps = 0.000001

! ******************* testspul **********
    goto 10
! comment goto statement uit om te testen ***
!*********************************

! test selector:
!    integer, dimension(10) :: list
!    integer, dimension(5)  :: selectednrs

!    double precision :: triange_area
!    double precision :: area
!    double precision :: x1,y1,z1, x2,y2,z2, x3,y3,z3


!    list = (/ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 /)

!    do i=1,20
!        call selectsublistfromlist(10, list, 5, selectednrs) 
!        write(*,*) 'selection completed: ', selectednrs
!    enddo

!    x1 = 1    y1 = 0
!    z1 = 0

!    x2 = 0
!    y2 = 1
!    z2 = 0

!    x3 = 0
!    y3 = 0
!    z3 = 1

!   area = triange_area(x1,y1,z1, x2,y2,z2, x3,y3,z3 )

!   write(*,*) 'area: ', area


!    stop
!******************* einde testspul *****


! commmand line argument variables are all in the globals module, hence no arguments
! (line numbered to possible skip or insert testing stuff)

10  call readcommandlineargs()
call validatecommandlineargs()
call readInput()


! depending on the choice of input paameter a,
! different algoritms can be used that work in various situations.

if (algorithm .eq. 1) then
     call algo1()
else if (algorithm .eq. 2) then 
    call algo2()
else if (algorithm .eq. 3) then 
    call algo3()
else
    ! revert to default algorithm 1
    call algo1()
endif

call freeglobals

end program
