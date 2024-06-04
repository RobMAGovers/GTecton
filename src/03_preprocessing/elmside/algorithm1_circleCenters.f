subroutine algo1()

!****** now that we have all the data in memory, we can compute the element sides.

! First of all we will create a list with all the information,
! Next we will print only what is needed, depending on the command
! line arguments.

! compute matching points and corresponding sides.

! Difficult case:

!         *--------*---------*
!         |      . .
!         |  1 .   .
!         | .   2  .
! *-------*. . . . .

! Let us have points collection marked with *
! marking a boundary marked with - and |
! Element 1 has multiple sides on this boundary, which must both be mentioned
! Element 2 has multiple points on the boundary, but no side. 
 
! a similar situation could be imagined in 3D   

! we check this in 2D: 
! If an element has two (and only two) vertices on the selection list.
! If those points have a common neighbour, and this neighbour is also in the selection list,
! but not on the element, it should not be included.

! Similarly in 3D, if an element has three (and only three) vertices in the 
! selection list, and they all have a common neighbour that is also in the 
! selection list, but not on the element, is not included.

! Moreover, of an element is adjacent to the boundary, with all points,
! we have to establish which sides are attached to the boundary. 

! Elmside is able to distinguish between side of a fraction, so that elements
! on one side can be distinguished from elements on the other side.
! The previous version of elmside did this by sorting the elements by location (polysort)
! and walking through this list. Because we do this also in 3D, it is not possible
! to construct a trivial sorted list that could help in this.
! In stead, both in 2D and in 3D, 2 or 3, respectively, are chosen to construct a basis.
! in 3D they form a triangle with the biggest surface area, in 2D they form the longest 
! line element. The orthogonal vector decides the direction. And side with a degree of
! less than 90 degrees gets a sign equal to that of the orthogonal, otherwise thse opposite.

! The sign of the orthogonal is determined by:
! if it is positive in x, it is +1, if negative -1
! if x = 0, then
! if it is positive in y, it is +1, if negative -1
! if both x and y are 0
! if it is positive in z, it is +1, if negative -1

! This should always yield a nonzero sign.


! find the points that control the basis, and determine the direction of the sign
! This can of course only happen when there is a coordinate file
! All variables used are from the globals module.


use globals
use stringHandling, only: appendtostring, &
                          cleanoutputstrings

implicit none

! internal stuff
integer, dimension(4)              :: selected    ! check, vertex of element is in selection list
                                                      ! in 2D, only the first 3 are used. in 3D, all are used.
integer, dimension(4)              :: isselected  ! same as selected, but holds 0/1 in stead of position

integer :: i,j, ierr, iLine, iOutput
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

! coords in 3D
double precision                   :: touchCoords(3)


!    write(*,*) 'endpoints before basis: ', endpoints

    if (iecho.eq.1) then
        write(0,*) "Building basis"
    endif

    call basis()

    if (iecho.eq.1) then
        write(0,*) "Finished basis"
        write(0,*) "n elems: ", nelems
    endif


! TODO: perhaps sort the vertexlist, to make searching through it easier.
! However, numinarray is seldom gigantic, so it takes not too much time.

    do i=1,nelems
        
        call cleanOutputStrings()

        ! check which nodes of the element occur in the cloud of selected points

        if (iecho.eq.1) then
            write(*,*) "-------------------------------------------------"
            write(*,*) "checking element with vertices", i, v(1,i),v(2,i),v(3,i),v(4,i)
        endif

        selected(1) = numinarray(v(1,i), vertices, size(vertices))
        selected(2) = numinarray(v(2,i), vertices, size(vertices))
        selected(3) = numinarray(v(3,i), vertices, size(vertices))

        if (ndimensions.eq.3) then
            ! tetrahedron has extra vertex
            selected(4) = numinarray(v(4,i), vertices, size(vertices))
        endif

        if (iecho.eq.1) then
            if (ndimensions.eq.3) then
                write(*,*) "checked ", v(1,i), v(2,i), v(3,i), v(4,i)
                write(*,*) "Has selected", selected(1:4)
             else
                write(*,*) "Has selected", selected(1:3)
            endif
        endif


        


        ! and count the number of vertices with which is borders on the cloud

        do j=1,3
            isselected(j) = 0
            if (selected(j).gt.0) then
                isselected(j) = 1
            endif
        enddo

        if (ndimensions.eq.3) then
            isselected(4) = 0
            if (selected(4).gt.0) then
                isselected(4) = 1
            endif
        endif         

        nmatch = isselected(1) + isselected(2) + isselected(3)
        if (ndimensions.eq.3) then
            nmatch = nmatch + isselected(4)
        endif


        if (iecho.eq.1 .and. nmatch.gt.0) then
            if (ndimensions.eq.3) then
                 write(*,*) "element ",i , &
                        " has ", nmatch, &
                        " point in selection ",selected(1:4)
            else
                 write(*,*) "element ",i , &
                        " has ", nmatch, &
                        " point in selection ",selected(1:3)
            endif
        endif

        ! seven possible cases to be treated separately
        ! they are:

        ! 2 dimension; 1 attached node -> only output if sideonly == false
        ! 2 dimension; 2 attached node -> check if it not a corner exclusion
        ! 2 dimension; 3 attached node -> check which two sides touch the boundary

        ! 3 dimension; 1 attached node -> only output if sideonly == false
        ! 3 dimension; 2 attached node -> only output if sideonly == false
        ! 3 dimension; 3 attached node -> check if is not an exclusion
        ! 3 dimension; 4 attached node -> check which two or three sides are attached
        !                                 This is a complicated case, but will not occur very often.




         ! theoretically, 0 attached nodes are also a case, 
         ! but we ignore those. ;)

         !******************  2D  ***********************
         if (ndimensions.eq.2) then

             if      (nmatch.eq.1) then

                 ! only include with nodal attachment                 
                 ! If it has only a single node, there is no risk
                 ! that the element is included twice

                 ! We compare the direction of the vector from the
                 ! attached point to the middle of the opposite
                 ! side to determine the sign, if this is required.

                 ! the elmsign routine computes the vector the other
                 ! way around, for when an element touches the
                 ! boundary with a side, in sted of with a single point,
                 ! hence a minus sign. 
                
                 if (.not.sideonly) then

                    ! all of the nodes of this element that are on the collection
                    ! get a separate entry


!                     write(*,*) 'coords of point', i, coords(1,i), coords(2,i)


                     ! when only a single point matches, there is nog adjacent side

                     ! in order to circumvent this problem, we extrapolate the basis vector through
                     ! the point that is on the curve.
                     ! The center of the element will be one side or the other,
                     ! determining the sign.

                     if (isselected(1).eq.1) then

                         if (printweights) then ! point 1 matches, opposite spans 2 and 3 -> side nr 2
                             thissign = elm_sign(v(1,i), v(2,i), v(3,i), 0, 2, nmatch)
!                             thissign = - thissign
!                             call output2D(i, v(1,i), 0, coords(1,i), coords(2,i), 0, thissign)
                             call output2D(i, v(1,i), 0, coords(1,v(1,i)), coords(2,v(1,i)), 0, thissign)
                         else
!                             call output2D(i, v(1,i), 0, coords(1,i), coords(2,i), 0, 0)
                             call output2D(i, v(1,i), 0, coords(1,v(1,i)), coords(2,v(1,i)), 0, 0)
                         endif

                     else if (isselected(2).eq.1) then
   
                         if (printweights) then ! point 2 matches, opposite spans 1 and 3 -> side nr 3
                             thissign = elm_sign(v(1,i), v(2,i), v(3,i), 0, 3, nmatch)
!                             thissign = - thissign
!                             call output2D(i, v(2,i), 0, coords(1,i), coords(2,i), 0, thissign)
                             call output2D(i, v(2,i), 0, coords(1,v(2,i)), coords(2,v(2,i)), 0, thissign)
                         else
!                             call output2D(i, v(2,i), 0, coords(1,i), coords(2,i), 0, 0)
                             call output2D(i, v(2,i), 0, coords(1,v(2,i)), coords(2,v(2,i)), 0, 0)
                         endif

                     else if (isselected(3).eq.1) then

                         if (printweights) then  ! point 3 matches, opposite spans 1 and 2 -> side nr 1
                             thissign = elm_sign(v(1,i), v(2,i), v(3,i), 0, 1, nmatch)
!                             thissign = - thissign
!                             call output2D(i, v(3,i), 0, coords(1,i), coords(2,i), 0, thissign)
                             call output2D(i, v(3,i), 0, coords(1,v(3,i)), coords(2,v(3,i)), 0, thissign)
                         else
!                             call output2D(i, v(3,i), 0, coords(1,i), coords(2,i), 0, 0)
                             call output2D(i, v(3,i), 0, coords(1,v(3,i)), coords(2,v(3,i)), 0, 0)
                        endif


                     endif ! isselected 1,2,3

                 endif ! .not. sideonly

             else if (nmatch.eq.2) then
                 ! with nodal attchament, the element is included always.
                 
                 ! check whether the two points have a common neighbour
                 ! that is on the vertex list, but not the element.

!                 write(*,*) '---------- elem: ', i, ' will get two entries  '

                 if (isselected(1).eq.1 .and. isselected(2).eq.1) then
                   ! side 1
                     if (iecho.eq.1) then
                         write(*,*) 'side 1 matches; with vertices', v(1,i), v(2,i)
                     endif
                     call mutualneighbours(v(1,i),v(2,i),0)
                 else if (isselected(1).eq.1 .and. isselected(3).eq.1) then
                     if (iecho.eq.1) then
                         write(*,*) 'side 2 matches; with vertices', v(1,i), v(3,i)
                     endif
                     call mutualneighbours(v(1,i),v(3,i),0)
                 else if (isselected(2).eq.1 .and. isselected(3).eq.1) then
                     if (iecho.eq.1) then
                         write(*,*) 'side 3 matches; with vertices', v(2,i), v(3,i)
                     endif
                     call mutualneighbours(v(2,i),v(3,i),0)
                 endif

                ! mutual neighbour IDs are now in 'neighbours'
                 publish = .true.
!                 write(*,*) 'neighbours: ', neighbours
                 do j=1,2  ! check for two possible neighbours
                    if (neighbours(j) .gt. 0) then
                        ! if point is in vertex list ...
                        if ((numinarray(neighbours(j),vertices,nvertices).gt.0).and. &
                        ! ... but not in element
                        (neighbours(j).ne.v(1,i)) .and. &
                        (neighbours(j).ne.v(2,i)) .and. &
                        (neighbours(j).ne.v(3,i)) .and. &
                        (neighbours(j).ne.v(4,i))  ) then
                            ! we have an imposter, do not publish the side!
                            publish = .false.

 !                           write(*,*) 'sadly, it is excluded'
                        else
                            ! nothing happens, this is OK.
                        endif

                    endif
                 enddo

                 ! if the element has passed this test, it can be published
                 if (publish) then
                     ! todo: write to output, we have a side

                     if (isselected(1).eq.1 .and. isselected(2).eq.1) then
!                         write(*,*) 'side 1: ', i, v(2,i), v(3,i)
                         xmidpoint = 0.5 * ( coords(1,v(1,i)) + coords(1,v(2,i)) )
                         ymidpoint = 0.5 * ( coords(2,v(1,i)) + coords(2,v(2,i)) )
                         thissign = elm_sign(v(1,i), v(2,i), v(3,i), 0, 1, nmatch)

                         call output2D(i, v(1,i), v(2,i) , xmidpoint, ymidpoint, 1, thissign)

                     else if (isselected(1).eq.1 .and. isselected(3).eq.1) then
!                         write(*,*) 'side 3: ', i, v(1,i), v(3,i)
                         xmidpoint = 0.5 * ( coords(1,v(1,i)) + coords(1,v(3,i)) )
                         ymidpoint = 0.5 * ( coords(2,v(1,i)) + coords(2,v(3,i)) )
                         thissign = elm_sign(v(1,i), v(2,i), v(3,i), 0, 3, nmatch)

                         call output2D(i, v(1,i), v(3,i) , xmidpoint, ymidpoint, 3, thissign)


                     else if (isselected(2).eq.1 .and. isselected(3).eq.1) then
!                         write(*,*) 'side 2: ', i, v(2,i), v(3,i))
                         xmidpoint = 0.5 * ( coords(1,v(2,i)) + coords(1,v(3,i)) )
                         ymidpoint = 0.5 * ( coords(2,v(2,i)) + coords(2,v(3,i)) )
                         thissign = elm_sign(v(1,i), v(2,i), v(3,i), 0, 2, nmatch)

                         call output2D(i, v(2,i), v(3,i) , xmidpoint, ymidpoint, 2, thissign)


                     endif
                 endif

             else if (nmatch.eq.3) then
                ! TODO: determine which side is not adjacent to the boundary, and print the other two to std out.
                ! this is done by checking for each... 



             endif
         !******************  3D  ***********************
         else if (ndimensions.eq.3) then

!            write(*,*) "element ", i, "has nmatch", nmatch

    ! 3D:
    !   face 1: vertices 1 2 4
    !        2:          1 3 2
    !        3:          1 4 3
    !        4:          2 3 4 

            iLine = 1

            if (.not.sideonly) then

                ! all possible output lines will be filled
                do iOutput=1,4
                    ! add element number i
                    call appendToString(output(iOutput), i, "(I12)")
                    ! add nodal point number
                    call appendToString(output(iOutput), v(iOutput,i), "(I12)")
                enddo

                if (isselected(1).eq.1) then

                    if (printweights) then ! point 1 matches, opposite spanned by 2,3 and 4 -> side nr 4
                          thissign = elm_sign(v(1,i), v(2,i), v(3,i), v(4,i), 4, nmatch)
!                         thissign = - thissign
                        call appendToString(output(iLine), thisSign, "(I5)")
                    endif

                    if (printcoords) then
                        ! add coordinates of points 1
                        call appendToString(output(iLine), coords(1,v(1,i)), "(f25.17)" )
                        call appendToString(output(iLine), coords(2,v(1,i)), "(f25.17)" )
                        call appendToString(output(iLine), coords(3,v(2,i)), "(f25.17)" )
                    endif

                    write(*,"(a)") trim(output(iLine))
                endif

                iLine = iLine + 1

                if (isselected(2).eq.1) then

                    if (printweights) then ! point 2 matches, opposite spanned by 1,3 and 4 -> side nr 3
                          thissign = elm_sign(v(1,i), v(2,i), v(3,i), v(4,i), 3, nmatch)
!                         thissign = - thissign
                        call appendToString(output(iLine), thisSign, "(I5)")
                    endif

                    if (printcoords) then
                        ! add coordinates of points 2
                        call appendToString(output(iLine), coords(1,v(2,i)), "(f25.17)" )
                        call appendToString(output(iLine), coords(2,v(2,i)), "(f25.17)" )
                        call appendToString(output(iLine), coords(3,v(2,i)), "(f25.17)" )
                    endif

                    write(*,"(a)") trim(output(iLine))
                endif

                iLine = iLine + 1

                if (isselected(3).eq.1) then

                    if (printweights) then ! point 3 matches, opposite spanned by 1,2 and 4 -> side nr 1
                          thissign = elm_sign(v(1,i), v(2,i), v(3,i), v(4,i), 1, nmatch)
!                         thissign = - thissign
                        call appendToString(output(iLine), thisSign, "(I5)")
                    endif

                    if (printcoords) then
                        ! add coordinates of points 3
                        call appendToString(output(iLine), coords(1,v(3,i)), "(f25.17)" )
                        call appendToString(output(iLine), coords(2,v(3,i)), "(f25.17)" )
                        call appendToString(output(iLine), coords(3,v(3,i)), "(f25.17)" )
                    endif

                    write(*,"(a)") trim(output(iLine))
                endif

                iLine = iLine + 1

                if (isselected(4).eq.1) then
                    if (printweights) then ! point 4 matches, opposite spanned by 1,2,3 -> side nr 2
                          thissign = elm_sign(v(1,i), v(2,i), v(3,i), v(4,i), 2, nmatch)
!                         thissign = - thissign
                        call appendToString(output(iLine), thisSign, "(I5)")
                    endif
                    if (printcoords) then
                        ! add coordinates of points 4
                        call appendToString(output(iLine), coords(1,v(4,i)), "(f25.17)" )
                        call appendToString(output(iLine), coords(2,v(4,i)), "(f25.17)" )
                        call appendToString(output(iLine), coords(3,v(4,i)), "(f25.17)" )
                    endif

                    write(*,"(a)") trim(output(iLine))
                endif

                iLine = iLine + 1

            else
                 !  a tetrahedron that has 4 point on the boundary, can have either 2 or 3 surfaces on the boundary
                 !  (assuming that the boundary is not a small enclosure of the single tetrahedron).
                 !  3 when it is a point of the boundary (punt van frietzak)
                 !  2 when it is in a fold (d4-tje midden in het adnd boek)

                if (nmatch.eq.3) then

                    call appendToString(output(1), i, "(I12)")

                    ! this will be the most common case. An element is adjacent with 

                    if     (isselected(2).eq.1 .and. &
                            isselected(3).eq.1 .and. &
                            isselected(4).eq.1 .and. &
                            isselected(1).eq.0) then
                        ! side 4
                        call appendToString(output(1), 4, "(I5)")

                        if (printweights) then 
                            thissign = elm_sign(v(1,i), v(2,i), v(3,i), v(4,i), 4, nmatch)
                            call appendToString(output(1), thisSign, "(I5)")
                        endif

                        if (printcoords) then
                            touchCoords = (coords(:,v(2,i)) + coords(:,v(3,i)) + coords(:,v(4,i))) / 3d0
                            call appendToString(output(1), touchCoords(1), "(f25.17)")
                            call appendToString(output(1), touchCoords(2), "(f25.17)")
                            call appendToString(output(1), touchCoords(3), "(f25.17)")
                        endif

                    else if (isselected(1).eq.1 .and. &
                             isselected(3).eq.1 .and. &
                             isselected(4).eq.1 .and. &
                             isselected(2).eq.0) then
                        ! side 3
                        call appendToString(output(1), 3, "(I5)")

                        if (printweights) then
                            thissign = elm_sign(v(1,i), v(2,i), v(3,i), v(4,i), 3, nmatch)
                            call appendToString(output(1), thisSign, "(I5)")
                        endif

                        if (printcoords) then
                            touchCoords = (coords(:,v(1,i)) + coords(:,v(3,i)) + coords(:,v(4,i))) / 3d0
                            call appendToString(output(1), touchCoords(1), "(f25.17)")
                            call appendToString(output(1), touchCoords(2), "(f25.17)")
                            call appendToString(output(1), touchCoords(3), "(f25.17)")
                        endif

                    else if (isselected(2).eq.1 .and. &
                             isselected(1).eq.1 .and. &
                             isselected(4).eq.1 .and. &
                             isselected(3).eq.0) then
                        ! side 1
                        call appendToString(output(1), 1, "(I5)")

                        if (printweights) then
                            thissign = elm_sign(v(1,i), v(2,i), v(3,i), v(4,i), 1, nmatch)
                            call appendToString(output(1), thisSign, "(I5)")
                        endif


                        if (printcoords) then
                            touchCoords = (coords(:,v(2,i)) + coords(:,v(1,i)) + coords(:,v(4,i))) / 3d0
                            call appendToString(output(1), touchCoords(1), "(f25.17)")
                            call appendToString(output(1), touchCoords(2), "(f25.17)")
                            call appendToString(output(1), touchCoords(3), "(f25.17)")
                        endif


                    else if (isselected(2).eq.1 .and. &
                             isselected(3).eq.1 .and. &
                             isselected(1).eq.1 .and. &
                             isselected(4).eq.0) then
                        ! side 2
                        call appendToString(output(1), 2, "(I5)")

                        if (printweights) then
                            thissign = elm_sign(v(1,i), v(2,i), v(3,i), v(4,i), 2, nmatch)
                            call appendToString(output(1), thisSign, "(I5)")
                        endif

                        if (printcoords) then
                            touchCoords = (coords(:,v(2,i)) + coords(:,v(3,i)) + coords(:,v(1,i))) / 3d0
                            call appendToString(output(1), touchCoords(1), "(f25.17)")
                            call appendToString(output(1), touchCoords(2), "(f25.17)")
                            call appendToString(output(1), touchCoords(3), "(f25.17)")
                        endif

                    else
                        write(*,*) "nmatch = 3, but the number of matched nodes is not 3."
                        write(*,*) "This is not supposed to happen. Contact model support."
                        write(*,*) "element: ", i
                        write(*,*) "selected: ", isselected(1:4)
                        stop "leaving elmside"
                    endif

                    write(*,"(a)") trim(output(1))


                else
                    ! matched no side
                endif
            endif


         endif
         !************************************************

    enddo

end subroutine

subroutine basis()

    use globals, only: coords, &
                       ux, uy, uz, &          ! vector to span up the base
                       vx, vy, vz, &          ! vector to span up the base, only for 3D
                       ox, oy, oz, normo, &   ! orthogonal vector to the base, and its norm
                       pm_sign, &             ! sign of the ortho vec +/- 1
                       eps, ndimensions, &
                       nvertices, vertices, & ! to check quality of anchor point
                       ax, ay, az, &          ! anchor point
                       iecho, &
                       xmax,ymax,zmax,xmin,ymin,zmin ! domain extrema
                       
    implicit none

    integer :: b1, b2, b3

    double precision :: cx, cy, cz  ! center of the fault. Average of outside base points

    double precision, external :: distanceRatio
    double precision :: plusRatio, minRatio, bestRatio
    double precision :: offset, bestOffset

    integer :: iVertex, iAnchor
    
    double precision :: domainRange

    logical :: faultIsPlane ! or almost a plane

    b1 = 0
    b2 = 0
    b3 = 0

    if (ndimensions.eq.2) then
        call findbasis2D(b1, b2)
        if (iecho.eq.1) then
            write(*,*) "Selected points", b1, b2
        endif
    else
        call findbasis3D(b1, b2, b3)
        if (iecho.eq.1) then
            write(*,*) "Selected points", b1, b2, b3
        endif

    endif

    ! from the basis, determine the center of the fault

    if (ndimensions.eq.2) then
        cx = 0.5 * (coords(1,b1) + coords(1,b2))
        cy = 0.5 * (coords(2,b1) + coords(2,b2))
    else
        cx = (coords(1,b1) + coords(1,b2) + coords(1,b3)) / 3d0
        cy = (coords(2,b1) + coords(2,b2) + coords(2,b3)) / 3d0
        cz = (coords(3,b1) + coords(3,b2) + coords(3,b3)) / 3d0
    endif

    ! and and orthogonal vector

    if (ndimensions.eq.2) then
        ox = coords(2,b1) - coords(2,b2)
        oy = coords(1,b2) - coords(1,b1)

        normo = sqrt(ox*ox + oy*oy)

    else ! dimensions = 3
        !    |  i   j   k   |
        !    |  ux  uy  uz  |
        !    |  vx  vy  vz  |
        !  with u = p2 - p1
        !  and  v = p3 - p1
        ux = coords(1,b2) - coords(1,b1)
        uy = coords(2,b2) - coords(2,b1)
        uz = coords(3,b2) - coords(3,b1)
        vx = coords(1,b3) - coords(1,b1)
        vy = coords(2,b3) - coords(2,b1)
        vz = coords(3,b3) - coords(3,b1)

!        write(*,*) "basis vector u: ", ux, uy, uz
!        write(*,*) "basis vector v: ", vx, vy, vz


        ox = uy * vz - uz * vy
        oy = uz * vx - ux * vz
        oz = ux * vy - uy * vx


        normo = sqrt(ox*ox + oy*oy + oz*oz)

    endif


!    write(*,*) "Basis says: ortho: ", ox, oy, oz

    if (normo.lt.eps .and. b1.ne.b2) then
        write(*,*) "Norm of orthogonal vector too small"
        write(*,*) "This is not supposed to happen"
        write(*,*) "Contact model support"
    endif

    ! find a point on the line through the c points, spanned by the
    ! orthogonal vector, for which the minimum and maximum distance are
    ! minimal.

    ! As this does not be done exactly, we do this by trial and error.

    ! set initial achor point to center
    ax = cx
    ay = cy
    az = 0d0
    if (ndimensions.eq.3) then
        az = cz
    endif


    bestRatio = abs(distanceRatio(ax,ay,az))
    offset = 0.01

    domainRange = xmax+ymax-xmin-ymin 
    if (ndimensions.eq.3) then
        domainRange = domainRange + zmax-zmin
    endif

    faultIsPlane = .false.

    do while (offset .lt. 10d0*domainRange)
        
        offset = offset * 2d0

        plusRatio = abs(distanceRatio(cx + offset*ox, &
                                      cy + offset*oy, &
                                      cz + offset*oz))

        minRatio  = abs(distanceRatio(cx - offset*ox, &
                                      cy - offset*oy, &
                                      cz - offset*oz))



!        write(*,*) "plusratio: ", plusRatio
!        write(*,*) "minratio:  ", minRatio


        if (plusRatio .lt. bestRatio) then
!            write(*,*) "plus is better. Adjusting best ratio to", plusRatio
            bestOffset = offset
            bestRatio = plusRatio
            if (offset .ge. 10d0*domainRange) then
                faultIsPlane = .true.
            endif
        endif

        if (minRatio .lt. bestRatio) then
!               write(*,*) "min is better. Switching signs and adjusting best ratio to", plusRatio
            bestOffset = -offset  
            bestRatio = minRatio
            if (offset .ge. 10d0*domainRange) then
                faultIsPlane = .true.
            endif
        endif

!    write(*,*) "Basis says: iteration: ", bestOffset

    enddo

    ax = cx + bestOffset*ox
    ay = cy + bestOffset*oy
    az = 0d0
    if (ndimensions.eq.3) then
        az = cz + bestOffset*oz
    endif


    if (faultIsPlane) then
        ! aim the vector in the proper direction, in case that is not the case.
        if (abs(ax) .gt. abs(ay) .and. abs(ax) .gt.    abs(az)) then
            ! ax is the dominant
            if (ax .lt. 0) then
                ax = -ax
                ay = -ay
                az = -az
            endif
        else if (abs(ay) .gt. abs(ax) .and. abs(ay) .gt. abs(az)) then
            ! ay is the dominant
            if (ay .lt. 0) then
                ax = -ax
                ay = -ay
                az = -az
            endif
        else if (abs(az) .gt. abs(ay) .and. abs(az) .gt. abs(ax)) then
            ! az is the dominant
            if (az .lt. 0) then
                ax = -ax
                ay = -ay
                az = -az
            endif
        endif

    endif

!    write(*,*) "Basis says: bestOffset: ", bestOffset
!    write(*,*) "Basis says: anchor point: ", ax, ay, az 

end subroutine


double precision function distanceRatio(ax,ay,az)
    ! helper routine for basis (above)
    ! compute ratio of distance to furthest point in collection
    ! to distance to closest point in collection
    use globals, only: coords, &
                       vertices, nvertices, &
                       ndimensions

    implicit none

    double precision :: ax,ay,az
    double precision :: distPlus, distMin

    double precision :: maxDistance, MinDistance
    double precision :: thisDistance    

    integer          :: iVertex

    maxDistance = 0
    minDistance = 9e20



    do iVertex = 1,nvertices
        thisDistance = (ax - coords(1,vertices(iVertex)))**2 + &
                       (ay - coords(2,vertices(iVertex)))**2

        if (ndimensions.eq.3) then
            thisDistance = thisDistance + (az - coords(3,vertices(iVertex)))**2
        endif

        thisDistance = sqrt(thisDistance)


        if (thisDistance.gt.maxDistance) then
            maxDistance = thisDistance
        endif
        if (thisDistance.lt.minDistance) then
            minDistance = thisDistance
        endif
        
    enddo

    distanceRatio = maxDistance - minDistance

!    write(*,*) "computed ratio: ", maxDistance, minDistance, distanceRatio


end function

!----------------------------------------------------------------



subroutine mutualneighbours(v1,v2,v3)

    use globals, only: ndimensions, &
                       nneighbours, &
                       neighbourIDs, &
                       neighbours, &
                       nallvertices, &
                       maxneighbours2D, &
                       maxneighbours3D


    implicit none
    ! subroutine returns a list of common neighbours of v1 and v2 when in 2D,
    ! and of v1, v2, and v3 when in 3D. There can be no more than 2.
    ! note there must be either one mutual neighbour, when on the edge of the domain,
    ! or two, when inside the domain, on a fault line.
    integer                                           :: maxneighbours
    integer                                           :: v1, v2, v3  ! , ndimensions
!    integer, dimension(nallvertices)                  :: nneighbours
!    integer, dimension(nallvertices, maxneighbours)   :: neighbourIDs
!    integer, dimension(2)                             :: neighbours

    integer                            :: i, pos, candID
    integer, dimension(:), allocatable :: candidates  

    integer :: numinarray  ! function below

    if (ndimensions.eq.3) then
        maxneighbours = maxneighbours3D
    else if (ndimensions.eq.2) then
        maxneighbours = maxneighbours2D
    endif

    neighbours = 0

    allocate(candidates(nneighbours(v1)))

    ! create an array with zeros, corresponding to the neighbours.
    ! at every other point, increase this number when a point is matched
    ! for two points, it should be 1 to have a mutual neighbour. for 3 points, it should be 2.

    ! noooo, can be more than 2:

    !                      B

    !                   .  * .
    !                .    / \  .   
    !             .      /   \   .
    !          .        /     \    .
    !         *--------*-------*----*
    !            .      \     /   .
    !               .    \   /  .
    !                  .  \ / .
    !                      *

    !                      A    

    ! for example, point A and B have four mutual neighbors.

    do i=1,nneighbours(v1)
         candidates(i) = 0
    enddo

    do i=1,nneighbours(v2)
         pos = numinarray(neighbourIDs(v2,i), neighbourIDs(v1,1:nneighbours(v1)), nneighbours(v1))
         if (pos.gt.0) then
             candidates(pos) = candidates(pos) + 1
         endif
    enddo

    if (ndimensions.eq.3) then
         pos = numinarray(neighbourIDs(v3,i), neighbourIDs(v1,1:nneighbours(v1)), nneighbours(v1))
         if (pos.gt.0) then
             candidates(pos) = candidates(pos) + 1
         endif
    endif

    ! extract the neighbour IDs and return them to the caller in 'neighbours'
    candID = 1
    do i=1,nneighbours(v1)
        if      (ndimensions .eq. 2 .and. candidates(i).eq.1) then
!            write(*,*) "checking neighbours: ", candID, "of",size(neighbours, 1)
            neighbours(candID) = neighbourIDs(v1,i)
            candID = candID + 1
        else if (ndimensions .eq. 3 .and. candidates(i).eq.2) then
            neighbours(candID) = neighbourIDs(v1,i)
            candID = candID + 1
        endif
    enddo
    
end subroutine

!****** subroutines to determine the sign.
! The sign is determined as follows:
! The two (in 2D) points or three (in 3D) that are on the edge of the domain,
! form a basis. Orthogonal vector to this basis determines the side.
! If the vector from the middle of the edge on the boundary to the other point
! has an angle less than 90 degrees, the sign is the same.
! Otherwise, it gets the opposite sign.

subroutine sidesign2D(x1, y1, x2, y2)
    implicit none
    
    double precision:: x1, y1, x2, y2 ! coordinates of points on the edge


end subroutine



! The findbasis routines take a 1000 random points (or less if there are less) in 2D
! or 100 in 3D
! from the selection list, and finds the points that are furthest apart.
! There will be a 0.5 million comparisons at most. This will take only a little time.
! in 2D probably less, because a line will rarely consist of as many as a 1000 points

! 2 in 2D, and 3 in 3D, to describe the plane.
! This allows the selection list to describe a non-straight edge sphere, although 
! there should not be too much differentiation.
! For example, a completely spherical edge is too much.

! the subs return the indices of the points that describethe basis p1 and p2, (and p3 in 3D)

subroutine findbasis2D(p1, p2)
    use globals, only: coords,nvertices,vertices

    implicit none

    integer, allocatable, dimension(:) ::  selectbasefrom
    integer :: i, j, nsearch
    integer :: p1, p2
    double precision :: maxdist, dist

    if (nvertices .lt. 1000) then
        ! select from the entire list
        allocate(selectbasefrom(nvertices))
        do i=1,nvertices
            selectbasefrom(i) = vertices(i)
        enddo
        nsearch = nvertices
    else
        ! randomly select a thousand points from the list
        allocate(selectbasefrom(1000))
        nsearch = 1000
        call selectsublistfromlist(nvertices, vertices, 1000, selectbasefrom)
    endif
    
    ! now search through these points to find the furthest apart pair
    maxdist = 0
    do i=1,nsearch-1
        do j=i+1, nsearch
            ! actually square of distance, but that is irrelevant
            dist = (coords(1,selectbasefrom(i)) - coords(1,selectbasefrom(j)))**2.0 + (coords(2,selectbasefrom(i)) - coords(2,selectbasefrom(j)))**2.0
!            write(*,*) 'comparing ', selectbasefrom(i), ' and ', selectbasefrom(j), ' distance: ', dist
            if (dist .gt. maxdist) then
!                write(*,*) 'a new record, yay: ', p1, p2, dist
                maxdist = dist
                p1 = selectbasefrom(i)
                p2 = selectbasefrom(j)
            endif
        enddo
    enddo



end subroutine


subroutine findbasis3D(p1, p2, p3)

    use globals, only: coords,nvertices,vertices

    implicit none

    double precision :: triangle_area

! findbasis3D is very much like the 2D variety.
! but in stead of using distance between two points to select,
! it uses max triangle size to select 3 points.
! (and it takes only 100 points)
! as 3D planes are likely to contain more than a 100 points,
! the subselection is likely to be used, here.

    integer, allocatable, dimension(:) ::  selectbasefrom
    integer :: i, j, k, nsearch
    integer :: p1, p2, p3
    double precision :: maxarea, area

    double precision :: x1, x2, x3, y1, y2, y3, z1, z2, z3

    if (nvertices .lt. 100) then
        ! select from the entire list
        allocate(selectbasefrom(nvertices))
        do i=1,nvertices
            selectbasefrom(i) = vertices(i)
        enddo
        nsearch = nvertices
    else
        ! randomly select a thousand points from the list
        allocate(selectbasefrom(100))
        nsearch = 100
        call selectsublistfromlist(nvertices, vertices, 100, selectbasefrom)
    endif
 
    ! now search through these points to find the furthest apart pair
    maxarea = 0
    do i=1,nsearch-2
        x1 = coords(1,selectbasefrom(i))
        y1 = coords(2,selectbasefrom(i))
        z1 = coords(3,selectbasefrom(i))
        do j=i+1, nsearch-1
            x2 = coords(1,selectbasefrom(j))
            y2 = coords(2,selectbasefrom(j))
            z2 = coords(3,selectbasefrom(j))
            do k=j+1, nsearch
                x3 = coords(1,selectbasefrom(k))
                y3 = coords(2,selectbasefrom(k))
                z3 = coords(3,selectbasefrom(k))
                ! actually square of distance, but that is irrelevant
                ! area = 1/2 det | u x v |
                ! where u is chosen to be p1 - p2 and v to be p1 - p3
                
                area = triangle_area(x1,y1,z1, x2,y2,z2, x3,y3,z3 )

                if (area .gt. maxarea) then
                    maxarea = area
                    p1 = selectbasefrom(i)
                    p2 = selectbasefrom(j)
                    p3 = selectbasefrom(k)
                endif
            enddo
        enddo
    enddo

end subroutine


subroutine  selectsublistfromlist(nlist, list, nselect, selected)
    ! returns a list of n elements
    implicit none

    integer :: i

    integer :: nlist, nselect
    integer, dimension(nlist) :: list
    integer, dimension(nselect) :: selected

    integer :: selectID

    integer, dimension(nlist) :: templist

    real    :: randomHarvest

    if (nselect .gt. nlist .or. nlist.lt.1) then
        ! prevent infinite loop
        write(*,*) 'IEPS. selectsublistfromlist says: list not correct.'
        write(*,*) 'This should not happen; contact model support'
        stop
    endif

    do i=1, nlist
        templist(i) = list(i)
    enddo

    do i=1,nselect
        ! select a number from the list
        call random_number(randomHarvest)
        selectID = ceiling((nlist-i+1) * randomHarvest)
        selected(i) = templist(selectID)
        ! and replace this number by the one on the end of the list
        if (selectID.lt.nlist-i+1) then
            templist(selectID) = templist(nlist-i+1)
        endif
    enddo

end subroutine


double precision function triangle_area(x1,y1,z1, x2,y2,z2, x3,y3,z3 )
    ! computes surface area of a triangle
    implicit none

    double precision :: x1, x2, x3, y1, y2, y3, z1, z2, z3
    double precision :: l1, l2, l3, s
    double precision :: det

    ! Using Heron's formula
    ! let length of the triangle sides be known, l1, l2, l3, then
    ! s = (l1 + l2 + l3) / 2
    ! area = sqrt( s * (s-l1) * (s-l2) * (s-l3) )

    l3 = sqrt((x1 - x2)**2.0 + (y1 - y2)**2.0 + (z1 - z2)**2.0)
    l2 = sqrt((x1 - x3)**2.0 + (y1 - y3)**2.0 + (z1 - z3)**2.0)
    l1 = sqrt((x2 - x3)**2.0 + (y2 - y3)**2.0 + (z2 - z3)**2.0)

    s = 0.5 * (l1 + l2 + l3)

    triangle_area = sqrt (s * (s-l1) * (s - l2) * (s - l3))

end function

integer function signOfElemWithOneNodeAdjacent2D(p1, p2, p3)
    ! this function checks the side of the elements with only
    ! a single node on the line/curve.
    ! The idea behind this is that elements at the ends of the
    ! line/curve are assigned the proper weights.

    ! so far only used/tested in 2D

    use globals, only: coords, & ! coordinates of all points
                       ox, oy, & ! orthogonal basis vectors
                       pm_sign ! sign of the orthogonal vector
    implicit none

    integer          :: p1, p2, p3

    double precision :: baseX, baseY
    double precision :: elemCenterX, elemCenterY    

    baseX = -oy
    baseY =  ox

    write(0,*) "basis: ", baseX, baseY

    elemCenterX = 0.33333 * (coords(1,p1) + coords(1,p2) + coords(1,p3))
    elemCenterY = 0.33333 * (coords(2,p1) + coords(2,p2) + coords(2,p3))

    
    signOfElemWithOneNodeAdjacent2D = 1

end function


integer function elm_sign(p1, p2, p3, p4, side_idx, nMatch)
    ! computes the sign of a side of an element.
    ! use to distinguish elements on different side of a fraction
    ! input is the vectr from the center of the edge that lies on
    ! th boundary, to the other point.
    use globals, only: ndimensions, & 
                       pm_sign, &
                       ox, oy, oz, &
                       ax, ay, az, &
                       normo, &
                       coords


    implicit none

    integer          :: p1, p2, p3, p4
    integer          :: side_idx
    double precision :: vx, vy, vz
    double precision :: toacos, theta
    integer          :: nMatch

    double precision  :: lineDist, pointDist
    integer, external :: setSign, setSign3D
    double precision  :: xm, ym
    ! Compute vector for middle of side that is adjacent to
    ! the boundary to the remaining vertex.

    ! again let us rephrase the convention for sides
    ! 2D:
    !   face 1: vertices 1 2
    !        2:          2 3
    !        3:          3 1 
    ! 3D:
    !   face 1: vertices 1 2 4
    !        2:          1 3 2
    !        3:          1 4 3
    !        4:          2 3 4    




    if (ndimensions.eq.2) then

        ! compute the vector to the point that 
!        write(*,*) "setting sign of element", side_idx

        if (side_idx.eq.1) then ! side 1, vertices 1, 2
            xm = 0.5 * (coords(1,p1) + coords(1,p2))
            ym = 0.5 * (coords(2,p1) + coords(2,p2))
!            lineDist  = sqrt( (0.5 * (coords(1,p1) + coords(1,p2)) - ax)**2 + &
!                              (0.5 * (coords(2,p1) + coords(2,p2)) - ay)**2    )
!            pointDist = sqrt( (coords(1,p3) - ax)**2 + &
!                              (coords(2,p3) - ay)**2 )

!            elm_sign = setSign(nMatch, lineDist, pointDist)
            elm_sign = setSign(nMatch, xm, ym, coords(1,p3), coords(2,p3))

        else if (side_idx.eq.2) then ! side 2, vertices 2, 3
!            lineDist  = sqrt( (0.5 * (coords(1,p2) + coords(1,p3)) - ax)**2 + &
!                              (0.5 * (coords(2,p2) + coords(2,p3)) - ay)**2    )
!            pointDist = sqrt( (coords(1,p1) - ax)**2 + &
!                              (coords(2,p1) - ay)**2 )

!            elm_sign = setSign(nMatch, lineDist, pointDist)

            xm = 0.5 * (coords(1,p3) + coords(1,p2))
            ym = 0.5 * (coords(2,p3) + coords(2,p2))
            elm_sign = setSign(nMatch, xm, ym, coords(1,p1), coords(2,p1))



        else if (side_idx.eq.3) then ! side 3, vertices 1, 3
!            lineDist  = sqrt( (0.5 * (coords(1,p1) + coords(1,p3)) - ax)**2 + &
!                              (0.5 * (coords(2,p1) + coords(2,p3)) - ay)**2    )
!            pointDist = sqrt( (coords(1,p2) - ax)**2 + &
!                              (coords(2,p2) - ay)**2 )

!            elm_sign = setSign(nMatch, lineDist, pointDist)

            xm = 0.5 * (coords(1,p1) + coords(1,p3))
            ym = 0.5 * (coords(2,p1) + coords(2,p3))
            elm_sign = setSign(nMatch, xm, ym, coords(1,p2), coords(2,p2))

        else
            write(*,*) "elm_sign says: illegal side nr. This should not happen; contact model support."
        endif

!        write(*,*) "element has middle: ", xm, ym, "found sign", elm_sign


    else if (ndimensions.eq.3) then

        if (side_idx.eq.1) then  ! side 1; vertices 1, 2, 4
            vx = coords(1,p3) - 0.333333 * (coords(1,p1) + coords(1,p2) + coords(1,p4))
            vy = coords(2,p3) - 0.333333 * (coords(2,p1) + coords(2,p2) + coords(2,p4))
            vz = coords(3,p3) - 0.333333 * (coords(3,p1) + coords(3,p2) + coords(3,p4))
            elm_sign = setSign3D(vx,vy,vz,coords(1,p3),coords(2,p3),coords(3,p3))

        else if (side_idx.eq.2) then ! side 2; vertices 1, 2, 3
            vx = coords(1,p4) - 0.333333 * (coords(1,p1) + coords(1,p2) + coords(1,p3))
            vy = coords(2,p4) - 0.333333 * (coords(2,p1) + coords(2,p2) + coords(2,p3))
            vz = coords(3,p4) - 0.333333 * (coords(3,p1) + coords(3,p2) + coords(3,p3))
            elm_sign = setSign3D(vx,vy,vz,coords(1,p4),coords(2,p4),coords(3,p4))

        else if (side_idx.eq.3) then ! side 3; vertices 1, 4, 3
            vx = coords(1,p2) - 0.333333 * (coords(1,p1) + coords(1,p4) + coords(1,p3))
            vy = coords(2,p2) - 0.333333 * (coords(2,p1) + coords(2,p4) + coords(2,p3))
            vz = coords(3,p2) - 0.333333 * (coords(3,p1) + coords(3,p4) + coords(3,p3))
            elm_sign = setSign3D(vx,vy,vz,coords(1,p2),coords(2,p2),coords(3,p2))

        else if (side_idx.eq.4) then ! side 4; vertices 2, 3, 4
            vx = coords(1,p1) - 0.333333 * (coords(1,p2) + coords(1,p3) + coords(1,p4))
            vy = coords(2,p1) - 0.333333 * (coords(2,p2) + coords(2,p3) + coords(2,p4))
            vz = coords(3,p1) - 0.333333 * (coords(3,p2) + coords(3,p3) + coords(3,p4))
            elm_sign = setSign3D(vx,vy,vz,coords(1,p1),coords(2,p1),coords(3,p1))

        else
            write(*,*) 'elm_sign says: illegal side nr. This should not happen; contact model support.'
        endif


    endif


end function

integer function setSign(nMatch, faceMidx, faceMidy, pointx, pointy)
! helper function for elm_sign

use globals, only: ax, ay, az          ! anchor point


implicit none

integer :: nMatch
double precision :: faceMidx, faceMidy, pointx, pointy

double precision :: anchorVec(2)
double precision :: elementVec(2)

double precision, external :: DOT, vectorLength
double precision :: takeAcosOf

double precision :: angle


if (nMatch .eq. 1) then
    anchorVec(1) = ax - pointx
    anchorVec(2) = ay - pointy

    elementVec(1) = faceMidx - pointx
    elementVec(2) = faceMidy - pointy

else if (nMatch .eq. 2) then
    anchorVec(1) = ax - faceMidx
    anchorVec(2) = ay - faceMidy

    elementVec(1) = pointx - faceMidx
    elementVec(2) = pointy - faceMidy

else
    write(0,*) "Element has three sides on fault. This is not well implemented"
    write(0,*) "Please refine mesh, or contact model support"
endif


takeAcosOf = ( dot(anchorVec,elementVec,2) / &
           (vectorLength(anchorVec,2) * vectorLength(elementVec,2)))

! make sure the number to be used for an arccos is not outside the [-1,1] domain
if(takeAcosOf .lt.-1) then
    takeAcosOf = -0.9999
else if (takeAcosOf .gt.1) then
    takeAcosOf = 0.99999
endif

angle = acos(takeAcosOf) * 180.0 / 3.1415926535

if  (abs(angle) .gt. 90.0) then
    ! inside of the curve
    setSign = -1
else
    ! outside of the curve
    setSign = 1
endif

end function


! this is done a bit wishy washy in 3D
integer function setSign3D(vx,vy,vz,pointx, pointy, pointz)


use globals, only: ax, ay, az          ! anchor point

implicit none

double precision :: vx,vy,vz
double precision :: pointx, pointy, pointz
double precision :: anchorVec(3)
double precision :: elementVec(3)
double precision :: angle

double precision, external :: DOT, vectorLength


anchorVec(1)=ax-pointx
anchorVec(2)=ay-pointy
anchorVec(3)=az-pointz

elementVec(1) = vx
elementVec(2) =    vy
elementVec(3) =    vz


angle = acos ( dot(anchorVec,elementVec,3) / &
           (vectorLength(anchorVec,3) * vectorLength(elementVec,3)))

angle =  angle * 180.0 / 3.1415926535

if  (angle .gt. 90.0) then
    ! inside of the curve
    setSign3D = -1
else
    ! outside of the curve
    setSign3D = 1
endif


end function


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

