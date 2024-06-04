!****** subroutines for printing the output data

! the output routines need:
! * variables about the element *
! - the element ID
! - the side ID with which is borders; 0 if bordering with point
! - the vertex ID with which it borders; could be one or two
! - the weight
! - the coordinates - of the point bordering. if borderdering with a point
!                   - of the middle of the side bordering, if bordering with a side
! * from the command line, to determine format; passed through the globals module *
! - globals % printcoords
! - globals % sideonly
! - whether opn-com


subroutine output2D(elemID, vv1, vv2, xx, yy, sideID, weight)
    use globals, only: sideonly, &
                       printcoords, &
                       opn, &
                       printweights, &
                       invertsign, &
                       v, &
                       endpoints, &
                       coords

    implicit none

    integer          :: elemID, vv1, vv2, sideID, weight
    double precision :: xx,yy



    if (invertsign) then
        weight = -weight
    endif


!    write(*,*) 'output2D called with: elemID, vv1, vv2, xx, yy, sideID, weight, endpoints: ', &
!                                      elemID, vv1, vv2, xx, yy, sideID, weight, endpoints


      if (sideonly) then
        ! when v2 = 0, this implies that only one point is included in the attachment,
        ! hence no side, hence exlcusion 
        if (vv2.ne.0) then
            if (endpoints) then
                if      (sideID.eq.1) then
                    write(*,'(I12,I5,2I12,4E14.6,I5)') elemID, sideID, v(1,elemID), v(2,elemID), &
                         coords(1,v(1,elemID)),coords(1,v(2,elemID)), coords(2,v(1,elemID)), coords(2,v(2,elemID)),weight
                else if (sideID.eq.2) then
                    write(*,'(I12,I5,2I12,4E14.6,I5)') elemID, sideID, v(2,elemID), v(3,elemID), &
                         coords(1,v(2,elemID)),coords(1,v(3,elemID)), coords(2,v(2,elemID)), coords(2,v(3,elemID)),weight
                else if (sideID.eq.3) then
                    write(*,'(I12,I5,2I12,4E14.6,I5)') elemID, sideID, v(1,elemID), v(3,elemID), &
                         coords(1,v(1,elemID)),coords(1,v(3,elemID)), coords(2,v(1,elemID)),coords(2,v(3,elemID)),weight
                else
                    write(*,*) 'sideID not equals 1, 2 or 3. This should not happen. Contact model support'
                endif
            else

                if (printweights) then
                    write(*,'(I12,3I5)') elemID, sideID, 0, weight
                else
                    write(*,'(I12,I5)') elemID, sideID
                endif
            endif

            if (printcoords) then
                if      (sideID.eq.1) then
                    write(*,'(I12,I5,2E14.6,I5)') elemID, sideID, coords(1,v(1,elemID))-coords(1,v(2,elemID)), coords(2,v(1,elemID))-coords(2,v(2,elemID)),weight
                else if (sideID.eq.2) then
                    write(*,'(I12,I5,2E14.6,I5)') elemID, sideID, coords(1,v(2,elemID))-coords(1,v(3,elemID)), coords(2,v(2,elemID))-coords(2,v(3,elemID)),weight
                else if (sideID.eq.3) then
                    write(*,'(I12,I5,2E14.6,I5)') elemID, sideID, coords(1,v(1,elemID))-coords(1,v(3,elemID)), coords(2,v(1,elemID))-coords(2,v(3,elemID)),weight
                else
                    write(*,*) 'sideID not equals 1, 2 or 3. This should not happen. Contact model support'
                endif

            endif

        endif

     else

        if (vv2.eq.0) then
            ! only a single point, so only one output line
            if (printweights) then
                if (printcoords) then
                    write(*,'(2I12,2I5,2f25.17)') elemID, vv1, 0, weight, coords(1,vv1), coords(2,vv1)
                else
                    write(*,'(2I12,2I5)') elemID, vv1, 0, weight
                endif
            else
                write(*,'(2I12)') elemID, vv1
            endif
        else
            ! one line for each of the points.
            if (printweights) then
                if (printcoords) then
                    write(*,'(2I12,2I5,2f25.17)') elemID, vv1, 0, weight, coords(1,vv1), coords(2,vv1)
                    write(*,'(2I12,2I5,2f25.17)') elemID, vv2, 0, weight, coords(1,vv2), coords(2,vv2)
                else
                    write(*,'(2I12,2I5)') elemID, vv1, 0, weight
                    write(*,'(2I12,2I5)') elemID, vv2, 0, weight
                endif
            else
                write(*,'(2I12)') elemID, vv1
                write(*,'(2I12)') elemID, vv2
            endif
        endif
    endif




end subroutine
