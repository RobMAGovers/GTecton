module conversions

public :: string2int
public :: string2double

contains


integer function nCommasInString(argumentString)

    implicit none

    character(len=256) :: argumentString
    integer            :: i

    nCommasInString = 0
    do i=1,len(trim(argumentString))
        if (argumentString(i:i).eq.",") then
            nCommasInString = nCommasInString + 1
        endif
    enddo

end function



logical function oneCharIsNumeric(oneChar)

    implicit none

    character(len=1) :: oneChar

    if (oneChar .eq."1"  .or. &
        oneChar .eq."2"  .or. &
        oneChar .eq."3"  .or. &
        oneChar .eq."4"  .or. &
        oneChar .eq."5"  .or. &
        oneChar .eq."6"  .or. &
        oneChar .eq."7"  .or. &
        oneChar .eq."8"  .or. &
        oneChar .eq."9"  .or. &
        oneChar .eq."0") then
        oneCharIsNumeric = .true.
    else
        oneCharIsNumeric = .false.
    endif

end function



integer function string2int (numberstring)
 
    implicit none

    character(len=*) :: numberstring

    read (numberstring, *) string2int

    return
end function


double precision function string2double (numberstring)

    implicit none
    
    character(len=*)    :: numberstring

    read (numberstring, *) string2double

    return
end function

subroutine resize_dblarray(array, oldsize, newsize)
 ! change an array to a new size.
    implicit none

    integer                  :: i, oldsize, newsize
    double precision, dimension(:), allocatable :: array    
    double precision, dimension(:), allocatable :: temparray

    ! copy the old array into the remporary array
    if (oldsize > newsize) then
        ! array is made smaller, backup only required data
        allocate(temparray(newsize))
        do i=1,newsize
            temparray(i) = array(i)
        enddo
    else
        ! array is made larger, backup all data
        allocate(temparray(oldsize))
        do i=1,oldsize
            temparray(i) = array(i)
        enddo
    endif

    ! resize the original array
    deallocate(array)
    allocate(array(newsize))

    ! copy the data back
    if (oldsize > newsize) then
        ! array is made smaller, copy only required data back
        do i=1,newsize
            array(i) = temparray(i)
        enddo
    else
        ! array is made larger, copy all the data back
        do i=1,oldsize
            array(i) = temparray(i)
        enddo
    endif

    ! remove the temporary data
    deallocate(temparray)

end subroutine


end module
