integer function numinarray(number, array, arraysize)
    ! appears brute force, but is actually quite fast.
    ! http://www.shocksolution.com/2011/03/finding-value-in-unordered-fortran-array/
    implicit none
    integer :: number, i, arraysize
    integer, dimension(arraysize) :: array

    do numinarray=1,size(array)
        if (array(numinarray).eq.number) then
            return
        endif
    enddo

    numinarray=0
    return

end function

