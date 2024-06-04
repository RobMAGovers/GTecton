module stringHandling

implicit none

public  :: appendToString
private :: appendFloatToString
private :: appendIntegerToString

interface appendToString
! main call to the module,
! depending on type and dimension of the argumented array
! a procedure is chosen.
module procedure appendFloatToString, &
                 appendIntegerToString

end interface

contains

subroutine appendFloatToString(string, float, format)

    use globals, only: stringLength

    implicit none

    character(len=stringLength) :: string
    double precision   :: float
    character(len=*)   :: format

    character(len=100) :: newString
    integer            :: currentLength
    integer            :: newLength

    call cleanString(newString, 100)
    currentLength = len(trim(string))
    write(newString, trim(format)) float
    newLength = len(trim(newString))
    string(currentLength+1:currentLength+1+newLength) = newString

end subroutine

subroutine appendIntegerToString(string, number, format)

    use globals, only: stringLength

    implicit none

    character(len=stringLength) :: string
    integer            :: number
    character(len=*)   :: format


    character(len=100) :: newString
    integer            :: currentLength
    integer            :: newLength

    call cleanString(newString, 100)
    currentLength = len(trim(string))
    write(newString, trim(format)) number
    newLength = len(trim(newString))
    string(currentLength+1:currentLength+1+newLength) = newString

end subroutine

subroutine cleanString(string, length)

    implicit none

    integer    :: length
    character(len=length) :: string

    integer :: i

    do i = 1,length
        string(i:i) = " "
    enddo

end subroutine

subroutine cleanOutputStrings()

    use globals, only: output, &
                       stringLength

    implicit none

    integer :: iOutputLine

    do iOutputLine = 1,4
        call cleanString(output(iOutputLine), stringLength)
    enddo

end subroutine


end module

