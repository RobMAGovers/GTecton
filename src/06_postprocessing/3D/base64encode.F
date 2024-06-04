module base64encode

!#############################################################################
!     below here subroutines to create the Base64 encoding for the double data
!   see http://en.wikipedia.org/wiki/Base64
!  VTK binary format wishes encoding in Base64, in stead of the regular base256,
!  because the XML standard does not support all the Base256 byte characters.
!  All the characters in the base64 standard are supported and hence safe to use.
!#############################################################################

! the first callee of the module EncodeArrayToFile.
! The arguments are a:
! - file pointer (assumed succesfully opened prior to call)
! - a 1d or 2d array of double precision numbers or integers
! - the precision in which the Base64 encoding is desired.
!   (1 for single precision; 2 for double)
!   (corresponding to types Float32 and Float64 for doubles, of Int32 an Int64 for integers)

! a secondary callee is RunSelfTest, which runs several tests.
! And the result is reported to the user. Program stops when test
! is completed

!**********************************************************************
! The algorithm is briefly described here.

! Every number read, is split up into its bits
! For unsigned integers, we use:
! 2^31 ......... 2^0 
! For single precision
! sign [power 2^7.... 2^0] [mantissa 1/2^23   1/2^1 ]
! For double 
! sign [power 2^11.... 2^0] [mantissa 1/2^52   1/2^1 ]

! after the split up, the bytes are reordered, but leaving the endianness of the individual bytes in tact
! single:
! 1 2 3 4 5 6 7 8            9 10 11 12 13 14 15 16   17 18 19 20 21 22 23 24   25 26 27 28 29 30 31 32
! to
! 25 26 27 28 29 30 31 32    17 18 19 20 21 22 23 34   9 10 11 12 13 14 15 16   1 2 3 4 5 6 7 8
! And for double precision a similar transformation, but with more bits.

! next we chop them up in blocks of six bits, 
! and compute the value of each block as with
! [2^5 ... 2^0]
! and we assign the character value of that position in the base64 encoding string
! ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/
! with A=0, B=1 ... +=62 and /=63

! Because the bits reprsenting the numbers, come in multiples of eight,
! and the base64 bytes in multiples of 6, there is overlap.
! A base64 character can consist of some bits of one number, and some bits of the next.
! In that case, we just keep on counting.
! Because of this it is possible that we have some bits left at the end,
! which are not a multiple of 6.
! We patch the string with an equals sign for every byte that is needed before 
! the bits reach a multiple of 8 again, according to the convention.
! This means we will only see 1 of 2 trailing equals signs, but in theory
! we could add any three-fold again to that.

! Finally, each encoding string gets a header, being an endoced 4-byte integer
! representing the size of the original array in 8-bit bytes.
! The 4 byte int consists of 32 bits, which is not a multiple of 6, so it is patched
! with equals signs just as described above. 

! A typical example would be

!  EAAAAA==AACAPwAAAED//z9AAACAQA==

! representing the four floating point numbers 1.0 2.0 3.0 4.0, encoded in single precision
! with the header '16' (4 x 4 bytes) -> 'EAAAAA=='

!-----------------------------------------------------------------------------
! The types of variable supported by the VTU format are:
! (from http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf page 15)
!
! Int8, Int16, Int32, Int64
! UInt8, UInt16, UInt32
! Float32, Float64
!
! . Note: the 64-bit integer types are only supported if
! VTK_USE_64BIT_IDS
! is on (a CMake variable?see ?CMake? on page 8) or the platform is 64-bit




implicit none


character(len=64), parameter :: Base64Scale = & 
   "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

integer :: precision
integer :: offset ! to initialize variables

!------------------------------------------------

!************************************************
public :: EncodeArrayToFile
!************************************************
private :: encodeDoubles1DarrayToBase64
private :: encodeDoubles2DarrayToBase64
private :: encodeDoubles3DarrayToBase64
private :: encodeUIntegers1DarrayToBase64
private :: encodeUIntegers2DarrayToBase64


!************************************************
public :: RunSelfTest
!************************************************

private :: NextEntry
private :: nextDoubleFrom1DArray
private :: nextDoubleFrom2DArray 
private :: nextDoubleFrom3DArray
private :: nextUIntegerFrom1DArray
private :: nextUIntegerFrom2DArray

private :: Int4ToBase64

private :: TestDoubleStreamer
private :: TestDoubleToBits
private :: TestIntToBase64
private :: TestRearrangeBitString
private :: TestFullStrings

private :: SplitIntoBits
private :: Double2Bits
private :: UInteger2Bits


interface EncodeArrayToFile
    ! main call to the module,
    ! depending on type and dimension of the argumented array
    ! a procedure is chosen.
    module procedure encodeDoubles1DarrayToBase64, &
                     encodeDoubles2DarrayToBase64, &
                     encodeDoubles3DarrayToBase64, &
                     encodeUIntegers1DarrayToBase64, &
                     encodeUIntegers2DarrayToBase64
end interface

interface NextEntry
    ! request the next double of the array.
    ! this can handle both 1D and 2D arrays, and it is easy to add higher order.
    ! It could be done easier with reshape, but that would require additional
    ! memory the size of the array to be written and for larger arrays
    ! that can become gigabytes, hence the stream approach.
    module procedure nextDoubleFrom1DArray, &
                     nextDoubleFrom2DArray, &
                     nextDoubleFrom3DArray, &
                     nextUIntegerFrom1DArray, &
                     nextUIntegerFrom2DArray
end interface

interface SplitIntoBits
    ! both procedures also depend on the module var 'precision',
    ! which was passed as 3rd argument of EncodeArrayToFile.
    module procedure Double2Bits, &
                     UInteger2Bits
end interface



contains

!----------------------------------------------------------------------------
!   First routine that stream the arrays
!----------------------------------------------------------------------------

subroutine nextDoubleFrom1DArray(doubles, currentDouble, firstCall)

    implicit none

    double precision :: doubles(:)
    integer          :: size1
    integer, save    :: currentPos1
    double precision :: currentDouble

    logical          :: firstCall

    if (firstCall) then
        currentPos1 = 0
        firstCall = .false.
    endif

    size1 = size(doubles, 1)

    if (currentPos1.eq.size1) then
        write(*,*) "Sorry, sahib. All out of doubles..."
        return
    else
        currentPos1 = currentPos1 + 1
        currentDouble = doubles(currentPos1)
    endif

    return
end subroutine

!-------------------------------------------------------------------------------------

subroutine nextDoubleFrom2DArray(doubles, currentDouble, firstCall)

    implicit none

    double precision :: doubles(:,:)
    integer          :: size1, size2
    integer, save    :: currentPos1, currentPos2
    double precision :: currentDouble

    logical          :: firstCall

    if (firstCall) then
        currentPos1 = 0
        currentPos2 = 1
        firstCall = .false.
    endif

    size1 = size(doubles, 1)
    size2 = size(doubles, 2)

    if ((currentPos1.eq.size1) .and. (currentPos2.eq.size2)) then
        write(*,*) "Sorry, sahib. All out of doubles..."
        return
    else
        if (currentPos1.eq.size1) then
            currentPos1 = 1
            currentPos2 = currentPos2 + 1
            currentDouble = doubles(currentPos1,currentPos2)
        else
            currentPos1 = currentPos1 + 1
            currentDouble = doubles(currentPos1,currentPos2)
        endif
    endif

    return
end subroutine

!-------------------------------------------------------------------------------------

subroutine nextDoubleFrom3DArray(doubles, currentDouble, firstCall)

    implicit none

    double precision :: doubles(:,:,:)
    integer          :: size1, size2, size3
    integer, save    :: currentPos1, currentPos2, currentPos3
    double precision :: currentDouble

    logical          :: firstCall

    if (firstCall) then
        currentPos1 = 0
        currentPos2 = 1
        currentPos3 = 1
        firstCall = .false.
    endif

    size1 = size(doubles, 1)
    size2 = size(doubles, 2)
    size3 = size(doubles, 3)

    if ((currentPos1.eq.size1) .and. \
        (currentPos2.eq.size2) .and. \
        (currentPos3.eq.size3)) then
        write(*,*) "Sorry, sahib. All out of doubles..."
        return
    else
        if (currentPos1.eq.size1) then
            if (currentPos2.eq.size2) then
                currentPos1 = 1
                currentPos2 = 1
                currentPos3 = currentPos3 + 1
                currentDouble = doubles(currentPos1,currentPos2,currentPos3)
            else
                currentPos1 = 1
                currentPos2 = currentPos2 + 1
                currentDouble = doubles(currentPos1,currentPos2,currentPos3)
            endif
        else
            currentPos1 = currentPos1 + 1
            currentDouble = doubles(currentPos1,currentPos2,currentPos3)
        endif
    endif

    return
end subroutine

!-------------------------------------------------------------------------------------


subroutine nextUIntegerFrom1DArray(integers, currentInteger, firstCall)

    implicit none

    integer          :: integers(:)
    integer          :: size1
    integer, save    :: currentPos1
    integer          :: currentInteger

    logical          :: firstCall

    if (firstCall) then
        currentPos1 = 0
        firstCall = .false.
    endif

    size1 = size(integers, 1)

    if (currentPos1.eq.size1) then
        write(*,*) "Sorry, sahib. All out of integers..."
        return
    else
        currentPos1 = currentPos1 + 1
        currentInteger = integers(currentPos1)
    endif

    return
end subroutine

!-------------------------------------------------------------------------------------


subroutine nextUIntegerFrom2DArray(integers, currentInteger, firstCall)

    implicit none

    integer          :: integers(:,:)
    integer          :: size1, size2
    integer, save    :: currentPos1, currentPos2
    integer          :: currentInteger

    logical          :: firstCall

    if (firstCall) then
        currentPos1 = 0
        currentPos2 = 1
        firstCall = .false.
    endif

    size1 = size(integers, 1)
    size2 = size(integers, 2)


    if ((currentPos1.eq.size1) .and. (currentPos2.eq.size2)) then
        write(*,*) "Sorry, sahib. All out of integers..."
        return
    else
        if (currentPos1.eq.size1) then
            currentPos1 = 1
            currentPos2 = currentPos2 + 1
            currentInteger = integers(currentPos1,currentPos2)
        else
            currentPos1 = currentPos1 + 1
            currentInteger = integers(currentPos1,currentPos2)
        endif
    endif

    return
end subroutine

!---------------------------------------------------------------------------












!-----------------------------------------------------------------

subroutine encodeDoubles1DarrayToBase64(fp, doubles, precisionFromArg)

    implicit none

    integer          :: fp
    double precision :: doubles(:)
    integer          :: size1
!    integer          :: precision

    integer          :: i
    double precision :: currentDouble
    integer          :: currentPosition1
    integer          :: precisionFromArg

    logical          :: isFirstEntry

    size1 = size(doubles,1)

    !import precision module-wide
    precision = precisionFromArg

    ! write the line header containing the size
    write(fp,'(a8)',advance='no') Int4ToBase64(4 * precision * size1)

    currentPosition1 = 0

    offset = 0

    isFirstEntry = .true.

    do i=1,size1    
        call nextEntry(doubles, currentDouble, isFirstEntry) 
        isFirstEntry = .false.

        if (size1 .eq. i) then
            ! last double; special treatment because of patching
            call encodeBitStringToBase64(fp, SplitIntoBits(currentDouble), .true.)
        else
            call encodeBitStringToBase64(fp, SplitIntoBits(currentDouble), .false.)
        endif

    enddo

    

end subroutine


!-------------------------------------------------------------------------------------

subroutine encodeDoubles2DarrayToBase64(fp, doubles, precisionFromArg)

    implicit none

    integer          :: fp
    double precision :: doubles(:,:)
    integer          :: size1, size2

    integer          :: i,j
    double precision :: currentDouble
    integer          :: precisionFromArg

    logical          :: isFirstEntry

    size1 = size(doubles,1)
    size2 = size(doubles,2)


    !import precision module-wide
    precision = precisionFromArg

    ! write the line header containing the size
    write(fp,'(a8)',advance='no') Int4ToBase64(4 * precision * size1 * size2)
    
    offset = 0

    isFirstEntry = .true.
    do i=1,size1
        do j=1,size2
            call nextEntry(doubles, currentDouble, isFirstEntry)
            isFirstEntry = .false.

            if ((size1 * size2) .eq. (i * j)) then
                ! last double; special treatment because of patching
                call encodeBitStringToBase64(fp, SplitIntoBits(currentDouble), .true.)
            else
                call encodeBitStringToBase64(fp, SplitIntoBits(currentDouble), .false.)
            endif
        enddo
    enddo

end subroutine

!---------------------------------------------------------------------------

subroutine encodeDoubles3DarrayToBase64(fp, doubles, precisionFromArg)

    implicit none

    integer          :: fp
    double precision :: doubles(:,:,:)
    integer          :: size1, size2, size3

    integer          :: i,j,k
    double precision :: currentDouble
    integer          :: precisionFromArg

    logical          :: isFirstEntry

    size1 = size(doubles,1)
    size2 = size(doubles,2)
    size3 = size(doubles,3)


    !import precision module-wide
    precision = precisionFromArg

    ! write the line header containing the size
    write(fp,'(a8)',advance='no') Int4ToBase64(4 * precision * size1 * size2 * size3)

    offset = 0

    isFirstEntry = .true.
    do i=1,size1
        do j=1,size2
            do k=1,size3

                call nextEntry(doubles, currentDouble, isFirstEntry)
                isFirstEntry = .false.

                if ((size1 * size2 * size3) .eq. (i * j * k)) then
                    ! last double; special treatment because of patching
                    call encodeBitStringToBase64(fp, SplitIntoBits(currentDouble), .true.)
                else
                    call encodeBitStringToBase64(fp, SplitIntoBits(currentDouble), .false.)
                endif
            enddo
        enddo
    enddo

end subroutine

!---------------------------------------------------------------------------


subroutine encodeUIntegers1DarrayToBase64(fp, integers, precisionFromArg)

    implicit none

    integer          :: fp
    integer          :: integers(:)
    integer          :: size1
!    integer          :: precision

    integer          :: i
    integer          :: currentInteger
    integer          :: currentPosition1
    integer          :: precisionFromArg

    logical          :: isFirstEntry

    currentInteger = 0
    size1 = size(integers,1)

    !import precision module-wide
    precision = precisionFromArg

    ! write the line header containing the size
    write(fp,'(a8)',advance='no') Int4ToBase64(4 * precision * size1)

    currentPosition1 = 0

    offset = 0

    isFirstEntry = .true.

    do i=1,size1
        call nextEntry(integers, currentInteger, isFirstEntry)
        isFirstEntry = .false.

        if (size1 .eq. i) then
            ! last double; special treatment because of patching
            call encodeBitStringToBase64(fp, SplitIntoBits(currentInteger), .true.)
        else
            call encodeBitStringToBase64(fp, SplitIntoBits(currentInteger), .false.)
        endif

    enddo



end subroutine


!-------------------------------------------------------------------------------------

subroutine encodeUIntegers2DarrayToBase64(fp, integers, precisionFromArg)

    implicit none

    integer          :: fp
    integer          :: integers(:,:)
    integer          :: size1, size2

    integer          :: i,j
    integer          :: currentInteger
    integer          :: precisionFromArg

    logical          :: isFirstEntry

    size1 = size(integers,1)
    size2 = size(integers,2)

    !import precision module-wide
    precision = precisionFromArg

    ! write the line header containing the size
    write(fp,'(a8)',advance='no') Int4ToBase64(4 * precision * size1 * size2)

    offset = 0

    isFirstEntry = .true.
    do i=1,size1
        do j=1,size2
            call nextEntry(integers, currentInteger, isFirstEntry)
            isFirstEntry = .false.

            if ((size1 * size2) .eq. (i * j)) then
                ! last double; special treatment because of patching
                call encodeBitStringToBase64(fp, SplitIntoBits(currentInteger), .true.)
            else
                call encodeBitStringToBase64(fp, SplitIntoBits(currentInteger), .false.)
            endif
        enddo
    enddo

end subroutine

!---------------------------------------------------------------------------


subroutine encodeBitStringToBase64(fp, newBitString, isLastEntry)

implicit none

integer          :: fp
logical          :: isLastEntry
integer          :: newBitString(64)

integer, save    :: bits(64)
integer, save    :: base64value
integer, save    :: bitsLeftFromPreviousBitString
logical, save    :: isSplitNumber
character(len=1) :: Base64Bit

integer          :: i

! 'offset' taken from module variables.
Base64Value = 0

if ((offset.eq.32*precision) .or. (offset.eq.0)) then
    ! nothing left from previous byte.
    ! chop the double mercilessly to bits
    offset = 0
    isSplitNumber = .false.
else
  
    ! some bits still need to be added to an overlapping Base64 character
    do i=1,bitsLeftFromPreviousBitString
        base64value = base64value + 2**(6 - i) * bits(32 * precision - bitsLeftFromPreviousBitString + i)
    enddo
    isSplitNumber = .true.

endif

offset = 0

! Chop the new double mercilessly to bits
bits = newBitString
call rearrangeBitString(bits)

! when necessary, finish the base64 char we started building with the previous bitstring
if (isSplitNumber) then
    do i=1,6-bitsLeftFromPreviousBitString
         base64value = base64value + 2**(6 - bitsLeftFromPreviousBitString-i) * bits(i)
    enddo

    ! write the recombined base64 chracter to file
    Base64Bit = Base64Scale(Base64Value+1:Base64Value+1)
    write(fp,'(a1)',advance='no') Base64Bit

    offset = 6 - bitsLeftFromPreviousBitString
endif

do while (offset.le.32*precision-6)
    ! we have enough data left to create a new encode character

    Base64Value = 0
    do i=1,6
        Base64Value = Base64Value + 2**(6 - i) * bits(offset +i)
    enddo

    Base64Bit = Base64Scale(Base64Value+1:Base64Value+1)
    write(fp,'(a1)',advance='no') Base64Bit

    offset = offset + 6

enddo



if (isLastEntry) then
    if ((offset.eq.32*precision - 2) .or. (offset.eq.32*precision - 4)) then
        base64value = 0
        ! we have some remaining bits of the last double.
        ! Encode it, and write it
        do i=1,32*precision - offset
            base64value = base64value + 2**(6 - i) * bits(offset + i)
        enddo

        ! work towards a patch
        if (offset.eq.32*precision - 2) then
            write(fp,'(a1)',advance="no") Base64Scale(Base64Value+1:Base64Value+1)
            write(fp,'(a2)') '=='
        else if (offset.eq.32*precision - 4) then
            write(fp,'(a1)',advance="no") Base64Scale(Base64Value+1:Base64Value+1)
            write(fp,'(a1)') '='
        else
            ! if offset 0, we need no patching
            write(fp,'(a1)') Base64Scale(Base64Value+1:Base64Value+1)
        endif

    else
        write(fp,'(a1)') " "
        ! write the space for the line break
        ! We are done!
    endif


else
    ! continue until out of double
    bitsLeftFromPreviousBitString = 32*precision-offset
endif



end subroutine


!--------------------------------------------------------------------------------

character(len=8) function Int4ToBase64(number)
! sets a 4 byte integer to a base 64 character string  

implicit none

integer          :: number, i,j
integer          :: thisDigit
character(len=8) :: ConstructionString
integer          :: localNumber

integer          :: bits(64)
integer          :: tempPrecision

! split the integer to bits
! this tag is always in single precision, even when the encoding is in double,
! hence we preserve the precision.


tempPrecision = precision
precision = 1
bits = SplitIntoBits(number)
precision = tempPrecision

!swap the bytes
call rearrangeBitString(bits)

! and encode to base 64

! Initialise to 0
ConstructionString = "AAAAAAAA"

! first digit outside loop, because of string concatenation not able to concatenate a zero length string
thisDigit = 0
do i = 1, 6
     thisDigit = thisDigit + 2**(6-i) * bits(i)
enddo

ConstructionString = base64Scale(thisDigit+1:thisDigit+1) // ConstructionString(2:8)

do i=2,5
     thisDigit = 0
     do j = 1, 6
         thisDigit = thisDigit + 2**(6-j) * bits(j+6*(i-1))
     enddo
     ConstructionString = ConstructionString(1:i-1) // base64Scale(thisDigit+1:thisDigit+1) // ConstructionString(i+1:8)
enddo

! do only the first four bits of the last digit. Then we are out of bits.
thisDigit = 0
do j = 1, 2
    thisDigit = thisDigit + 2**(6-j) * bits(j+30)
enddo
ConstructionString = ConstructionString(1:i-1) // base64Scale(thisDigit+1:thisDigit+1) // ConstructionString(i+1:8)

! padding (because of the fixed size of input and output, 
! we know the padding by hard.
Int4ToBase64 = ConstructionString(1:6) // '=='



return
end function

!-------------------------------------------------------------------------------------
function Double2bits(number)
! split one double precision number into its bits
! return answer is always 64 wide, for double precision
! If single precision is requested, only the first 32 are used.

implicit none

integer          :: Double2bits(64)
double precision :: number
double precision :: localNumber

integer          :: nExponentBits
integer          :: nSignificantBits
integer          :: baseExponent

integer          :: biggestPowerOfTwo ! biggest power of two that fits withing the double
integer          :: power
integer          :: i

if (precision.eq.1) then
    nExponentBits    = 8
    nSignificantBits = 23
else if (precision.eq.2) then
    nExponentBits    = 11
    nSignificantBits = 52
else
    write(*,*) "Error: precision should be either 1 (single) or 2 (double); assuming double precision"
    nExponentBits    = 11
    nSignificantBits = 52
endif

baseExponent = 2**(nExponentBits-1) - 1

! cannot edit argument of a function, so copy to local variable
localNumber = number

! Initiliase
Double2bits = 0

! determine minus sign bit
if (localNumber.lt.0) then
    Double2bits(1) = 1
    localNumber = -localNumber
endif

if (localNumber.eq.0.0d0) then
! special case, cannot determine a power.
    Double2bits = 0
    return
endif

! determine exponent, and scale it to withing the [ 1,2 > interval
biggestPowerOfTwo = floor(log(localNumber) / log(2.0d0))

localNumber = localNumber / (dble(2)**dble(biggestPowerOfTwo))

power = baseExponent + biggestPowerOfTwo

! ... and convert the power to bits
if ((power.gt.2**(nExponentBits) - 1) .or. (power.lt.0)) then
    ! overflow error. Double too large or small
    write(*,*) 'power: ', power
    write(*,*) 'baseExponent ', baseExponent
    write(*,*) 'precision overflow, trying to convert number ', number, 'to binary. Did not fit'
endif

i = nExponentBits
do while (power.gt.0)
    if (power.ge.2**i) then
        Double2bits(nExponentBits - i + 1) = 1

        power = power - 2**i
    endif
    i = i - 1
enddo

! and do the significant digits

! bit notation = 1 + sum over i   bit * 2^-i
! hence we remove the 1. 
localNumber = localNumber - 1

do i=1,nSignificantBits
    if (localNumber.gt.0.5) then
        
        Double2bits(1 + nExponentBits + i) = 1
        localNumber = localNumber - 0.5d0
    endif

    localNumber = localNumber * 2.0d0
enddo

return



end function


!-------------------------------------------------------------------------------------

function UInteger2bits(number)
! split one double precision number into its bits
! return answer is always 64 wide, for double precision
! If single precision is requested, only the first 32 are used.

implicit none

integer          :: UInteger2bits(64)
integer          :: number, localNumber
integer          :: i

localNumber = number

UInteger2bits = 0
! fill the bits
i = 1
do while (localNumber.gt.0)
    if (mod(localNumber,2).eq.1) then
        UInteger2bits(32*precision+1-i) = 1
        localNumber = localNumber - 1
    endif
    localNumber = localNumber / 2
    i = i + 1
enddo

end function


!----------------------------------------------------------

subroutine rearrangeBitString(bits)

implicit none

integer :: bits(64), temp
integer :: i

! reverse the order of the bytes, but leave the endianness in tact
! single:
! 1 2 3 4 5 6 7 8            9 10 11 12 13 14 15 16   17 18 19 20 21 22 23 24   25 26 27 28 29 30 31 32
! to
! 25 26 27 28 29 30 31 32    17 18 19 20 21 22 23 34   9 10 11 12 13 14 15 16   1 2 3 4 5 6 7 8
! And for double precision a similar transformation, but with more bits.

if (precision.eq.1) then
    call SwapBytesDouble(bits,1,4)
    call SwapBytesDouble(bits,2,3)
else if (precision.eq.2) then
    call SwapBytesDouble(bits,1,8)
    call SwapBytesDouble(bits,2,7)
    call SwapBytesDouble(bits,3,6)
    call SwapBytesDouble(bits,4,5)
endif

end subroutine

!----------------------------------------------------------


subroutine SwapBytesDouble(bits, b1, b2)
    implicit none
    integer :: b1, b2, bits(64)
    integer :: temp, i

    do i=1,8
        temp = bits(8*(b1-1) + i)
        bits(8*(b1-1) + i) = bits(8*(b2-1) + i)
        bits(8*(b2-1) + i) = temp
    enddo

end subroutine

subroutine SwapBytesSingle(bits, b1, b2)
    implicit none
    integer :: b1, b2, bits(32)
    integer :: temp, i

    do i=1,8
        temp = bits(8*(b1-1) + i)
        bits(8*(b1-1) + i) = bits(8*(b2-1) + i)
        bits(8*(b2-1) + i) = temp
    enddo

end subroutine


!##########################################################################################
!        below here are the test routines
!##########################################################################################

subroutine RunSelfTest()
    implicit none
    call TestDoubleStreamer()
    call TestDoubleToBits()
    call TestRearrangeBitString()
    call TestIntToBase64()
    call TestRearrangeBitString()
    call TestFullStrings() 
    stop
end subroutine

!-------------------------------------------------------------------------------------

subroutine TestDoubleStreamer()
    ! the function NextEntry should scroll through the array giving all of the doubles
    implicit none
  
    logical          :: OK

    integer          :: i, j
    double precision :: eps = 1e-10
    double precision :: TestArray1D(5) = &
                         (/6.0d0, 8.0d0, 123.0d0, -4.0d0, 3.1415926d0/)

    double precision :: TestArray2D(3,4)

    double precision :: currentDouble
   
    logical          :: firstCall

! test with 1D array
    OK = .true.

    firstCall = .true.
    do i=1,5
        call nextEntry(TestArray1D, currentDouble, firstCall)
        firstCall = .false.
        if      (i.eq.1) then
            if (currentDouble - 6.0d0 .gt.eps)     OK = .false.
        else if (i.eq.2) then
            if (currentDouble - 8.0d0 .gt.eps)     OK = .false.
        else if (i.eq.3) then
            if (currentDouble - 123.0d0 .gt.eps)   OK = .false.
        else if (i.eq.4) then
            if (currentDouble + 4.0d0 .gt.eps)     OK = .false.
        else if (i.eq.5) then    
            if (currentDouble - 3.1415926d0 .gt.eps) OK = .false.
        endif
    enddo

    if (OK .eqv. .false.) then
        write(*,*) 'Panic! BitStreamer in 1D is broken!'
    else
        write(*,*) 'BitStreamer in 1D is in tiptop condition'
    endif

! test with 2D array
    OK = .true.

! to make sure we have the correct index sequence, 
! we fill the test array by hand, in stead of letting 'reshape' do it.
    TestArray2D(1,1) = 6.0d0
    TestArray2D(2,1) = 8.0d0       
    TestArray2D(3,1) = 123.0d0
    TestArray2D(1,2) = -4.0d0
    TestArray2D(2,2) = 3.1415926d0
    TestArray2D(3,2) = 16.0d0
    TestArray2D(1,3) = 18.0d0
    TestArray2D(2,3) = 1234.0d0
    TestArray2D(3,3) = -41.0d0
    TestArray2D(1,4) = 38.1415926d0
    TestArray2D(2,4) = 88.0d0
    TestArray2D(3,4) = -123.0d0

    firstcall = .true.

    do i=1,3
        do j=1,4
            call nextEntry(TestArray2D, currentDouble, firstcall)
            firstcall = .false.
            if (i.eq.1) then
                if      (j.eq.1) then
                   if (currentDouble - 6.0d0 .gt.eps)         OK = .false.
                else if (j.eq.2) then
                   if (currentDouble - 8.0d0 .gt.eps)         OK = .false.
                else if (j.eq.3) then
                   if (currentDouble - 123.0d0 .gt.eps)       OK = .false.
                else if (j.eq.4) then
                   if (currentDouble + 4.0d0 .gt.eps)         OK = .false.
                endif
            else if (i.eq.2) then
                if      (j.eq.1) then
                   if (currentDouble - 3.1415926d0 .gt.eps)   OK = .false.
                else if (j.eq.2) then
                   if (currentDouble - 16.0d0 .gt.eps)        OK = .false.
                else if (j.eq.3) then
                   if (currentDouble - 18.0d0 .gt.eps)        OK = .false.
                else if (j.eq.4) then
                   if (currentDouble - 1234.0d0 .gt.eps)      OK = .false.
                endif
            else if (i.eq.3) then
                if      (j.eq.1) then
                   if (currentDouble + 41.0d0 .gt.eps)        OK = .false.
                else if (j.eq.2) then
                   if (currentDouble - 38.1415926d0 .gt.eps)  OK = .false.
                else if (j.eq.3) then
                   if (currentDouble - 88.0d0 .gt.eps)        OK = .false.
                else if (j.eq.4) then
                   if (currentDouble + 123.0d0 .gt.eps)       OK = .false.
                endif
            endif

        enddo
    enddo

    if (OK .eqv. .false.) then
        write(*,*) 'Panic! BitStreamer in 2D is broken!'
    else
        write(*,*) 'BitStreamer in 2D is in tiptop condition'
    endif

    return
end subroutine

!-------------------------------------------------------------------------------------

subroutine TestDoubleToBits()

    implicit none

    logical :: OK

    double precision, parameter :: testnum1 = 12345.6789d0
    double precision, parameter :: testnum2 = -0.00000031415926d0

    integer :: i
    integer :: bits(64)

    ! single bits: 1 minus sign bit; 8 power bits; 23 mantissa bits
    ! double bits: 1 minus sign bit; 11 power bits; 52 mantissa bits

    ! test value 1: 12345.6789
    ! sign        : plus             = 0
    ! power       : 
    !   - single  : 127 + 13  =  140 = 1 0 0 0 1 1 0 0
    !   - double  : 1023 + 13 = 1036 = 1 0 0 0 0 0 0 1 1 0 0
    ! mantissa    : 12345.6789 / 2^13 - 1 = 0.507040881347656
    !               
    integer :: referencebits1s(64) = &
      (/0, & ! sign
        1,0,0,0,1,1,0,0, & ! power
        1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,1,1,0,1,1,1, & ! mantissa
        ! + 32 unused bits in single precision
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
    integer :: referencebits1d(64) = &
      (/0,& ! sign
        1,0,0,0,0,0,0,1,1,0,0,& ! power
        1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,1,1,0,1,1,1,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,1,0,0,0,1,0,1,0,0,0,0,1 /) ! mantissa

    ! test value 2: -0.00000031415926 
    ! sign        : minus            = 1
    ! power       : 
    !   - single  : 127 -  21 =  106 = 0 1 1 0 1 0 0 1
    !   - double  : 1023 - 21 = 1002 = 0 1 1 1 1 1 0 1 0 0 1
    ! mantissa    : 0.00000031415926 / 2^-22 - 1 = 0.31767944085504
    !

    integer :: referencebits2s(64) = &
      (/1, & ! sign
        0,1,1,0,1,0,0,1, & ! power
        0,1,0,1,0,0,0,1,0,1,0,1,0,0,1,1,0,1,1,1,0,0,0, & ! mantissa
        ! + 32 unused bits in single precision
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, & 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
    integer :: referencebits2d(64) = &
      (/1,& ! sign
        0,1,1,1,1,1,0,1,0,0,1,& ! power
        0,1,0,1,0,0,0,1,0,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,0,1,0,1,0 /) ! mantissa


    OK = .true.

! single precision
! only the first 32 bits of the string are checked.
    precision = 1
    bits = SplitIntoBits(testnum1)

    do i=1,30 ! omit last two bits to prevent test failure due to rounding errors
        if (bits(i) .ne. referencebits1s(i)) then
            OK = .false.
            write(*,*) "Single precision does not convert correctly to digits", i
        endif
    enddo
    
    bits = SplitIntoBits(testnum2)

    do i=1,30 ! omit last two bits to prevent test failure due to rounding errors
        if (bits(i) .ne. referencebits2s(i)) then
            OK = .false.
            write(*,*) "Single precision does not convert correctly to digits", i
        endif
    enddo


! double precision    

    precision = 2
    bits = SplitIntoBits(testnum1)

    do i=1,62 ! omit last two bits to prevent test failure due to rounding errors
        if (bits(i) .ne. referencebits1d(i)) then
            OK = .false.
            write(*,*) "Double precision does not convert correctly to digits", i
        endif
    enddo 

    bits = SplitIntoBits(testnum2)

    do i=1,62 ! omit last two bits to prevent test failure due to rounding errors
        if (bits(i) .ne. referencebits2d(i)) then
            OK = .false.
            write(*,*) "Double precision does not convert correctly to digits", i
        endif
    enddo 

    if (OK) then
        write(*,*) "Double 2 bits convertor fully operational"
    else
        write(*,*) "Please tighten the loose screws on the double 2 bits convertor."
    endif


end subroutine

!-------------------------------------------------------------------------------------

subroutine TestIntToBase64()
    ! the integer is 4 bytes by definition of the VTK standard.
    implicit none

    logical :: OK

    integer :: testint1 = 0
    integer :: testint2 = 100000
    integer :: testint3 = 87654321

    character(len=8) :: base64String
     
    OK = .true.

    base64String = Int4ToBase64(testint1)

    if (base64String .ne. "AAAAAA==") then
!       A = 0
        OK = .false.
    endif

    base64String = Int4ToBase64(testint2)

    if (base64String .ne. "oIYBAA==") then
        ! 100000 = 2^16 + 2^15 + 2^10 + 2^9 + 2^7 + 2^5
        ! in bits:  0000 0000   0000 0001   1000 0110   1010 0000
        ! swaps to: 1010 0000   1000 0110   1000 0001   0000 0000
        ! in hexes: 101000  001000 011010 000001 000000 00~~~~
        ! adds to:  40       8      24       1      0     0
        ! in base64:  o      I      Y        B      A     A  == 
        OK = .false.
    endif

    base64String = Int4ToBase64(testint3)

    if (base64String .ne. "sX85BQ==") then
        ! 87654321 = 2^26 + 2^24 + 2^21 + 2^20 + 2^19 + 2^16 + 2^14 + 2^13 + 2^12 + 2^11 + 2^10 + 2^9 + 2^8 + 2^7 + 2^5 + 2^4 + 2^0
        ! in bits:  0000 0101   0011 1001   0111 1111   1011 0001
        ! swaps to: 1011 0001   0111 1111   0011 1001   0000 0101
        ! in hexes: 101100  010111  111100  111001  000001  01~~~~
        ! adds to:  44       23      60       57      1     16    
        ! in base64: s       X       8        5       B     Q ==
        OK = .false.
    endif


    if (OK) then
        write(*,*) 'Int2bits convertor in perfect working condition'
    else
        write(*,*) 'Int2bits convertor is broken'
    endif

end subroutine

!-------------------------------------------------------------------------------------

subroutine TestRearrangeBitString()
! due the integer nature of the values, we can use non-binary integers and track them thoroughly
implicit none

integer :: bits(64)
integer :: refbitsSingle(64), refbitsDouble(64)
logical :: OK
integer :: i

bits = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64/)

! single; only first four bytes swapped
refbitsSingle=(/25,26,27,28,29,30,31,32,17,18,19,20,21,22,23,24,9,10,11,12,13,14,15,16,1,2,3,4,5,6,7,8,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64/)

! double: all bytes swapped
refbitsDouble=(/57,58,59,60,61,62,63,64,49,50,51,52,53,54,55,56,41,42,43,44,45,46,47,48,33,34,35,36,37,38,39,40,25,26,27,28,29,30,31,32,17,18,19,20,21,22,23,24,9,10,11,12,13,14,15,16,1,2,3,4,5,6,7,8/)

OK = .true.

precision = 1
call RearrangeBitString(bits)

do i=1,64
    if (bits(i).ne.refbitsSingle(i)) then
        OK = .false.
    endif
enddo

!and return them to their original place
call RearrangeBitString(bits)

precision = 2

call RearrangeBitString(bits)

do i=1,64
    if (bits(i).ne.refbitsDouble(i)) then
        OK = .false.
    endif
enddo

if (OK) then
    write(*,*) 'RearrangeBitString = OK'
else
    write(*,*) 'RearrangeBitString = FAIL'
endif

end subroutine

!-------------------------------------------------------------------------------------

subroutine TestFullStrings()
! test taken from a VTU file 'vox8', known to be working correctly
! The file exists both in binary and in ASCII (and in several other varieties).
! A choice selection has been made from its data.
! These files are created by Prof. Hans-Peter Bunge
! and is available from:
! http://www.earthmodels.org/software/vtk-and-paraview/vtk-file-formats

! due to rounding differences, the test result and reference encoding can differ significantly,
! while the numbers they represent are nearly identical.
! To test this for real, we would actually need a base 64 decoder as well.
! But that will come in the future

implicit none




double precision, dimension(27)   :: doubles1 = &
       (/ -0.4, -0.4, -0.4, -0.4, -1.0, -0.4, -0.4, -0.4, -0.4,  &
           0.4,  1.0,  0.4,  1.0,  0.4,  1.0,  0.4,  1.0,  0.4, &
           0.4,  0.4,  0.4,  0.4,  1.0,  0.4,  0.4,  0.4,  0.4 /)

                                ! strings are chopped off into pieces of a hundred characters
character(len=152) :: string1 = "bAAAAA==zczMvs3MzL7NzMy+zczMvgAAgL/NzMy+zczMvs3MzL7NzMy+zczMPgAAgD/NzMw+AACAP83MzD4AAIA/zczMPgAAgD/N" // &
                                "zMw+zczMPs3MzD7NzMw+zczMPgAAgD/NzMw+zczMPs3MzD7NzMw+"



double precision, dimension(27)   :: doubles2 = &
       (/ -13.0, -12.0, -11.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, &
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0 /)

character(len=152) :: string2 = "bAAAAA==AABQwQAAQMEAADDBAAAgwQAAEMEAAADBAADgwAAAwMAAAKDAAACAwAAAQMAAAADAAACAvwAAAAAAAIA/AAAAQAAAQEAA" // &
                                "AIBAAACgQAAAwEAAAOBAAAAAQQAAEEEAACBBAAAwQQAAQEEAAFBB"


double precision, dimension(27)   :: doubles3 = &
       (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
          0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, &
          0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 /)

character(len=152) :: string3 = "bAAAAA==AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAzczMPs3MzD7NzMw+zczMPs3MzD7NzMw+zczMPs3MzD7N" // &
                                "zMw+zcxMP83MTD/NzEw/zcxMP83MTD/NzEw/zcxMP83MTD/NzEw/"

double precision, dimension(81)   :: doubles4 = &
       (/ -1, -1, -1, &
          -1, -1,  0, &
          -1, -1,  1, &
          -1,  0, -1, &
          -1,  0,  0, &
          -1,  0,  1, &
          -1,  1, -1, &
          -1,  1,  0, &
          -1,  1,  1, &
           0, -1, -1, &
           0, -1,  0, &
           0, -1,  1, &
           0,  0, -1, &
           0,  0,  0, &
           0,  0,  1, &
           0,  1, -1, &
           0,  1,  0, &
           0,  1,  1, &
           1, -1, -1, &
           1, -1,  0, &
           1, -1,  1, &
           1,  0, -1, &
           1,  0,  0, &
           1,  0,  1, &
           1,  1, -1, &
           1,  1,  0, &
           1,  1,  1 /)

character(len=440) :: string4 = "RAEAAA==AACAvwAAgL8AAIC/AACAvwAAgL8AAAAAAACAvwAAgL8AAIA/AACAvwAAAAAAAIC/AACAvwAAAAAAAAAAAACAvwAAAAAA" // &
                                "AIA/AACAvwAAgD8AAIC/AACAvwAAgD8AAAAAAACAvwAAgD8AAIA/AAAAAAAAgL8AAIC/AAAAAAAAgL8AAAAAAAAAAAAAgL8AAIA/" // &
                                "AAAAAAAAAAAAAIC/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAIA/AAAAAAAAgD8AAIC/AAAAAAAAgD8AAAAAAAAAAAAAgD8AAIA/AACA" // &
                                "PwAAgL8AAIC/AACAPwAAgL8AAAAAAACAPwAAgL8AAIA/AACAPwAAAAAAAIC/AACAPwAAAAAAAAAAAACAPwAAAAAAAIA/AACAPwAA" // &
                                "gD8AAIC/AACAPwAAgD8AAAAAAACAPwAAgD8AAIA/"




double precision, dimension(64)   :: doubles5 = &
     (/    0,  9,  3, 12,  1, 10,  4, 13, &
           1, 10,  4, 13,  2, 11,  5, 14, &
           3, 12,  6, 15,  4, 13,  7, 16, &
           4, 13,  7, 16,  5, 14,  8, 17, &
           9, 18, 12, 21, 10, 19, 13, 22, &
          10, 19, 13, 22, 11, 20, 14, 23, &
          12, 21, 15, 24, 13, 22, 16, 25, &
          13, 22, 16, 25, 14, 23, 17, 26 /)


character(len=352) :: string5 = "AAEAAA==AAAAAAkAAAADAAAADAAAAAEAAAAKAAAABAAAAA0AAAABAAAACgAAAAQAAAANAAAAAgAAAAsAAAAFAAAADgAAAAMAAAAM" // &
                                "AAAABgAAAA8AAAAEAAAADQAAAAcAAAAQAAAABAAAAA0AAAAHAAAAEAAAAAUAAAAOAAAACAAAABEAAAAJAAAAEgAAAAwAAAAVAAAA" // &
                                "CgAAABMAAAANAAAAFgAAAAoAAAATAAAADQAAABYAAAALAAAAFAAAAA4AAAAXAAAADAAAABUAAAAPAAAAGAAAAA0AAAAWAAAAEAAA" // &
                                "ABkAAAANAAAAFgAAABAAAAAZAAAADgAAABcAAAARAAAAGgAAAA=="


double precision, dimension(4)   :: doublesPyra = (/1.0d0, -2.0d0, 3.0d0 , 4.0d0/)



! run test 1    fp, doubles, size1, precisionFromArg
!open(unit=61, file = "test1.dat")
!call EncodeDoublesToFile(61, doubles1, 27, 1)
!close(61)

open(unit=62, file = "test2.dat")
call EncodeArrayToFile(62, doublesPyra, 1)
close(62)


open(unit=63, file = "teststring1.base64")
call EncodeArrayToFile(63, doubles1, 1)
close(63)

open(unit=64, file = "teststring2.base64")
call EncodeArrayToFile(64, doubles2, 1)
close(64)

open(unit=65, file = "teststring3.base64")
call EncodeArrayToFile(65, doubles3, 1)
close(65)

open(unit=66, file = "teststring4.base64")
call EncodeArrayToFile(66, doubles4, 1)
close(66)

open(unit=67, file = "teststring5.base64")
call EncodeArrayToFile(67, doubles5, 1)
close(67)




end subroutine

end module
