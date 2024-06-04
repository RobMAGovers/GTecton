integer function IPEVEN(IP, wordSize)

!
! Promotes doubles to even (actually odd) word boundaries
! A(1) is even, A(3) is even and so on.
!
implicit none

integer :: IP, wordSize

if (wordSize .ne. 1 .and. wordSize .ne. 2) then
    write(0,*) 'IPEVEN says: unknown precision wordSize.'
    write(0,*) 'Must be 1 of 2 but it is: ', wordSize
    write(0,*) 'Please contact model support'
    ipeven = -1
    return
endif

IPEVEN = 1 + INT((IP-1)/wordSize)*wordSize

if (MOD(IP-1,wordSize).ne.0) then
    IPEVEN = IPEVEN + wordSize
endif

return
end

