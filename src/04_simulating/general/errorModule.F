module errors

implicit none

! flags for error handling
integer,          parameter :: OK                   = 0
integer,          parameter :: EOF                  = -1
integer,          parameter :: EXEC                 = 1
integer,          parameter :: SYNTAX               = 2
integer,          parameter :: PRGERR               = 3
integer,          parameter :: STOP                 = 4

end module
