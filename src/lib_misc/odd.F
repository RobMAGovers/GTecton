logical function odd(i)

implicit none

integer :: i

odd = (MOD(i,2).ne.0)

return

end function
