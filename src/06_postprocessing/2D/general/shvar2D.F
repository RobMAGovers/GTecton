subroutine shvar (input)

! substitutes shell variables contents

implicit none
!-pass
character              :: input*(*)
!-locl
#include "lenchar.i"
character(len=LENCHAR) :: variable,output
integer, external      :: fnblnk,lnblk,indx
integer                :: l,n,m,n1,n2


l = LEN(input)

if (l.gt.0) then
    l = fnblnk(input)
    100        m = lnblk(input)
    if (l.le.m) then
        if (((indx(input,'~')).eq.l) .and. &
             (l.eq.fnblnk(input))) then
!           string starts with a tilde
            if (indx(input,'/').eq.l+1) then
!               substitute HOME variable
                call get_environment_variable('HOME',variable)
                output = variable(1:lnblk(variable)) // &
                 input(l+1:m)
                l = l + 2
                input = output
                goto 100
            else
!               substitute USER home directory
                n = indx(input(l+1:m),'/')

                if (n.le.0) then
                    n=m
                endif

                call get_environment_variable('HOME',variable)

                if (n.lt.m) then
                    output = variable(1:lnblk(variable)) // &
                             input(n+1:m)
                else
                    output = variable
                endif

                l = n+1
                input = output
                goto 100
            endif

        else if (indx(input,'$').gt.0) then
            n = indx(input,'$')
            if (n.lt.m) then
                call varfit (input(n+1:m),variable,n1,n2)
                if (lnblk(variable).gt.0) then
                    n2 = n2 + n
                    if (n2.lt.m) then
                        output = variable(1:lnblk(variable)) // &
                                 input(n2+1:m)
                    else
                        output = variable
                    endif
                    input = output
                endif
            endif
        endif
    endif
endif
return
end subroutine
