subroutine varfit (string,value,k,l)

implicit none

character         :: string*(*),value*(*)
integer, external :: lnblk, fnblnk

integer :: k, l, m

k = fnblnk(string)
m = lnblk(string)
l = m

if (k.le.m) then
100 if (l.ge.k) then
        call get_environment_variable(string(k:l),value)
        if (lnblk(value).gt.0) then
            goto 200
        endif
        l = l - 1
        goto 100
    else
        value = ' '
    endif
else
    value = ' '
endif
200    return
end
