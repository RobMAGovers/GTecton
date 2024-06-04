subroutine exitp(n)

implicit none

!-pass
integer           :: n
!-local
integer, external :: iflu
integer           :: i, ios, lu, lue

logical open
lue=iflu('stderr')

do i=0,99

    lu=i

    if (lu.ne.lue) then

        inquire(unit=lu,opened=open)
        if (open) then
            close(unit=lu,err=200,iostat=ios)
        endif

    endif

enddo

call exitpc(n)
return

200    write(lue,1) ios,lu
1    format(1x,'exitp: error ',I4,' on closing unit ',I2)
call exitpc(0)

end subroutine
