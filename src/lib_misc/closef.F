subroutine closef(lu)

implicit none

!-pass
integer           :: lu
!-local
integer, external :: iflu
logical           :: od
integer           :: stderr,stdout,stdin
integer           :: ios

stderr=iflu('stderr')
stdout=iflu('stdout')
stdin =iflu('stdin')
!
if (lu.ne.stderr .and. &
    lu.ne.stdin  .and. &
    lu.ne.stdout) then
    inquire(unit=lu,opened=od)
    if (od) then
        close(unit=lu,err=20,iostat=ios)
    endif
else
    write(stderr,1)
1        format(1x,'closef: unit number error')
endif

10    return
20    write(stderr,2) ios
2    format(1x,'closef: close error ',I3)
call prerr('closef')
goto 10
end subroutine
