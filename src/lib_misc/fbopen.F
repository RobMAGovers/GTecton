subroutine fbopen(lu,f,stat)

implicit none

!-pass
character         :: f*(*),stat*(*)
integer           :: lu
!-local
integer, external :: lnblk, iflu
integer           :: stderr,stdin,stdout
logical           :: AlOpen
integer           :: l, m, ios

stderr=iflu('stderr')
stdout=iflu('stdout')
stdin =iflu('stdin')

l = LEN(f)
if (stat.eq.'old' .or. stat.eq.'OLD' .or. &
    stat.eq.'new' .or. stat.eq.'NEW' .or. &
    stat.eq.'unknown' .or. stat.eq.'UNKNOWN') then
    if (l.gt.0) then
        m=lnblk(f)
        if (m.gt.0) then
            if (lu.ne.stderr.and.lu.ne.stdin.and.lu.ne.stdout) then
            inquire(unit=lu,opened=AlOpen)
                if (AlOpen) then
                    close (lu,err=10,iostat=ios)
                endif
                open(unit=lu,file=f(1:m),status=stat, &
                  form='unformatted',err=10,iostat=ios)
                rewind (unit=lu,err=20,iostat=ios)
            else
                write(stderr,1)
1                    format(1x,'fbopen: fatal unit number error')
                call exitp(1)
            endif
        else
            write(stderr,2)
2                format(1x,'fbopen: no file')
            call exitp(1)
        endif
    else
        write(stderr,2)
        call exitp(1)
    endif
else
    write(stderr,3)
3        format(1x,'fbopen: fatal status error')
endif

return

10    write(stderr,4) ios,'"',f(1:m),'"'
4    format(1x,'fbopen: fatal open error ',I3,' on ',80a)

call prerr('fbopen')
call exitp(1)

20    write(stderr,5) ios,'"',f(1:m),'"'
5    format(1x,'fbopen: rewind error ',I3,' on ',80a)
call prerr('fbopen')

end subroutine
