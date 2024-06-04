integer function ifbopen(lu,f,stat)

implicit none

!-pass
integer           :: lu     
character         :: f*(*),stat*(*)
!-local
integer, external :: lnblk, iflu
integer                 :: lui, lue, luo
integer           :: ierr
integer           :: l, m
integer           :: ios

lue  = iflu('stderr')
luo  = iflu('stdout')
lui  = iflu('stdin')
ierr = 0
l    = LEN(f)

if (stat.eq.'old' .or. stat.eq.'OLD' .or. &
    stat.eq.'new' .or. stat.eq.'NEW' .or. &
    stat.eq.'unknown' .or. stat.eq.'UNKNOWN') then
    if (l.gt.0) then
        m=lnblk(f)
        if (m.gt.0) then
            if (lu.ne.lue .and. lu.ne.lui .and. lu.ne.luo) then
                open(unit=lu,file=f(1:m),status=stat, &
              form='unformatted',err=10,iostat=ios)
                rewind (unit=lu,err=20,iostat=ios)
            else
                write(lue,1)
1                format(1x,'ifbopen: unit number error')
                ierr = 1000
            endif
        else
            write(lue,2)
2                format(1x,'ifbopen: no file')
            ierr = 1000
        endif
    else
        write(lue,2)
        ierr = 1000
    endif
else
    write(lue,3)
3        format(1x,'ifbopen: status error')
    ierr = 1000
endif
100    ifbopen = ierr
return

10    write(lue,4) ios,'"',f(1:m),'"'
4    format(1x,'ifbopen: open error ',I3,' on ',80a)
call prerr('ifbopen')
ierr = ios
goto 100

20    write(lue,5) ios,'"',f(1:m),'"'
5    format(1x,'ifbopen: rewind error ',I3,' on ',80a)
call prerr('ifbopen')
ierr = ios
goto 100

end function
