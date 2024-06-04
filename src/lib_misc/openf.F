subroutine openf(lu,f,stat)

implicit none

!-pass
integer           :: lu
character         :: f*(*),stat*(*)
!-local
logical           :: AlOpen,goodlu
integer, external :: lnblk, iflu
integer                 :: lui, luo, lue
integer           :: l, m
integer           :: ios

lui=iflu('stdin')
luo=iflu('stdout')
lue=iflu('stderr')

goodlu = (lu.ge.0.and.lu.le.99)
goodlu = (goodlu.and.lu.ne.lue.and.lu.ne.luo.and.lu.ne.lui)

l = LEN(f)

if (stat.eq.'old' .or. stat.eq.'OLD' .or. &
    stat.eq.'new' .or. stat.eq.'NEW' .or. &
    stat.eq.'unknown' .or. stat.eq.'UNKNOWN') then
    if (l.gt.0) then
        m=lnblk(f)
        if (m.gt.0) then
            if (goodlu) then
                inquire(unit=lu,opened=AlOpen)
                if (AlOpen) then
                    close (lu,err=10,iostat=ios)
                endif

                open(unit=lu,file=f(1:m),status=stat,&
                      err=10,iostat=ios)

                rewind (unit=lu,err=20,iostat=ios)
            else
                write(lue,1) lu
1                    format(1x,'openf: unit number error: lu=',I3)
                call exitp(1)
            endif
        else
            write(lue,2)
2               format(1x,'openf: no file')
            call exitp(1)
        endif
    else
        write(lue,2)
        call exitp(1)
    endif
else
    write(lue,3)
3        format(1x,'openf: fatal status error')
    stop "Leaving GTECTON..."
endif
return

10    write(lue,4) ios,'"',f(1:m),'"'
4    format(1x,'openf: fatal open error ',I3,' on ',80a)
call prerr('openf')
call exitp(1)

20    write(lue,5) ios,'"',f(1:m),'"'
5    format(1x,'openf: rewind error ',I3,' on ',80a)
call prerr('openf')

end subroutine
