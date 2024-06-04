integer function ifopen(lu,f,stat)

implicit none

character          :: f*(*),stat*(*)
integer, external  :: lnblk, iflu
logical            :: goodlu,isOpen
integer            :: lu, lue, luo, lui

integer            :: l,m
integer            :: ierr, ios

lue=iflu('stderr')
luo=iflu('stdout')
lui=iflu('stdin')

! determining whether lthe file handler lu qualifies
goodlu = (lu.ge.0.and.lu.le.99)

goodlu = (goodlu    .and. &
          lu.ne.lue .and. &
          lu.ne.luo .and. &
          lu.ne.lui)

ierr = 0
l = LEN(f)

if (stat.eq.'old' .or. stat.eq.'OLD' .or. &
    stat.eq.'new' .or. stat.eq.'NEW' .or. &
    stat.eq.'unknown' .or. stat.eq.'UNKNOWN') then
    if (l.gt.0) then
        m=lnblk(f)
        if (m.gt.0) then
            if (goodlu) then

                inquire(unit=lu,opened=isOpen)

                if (isOpen) then
                    close (lu,err=10,iostat=ios)
                endif

                open(unit=lu,file=f(1:m),status=stat, &
                  form='formatted',err=15,iostat=ios)
                rewind (unit=lu,err=20,iostat=ios)

            else
                write(lue,1) lu
    1                format(1x,'ifopen: unit number error: lu=',I3)
                ierr = 1000
            endif
        else
            write(lue,2)
    2            format(1x,'ifopen: no file')
            ierr = 1000
        endif
    else
        write(lue,2)
        ierr = 1000
    endif
else
    write(lue,3)
3        format(1x,'ifopen: status error')
    ierr = 1000
endif


100    ifopen = ierr
    return

! error handling

10    write(lue,4) ios,'"',f(1:m),'"'
    4    format(1x,'ifopen: close error ',I3,' on ',80a)
    call prerr('ifopen')
    ierr = ios
    goto 100

15    write(lue,5) ios,'"',f(1:m),'"'
    5    format(1x,'ifopen: open error ',I3,' on ',80a)
    call prerr('ifopen')
    ierr = ios
    goto 100

20    write(lue,6) ios,'"',f(1:m),'"'
    6    format(1x,'ifopen: rewind error ',I3,' on ',80a)
    call prerr('ifopen')
    ierr = ios
    goto 100

    end
