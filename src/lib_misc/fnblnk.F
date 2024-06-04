integer function fnblnk(a)

implicit none

!-pass
character         :: a*(*)
!-local
integer           :: lue
integer, external :: lnblk, iflu
integer           :: i, j, l
integer           :: ispace, itab

lue = iflu('stderr')
ISPACE = ICHAR(' ')
ITAB   = ICHAR('	')
l = LEN(a)

if (l.le.0) then
    goto 400
endif

if (a.eq.' ') then
    fnblnk = 0
    return
endif

l = lnblk(a)
i = 1

if (l.gt.0) then
100 if (i.gt.l) then
        goto 400
    endif
    j = ICHAR(a(i:i))
    if (j.ne.ISPACE.and.j.ne.ITAB) then
        goto 200
    endif
    i = i + 1
    goto 100
endif
200    j = ICHAR(a(1:1))
fnblnk = i
300    return

400    write(lue,1) 'fnblnk: error on "',a,'"'
1    format(1x,160a)
fnblnk = -1
goto 300
end function
