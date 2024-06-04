subroutine HowToUse (string,usage)

implicit none

!-pass
character*(*)     :: string, usage
!-local
integer, external :: fnblnk, lnblk, iflu
character(len=80) :: prgnam
integer           :: lue, k, l, ll, m, n, length


lue=iflu('stderr')
length = LEN(usage)
usage  = ' '
ll     = 10

if (length.le.ll) then
    write(lue,3) ll
    return
endif

k = fnblnk(string)
l = lnblk(string)

if (k.gt.l .or. l.le.0) then
    write(lue,1)
1        format(1x,'HowToUse: zero usage string error')
    return
endif

call get_command_argument(0, prgnam)

m = fnblnk(prgnam)
n = lnblk(prgnam)

if (m.gt.n .or. n.le.0) then
    write(lue,2)
2        format(1x,'HowToUse: zero prgname error')
    return
endif

ll = 7 + (n-m+1) + 1 + (l-k+1) 

if (ll.gt.length) then
    write(lue,3) ll
3        format(1x,'HowToUse: string length of "usage" should be',1x,'>= ',I3)
    return
endif

usage = 'Usage: '//prgnam(m:n)//' '//string(k:l)

return

end subroutine
