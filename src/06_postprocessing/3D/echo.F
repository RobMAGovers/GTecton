!-------------------------------------------------------------------------------
subroutine echo (cmdln,cmd,ncmd)

use iomodule, only: stderr, stdin, stdout

implicit none


!-pass
character cmd*(*),cmdln*(*)
dimension cmd(*)
integer   :: ncmd

!-local
integer :: k,n,m
integer :: lnb

!-external
integer, external :: lnblk ! from tools
integer, external :: indx

if (ncmd.gt.1) then
    k=2
    if (cmd(2)(1:3).eq.'-n ') k=3
    if (ncmd.ge.k) then
        lnb = lnblk(cmd(k))
        if (lnb.le.0) then
            write(stderr,*) ' echo: parser error'
            call exitp(1)
        endif
        n = indx(cmdln,cmd(k)(1:lnb))
        if (n.gt.0) then
            m = lnblk(cmdln)
            if (cmdln(n:n).eq.'"' .or. cmdln(n:n).eq.'''') n=n+1
            if (cmdln(m:m).eq.'"' .or. cmdln(m:m).eq.'''') m=m-1
            if (n.le.m) then
                if (k.eq.2) then
                    write(stdout,1) cmdln(n:m)
 1                format(1x,256A)
                else
                    write(stdout,2,advance='no') cmdln(n:m)
 2                format(1x,A)
                endif
            else
                write(stdout,4)
            endif
        else
            write(stderr,3)
 3            format(1x,'echo: cmdln indx error')
        endif
    endif
else
    write(stdout,4)
 4        format(1X)
endif

return
end
