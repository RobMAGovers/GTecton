    subroutine flushf(lu)

! Flushes logical unit LU.
! Note: files, opened with STATUS='SCRATCH', are not flushed
!    direct access files are not flushed either.

    implicit none

!-pass
    integer lu
!-locl
    logical opn,name
    integer stderr,stdin,stdout,flushc
    character(len=256) :: fname
    character(len=11)  :: fm
    character(len=4)   :: blnk
    character(len=10)  :: acc
    character(len=1)   :: oneChar

integer :: i, iflu, ios

    integer, external :: lnblk

!-init
    fname = ' '
    acc   = 'SEQUENTIAL'


    stdin  = iflu('stdin')
    stdout = iflu('stdout')
    stderr = iflu('stderr')

    if (lu.eq.stderr) then
        i = 0
        ios = flushc(i)
        if (ios.ne.0) goto 1000
        return
    else if (lu.eq.stdin) then
        i = 1
        ios = flushc(i)
        if (ios.ne.0) goto 1000
        return
    else if (lu.eq.stdout) then
        i = 1
        ios = flushc(i)
        if (ios.ne.0) goto 1000
        return
    endif

    if (lu.lt.0) then
        write(stderr,1)
    1        format(1x,'flush: called with negative logical unit')
        return
    endif

    inquire (unit=lu,opened=opn,err=1000,iostat=ios)
    if (.not.opn) then
        write(stderr,2) lu
    2        format(1x,'Unit ',I3,' not open')
        return
    endif

!    Scratch file?
    inquire (unit=lu,named=name,err=1000,iostat=ios)
    if (.not.name) then
        return
    endif

    inquire (unit=lu,name=fname,err=1000,iostat=ios)
    if (lnblk(fname).le.0) then
        write(stderr,3)
    3        format(1x,'flush: cannot retrieve file name')
        return
    endif

    inquire (unit=lu,access=acc,err=1000,iostat=ios)
    if (acc.eq.'DIRECT' .or. acc.eq.'direct') return

    inquire (unit=lu,form=fm,err=1000,iostat=ios)
    if (fm.eq.'FORMATTED' .or. fm.eq.'formatted') then
          inquire (unit=lu,blank=blnk,err=1000,iostat=ios)
    endif

    close(unit=lu,status='keep',err=1000,iostat=ios)

    if (fm.eq.'FORMATTED' .or. fm.eq.'formatted') then
        open(unit=lu,blank=blnk,form=fm,file=fname,status='old',err=1000,iostat=ios)
100        read(unit=lu,fmt='(a)',err=1000,end=300) oneChar
        goto 100
    else
        open(unit=lu,form=fm,file=fname,status='old',err=1000,iostat=ios)
200        read(unit=lu,err=1000,end=300) i
        goto 200
    endif

300    return

1000    write(stderr,4) ios
    4    format(1x,'flushf: I/O error ',I3)
    call prerr('flushf')
    return
    end
