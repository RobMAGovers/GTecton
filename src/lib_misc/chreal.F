    double precision function chreal (number)

    implicit none
    double precision, parameter :: ZERO=0D0

!-pass
    character*(*) number
!-locl
    character(len=10) :: cnumber
    character(len=15) :: numberFormat
    integer fnblnk,lnblk,lue,iflu
    external lnblk,fnblnk,iflu

    integer :: l,n,m
    integer :: ierr, ios

!-init
    lue=iflu('stderr')

    chreal = ZERO
    l = lnblk(number)
    if (l.le.0) then
        return
    endif

    write(cnumber,'(I10)',err=20,iostat=ierr) l
    n = fnblnk(cnumber)
    m = lnblk(cnumber)
    numberFormat = '(F' // cnumber(n:m) // '.0)'
    n = fnblnk(numberFormat)
    m = lnblk(numberFormat)
    read (number,fmt=numberFormat(n:m),err=10,end=20,iostat=ios) chreal

    return

10    write(lue,11) ios,'"',number(1:l),'"'
11    format(1x,'chreal: I/O error: ',I3,' on ',80a)
    call exitp(1)

20    write(lue,21) ios,'"',number(1:l),'"'
21    format(1x,'chreal: EOF error: ',I3,' on ',80a)
    call exitp(1)
    end
