    character*(*) function tmpnm (init)

!-pass
    character*(*) init
!-locl
    character(len=19)   :: time
    character(len=24)   :: user
    character(len=1024) :: linit, tmpdir, lpath

    integer :: values(8)
    character(len=4) :: year
    character(len=2) :: month
    character(len=2) :: day
    character(len=2) :: hour
    character(len=2) :: minute
    character(len=2) :: second

    integer :: i

    integer fnblnk,lnblk,k,l,m,n
    external fnblnk,lnblk

    lpath = tmpdir()

    k = fnblnk(lpath)
    l = lnblk(lpath) + 1

    lpath(l:l) = '/'

    m = fnblnk(init)
    n = lnblk(init)

    if (m.gt.n .or. n.le.0) then
        linit = 'tmp'
        m = 1
        n = 3
    else
        linit = init
    endif
 
    call get_environment_variable('USER',USER)
!    call date_and_time(DATE=date)
!    date = fdate()
    call date_and_time(VALUES=values)

    write(year,"(i4)") values(1)
    write(month,"(i2)") values(2)
    write(day,"(i2)") values(3)
    write(hour,"(i2)") values(5)
    write(minute,"(i2)") values(6)
    write(second,"(i2)") values(7)

    time = year//"-"//month//"-"//day//"_"//hour//"-"//minute//"-"//second

    do i=1,19
        if (time(i:i).eq.' ') then
            time(i:i)='0'
        endif
    enddo

    tmpnm = linit(m:n) // '.' // time // '.' // USER(1:lnblk(USER))

    return
    end
