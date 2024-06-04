logical function numeric (str)

implicit none

!-pass
character*(*) str
!-locl
character(len=1)    :: zero,nine,dot,exps,signs,letter
character(len=1024) :: string
logical             :: dotOK,signOK,expOK,num,echo
integer, external   :: lnblk,fnblnk, iflu
dimension           :: exps(4),signs(2)
integer             :: i, n, m, nexps
integer             :: length, lue

save zero,nine,dot,exps,signs
!-init
data zero,nine,dot,exps,signs /'0','9','.','e','E','d','D','-','+'/

lue     = iflu('stderr')
dotOK   = .true.    
signOK  = .true.
expOK   = .false.
num     = .false.
nexps   = 0
echo    = .false.
letter  = ' '
i       = 1
numeric = .false.

n = fnblnk(str)
m = lnblk(str)

if (n.gt.m .or. m.le.0) then
    return
endif

if (echo) then
    write(lue,1) str(n:m)
1    format(1x,'numeric: input string = ',80a)
endif

length = m - n + 1

if (echo) then
    write(lue,2) length
2    format(1x,'numeric: length of input string = ',i3)
endif

string = str(n:m)

if (echo) then
    write(lue,3) string(1:length)
3        format(1x,'numeric: left-adjusted string = ',80a)
    write(lue,4)
4        format(/5x,'|   | dotOK | signOK | expOK | num | nexps')
    write(lue,5) i,letter,dotOK,signOK,expOK,num,nexps
5        format(1x,i3,' | ',a,' |   ',l1,'   |   ',l1,'    |   ',l1,'   |  ',l1,'  |   ',i1)
    write(lue,6)
6        format(1x,45('-'))
endif

100    if ( i.le.length ) then
    letter = string(i:i)
    if ( letter.ge.zero .and. letter.le.nine ) then

        if (nexps.lt.1) then
            expOK = .true.
        endif

        signOK = .false.
        num    = .true.

    elseif ( letter.eq.dot ) then

        if (dotOK) then
            dotOK = .false.
        else
            num = .false.
            goto 200
        endif

    elseif ( letter.eq.exps(1) .or. &
             letter.eq.exps(2) .or. &
             letter.eq.exps(3) .or. &
             letter.eq.exps(4) ) then

        num = .false.
        if (expOK) then
            dotOK  = .false.
            signOK = .true.
            expOK  = .false.
            nexps  = nexps + 1
        else
            goto 200
        endif

    elseif ( letter.eq.signs(1).or.letter.eq.signs(2) ) then

        if ( signOK ) then
            signOK = .false.
        else
            num = .false.
            goto 200
        endif

    else
        num = .false.
        goto 200
    endif

    if (echo) then
        write(lue,5) i,letter,dotOK,signOK,expOK,num,nexps
    endif

    i = i + 1
    goto 100
endif
200    numeric = num
return
end
