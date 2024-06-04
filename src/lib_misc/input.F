    integer function input (word,rnumber,lu)

! subroutine to read an input item until a backslash" is encountered
! and, depending whether input
! is a character, integer or float, to return it into WORD
! or RNUMBER.
! INPUT returns  1: character string
!         2: real
!        -1: EOF
!        -2: input error
!
    implicit none

    character word*(*)
    character(len=13) :: varifmt, frmt
    character(len=2)  :: special
    integer lu,lue
    logical numeric

    double precision :: rnumber,chreal

    integer, external :: fnblnk,lnblk, iflu, indx
    integer :: i,j,l
    integer :: nf, mf

    input = 0

    lue=iflu('stderr')
    special(1:1) = CHAR(92)
    special(2:2) = '"'

    l = LEN(word)
    if (l.gt.0) then
!        write variable format "(AL)"
        frmt=varifmt(l)
        nf=fnblnk(frmt)
        mf=lnblk(frmt)
        if (nf.gt.mf .or. mf.le.0) then
        write(lue,*) ' input: varifmt error'
        input = -2
        return
        endif
        read(unit=lu,fmt=frmt(nf:mf),err=20,end=30) word
!        test for \" (ie start of comment)
        i = indx(word,special)
        if (i.gt.0) then

        do j=i,l
            word(j:j)=CHAR(32)
        enddo

        endif
        if (numeric(word)) then
        rnumber  = chreal(word)
        input = 2
        else
        input = 1
        endif
    else
        write(lue,1)
    1        format(1x,'input: zero length input word')
        call exitp(1)
    endif

10    return

20    input=-2
    call prerr('input')
    goto 10
30    input=-1
    goto 10
    end
