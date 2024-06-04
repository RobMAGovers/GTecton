integer function parseCommandLineArguments(input,output,NF,maxpar)

use iomodule, only: stderr

! subroutine parses an input line INPUT into OUTPUT. Parseing characters
! are (unbackslashed) spaces and tabs. Unbackslashed quotes or apostrophes
! may be used to prevent from parseing. Backslashes are used to protect
! single characters.

implicit none

!-pass
integer           :: nf, maxpar
character(len=*)  :: input,output(maxpar)
!-local
character(len=1)  :: letter,TAB,SPACE,QUOTES,QUOTE,BACKSLASH
integer           :: base,length,i,j, parsln
logical           :: skipspaces,backsl,quoted,inputok,OK
integer, external :: lnblk


QUOTE      = CHAR(39)
QUOTES     = CHAR(34)
SPACE      = CHAR(32)
BACKSLASH  = CHAR(92)
TAB        = CHAR(9)
NF         = 0
base       = 1
quoted     = .false.
backsl     = .false.
skipspaces = .false.
OK         = .true.

parseCommandLineArguments = 0


length = lnblk(input)
if ( length.gt.0 ) then
    OK = inputok(input,length)
    if ( OK ) then
        do i=base,length
            letter = input(i:i)
            if (quoted) then
                if ( (letter.eq.QUOTE .or. letter.eq.QUOTES) &
                .and. (.not.backsl) ) then
                    if (i.gt.base) then
                        NF = NF + 1
                        call cpystr(input,base,i-1,output,NF, maxpar)
                    endif
                    base       = i + 1
                    quoted     = .false.
                    backsl     = .false.
                    skipspaces = .false.
                else if (letter.eq.BACKSLASH) then
                    backsl = .true.
                    base = i + 1
                else
                    backsl = .false.
                endif
            else
!               toggle BACKSL logical
                if (backsl) then
                    backsl = .false.
                else if (letter.eq.BACKSLASH) then
                    backsl = .true.
                    base   = i + 1
                elseif (letter.eq.SPACE .or. letter.eq.TAB) then
                    if (skipspaces) then
                        base = i + 1
                    else
                        if (i.gt.base) then
                            NF = NF + 1
                            call cpystr (input,base,i-1,output, &
                                        NF,maxpar)
                        endif
                        base       = i + 1
                        skipspaces = .true.
                    endif
                else if (letter.eq.QUOTE .or. letter.eq.QUOTES) then
                    if (i.gt.base) then
                        NF = NF + 1
                        call cpystr(input,base,i-1,output,NF, &
                                    maxpar)
                    endif
                    quoted     = .true.
                    skipspaces = .false.
                    base       = i + 1
                else
                    skipspaces = .false.
                endif
            endif
        enddo

!        write last word
        if (skipspaces) then
            write(stderr,1)
 1            format(1x,'parsln: lnblk malfunction')
            call exitp(1)
        else
            if (i.gt.base) then
                NF = NF + 1
                call cpystr(input,base,i-1,output,NF,maxpar)
            endif
        endif
    endif
endif

if ( OK ) then
    parsln = 0
else
    parsln = 1
endif

return
end

!-------------------------------------------------------------------------------
logical function inputok(input,length)

use iomodule,           only: stderr

! Checks for unbalanced quotes and apostrophes. Return values: true if OK,
! false if non-printable ASCII or unbalanced input.

implicit none

character(len=*) input
logical backsl,OK,quoting
integer space,tilde,tab,backslash,quotes,apostrophe
integer :: length,nonasc,i,letter

 data space,tilde,tab,backslash,quotes,apostrophe &
/32,126,9,92,34,39/

nonasc = 0
OK     = .true.
backsl = .false.

if (length.gt.0) then
!        check if ASCII based normal input
    do i=1,length
        letter = ICHAR(input(i:i))

        if ( (letter .lt. space .or. &
              letter .gt. tilde) .and. &
              letter .ne. tab ) then
            nonasc = nonasc + 1
            OK = .false.
        endif

    enddo

    if ( .not. OK ) then
        write(stderr,1) nonasc
 1        format(1x,i3,' non-printable ASCII character(s)')
    endif
!        check for unbalanced quotes and apostrophes
    i = 1
20  if ( i.le.length ) then
        if (OK) then
            letter = ICHAR(input(i:i))
            if (backsl) then
                backsl = .false.
            else
                if ( letter.eq.backslash ) then
                    backsl = .true.
                else if ( letter.eq.apostrophe ) then
                    OK = quoting (i,input,1,length)
                    if (.not. OK) then
                        write(stderr,2)
 2                        format(1x,'unbalanced apostrophes')
                    endif
                else if ( letter.eq.quotes ) then
                    OK = quoting (i,input,0,length)
                    if (.not. OK) then
                        write(stderr,3)
 3                    format(1x,'unbalanced quotes')
                    endif
                endif
            endif
        endif
        i = i + 1
        goto 20
    endif
endif

inputok = OK
return
end function

!-------------------------------------------------------------------------------

logical function quoting (i,input,mode,length)

use iomodule,           only: stderr

implicit none

character(len=*) input
integer i,j,k,mode,length,quast,backslash
logical OK
data backslash/92/

OK = .true.

if (mode.eq.0) then
!        check for quotes
    quast = 34
else if (mode.eq.1) then
!        check for apostrophes
    quast = 39
else
    write(stderr,1)
 1        format(1x,'quoting: fatal mode error')
    call exitp(1)
endif

i = i + 1
10    if ( i.le.length ) then
    j = ICHAR( input(i:i) )
    k = ICHAR( input(i-1:i-1) )
    if (.not.(j.eq.quast .and. k.ne.backslash)) then
        i = i + 1
        goto 10
    endif
else
    OK = .false.
endif

quoting = OK
return
end function

!-------------------------------------------------------------------------------
subroutine cpystr (input,i,j,output,NF,maxpar)

use iomodule,only: stderr

implicit none

!-pass
integer          :: maxpar, i, j, nf
character(len=*) :: input,output(maxpar)
!-local
integer          :: l

l = LEN(input)
if ( i.ge.1 .and. j.le.l .and. i.le.j ) then
    if (NF.le.maxpar) then
        output(NF) = input(i:j)
    else
        write(stderr,*) 'More than MAXPAR=',maxpar,' command line arguments'
        write(stderr,*) 'Most likely you are monitoring too many nodes or elements'
        call exitp(1)
    endif
else
    write(stderr,2)
2   format(1x,'cpystr: fatal string index error')
    call exitp(1)
endif
return
end subroutine
!-------------------------------------------------------------------------------

subroutine parseSampleRange(cmd,min,max,inc,ierr)

! this gets a number range of:
!
! x
! or
! xmin, xmax
! or
! xmin,xmax,xstep
! or
! xmin,xmax;nSamplePoints
!

use iomodule,    only: stderr
use errors,      only: prgerr

implicit none


!    parses real(s) from cmd
!
!-pass
integer ierr
character(len=*) cmd
double precision min,max,inc
!-locl
logical, external :: numeric
integer :: n,m,iString
double precision, external :: chreal
integer, external :: lnblk,fnblnk

integer :: nSamplePoints

integer :: subStart, subEnd

integer :: nCommas
integer :: nSemiColons
integer :: commaPositions(2)
integer :: semiColonsPositions(2)

logical :: minSet, maxSet, incSet

! initialize local data
nCommas    = 0
nSemiColons = 0
commaPositions = 0
semiColonsPositions = 0

! initialize output values
min = 0d0
max = 0d0
inc = 0d0

minSet = .False.
maxSet = .False.
incSet = .False.

! string in CMD(N:M)
n = fnblnk(cmd)  ! first non-blank
m = lnblk(cmd)   ! last  non-blank

if (n.gt.m) then
    write(stderr,1)
 1        format(1x,'sample: n>m error')
    ierr = PRGERR
    return
endif

do iString = n, m
    if (cmd(iString:iString) .eq. ",") then
        nCommas = nCommas + 1
        commaPositions(nCommas) = iString
    else if (cmd(iString:iString) .eq. ";") then
        nSemiColons = nSemiColons + 1
        semiColonsPositions(nSemiColons) = iString
    else
        ! the character is a number
    endif
enddo


! the number of commas and semisolons tells us the format 
if    (nCommas.eq.0 .and. nSemiColons.eq.0) then
    ! only one coordinate given

    read(cmd(n:m),*) min
    max = min
    inc = 1d0

else if (nCommas.eq.1 .and. nSemiColons.eq.0) then
    ! only a min and max coordinate given

    read(cmd(n:commaPositions(1)-1),*) min
    read(cmd(commaPositions(1)+1:m),*) max
    ! set default to 20 points
    inc    = (max - min) / 19d0

else if    (nCommas.eq.2 .and. nSemiColons.eq.0) then
    ! min and max coordinate and increment given

    read(cmd(n:commaPositions(1)-1),*) min
    read(cmd(commaPositions(1)+1:commaPositions(2)-1),*) max
    read(cmd(commaPositions(2)+1:m),*) inc

else if    (nCommas.eq.1 .and. nSemiColons.eq.1) then
    ! min and max coordinate and number of sample points given

    read(cmd(n:commaPositions(1)-1),*) min
    read(cmd(commaPositions(1)+1:semiColonsPositions(1)-1),*) max
    read(cmd(semiColonsPositions(1)+1:m),*) nSamplePoints
    inc = (max - min) / dble(nSamplePoints-1)

else
    write(stderr,*) "Sample argument with", nCommas, &
                    "commas and", nSemiColons, &
                    "semicolons not recognized."
    write(stderr,*) "Allowed formats are:"
    write(stderr,*) "min"
    write(stderr,*) "min,max"
    write(stderr,*) "min,max,increment"
    write(stderr,*) "min,max;nPoints"

endif


return
end subroutine

