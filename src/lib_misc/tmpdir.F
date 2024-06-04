character*(*) function tmpdir()

! Function for setting temporary files directory
! If environment variable MYTMP has been set, this directory is used
! TMPDIR returns a stripped directory name, i.e. without a '/' at its tail

implicit none

integer,  parameter :: NAMLEN=256
!-locl
character(len=NAMLEN) :: DEFDIR,MYTMP
integer, external :: fnblnk,lnblk, iflu
integer :: lue
integer :: ld, lf, lm, md, mm
integer :: ierr
logical :: there
integer, external :: csystem

lue = iflu('stderr')

!    Default directory
DEFDIR = '/tmp'



lf = LEN(tmpdir)
!    Remove trailing "/"'s from DEFDIR
call strpit(DEFDIR)
ld = fnblnk(DEFDIR)
md = lnblk(DEFDIR)

MYTMP = ' '
call get_environment_variable('MYTMP',MYTMP)

if (MYTMP(1:1).ne.' ') then
    inquire (file=MYTMP,exist=there)

    if (.not. there) then
        ! see if we can make it
        ierr = csystem('mkdir -p '//trim(MYTMP))
        if (ierr.ne.0) then
            write(0,*) "Error: Not found temp directory MYTMP: ", trim(MYTMP)
            write(0,*) "Could not make directory ", trim(MYTMP), ". Error:", ierr
            write(0,*) "Possible solution: "
            write(0,*) "- create directory ", trim(MYTMP), "manually."
            write(0,*) "- change environement variable ", trim(MYTMP)
            write(0,*) "- change permissions of parent directories of", trim(MYTMP)
            stop "Leaving"
        endif
    endif
else
    MYTMP = "."
endif


!    Remove trailing "/"'s from MYTMP
call strpit(MYTMP)
lm = fnblnk(MYTMP)
mm = lnblk(MYTMP)
if (mm.le.lm .or. mm.le.0) then
!        set to default
    if ((md-ld).gt.lf) then
    write(lue,1)
1        format(1x,'tmpdir: character length not sufficient')
    call exitp(1)
    endif
!        do a write test on scratch file
    call wperm(DEFDIR)
    tmpdir = DEFDIR(ld:md)
else
!        set to MYTMP
    if ((mm-lm).gt.lf) then
    write(lue,1)
    call exitp(1)
    endif
!        do a write test on scratch file
    call wperm(MYTMP)
    tmpdir = MYTMP(lm:mm)
endif
return
end
!-----------------------------------------------------------------------
subroutine strpit (DIR)

implicit none

!-pass
character*(*) DIR
!-local
integer, external :: lnblk,fnblnk
integer           :: i, istart, iend

istart=fnblnk(DIR)
iend=lnblk(DIR)
i = iend

100    if (i.lt.istart .or. i.le.0) then
    goto 200
endif

if (DIR(i:i).ne.'/') then
    goto 200
endif

DIR(i:i) = ' '
i = i - 1

goto 100
200    return
end subroutine
!-----------------------------------------------------------------------
subroutine wperm(DIR)

implicit none

!-pass
character*(*)      :: DIR
!-locl
character(len=512) :: f
integer, external  :: nextlu, fnblnk, lnblk, unlink, iflu
integer            :: i, j, lu, lue


lue = iflu('stderr')
i   = fnblnk(DIR)
j   = lnblk(DIR)

if (j.lt.i .or. j.le.0) then
    write(lue,1)
1        format(1x,'tmpdir: directory passed to write test missing')
    call exitp(1)
endif

f  = DIR(i:j) // '/Scratch'
lu = nextlu(0)

open(unit=lu,file=f,form='UNFORMATTED',err=100)
rewind(unit=lu,err=100)
write(lu,err=100) 1
close(unit=lu,err=100)
!todo test whether this tm call actually gives the errors in the run
call rm(f)
return
100    call prerr('tmpdir')
call exitp(1)
end
