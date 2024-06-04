subroutine wrvect (FILE,NDIM,XPT,IERR)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: crdfil, elmfil, ascout, icoord, ifelm, stderr
use monitorVars,     only: elmon
use plotControl,     only: npbased
use errors,          only: exec
!
! Writes connection-dependent vector data to a DX file
!
implicit none
!-pass
character(len=*) FILE
integer NDIM,IERR
real XPT(NDIM,*)


!-locl
#include "lenchar.i"
integer, external :: lnblk,wrdxFloat,nextlu,ifopen
integer :: i, n, lu, ios, num

!
 if (.not.elmon) then
!        Write header data and, in case of ASCII output, data itself
lu = nextlu(0)
ierr = ifopen(lu,file,'unknown')
if (ierr.ne.0) then
ierr = EXEC
return
endif
    if (npbased) then
        num = NUMNP
    else
        num = NUMEL*NEN
    endif
if (icoord.eq.1) then
write(lu,2) NDOF,num,'binary',crdfil(1:lnblk(crdfil))
else
write(lu,2) NDOF,num,'ascii ',crdfil(1:lnblk(crdfil))
endif
 2       format('# nodal points:'/'object 1 class array type',1x, &
     'float rank 1 shape ',I1,' items ',i8,1X,A6/'data file ', &
     256A)
if (ifelm.eq.1) then
write(lu,3) NEN,NUMEL,'binary',elmfil(1:lnblk(elmfil))
else
write(lu,3) NEN,NUMEL,'ascii ',elmfil(1:lnblk(elmfil))
endif
 3       format(/'# connectivity:'/'object 2 class array type',1x, &
     'int rank 1 shape ',I1,' items ',i8,1X,A6/'data file ', &
     256A)
if (NTYPE.eq.5) then
    write(lu,4)
 4           format('attribute "element type" string "cubes"'/ &
         'attribute "ref" string "positions"'/)
else
    write(lu,5)
 5           format('attribute "element type" string "tetrahedra"'/ &
         'attribute "ref" string "positions"'/)
endif
if (ascout) then
write(lu,6) NDIM,NUMEL
 6           format('# vector data:'/ &
         'object 3 class array type float rank 1',1x, &
         'shape ',I1,' items ',i8,' data follows')
write(lu,7,err=1000,iostat=ios) ((XPT(i,n),i=1,NDIM), &
       n=1,NUMEL)
 7        format(3(1x,1PG14.6))
else
write(lu,8) NDIM,NUMEL
 8           format('# vector data:'/ &
         'object 3 class array type float rank 1',1x, &
         'shape ',I1,' items ',i8,' ieee data 0')
endif
write(lu,9)
 9       format('attribute "dep" string "connections"'// &
     'object "irregular positions irregular',1x, &
     'connections" class field'/'component "positions"',1x, &
     'value 1'/'component "connections" value 2'/'component',1x, &
     '"data" value 3'/'end')
call closef(lu)
!        write binary data to file
if (.not.ascout) then
ierr = wrdxFloat (file,'a',XPT,NUMEL*NDIM)
if (ierr.ne.0) ierr = EXEC
endif
 endif
!
return
!
1000    write(stderr,11) ios,'"',file(1:lnblk(file)),'"'
11   format(1x,'wrvect: write error ',I3,' on file ',80A)
call prerr('wrvect')
ierr = EXEC
return
end
