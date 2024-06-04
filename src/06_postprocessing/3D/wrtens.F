subroutine wrtens (file,T,XPT,name,scale,tecton,ierr)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: crdfil, elmfil, ascout, icoord, ifelm, &
                           lumon, stderr, verbose
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: npbased
use errors,          only: exec

!
! Writes tensor symmetric T (TECTON format, rank 2, dim 3x3) data to DX file
!
implicit none
!-pass
character(len=*) file,name
real XPT
logical tecton
double precision :: T
dimension T(NSTR,*),XPT(3,3,*)
double precision :: scale
integer :: ierr


!-locl
#include "lenchar.i"
integer, external :: lnblk,wrdxFloat,nextlu,ifopen
integer :: i, j, n, nn, num
integer :: lu, ios
!
 if ((.not.elmon).and.(.not.tecton)) then
!        Set tensor
do n=1,NUMEL
XPT(1,1,n) = REAL(T(1,n)*scale)
XPT(1,2,n) = REAL(T(4,n)*scale)
XPT(1,3,n) = REAL(T(5,n)*scale)
XPT(2,1,n) = REAL(T(4,n)*scale)
XPT(2,2,n) = REAL(T(2,n)*scale)
XPT(2,3,n) = REAL(T(6,n)*scale)
XPT(3,1,n) = REAL(T(5,n)*scale)
XPT(3,2,n) = REAL(T(6,n)*scale)
XPT(3,3,n) = REAL(T(3,n)*scale)
enddo
!
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
write(lu,1) NDOF,num,'binary',crdfil(1:lnblk(crdfil))
else
write(lu,1) NDOF,num,'ascii ',crdfil(1:lnblk(crdfil))
endif
 1       format('# nodal points:'/'object 1 class array type',1x, &
     'float rank 1 shape ',I1,' items ',i8,1X,A6/'data file ', &
     256A)
if (ifelm.eq.1) then
write(lu,2) NEN,NUMEL,'binary',elmfil(1:lnblk(elmfil))
else
write(lu,2) NEN,NUMEL,'ascii ',elmfil(1:lnblk(elmfil))
endif
 2       format(/'# connectivity:'/'object 2 class array type',1x, &
     'int rank 1 shape ',I1,' items ',i8,1X,A6/'data file ', &
     256A)
if (NTYPE.eq.5) then
    write(lu,3)
 3           format('attribute "element type" string "cubes"'/ &
         'attribute "ref" string "positions"'/)
else
    write(lu,4)
 4           format('attribute "element type" string "tetrahedra"'/ &
         'attribute "ref" string "positions"'/)
endif
if (ascout) then
write(lu,5) name,NUMEL
 5           format('# ',A6,' tensor data:'/ &
         'object 3 class array type float rank 2',1x, &
         'shape 3 3 items ',i8,' data follows')
write(lu,6,err=1000,iostat=ios) &
       (((XPT(i,j,n),j=1,3),i=1,3),n=1,NUMEL)
 6        format(3(1x,1PG14.6))
else
write(lu,7) name,NUMEL
 7           format('# ',A6,' tensor data:'/ &
         'object 3 class array type float rank 2',1x, &
         'shape 3 3 items ',i8,' ieee data 0')
endif
write(lu,8)
 8       format('attribute "dep" string "connections"'// &
     'object "irregular positions irregular',1x, &
     'connections" class field'/'component "positions"',1x, &
     'value 1'/'component "connections" value 2'/'component',1x, &
     '"data" value 3'/'end')
call closef(lu)
!        write binary data to file
if (.not.ascout) then
ierr = wrdxFloat (file,'a',XPT,NUMEL*3*3)
if (ierr.ne.0) ierr = EXEC
endif
 else if (.not.tecton) then
do nn=1,nnel
n = monel(nn)
if (VERBOSE) then
    write(luMon,9,err=1100,iostat=ios) name,n, &
           (T(i,n),i=1,6)
 9            format(1x,A6,' in element ',I5,' =',6(1X,1PG14.6))
else
    write(luMon,14,err=1100,iostat=ios) (T(i,n),i=1,6)
14            format(6(1X,1PG14.6))
endif
enddo
 else
    do n=1,NUMEL
        write(luMon,10,err=1100,iostat=ios) n,0,1,(T(i,n),i=1,6)
10           format(2i5,i6,6(1PG14.6))
    enddo
    write(luMon,11)
11       format('end stresses')
 endif
!
return
!
1000    write(stderr,12) ios,'"',file(1:lnblk(file)),'"'
12   format(1x,'wrtens: write error ',I3,' on file ',80A)
call prerr('wrtens')
ierr = EXEC
return
!
1100    write(stderr,13) ios,luMon
13   format(1x,'wrtens: write err ',I3,' on monitor logical unit ', &
   I3)
call prerr('wrtens')
ierr = EXEC
return
end
