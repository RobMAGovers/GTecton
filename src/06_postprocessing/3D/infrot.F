subroutine INFROT (file,XPT,IEN,X,D,V,LMF,TFAULT, &
                   SKEW,VX,RATE,ierr)

use modeltopology,   only: ndof, nsd, nen
use constants,       only: useDXE, half, one
USE MODELDEFINITION
use iomodule,        only: crdfil, elmfil, ascout, icoord, ifelm, &
                           lumon, verbose, stderr, ifelm
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: npbased
!
! displacement gradients or infinitesimal rotation
!

! VX is actually modeldatact%DELX...


implicit none
!-pass
logical RATE
integer ierr
character(len=*) file
real XPT
integer :: ien, lmf
double precision :: x, d, v, tfault, vx, skew
dimension XPT(3,3,*), &
   IEN(NEN,*),X(NSD,*),D(NDOF,*),V(NDOF,*), &
   LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   VX(NDOF,NEN,*),SKEW(2,*)

#include "lenchar.i"
!-locl
character(len=lenchar) :: name
integer           :: num,nn,n,lu
double precision  :: BP(NDOF,NSD),XL(NSD,NEN),DL(NDOF,NEN)
integer, external :: nextlu,ifopen,wrdxFloat, lnblk
integer           :: i, j, ios
!-init
!
if (elmon) then
    num = nnel
else
    num = NUMEL
endif

do nn=1,num
    n = nn
    if (elmon) then
        n = monel(nn)
    endif

    call LCOORD (X,XL,IEN(1,n))
    call LDISP (DL,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,102)
    call ADDSNE (DL,n,NDOF,NEN,0)

    call REZONE (XL,DL,ONE)

    if (RATE) then
        call LDISP (DL,V,IEN(1,n),NDOF,NEN)
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
        call BPMATRIX (XL,DL,.false.,BP,ierr)
    else
        call BPMATRIX (XL,DL,.false.,BP,ierr)
    endif

    do i=1,3
        do j=1,3
            XPT(i,j,nn) = HALF * (BP(j,i)-BP(i,j))
        enddo
    enddo

    if (elmon) then
        if (verbose) then
            write(luMon,14) n,XPT(1,2,nn),XPT(1,3,nn),XPT(2,3,nn)
14               format(1x,'element ',I5,': Rxy, Rxz, Ryz = ', &
             3(1X,1PG16.6))
        else
            write(luMon,*) XPT(1,2,nn),XPT(1,3,nn),XPT(2,3,nn)
        endif
    endif

enddo

if (elmon) then
    return
endif
!
!    Write header data and, in case of ASCII output, data itself
lu = nextlu(0)
ierr = ifopen(lu,file,'unknown')

if (ierr.ne.0) then
    write(*,*) "Could not open file", trim(file), "error: ", ierr
    stop "Leaving Gtecton"
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

1   format('# nodal points:'/'object 1 class array type',1x, &
'float rank 1 shape ',I1,' items ',i8,1X,A6/ &
'data file ',256A)

if (ifelm.eq.1) then
    write(lu,2) NEN,NUMEL,'binary',elmfil(1:lnblk(elmfil))
else
    write(lu,2) NEN,NUMEL,'ascii ',elmfil(1:lnblk(elmfil))
endif

2   format(/'# connectivity:'/'object 2 class array type',1x, &
'int rank 1 shape ',I1,' items ',i8,1X,A6/'data file ', &
256A)

if (NTYPE.eq.5) then
    write(lu,3)
 3        format('attribute "element type" string "cubes"'/ &
      'attribute "ref" string "positions"'/)
else
    write(lu,4)
 4       format('attribute "element type" string "tetrahedra"'/ &
     'attribute "ref" string "positions"'/)
endif

if (ascout) then
    write(lu,5) name,NUMEL
 5       format('# ',A6,' tensor data:'/ &
     'object 3 class array type float rank 2',1x, &
     'shape 3 3 items ',i8,' data follows')
    write(lu,6,err=1000,iostat=ios) &
       (((XPT(i,j,n),j=1,3),i=1,3),n=1,NUMEL)
 6        format(3(1x,1PG14.6))
else
    write(lu,7) name,NUMEL
 7       format('# ',A6,' tensor data:'/ &
     'object 3 class array type float rank 2',1x, &
     'shape 3 3 items ',i8,' ieee data 0')
endif

write(lu,8)
 8   format('attribute "dep" string "connections"'// &
 'object "irregular positions irregular',1x, &
 'connections" class field'/'component "positions"',1x, &
 'value 1'/'component "connections" value 2'/'component',1x, &
 '"data" value 3'/'end')
 call closef(lu)
!    write binary data to file

if (.not.ascout) then
    ierr = wrdxFloat (file,'a',XPT,NUMEL*3*3)
    if (ierr.ne.0) then
        write(*,*) "Could not open file", trim(file), " Error:", ierr
        stop "Leaving GTecton"
    endif
endif
!
return
!
1000    write(stderr,12) ios,'"',file(1:lnblk(file)),'"'
12   format(1x,'infrot: write error ',I3,' on file ',80A)
call prerr('infrot')
!ierr = EXEC
return
!
1100    write(stderr,13) ios,luMon
13   format(1x,'infrot: write err ',I3,' on monitor logical unit ', &
    I3)
call prerr('infrot')
!ierr = 
return

end
