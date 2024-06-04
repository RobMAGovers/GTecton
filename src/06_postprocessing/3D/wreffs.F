subroutine wreffs (file,T,XPT,name,flag,ierr)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: ascout, icoord, stderr, ifelm, verbose, &
                           crdfil, elmfil, lumon
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: npbased
use errors,          only: exec
!
! Writes second invariants or components of the deviator of tensor T
! to a DX file
!
implicit none
!-pass
character(len=*) file, name, flag
double precision :: XPT
double precision :: T
dimension T(NSTR,*),XPT(*)
integer :: ierr


!-locl
#include "lenchar.i"
character(len=15)component
double precision  :: sdev(6)
integer, external :: effs,lnblk,wrdxFloat,nextlu,ifopen
integer :: lu
integer :: ios
integer :: n, nn, num
double precision :: SE
!
if (flag.eq.'e') then
    component = 'effective'
else if (flag.eq.'xx') then
    component = 'xx-component of'
else if (flag.eq.'yy') then
    component = 'yy-component of'
else if (flag.eq.'zz') then
    component = 'zz-component of'
else if (flag.eq.'xy'.or.flag.eq.'yx') then
    component = 'xy-component of'
else if (flag.eq.'xz'.or.flag.eq.'zx') then
    component = 'xz-component of'
else if (flag.eq.'yz'.or.flag.eq.'zy') then
    component = 'zz-component of'
else
    write(stderr,*) 'wreffs: illegal flag argument'
    call exitp(1)
endif
!
if (.not.elmon) then
    do n=1,NUMEL
        if (flag.eq.'e') then
            XPT(n) = REAL(effs(T(1,n)))
        else if (flag.eq.'xx') then
            call DEVSTN(T(1,n),sdev)
            XPT(n) = REAL(sdev(1))
        else if (flag.eq.'yy') then
            call DEVSTN(T(1,n),sdev)
            XPT(n) = REAL(sdev(2))
        else if (flag.eq.'zz') then
            call DEVSTN(T(1,n),sdev)
            XPT(n) = REAL(sdev(3))
        else if (flag.eq.'xy'.or.flag.eq.'yx') then
            call DEVSTN(T(1,n),sdev)
            XPT(n) = REAL(sdev(4))
        else if (flag.eq.'xz'.or.flag.eq.'zx') then
            call DEVSTN(T(1,n),sdev)
            XPT(n) = REAL(sdev(5))
        else if (flag.eq.'yz'.or.flag.eq.'zy') then
            call DEVSTN(T(1,n),sdev)
            XPT(n) = REAL(sdev(6))
        endif
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
        write(lu,5) component,name,NUMEL
 5           format('# ',A15,1X,A20,' data:'/ &
         'object 3 class array type float rank 0',1x, &
         'items ',i8,' data follows')
        write(lu,6,err=1000,iostat=ios) (XPT(n),n=1,NUMEL)
 6        format(4(1x,1PG14.6))
    else
        write(lu,7) component,name,NUMEL
 7           format('# ',A15,1X,A20,' data:'/ &
         'object 3 class array type float rank 0',1x, &
         'items ',i8,' ieee data 0')
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
        ierr = wrdxFloat (file,'a',XPT,NUMEL)
        if (ierr.ne.0) then
            ierr = EXEC
        endif
    endif

else
    do nn=1,nnel
        n = monel(nn)
        if (flag.eq.'e') then
            SE = effs(T(1,n))
        else if (flag.eq.'xx') then
            call DEVSTN(T(1,n),sdev)
            SE = sdev(1)
        else if (flag.eq.'yy') then
            call DEVSTN(T(1,n),sdev)
            SE = sdev(2)
        else if (flag.eq.'zz') then
            call DEVSTN(T(1,n),sdev)
            SE = sdev(3)
        else if (flag.eq.'xy'.or.flag.eq.'yx') then
            call DEVSTN(T(1,n),sdev)
            SE = sdev(4)
        else if (flag.eq.'xz'.or.flag.eq.'zx') then
            call DEVSTN(T(1,n),sdev)
            SE = sdev(5)
        else if (flag.eq.'yz'.or.flag.eq.'zy') then
            call DEVSTN(T(1,n),sdev)
            SE = sdev(6)
        endif

        if (VERBOSE) then
            write(luMon,9,err=1100,iostat=ios) component,name,n, &
           SE
 9               format(1x,A15,1x,A20,' in element ',i8, &
             ' = ',1PG15.6)
        else
            write(luMon,*,err=1100,iostat=ios) SE
        endif
    enddo
endif
!
return
!
1000    write(stderr,10) ios,'"',file(1:lnblk(file)),'"'
10   format(1x,'wreffs: write error ',I3,' on file ',80A)
call prerr('wreffs')
ierr = EXEC
return
!
1100    write(stderr,11) ios,luMon
11   format(1x,'wreffs: write err ',I3,' on monitor logical unit ', &
   I3)
call prerr('wreffs')
ierr = EXEC
return
end
!-------------------------------------------------------------------------------
subroutine wrcval (file,XPT,name,ierr)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: crdfil, elmfil, ascout, lumon, icoord, ifelm,&
                           stderr, verbose
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: npbased
use errors,          only: exec

!
! Writes body centered data XPT to a DX file
!
implicit none
!-pass
character(len=*) file, name
double precision :: XPT
dimension XPT(*)
integer :: ierr


!-locl
#include "lenchar.i"
integer, external :: lnblk,wrdxFloat,nextlu,ifopen
integer :: ios, lu, n, nn, num

!
if (.not.elmon) then
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
 5           format('# ',A20,' data:'/ &
         'object 3 class array type float rank 0',1x, &
         'items ',i8,' data follows')
        write(lu,6,err=1000,iostat=ios) (XPT(n),n=1,NUMEL)
 6        format(4(1x,1PG14.6))
    else
        write(lu,7) name,NUMEL
 7           format('# ',A20,' data:'/ &
         'object 3 class array type float rank 0',1x, &
         'items ',i8,' ieee data 0')
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
        ierr = wrdxFloat (file,'a',XPT,NUMEL)
        if (ierr.ne.0) then
            ierr = EXEC
        endif
    endif

else
    ! loop on elements
    do nn=1,nnel
        ! find monitored element
        n = monel(nn)
        if (VERBOSE) then
            write(luMon,9,err=1100,iostat=ios) name,n, XPT(nn)
 9               format(1x,A20,' in element ',i8,' = ',1PG15.6)
        else
            write(luMon,*,err=1100,iostat=ios) XPT(nn)
        endif
    enddo
endif
!
return
!
1000    write(stderr,10) ios,'"',file(1:lnblk(file)),'"'
10   format(1x,'wrcval: write error ',I3,' on file ',80A)
call prerr('wrcval')
ierr = EXEC
return
!
1100    write(stderr,11) ios,luMon
11   format(1x,'wrcval: write err ',I3,' on monitor logical unit ', &
       I3)
call prerr('wrcval')
ierr = EXEC
return

end
