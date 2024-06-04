subroutine BCWNKP (IWELM,IWSIDE,IWTIME,WPRES,IEN,X,D,LMF, &
                   TFAULT,SKEW,FC,FILE)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK, refact
use errors,          only: exec

!
implicit none
!-pass
character FILE*(*)

integer :: iwelm, iwside, iwtime, ien, lmf
double precision :: wpres, x, d, tfault, skew, fc
dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
 X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),FC(NSD,4,*)
!-locl
integer           :: nwnk
integer, external :: nextlu, ifopen
double precision  :: xn(3)
integer :: lu
integer :: i, j, k, n, nn
integer :: ierr, ios


if (NUMWNK.le.0) then
    return
endif

nwnk = 0
do nn=1,NUMWNK
    n = IWTIME(nn)
    if (n.lt.0) then
        if (n.ne.-1 .and. NSTEP.ge.-n) then
            cycle
        endif
    else if (n.eq.0) then
        cycle
    else
        if (NSTEP.lt.(n-1)) then
            cycle
        endif
    endif
    n = IWELM(nn)
    nwnk = nwnk + 1
    call FACECRD (IEN(1,n),X,D,LMF(1,1,n),TFAULT, &
                   SKEW,n,IWSIDE(nn), &
                   FC(1,1,nwnk),xn)
enddo

!
lu = nextlu(0)
ierr = ifopen(lu,FILE,'unknown')
if (ierr.ne.0) then
    ierr = EXEC
    return
endif

write(lu,1,err=1000,iostat=ios) nwnk*4
 1   format('# face coordinates:'/'object 1 class array type',1x, &
 'float rank 1 shape 3 items ',i8,' ASCII data follows')

do k=1,nwnk
    write(lu,2,err=1000,iostat=ios) ((FC(i,j,k),i=1,NSD),j=1,4)
 2        format(1PG14.6,1X,1PG14.6,1X,1PG14.6)
enddo

write(lu,3,err=1000,iostat=ios) 4,nwnk
 3   format(/'# connectivity:'/'object 2 class array type',1x, &
 'int rank 1 shape ',I1,' items ',i8,' ASCII data follows')

do i=1,nwnk
    write(lu,4,err=1000,iostat=ios) ((i-1)*4+k, k=0,3)
 4        format(I5,1X,I5,1X,I5,1X,I5)
enddo

write(lu,5,err=1000,iostat=ios)
 5   format('attribute "element type" string "quads"'/ &
 'attribute "ref" string "positions"'/)

write(lu,6,err=1000,iostat=ios) nwnk
 6   format('object 3 class array type float rank 1 shape 3 ', &
 'items ',i8,' ASCII data follows')

nwnk = 1

do nn=1,NUMWNK
n = IWTIME(nn)
    if (n.lt.0) then
        if (n.ne.-1 .and. NSTEP.ge.-n) then
            cycle
        endif
    else if (n.eq.0) then
        cycle
    else
        if (NSTEP.lt.(n-1)) then
            cycle
        endif
    endif
    n = IWELM(nn)

    call FACECRD (IEN(1,n),X,D,LMF(1,1,n),TFAULT, &
                   SKEW,n,IWSIDE(nn), &
                   FC(1,1,nwnk),xn)

    write(lu,7,err=1000,iostat=ios) xn(1)*WPRES(nn), &
           xn(2)*WPRES(nn),xn(3)*WPRES(nn)
 7        format(1PG14.6,1X,1PG14.6,1X,1PG14.6)
enddo

write(lu,8,err=1000,iostat=ios)
 8   format('attribute "dep" string "connections"'// &
'object 4 class field'/'component "positions"',1x, &
 'value 1'/'component "connections" value 2'/'component',1x, &
 '"data" value 3'/'end')

 call closef(lu)

return

1000    write(stderr,1001) ios
1001   format(1x,'bc: write error ',I4,' while writing Winkler ', &
 'pressure data')

call prerr('bc')
ierr = EXEC
return
end subroutine
!-----------------------------------------------------------------------
subroutine BCSTRS (IELSTR,ISSIDE,ISTR,STRS,IEN,X,D,LMF, &
                   TFAULT,SKEW,FILE)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact
use errors,          only: exec

!
implicit none
!-pass
character FILE*(*)
integer          :: ielstr, isside, istr, ien, lmf
double precision :: strs, x, d, tfault, skew

dimension STRS(NSTR,*),IELSTR(*),ISSIDE(*),ISTR(2,*),IEN(NEN,*), &
 X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*)

double precision, allocatable :: faceCoords(:,:,:),TRACT(:,:)

!-locl
integer           :: lu
integer, external :: nextlu, ifopen
double precision  :: xn(3)
integer           :: it0, it1
integer           :: i, j, k, n, nn
integer           :: nstrs
integer           :: ierr, ios


if (NUMSTR.le.0) then
    return
endif

allocate(faceCoords(nsd,4,NUMSTR))
allocate(tract(nsd,NUMSTR))

nstrs = 0
do nn = 1,NUMSTR
    it0 = ISTR(1,nn)
    it1 = ISTR(2,nn)
!        if (DBLE(NSTEP-it0)*DBLE(NSTEP-it1).le.0d0) then
    nstrs = nstrs + 1
    n = IELSTR(nn)
    call FACECRD (IEN(1,n),X,D,LMF(1,1,n),TFAULT, &
                  SKEW,n, ISSIDE(nn),faceCoords(1,1,nstrs),xn)

    TRACT(1,nstrs) = STRS(1,nn)*xn(1) + &
                     STRS(4,nn)*xn(2) + &
                     STRS(5,nn)*xn(3)

    TRACT(2,nstrs) = STRS(4,nn)*xn(1) + &
                     STRS(2,nn)*xn(2) + &
                     STRS(6,nn)*xn(3)

    TRACT(3,nstrs) = STRS(5,nn)*xn(1) + &
                     STRS(6,nn)*xn(2) + &
                     STRS(3,nn)*xn(3)
!        endif
enddo
!
lu = nextlu(0)
ierr = ifopen(lu,FILE,'unknown')

if (ierr.ne.0) then
    ierr = EXEC
    return
endif

write(lu,1,err=1300,iostat=ios) nstrs*4
 1   format('# face coordinates:'/'object 1 class array type',1x, &
 'float rank 1 shape 3 items ',i8,' ASCII data follows')

do k=1,nstrs
    write(lu,2,err=1300,iostat=ios) ((faceCoords(i,j,k),i=1,NSD),j=1,4)
 2        format(1PG14.6,1X,1PG14.6,1X,1PG14.6)
enddo

write(lu,3,err=1300,iostat=ios) 4,nstrs
 3   format(/'# connectivity:'/'object 2 class array type',1x, &
 'int rank 1 shape ',I1,' items ',i8,' ASCII data follows')

do i=1,nstrs
    write(lu,4,err=1300,iostat=ios) ((i-1)*4+k, k=0,3)
 4        format(I5,1X,I5,1X,I5,1X,I5)
enddo

write(lu,5,err=1300,iostat=ios)
 5   format('attribute "element type" string "quads"'/ &
 'attribute "ref" string "positions"'/)

write(lu,6,err=1300,iostat=ios) nstrs
 6   format('object 3 class array type float rank 1 shape 3 ', &
 'items ',i8,' ASCII data follows')

do nn = 1,NUMSTR
    it0 = ISTR(1,nn)
    it1 = ISTR(2,nn)
!        if (DBLE(NSTEP-it0)*DBLE(NSTEP-it1).le.0d0) then
        write(lu,7,err=1300,iostat=ios) (TRACT(i,nn),i=1,NSD)
 7            format(1PG14.6,1X,1PG14.6,1X,1PG14.6)
!        endif
enddo

write(lu,8,err=1300,iostat=ios)
 8   format('attribute "dep" string "connections"'// &
'object 4 class field'/'component "positions"',1x, &
 'value 1'/'component "connections" value 2'/'component',1x, &
 '"data" value 3'/'end')

call closef(lu)
!
return
!
1300    write(stderr,1301) ios
 1301   format(1x,'bc: write error ',I4,' while writing Winkler ', &
 'pressure data')

call prerr('bc')
ierr = EXEC
return
end subroutine
!-----------------------------------------------------------------------
subroutine FACECRD (IEN,X,D,LMF,TFAULT, &
                    SKEW,elemID,ISIDE,FC,XN)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use postProcessData, only: deffac,setria,triang,deformed,CST,opn

!
! returns face coordinates and normal vector
!
implicit none

!-pass
integer IEN,LMF,ISIDE, elemID
double precision X,D,TFAULT,SKEW,FC,XN
dimension IEN(NEN),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN), &
   TFAULT(NDOF,*), &
   SKEW(2,*),FC(NSD,4),XN(NSD)
!-locl
logical quad
integer iien
double precision xl,dl,area
dimension xl(3,8),dl(3,8),iien(4)
!
call LCOORD (X,xl,IEN)
if (deformed) then
    call LDISP  (dl,D,IEN,NDOF,NEN)
    call AddFaultDisplacement (dl,LMF,TFAULT,NDOF,NEN,112)
    call ADDSNE (dl,elemID,NDOF,NEN,0)
    call REZONE (xl,dl,deffac)
endif

if (NTYPE.eq.5) then
    quad = (IEN(3).ne.IEN(4))
else if (NTYPE.eq.6) then
    quad = .false.
endif

call SIDENP (ISIDE,quad,iien)
call FCGEOM (xl,iien,area,XN)

! todo
! really? 1,2,3 to 2,1,3 ?
! Rob, what is happening here?

FC(1,2) = xl(1,iien(1))
FC(1,1) = xl(1,iien(2))
FC(1,3) = xl(1,iien(3))
if (quad) then
    FC(1,4) = xl(1,iien(4))
else
    FC(1,4) = xl(1,iien(3))
endif

FC(2,2) = xl(2,iien(1))
FC(2,1) = xl(2,iien(2))
FC(2,3) = xl(2,iien(3))

if (quad) then
    FC(2,4) = xl(2,iien(4))
else
    FC(2,4) = xl(2,iien(3))
endif

FC(3,2) = xl(3,iien(1))
FC(3,1) = xl(3,iien(2))
FC(3,3) = xl(3,iien(3))

if (quad) then
    FC(3,4) = xl(3,iien(4))
else
    FC(3,4) = xl(3,iien(3))
endif

return
end subroutine
