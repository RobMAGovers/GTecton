subroutine centrs (file,IEN,X,D,LMF,TFAULT,DXE, &
                   SKEW,XPT,ierr)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: verbose, stderr, lumon, ascout, append
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: elmon, nnel, monel
use errors,          only: ok, exec

!
! Routine writes element center coordinates to a DX file
!
implicit none
!-pass
real XPT(NDOF,*)
character(len=*) file

integer :: ien, lmf
double precision :: x, d, tfault, dxe, skew

dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),DXE(NDOF,NEN,*),SKEW(2,*)

!-locl
integer, external :: ifopen,fnblnk,lnblk,lueof, nextlu
integer           :: i,j,n,lu
integer, external :: wrdxFloat
double precision  :: CENTR(3)
integer           :: ierr, ios
!-init
ierr = OK

if (.not.elmon) then


do n=1,NUMEL

call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
         SKEW, &
       deffac,ierr)

if (ierr.ne.0) then
    write(stderr,2) n
 2            format(1x,'centrs: error in element ',I8)
    ierr = EXEC
    return
endif

XPT(1,n) = REAL(CENTR(1))
XPT(2,n) = REAL(CENTR(2))
XPT(3,n) = REAL(CENTR(3))
enddo

if (ascout) then
        lu = nextlu(0)
        ierr = ifopen(lu,file,'unknown')
        if (ierr.ne.0) then
            ierr = EXEC
            return
        endif
        if (APPEND) then
            ierr = lueof(lu)
            if (ierr.ne.0) then
                ierr = EXEC
                return
            endif
        endif
        do j=1,NUMEL
            write(lu,3,err=1000,iostat=ios) (XPT(i,j),i=1,NDOF)
 3               format(3(1X,1PG14.6))
enddo
        call closef(lu)
    else
        if (APPEND) then
            ierr = wrdxFloat (file,'a',XPT,NUMEL*NDOF)
        else
            ierr = wrdxFloat (file,'w',XPT,NUMEL*NDOF)
        endif
        if (ierr.ne.0) ierr = EXEC
    endif
 else
do i=1,nnel
        n = monel(i)
call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
         SKEW, deffac,ierr)
if (ierr.ne.0) then
    write(stderr,2) n
    ierr = EXEC
    return
endif
        if (VERBOSE) then
            write(luMon,4,err=1100,iostat=ios) n, &
             (CENTR(j),j=1,NSD)
 4               format(1x,'element ',I8,': ',3(1PG14.5,1X))
        else
            write(luMon,5,err=1100,iostat=ios) &
             (CENTR(j),j=1,NSD)
 5               format(3(1x,F15.6))
        endif
enddo
endif
!
return
!
1000    write(stderr,6) ios,'"',file(fnblnk(file):lnblk(file)),'"'
 6   format(1x,': write error ',I3,' on file ',80A)
call prerr('centr')
ierr = EXEC
return
!
1100    write(stderr,7) ios,luMon
 7   format(1x,'centr: write err ',I3,' on monitor logical unit ',I3)
call prerr('centr')
ierr = EXEC
return

end
!-------------------------------------------------------------------------------
subroutine nodals(file,IEN,X,D,LMF,TFAULT,DXE,SKEW,XPT,iflg,ierr)

USE ALGEBRA
USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use meshdatamodule,  only: meshdatactx
use iomodule,        only: stderr, stdout, verbose, append, ascout, lumon
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use monitorVars,     only: fmon, elmon, nnel, nnp, npmon, monnp, monel, monf
use plotControl,     only: npbased, shrink, shrinkFactor
use errors,          only: ok, exec

!
! Subroutine to write nodal coordinates to DX output file.
! In case deformed coordinates have been requested, IFLG determines the
! update type:
!    0=total displacements,
!    1=displacements,
!    2=slippery node displacements,
!    3=faulted node displacements.
!
 implicit none
!-pass

#ifdef EXPERIMENTAL_ReducedSlipAlloc
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

 integer ierr, iflg
 character(len=*) file
 real XPT(NDOF,*)

integer          :: ien,lmf
double precision :: x,d,tfault,skew
 dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),SKEW(2,*)

!-locl
integer i,j,lu,lueof,ios
integer, external :: wrdxFloat,ifopen,fnblnk,lnblk
logical QUAD
double precision :: XL(3,8),DL(3,8),xc(3)
integer          :: IIEN(4)
double precision :: def
integer          :: k,m,n,nn,node,nxpt
integer, external :: nextlu

!-init
 ierr = OK
 def = 0d0
if (DEFORMED.and.(iflg.eq.0.or.iflg.eq.1)) then
    def=deffac
endif
!
 if (fmon) goto 100
!
 if (npbased) then
!        write straight nodal point coordinates
if (.not.npmon) then

    do i=1,NUMNP
        XPT(1,i) = REAL(X(1,i)+def*D(1,i))
        XPT(2,i) = REAL(X(2,i)+def*D(2,i))
        XPT(3,i) = REAL(X(3,i)+def*D(3,i))
    enddo

    if (ascout) then
    lu = nextlu(0)
    ierr = ifopen(lu,file,'unknown')

    if (ierr.ne.0) then
        ierr = EXEC
        return
    endif

    if (APPEND) then
        ierr = lueof(lu)
        if (ierr.ne.0) then
        ierr = EXEC
        return
        endif
    endif
    write(lu,2,err=1000,iostat=ios) &
           ((XPT(i,j),i=1,NDOF),j=1,NUMNP)
 2            format(3(1X,1PG14.6))
    call closef(lu)
    else
    if (APPEND) then
        ierr = wrdxFloat (file,'a',XPT,NUMNP*NDOF)
    else
        ierr = wrdxFloat (file,'w',XPT,NUMNP*NDOF)
    endif
    if (ierr.ne.0) ierr = EXEC
    endif
else
    do i=1,nnp
    n = monnp(i)
    XPT(1,1) = REAL(X(1,n)+def*D(1,n))
    XPT(2,1) = REAL(X(2,n)+def*D(2,n))
    XPT(3,1) = REAL(X(3,n)+def*D(3,n))
    if (VERBOSE) then
        write(luMon,3,err=1100,iostat=ios) n, &
               (XPT(j,1),j=1,NSD)
 3                format(1x,'nodal point ',I8,': ',3(1PG14.5,1X))
    else
        write(luMon,4,err=1100,iostat=ios) &
               (XPT(j,1),j=1,NSD)
 4                format(3(1x,F15.6))
    endif
enddo
endif
 else
!        Write on element by element basis.
if (.not.elmon) then

nxpt = 0
do n=1,NUMEL
!                   Localize coordinates
            call LCOORD (X,XL,IEN(1,n))
    if (DEFORMED.and.(iflg.eq.0.or.iflg.eq.1)) then
!                       Localize displacements
                call LDISP (DL,D,IEN(1,n),NDOF,NEN)
    else
    call CLEAR(DL,NDOF*NEN,"DL")
    endif
!            Add slippery node displacements
    if (DEFORMED.and.NUMSLP.gt.0.and. &
           (iflg.eq.0.or.iflg.eq.2)) then
        call ADDSNE (DL,n,NDOF,NEN,0)

    endif
!            Add faulted node displacements
    if (DEFORMED.and.(iflg.eq.0.or.iflg.eq.3)) &
           call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,105)
!            Update coordinates with displacements
    call REZONE (XL,DL,def)
    if (shrink) then

    call ELCTR (xc,IEN(1,n),X,D,n, &
           TFAULT,SKEW,def,ierr)

    if (ierr.ne.0) return
    do i=1,NEN
        XL(1,i) = shrinkFactor*XL(1,i)+(1d0-shrinkFactor)*xc(1)
        XL(2,i) = shrinkFactor*XL(2,i)+(1d0-shrinkFactor)*xc(2)
        XL(3,i) = shrinkFactor*XL(3,i)+(1d0-shrinkFactor)*xc(3)
    enddo
    endif
    do i=1,NEN
    nxpt = nxpt + 1
    XPT(1,nxpt) = REAL(XL(1,i))
    XPT(2,nxpt) = REAL(XL(2,i))
    XPT(3,nxpt) = REAL(XL(3,i))
    enddo
enddo
    if (ascout) then
    lu = nextlu(0)
    ierr = ifopen(lu,file,'unknown')
    if (ierr.ne.0) then
        ierr = EXEC
        return
    endif
    if (APPEND) then
        ierr = lueof(lu)
        if (ierr.ne.0) then
        ierr = EXEC
        return
        endif
    endif
    write(lu,2,err=1000,iostat=ios) &
           ((XPT(i,j),i=1,NDOF),j=1,nxpt)
    call closef(lu)
    else
    if (APPEND) then
        ierr = wrdxFloat (file,'a',XPT,nxpt*NDOF)
    else
        ierr = wrdxFloat (file,'w',XPT,nxpt*NDOF)
    endif
    if (ierr.ne.0) ierr = EXEC
    endif
else
!        Monitor elements
    do i=1,nnel
    n = monel(i)
!                   Localize coordinates
            call LCOORD (X,XL,IEN(1,n))
    if (DEFORMED.and.(iflg.eq.0.or.iflg.eq.1)) then
!                       Localize displacements
                call LDISP (DL,D,IEN(1,n),NDOF,NEN)
    else
    call CLEAR(DL,NDOF*NEN,"DL")
    endif
!            Add slippery node displacements
    if (DEFORMED.and.NUMSLP.gt.0.and. &
           (iflg.eq.0.or.iflg.eq.2)) then

        call ADDSNE (DL,n,NDOF,NEN,0)

    endif
!            Add faulted node displacements
    if (DEFORMED.and.(iflg.eq.0.or.iflg.eq.3)) &
           call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,106)
!            Update coordinates with displacements
    call REZONE (XL,DL,def)
    do k=1,NEN
    XPT(1,k) = REAL(XL(1,k))
    XPT(2,k) = REAL(XL(2,k))
    XPT(3,k) = REAL(XL(3,k))
    enddo
    if (VERBOSE) then
        write(luMon,5,err=1100,iostat=ios) n, &
               (j,(XPT(k,j),k=1,NSD),j=1,NEN)
 5                   format(1X,'element ',I8,':', &
                 8(/1X,I1,2X,3(1X,1PG14.5)))
    else
        write(luMon,6,err=1100,iostat=ios) &
               ((XPT(k,j),k=1,NSD),j=1,NEN)
 6                format(3(1X,F15.6))
    endif
enddo
endif
 endif
 goto 200

!    monitor faces

100 nxpt = 0
 do nn=1,nnel
n = monel(nn)
m = monf(nn)
call LCOORD (X,XL,IEN(1,n))
if (DEFORMED.and.(iflg.eq.0.or.iflg.eq.1)) then
!        Localize displacements
call LDISP (DL,D,IEN(1,n),NDOF,NEN)
else
call CLEAR(DL,NDOF*NEN,"DL")
endif
!        Add slippery node displacements
if (DEFORMED.and.NUMSLP.gt.0.and.(iflg.eq.0.or.iflg.eq.2)) then

    call ADDSNE (DL,n,NDOF,NEN,0)

endif
!        Add faulted node displacements
if (DEFORMED.and.(iflg.eq.0.or.iflg.eq.3)) &
       call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,107)
quad = (NTYPE.eq.5 .and. IEN(3,n).ne.IEN(4,n))
!        Update coordinates by displacements
call REZONE (XL,DL,def)
call sidenp(m,quad,IIEN)
do i=1,3
node = IIEN(i)
nxpt = nxpt + 1
XPT(1,nxpt) = REAL(XL(1,node))
XPT(2,nxpt) = REAL(XL(2,node))
XPT(3,nxpt) = REAL(XL(3,node))
enddo
if (quad) then
    node = IIEN(4)
else
node = IIEN(3)
endif
nxpt = nxpt + 1
XPT(1,nxpt) = REAL(XL(1,node))
XPT(2,nxpt) = REAL(XL(2,node))
XPT(3,nxpt) = REAL(XL(3,node))
 enddo
 if (ascout) then
lu = nextlu(0)
ierr = ifopen(lu,file,'unknown')
if (ierr.ne.0) then
ierr = EXEC
return
endif
if (APPEND) then
ierr = lueof(lu)
if (ierr.ne.0) then
    ierr = EXEC
    return
endif
endif
write(lu,2,err=1000,iostat=ios) &
       ((XPT(i,j),i=1,NDOF),j=1,nxpt)
call closef(lu)
 else
if (APPEND) then
ierr = wrdxFloat (file,'a',XPT,nxpt*NDOF)
else
ierr = wrdxFloat (file,'w',XPT,nxpt*NDOF)
endif
if (ierr.ne.0) ierr = EXEC
 endif
!
200    return
!
1000    write(stderr,7) ios,'"',file(fnblnk(file):lnblk(file)),'"'
 7    format(1x,'np: write error ',I3,' on file ',80A)
 call prerr('np')
 ierr = EXEC
 return
!
1100    write(stderr,8) ios,luMon
 8    format(1x,'np: write error ',I3,' on monitor logical unit ',I3)
 call prerr('np')
 ierr = EXEC
 return
 end
!-------------------------------------------------------------------------------
subroutine temprs(file,IEN,T,TLST,TX,TLSTX,LMTX,XPT,iflag,ierr)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr, verbose, ifelm, ascout, lumon, &
                           icoord, crdfil, elmfil
use monitorVars,     only: elmon, npmon, nnel, nnp, monnp, monel
use plotControl,     only: irelt, npbased
use errors,          only: ok, exec

!
! Subroutine to write nodal temperatures to DX output file.
!
implicit none
!-pass
integer :: ien, lmtx, iflag
double precision :: t, tlst, tx, tlstx

 integer ierr
 character(len=*) file
 real XPT(*)
 dimension T(*),TLST(*),TX(*),TLSTX(*),IEN(NEN,*),LMTX(NEN,*)



!-locl
#include "lenchar.i"
 integer i,lu
integer, external :: ifopen,fnblnk,lnblk, nextlu, wrdxFloat
 logical monitor
 character(len=32) type
 character(len=6) ftype
double precision :: TL(8)
integer :: ios
integer :: j, k, n, num
double precision :: Te
!-init
 ierr = OK
 if (npbased) then
    monitor = npmon
else
    monitor = elmon
endif
 if (iflag.eq.0) then
type = 'Temperatures'
 else if (iflag.eq.1) then
type = 'Temperature increments'
 else
write(type,1) irelt
 1        format('Temperature change since t=',i5)
 endif
!
 if (.not.monitor) then
!
!        Prepare temperatures output
if (npbased) then
    do i=1,NUMNP
    if (iflag.eq.0) then
        XPT(i) = REAL(T(i))
    else
        XPT(i) = REAL(T(i)-TLST(i))
    endif
enddo
num = NUMNP
else
do n=1,NUMEL
    call LTEMP (IEN(1,n),T,TL,NEN)
    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) &
           call LTEMPX (LMTX(1,n),IEN(1,n),TX,TL,NEN)
    do j=1,NEN
    k = NEN*(n-1)+j
    XPT(k) = REAL(TL(j))
    enddo
    if (iflag.ne.0) then
        call LTEMP (IEN(1,n),TLST,TL,NEN)
    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) &
           call LTEMPX (LMTX(1,n),IEN(1,n),TLSTX,TL,NEN)
        do j=1,NEN
        k = NEN*(n-1)+j
        XPT(k) = XPT(k)-REAL(TL(j))
    enddo
    endif
enddo
num = NUMEL*NEN
endif
!
!        Write file header
    lu = nextlu(0)
    ierr = ifopen(lu,file,'unknown')
    if (ierr.ne.0) then
        ierr = EXEC
        return
    endif
    if (icoord.eq.1) then
        ftype = 'binary'
    else
        ftype = 'ascii '
    endif
    write(lu,3) NDOF,num,ftype,crdfil(1:lnblk(crdfil))
 3       format('# nodal points:'/'object 1 class array type',1x, &
     'float rank 1 shape ',I1,' items ',I8,1X,A6/ &
     'data file ',256A)
    if (ifelm.eq.1) then
        ftype = 'binary'
    else
        ftype = 'ascii '
    endif
    write(lu,4) NEN,NUMEL,ftype,elmfil(1:lnblk(elmfil))
 4       format(/'# connectivity:'/'object 2 class array type',1x, &
     'int rank 1 shape ',I1,' items ',I8,1X,A6/'data file ', &
   256A)
    if (NEN.eq.8) then
write(lu,5)
 5           format('attribute "element type" string "cubes"'/ &
         'attribute "ref" string "positions"'/)
else
write(lu,6)
 6           format('attribute "element type" string "tetrahedra"'/ &
         'attribute "ref" string "positions"'/)
endif

!        if (ascout) then
!                write(lu,7) type,num
!    7           format('# ',A32,' data:'/ &
!                 'object 3 class array type float rank 0',1x, &
!                 ' items ',I8,' data follows')
!                write(lu,8,err=1000,iostat=ios) (XPT(n),n=1,num)
!    8           format(4(1x,1PG14.6))
!            else
!                write(lu,9) type,num
!    9           format('# ',A32,' data:'/ &
!                 'object 3 class array type float rank 0',1x, &
!                 ' items ',I8,' ieee data 0')
!            endif

!            write(lu,10)
!   10       format('attribute "dep" string "positions"'// &
!             'object "irregular positions irregular',1x, &
!             'connections" class field'/'component "positions"',1x, &
!             'value 1'/'component "connections" value 2'/'component',1x, &
!             '"data" value 3'/'end')
call closef(lu)
if (.not.ascout) then
ierr = wrdxFloat (file,'a',XPT,num)
if (ierr.ne.0) ierr = EXEC
endif
 else
if (npbased) then
    do i=1,nnp
    n = monnp(i)
    if (iflag.eq.0) then
        Te = T(n)
    else
        Te = T(n)- TLST(n)
    endif
    if (VERBOSE) then
                if (iflag.eq.0) write(luMon,11,err=1100, &
           iostat=ios) n,Te
11                   format(1x,'temperature at node ',I8, &
                 ' = ',1PG15.6)
                if (iflag.eq.1) write(luMon,12,err=1100, &
           iostat=ios) n,Te
12                   format(1x,'temperature increment at node ',I8, &
                 ' = ',1PG15.6)
                if (iflag.eq.2) write(luMon,13,err=1100, &
           iostat=ios) irelt,n,Te
13                   format(1x,'temperature change since t=',i5,1x, &
                 'at node ',I8,' = ',1PG15.6)
            else
                write(luMon,*,err=1100,iostat=ios) Te
            endif
enddo
else
do i=1,nnel
    n = monel(i)
            call LTEMP (IEN(1,n),T,TL,NEN)
    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) &
           call LTEMPX (LMTX(1,n),IEN(1,n),TX,TL,NEN)
            do j=1,NEN
                XPT(j) = REAL(TL(j))
    enddo
            if (iflag.ne.0) then
                call LTEMP (IEN(1,n),TLST,TL,NEN)
    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) &
           call LTEMPX (LMTX(1,n),IEN(1,n),TLSTX,TL,NEN)
                do j=1,NEN
                    XPT(j) = XPT(j)-REAL(TL(j))
    enddo
    endif
enddo
if (VERBOSE) then
            if (iflag.eq.0) write(luMon,14,err=1100,iostat=ios) &
           n,(j,XPT(j),j=1,NEN)
14               format(1X,'Temperatures in element ',I8,':', &
             4(/1X,I1,3X,1PG14.5),4(/1X,I1,3X,1PG14.5))
            if (iflag.eq.1) write(luMon,15,err=1100,iostat=ios) &
           n,(j,XPT(j),j=1,NEN)
15               format(1X,'Temperature increments in element ',I8, &
             ':',4(/1X,I1,3X,1PG14.5),4(/1X,I1,3X,1PG14.5))
            if (iflag.eq.2) write(luMon,16,err=1100,iostat=ios) &
           irelt,n,(j,XPT(j),j=1,NEN)
16               format(1X,'Temperature change since t=',i5,1x, &
             'in element ',I8,':',4(/1X,I1,3X,1PG14.5), &
             4(/1X,I1,3X,1PG14.5))
else
    write(luMon,*,err=1100,iostat=ios) (XPT(j),j=1,NEN)
endif
endif
 endif
!
 return
!
1000    write(stderr,17) ios,'"',file(fnblnk(file):lnblk(file)),'"'
17    format(1x,'te: write error ',I3,' on file ',80A)
 call prerr('te')
 ierr = EXEC
 return
!
1100    write(stderr,18) ios,luMon
18    format(1x,'te: write error ',I3,' on monitor logical unit ',I3)
 call prerr('te')
 ierr = EXEC
 return
 end
!-------------------------------------------------------------------------------
subroutine wrconn (file,IEN,IPT,ierr)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: verbose, stderr, lumon, ascout, append
use monitorVars,     only: elmon, fmon, nnel, monel
use plotControl,     only: npbased
use errors,          only: ok, exec        

!
! writes nodal connectivity arrays to a DX file
!
implicit none

!-pass
 integer ierr
 character(len=*) file
integer :: IEN, ipt(nen,numel)
 dimension IEN(NEN,*)
!-locl
logical ODD
integer wrdxInt,lu,i,j,n, nn
integer, external :: lnblk,fnblnk,ifopen,lueof, nextlu
integer :: LIEN(nen)
integer :: ios
integer :: nipt
!-init
 ierr = OK
!
 if (NUMEL.le.0) return
 if (fmon) goto 300
!
 if (.not.elmon) then
!        Check if sufficient work space

if (npbased) then
!            Reorder (TECTON) ordering to DataExplorer model
    do i=1,NUMEL
        call TtoDX(IEN(:,i),IPT(:,i))
enddo
else
!        Element model: local coordinates
do i=1,NUMEL
    do j=1,NEN
    LIEN(j) = NEN*(i-1)+j
    enddo
    call TtoDX(LIEN,IPT(:,i))
enddo
endif
if (ascout) then
    lu = nextlu(0)
    ierr = ifopen(lu,file,'unknown')
if (ierr.ne.0) then
            ierr = EXEC
            goto 200
        endif
if (APPEND) then
            ierr = lueof(lu)
            if (ierr.ne.0) then
                ierr = EXEC
                goto 200
            endif
        endif
    do j=1,NUMEL
    write(lu,2,err=100,iostat=ios) (IPT(i,j),i=1,NEN)
 2            format(8(1X,I6))
enddo
goto 200
!
100        write(stderr,3) ios,'"', &
       file(fnblnk(file):lnblk(file)),'"'
 3        format(1x,'conn: write error ',I3,' on file ',256A)
call prerr('conn')
ierr = EXEC
!
200            call closef(lu)
else
if (APPEND) then
        ierr = wrdxInt(file,'a',IPT,NEN*NUMEL)
else
        ierr = wrdxInt(file,'w',IPT,NEN*NUMEL)
endif
    if (ierr.ne.0) ierr = EXEC
endif
 else
do i=1,nnel
        n = monel(i)
        if (VERBOSE) then
            write(luMon,4,err=700,iostat=ios) n, &
             (IEN(j,n),j=1,NEN)
 4               format(1x,'element ',I8,':',8(1X,I8))
        else
            write(luMon,5,err=700,iostat=ios) &
             (IEN(j,n),j=1,NEN)
 5               format(8(1X,I8))
        endif
enddo
endif
 goto 600
!
! NEXT PART NOT ADJUSTED FOR TETS
!    Face monitor data
!    In order to pack data, store two elements per loop
300    do nn=1,nnel/2
n = 2*nn - 1
IPT(1,nn) = (n-1)*4
IPT(2,nn) = (n-1)*4+1
IPT(3,nn) = (n-1)*4+3
IPT(4,nn) = (n-1)*4+2
n = 2*nn
IPT(5,nn) = (n-1)*4
IPT(6,nn) = (n-1)*4+1
IPT(7,nn) = (n-1)*4+3
IPT(8,nn) = (n-1)*4+2
 enddo

 nipt = (nnel/2)*8

 if (ODD(nnel)) then
nipt = nipt + 4
IPT(1,nnel/2+1) = (nnel-1)*4
IPT(2,nnel/2+1) = (nnel-1)*4+1
IPT(3,nnel/2+1) = (nnel-1)*4+3
IPT(4,nnel/2+1) = (nnel-1)*4+2
 endif

 if (ascout) then
lu = nextlu(10)
ierr = ifopen(lu,file,'unknown')
rewind(lu)
if (ierr.ne.0) then
        ierr = EXEC
        goto 500
    endif
if (APPEND) then
        ierr = lueof(lu)
        if (ierr.ne.0) then
            ierr = EXEC
            goto 500
        endif
    endif
do j=1,nnel/2
write(lu,2,err=400,iostat=ios) (IPT(i,j),i=1,8)
enddo
if (ODD(nnel)) then
j=nnel/2+1
write(lu,2,err=400,iostat=ios) (IPT(i,j),i=1,4)
endif
goto 500
!
400        write(stderr,3) ios,'"',file(fnblnk(file):lnblk(file)),'"'
call prerr('conn')
ierr = EXEC
!
500        call closef(lu)
 else
if (APPEND) then
ierr = wrdxInt(file,'a',IPT,nipt)
else
ierr = wrdxInt(file,'w',IPT,nipt)
endif
if (ierr.ne.0) ierr = EXEC
 endif
!
600    return
!
700    write(stderr,6) ios,luMon
 6    format(1x,'conn: write err ',I3,' on monitor logical unit ',I3)
 ierr = EXEC
 return
 end
!-------------------------------------------------------------------------------
 subroutine TtoDX (IEN,LIEN)
use modeltopology,   only: NSD,NDOF,NSTR,NEN
!
! permutes TECTON node order IEN to DX order LIEN
!
!                   T E C T O N                         D X
!
!                  8             7                 7             8
!                  o-------------o                 o-------------o
!                 /|            /|                /|            /|
!                / |           / |               / |           / |
!               /  |          /  |              /  |          /  |
!              /   |         /   |             /   |         /   |
!            5/    |       6/    |           3/    |       4/    |
!            o-------------o     |           o-------------o     |
!            |     |       |     |           |     |       |     |
!            |     o-------|-----o           |     o-------|-----o
!            |    /4       |    /3           |    /5       |    /6
!            |   /         |   /             |   /         |   /
!            |  /          |  /              |  /          |  /
!            | /           | /               | /           | /
!            |/            |/                |/            |/
!            o-------------o                 o-------------o
!            1             2                 1             2
!
!
 implicit none
!-pass
 integer IEN,LIEN
 dimension IEN(NEN),LIEN(NEN)
!
 if (NEN.eq.8) then
!        Bricks and wedges
LIEN(1) = IEN(1) - 1
LIEN(2) = IEN(2) - 1
LIEN(3) = IEN(5) - 1
LIEN(4) = IEN(6) - 1
LIEN(5) = IEN(4) - 1
LIEN(6) = IEN(3) - 1
LIEN(7) = IEN(8) - 1
LIEN(8) = IEN(7) - 1
 else
!        Tetrahedra
LIEN(1) = IEN(1) - 1
LIEN(2) = IEN(2) - 1
LIEN(3) = IEN(3) - 1
LIEN(4) = IEN(4) - 1
 endif
!
return
end
!-------------------------------------------------------------------------------
subroutine wrdisp (file,IEN,D,LMF,TFAULT,DXE,XPT,ierr)

! output displacements or velocities, depending on the arguments and flag iDisp.
! As a consequence, some of the arguments cannot be replaced by direct references to the modeldata structure in this subroutine.
! XPT is called with B

USE ALGEBRA
USE MODELDEFINITION
use meshdatamodule,  only: meshdatactx
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: modeldatactx
#endif
use modeltopology,   only: NDOF,NEN
use iomodule,        only: verbose, stderr, lumon, ifelm, icoord, crdfil
use constants,       only: useDXE,useDELDXE
use monitorVars,     only: elmon, npmon, nnel, nnp, monnp, monel
use plotControl,     only: iDisp, npbased
use errors,          only: exec

!
! Writes displacement or velocity vectors to a DX file.
! IDISP:first digit     1=displacement
!                       2=velocity
!      second digit     0=total
!                       1=average
!                       2=slippery
!                       3=faulted
!
implicit none
!-pass

#ifdef EXPERIMENTAL_ReducedSlipAlloc
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

character(len=*) file
real XPT(NDOF,*)
integer          :: ierr
integer          :: ien, lmf
double precision :: d, tfault
dimension IEN(NEN,*),D(NDOF,*),LMF(NDOF,NEN,*), TFAULT(NDOF,*)

!-locl
#include "lenchar.i"
integer           :: i,j,lu,n,num,ivelo,isort
integer           :: io,iel,il
character(len=26) :: type
character(len=6)  :: ftype
logical monitor
double precision  :: DL(3,8), XL(3,8)
integer, external :: lnblk, ifopen, nextlu
real, allocatable :: XLT(:,:)
integer           :: slipperySequenceNr, slipperytype
save type

!-init
if (npbased) then
    monitor = npmon
else
    monitor = elmon
endif

ivelo = INT(iDisp/10)
isort = MOD(iDisp,10)

if (isort.eq.0) then
    type = 'total '
    i = 7
else if (isort.eq.1) then
    type = 'average '
    i = 9
else if (isort.eq.2) then
    type = 'slippery '
    i = 10
else if (isort.eq.3) then
    type = 'faulted '
    i = 9
endif

if (ivelo.eq.1) then
    type(i:i+11) = 'displacement'
    SlipperyType = UseDXE
else if (ivelo.eq.2) then
    type(i:i+7) = 'velocity'
    SlipperyType = UseDELDXE
endif

! monitor specific elements or nodes?
if ((.not.npmon) .and. (.not.elmon)) then

    ! Write file header
    lu = nextlu(0)
    ierr = ifopen(lu,file,'unknown')
    if (ierr.ne.0) then
        ierr = EXEC
        return
    endif
    if (icoord.eq.1) then
        ftype = 'binary'
    else
        ftype = 'ascii '
    endif
    if (npbased) then
        num=meshdatactx%nvglobal
    else
        num=meshdatactx%neglobal*NEN
    endif
    if (ifelm.eq.1) then
        ftype = 'binary'
    else
        ftype = 'ascii '
    endif

    if (npbased) then
        ! remember: npbased means that there are no slippery or split nodes
        do n=1,meshdatactx%nvglobal
            if (isort.eq.0 .or. isort.eq.1) then
                XPT(1,n) = REAL(D(1,n))
                XPT(2,n) = REAL(D(2,n))
                XPT(3,n) = REAL(D(3,n))
            else
                XPT(1,n) = 0.
                XPT(2,n) = 0.
                XPT(3,n) = 0.
            endif
        enddo
        num = meshdatactx%nvglobal
    else 
        allocate(XLT(NDOF,NEN*meshdatactx%neglobal))
        num = 0
        do n=1,meshdatactx%neglobal
            if (isort.eq.0. .or. isort.eq.1) then
                call LDISP (DL,D,IEN(1,n),NDOF,NEN)
            else
                call CLEAR (DL,NDOF*NEN,"DL")
            endif
            if (NUMSLP.gt.0.and. isort.eq.0 .or. isort.eq.2) then
                call ADDSNE (DL,n,NDOF,NEN,slipperytype)
            endif
            if (isort.eq.0 .or. isort.eq.3) then
                call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,108)
            endif
            call LCOORD (meshdatactx%X,XL,IEN(1,n))
            do i=1,NEN
                num = num + 1
                ! array with displacements of nodes per element
                XPT(1,num) = REAL(DL(1,i))
                XPT(2,num) = REAL(DL(2,i))
                XPT(3,num) = REAL(DL(3,i))
                ! array with original node coordinates per element
                XLT(1,num) = XL(1,i)
                XLT(2,num) = XL(2,i)
                XLT(3,num) = XL(3,i)
            enddo
        enddo
    endif

    if (npbased) then
        do i=1,meshdatactx%nvglobal
            write(lu,*,err=1000,iostat=io) modeldatactx%D(1,i),modeldatactx%D(2,i),modeldatactx%D(3,i)
        enddo
    else
        do i=1,meshdatactx%neglobal*NEN
            write(lu,*,err=1000,iostat=io) XLT(1,i),XLT(2,i),XLT(3,i),XPT(1,i), XPT(2,i), XPT(3,i)
        enddo
        deallocate(XLT)
    endif
    call closef(lu)

else

    !   Monitoring nodes or elements
    if (npbased) then
        ! remember: npbased means that there are no slippery or split nodes
        if (npmon) then
            do i=1,nnp
                n = monnp(i)
                iel = 1 ! find the first element in which this node appears
                loop_iel: do while (iel.le.meshdatactx%neglobal) 
                    il = 1 ! identify element-local node number in the element
                    do while (il.le.NEN) 
                        if (IEN(il,iel).eq.n) exit loop_iel
                        il = il + 1
                    enddo
                    iel = iel + 1
                enddo loop_iel
                if (isort.eq.0. .or. isort.eq.1) then
                    call LDISP (DL,D,IEN(1,iel),NDOF,NEN)
                else
                    call CLEAR (DL,NDOF*NEN,"DL")
                endif
                if (VERBOSE) then
                    write(luMon,12,err=1100,iostat=io,advance='no') type(1:lnblk(type))
12                  format(1x,A)
                    write(luMon,13,err=1100,iostat=io) n,(DL(j,il),j=1,NDOF)
13                  format(' of nodal point ',I12,': ',3(1PG14.5,1X))
                else
                    write(luMon,14,err=1100,iostat=io) (DL(j,il),j=1,NDOF)
14                  format(3(1x,1PG14.6))
                endif
            enddo
        endif
        if (elmon) then
            do i=1,nnel
                iel = monel(i)
                if (isort.eq.0. .or. isort.eq.1) then
                    call LDISP (DL,D,IEN(1,iel),NDOF,NEN)
                else
                    call CLEAR(DL,NDOF*NEN,"DL")
                endif
                if (VERBOSE) then
                    write(luMon,12,err=1100,iostat=io,advance='no') type(1:lnblk(type))
                    write(luMon,15,err=1100,iostat=io) iel,((DL(j,il),j=1,NDOF),il=1,NEN)
15                  format(' of element ',I12,':',8(/1PG14.5,1X,1PG14.5,1X,1PG14.5))
                else
                    write(luMon,16,err=1100,iostat=io) ((DL(j,il),j=1,NDOF),il=1,NEN)
16                  format(8(1PG14.5,1X,1PG14.5,1X,1PG14.5/))
                endif
            enddo
        endif
    else ! there are slippery or split nodes
        if (npmon) then
            ! monitor displacement of a node as it appears in all elements
            do i=1,nnp
                n = monnp(i)
                do iel = 1,meshdatactx%neglobal
                    il = 1 ! identify element-local node number in the element
                    do while (il.le.NEN) 
                        if (IEN(il,iel).eq.n) exit
                        il = il + 1
                    enddo
                    if (il.le.NEN) then
                        if (isort.eq.0. .or. isort.eq.1) then
                            call LDISP (DL,D,IEN(1,iel),NDOF,NEN)
                        else
                            call CLEAR (DL,NDOF*NEN,"DL")
                        endif
                        if (isort.eq.0 .or. isort.eq.2) then
                            call ADDSNE (DL,iel,NDOF,NEN,slipperytype)
                        endif
                        if (isort.eq.0 .or. isort.eq.3) then
                            call AddFaultDisplacement (DL,iel,TFAULT,NDOF,NEN,109)
                        endif
                        if (VERBOSE) then
                            write(luMon,12,err=1100,iostat=io,advance='no') type(1:lnblk(type))
                            write(luMon,17,err=1100,iostat=io) n,iel,(DL(j,il),j=1,NDOF)
17                          format(' of nodal point ',I12,' in element ',I12,': ',3(1PG14.5,1X))
                        else
                            write(luMon,14,err=1100,iostat=io) (DL(j,il),j=1,NDOF)
                        endif
                    endif
                enddo
            enddo
        endif
        if (elmon) then
            if (VERBOSE) write(luMon,*) 'elno ((DL(i,j),i=1,NDOF),j=1,NEN)'
            do i=1,nnel
                iel = monel(i)
                if (isort.eq.0. .or. isort.eq.1) then
                    call LDISP (DL,D,IEN(1,iel),NDOF,NEN)
                else
                    call CLEAR(DL,NDOF*NEN,"DL")
                endif
                if (isort.eq.0 .or. isort.eq.2) then
                    call ADDSNE (DL,iel,NDOF,NEN,slipperytype)
                endif
                if (isort.eq.0 .or. isort.eq.3) then
                    call AddFaultDisplacement (DL,iel,TFAULT,NDOF,NEN,109)
                endif
                if (VERBOSE) then
                    write(luMon,12,err=1100,iostat=io,advance='no') type(1:lnblk(type))
                    write(luMon,15,err=1100,iostat=io) iel,((DL(j,il),j=1,NDOF),il=1,NEN)
                else
                    write(luMon,16,err=1100,iostat=io) ((DL(j,il),j=1,NDOF),il=1,NEN)
                endif
            enddo
        endif
    endif

endif ! monitor

return

1000 write(stderr,18) io,'"',file(1:lnblk(file)),'"'
18   format(1x,'wrdisp: write error ',I3,' on file ',80A)
call prerr('wrdisp')
ierr = EXEC
return

1100 write(stderr,19) io,luMon
19   format(1x,'wrdisp: write err ',I3,' on monitor logical unit ', &
     I3)
call prerr('wrdisp')
ierr = EXEC
return

end

!-----------------------------------------------------------------------
subroutine tempc (IEN,X,D,LMF,TFAULT,DXE,T,LMTX,TX,XPT,ierr)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDEFINITION, only: nElemsWithSlipperyNodes, NUMEL
#else
USE MODELDEFINITION, only: NUMEL
#endif
use modeltopology,   only: NSD,NDOF,NEN
use constants,       only: one, useDXE
use meshdatamodule,  only: meshdatactx
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel
use errors,          only: exec, ok
!
! Routine loads element center temperatures into XPT
!
implicit none
!-pass
real XPT(*)
integer          :: ien,lmf,lmtx
double precision :: x,d,tfault,t,tx
dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
          TFAULT(NDOF,*),T(*),TX(*),LMTX(NEN,*)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif




!-locl
integer i
double precision :: XL(24),DL(24)
integer :: ierr, n
double precision :: Tc
double precision, external :: Tcentr
!-init
 ierr = OK
!
 if (.not.elmon) then

do n=1,NUMEL
call LCOORD (X,XL,IEN(1,n))
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,110)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)
        call REZONE (XL,DL,ONE)
        Tc = Tcentr (T,TX,XL,IEN(1,n),LMTX(1,n),ierr)
if (ierr.ne.0) then
    write(stderr,2) n
 2            format(1x,'tempc: error in element ',I8)
    ierr = EXEC
    return
endif
XPT(n) = REAL(Tc)
enddo
 else
do i=1,nnel
        n = monel(i)
call LCOORD (X,XL,IEN(1,n))
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,111)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)
        call REZONE (XL,DL,ONE)
        Tc = Tcentr (T,TX,XL,IEN(1,n),LMTX(1,n),ierr)
if (ierr.ne.0) then
    write(stderr,2) n
    ierr = EXEC
    return
endif
XPT(i) = REAL(Tc)
enddo
endif
!
return
end
!-----------------------------------------------------------------------
subroutine matc (MAT,XPT)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel

!
! Routine loads element center material data into XPT
!
implicit none

!-pass
real XPT(*)
integer MAT(*),ierr
!-locl
integer :: n,i
!-init
!
if (.not.elmon) then
    do n=1,NUMEL
        XPT(n) = REAL(MAT(n))
    enddo
else
    do i=1,nnel
        n = monel(i)
        XPT(i) = REAL(MAT(n))
    enddo
endif
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine wrnbc (file,IBOND,BOND,SKEW,XPT,ibcflg,ierr)

USE ALGEBRA
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use meshdatamodule,  only: meshdatactx
use iomodule,        only: verbose, stderr, lumon, icoord, ascout, &
                           crdfil
use monitorVars,     only: npmon, nnp, monnp
use plotControl,     only: hatangle,hatsize,sysi,iarcont, &
                           readsy,sbhor, sysc
use errors,          only: exec
!
! Writes nodal boundary condition arrows.
! ibcflg = 1: write displacement bc's only
! ibcflg = 2: write velocity bc's only
! ibcflg = 3: write force bc's only
! ibcflg = 4: write initial displacement bc's only
! ibcflg = 5: write strain rate bc's only
! ibcflg > 5: all bcs will be written.
!


implicit none
!-pass
 character(len=*) file

integer :: ibcflg, ibond
real :: XPT(NDOF,meshdatactx%nvglobal)
double precision :: bond, skew
 dimension IBOND(NDOF,meshdatactx%nvglobal),BOND(NDOF,meshdatactx%nvglobal),SKEW(2,*)
integer :: ierr




!-locl

integer :: i, ios, lu

double precision, parameter :: ZERO=0d0, ONE=1d0
#include "lenchar.i"
 logical ALL,isbc
 character(len=4) KEY(0:5)
 save KEY
integer, external :: lnblk,wrdxFloat,nextlu,ifopen
double precision :: v(3),w(3),rot(3,3)
 data KEY/'FREE','DISP','VELO','INIT','FORC','RATE'/
integer :: n, nn


 ALL = (ibcflg.lt.1.or.ibcflg.gt.5)

 if (.not.npmon) then
call CLEAR(XPT,meshdatactx%nvglobal*NDOF,"XPT")
do n=1,meshdatactx%nvglobal
do i=1,NDOF
    v(i) = ZERO
    if ((ALL.and.IBOND(i,n).ne.0).or. &
           IBOND(i,n).eq.ibcflg) then
    if (sysc) then
        v(i) = BOND(i,n)*sysi
    else
        v(i) = DSIGN(ONE,BOND(i,n))*sysi
    endif
    endif
enddo
if (SKEW(1,n).ne.ZERO .or. SKEW(2,n).ne.ZERO) then
    call FORMRT(SKEW(:,n),rot,NDOF)
    call EQUATE(w,v,NDOF)
    call VMPRD (rot,w,v,NDOF,3)
endif
do i=1,NDOF
    XPT(i,n) = REAL(v(i))
enddo
enddo
!
!        Write header data and, in case of ASCII output, data itself
lu = nextlu(0)
ierr = ifopen(lu,file,'unknown')
if (ierr.ne.0) then
ierr = EXEC
return
endif
if (icoord.eq.1) then
write(lu,1) NDOF,meshdatactx%nvglobal,'binary', crdfil(1:lnblk(crdfil))
else
write(lu,1) NDOF,meshdatactx%nvglobal,'ascii ',crdfil(1:lnblk(crdfil))
endif
 1       format('# nodal points:'/'object 1 class array type',1x, &
     'float rank 1 shape ',I1,' items ',I8,1X,A6/'data file ', &
     256A)
if (ascout) then
write(lu,2) meshdatactx%nvglobal
 2           format('object 2 class array type float rank 1 shape 3', &
         1x,'items ',I8,' data follows')
        write(lu,3,err=1000,iostat=ios) &
         ((XPT(i,n),i=1,NDOF),n=1,meshdatactx%nvglobal)
 3        format(3(1x,1PG14.6))
else
write(lu,4) meshdatactx%nvglobal
 4           format('object 2 class array type float rank 1 shape 3', &
         1x,'items ',i8,' ieee data 0')
endif

write(lu,5)
 5       format('attribute "dep" string "positions"'// &
     'object "BcArrows" class field'/ &
     'component "positions"',1x,'value 1'/'component',1x, &
     '"data" value 2'/'end')
    call closef(lu)
!        write binary data to file
if (.not.ascout) then
ierr = wrdxFloat (file,'a',XPT,meshdatactx%nvglobal*NDOF)
if (ierr.ne.0) ierr = EXEC
endif
 else
do nn=1,nnp
n = monnp(nn)
isbc = .false.
do i=1,NDOF
    if (ALL.or.IBOND(i,n).eq.ibcflg) isbc = .true.
enddo
if (isbc) then
    if (VERBOSE) then
    write(luMon,6,err=1100,iostat=ios) n,(BOND(i,n), &
           KEY(IBOND(i,n)),i=1,NDOF)
 6                   format(1x,'Nodal bc on node ',I6,':', &
                 3(2X,1PG12.4,1X,A4))
    else
    write(luMon,*,err=1100,iostat=ios) (BOND(i,n), &
           i=1,NDOF)
    endif
endif
enddo
 endif
!
return
!
1000    write(stderr,9) ios,'"',file(1:lnblk(file)),'"'
 9   format(1x,'wrnbc: write error ',I3,' on file ',80A)
call prerr('wrnbc')
ierr = EXEC
return
!
1100    write(stderr,10) ios,luMon
10   format(1x,'wrnbc: write err ',I3,' on monitor logical unit ', &
    I3)
call prerr('wrnbc')
ierr = EXEC
return
end
!-----------------------------------------------------------------------
subroutine coulomb (STN,FRMHU,XPT,IERR)

USE ALGEBRA
USE MODELDEFINITION
USE CONSTANTS,       only: PI, ONE, TWO, DEG2RAD
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: verbose, stderr, lumon
use monitorVars,     only: elmon, monel, nnel
use errors,          only: exec, ok
!
! computes direction perpendicular to optimal fault plane. The magnitude of the
! resulting vector is the resolved shear stress
!
implicit none

!-pass
integer IERR
real XPT
double precision STN,FRMHU
dimension STN(NSTR,*),XPT(3,*)




!-locl
integer nn,n
double precision alpha,val,dir,vec,f,R,phi,theta
dimension val(3),dir(3,3),vec(3)
!-init
IERR = OK

!    orientation of optimal plane
if (FRMHU.lt.1d-6) then
    alpha = PI*0.25d0
else
    alpha = 5D-1*ATAN(ONE/FRMHU)
endif

if (elmon) then
    do nn=1,nnel
        n = monel(nn)

!       compute eigenvalues and eigenvectors
        call eigen (STN(1,n),val,dir)
!       sort in descending order s1 >= s2 => s3
        call eigensort (val,dir)
!       determine polar coordinates of intermediate eigenvector
        call cartopolar (dir(1,2),dir(2,2),dir(3,2),R,phi,theta)
!       rotate s1-vector over angle alpha around s2-direction
        call arbitrot (dir,phi,theta,alpha,vec)

        f = (val(1)-val(3))/(TWO*SQRT(ONE+FRMHU*FRMHU))
        call BMULT (vec,3,f)
        call cartopolar (vec(1),vec(2),vec(3),R,phi,theta)
        if (VERBOSE) then
            write(luMon,10,err=1000) n,R,phi/DEG2RAD,theta/DEG2RAD
10               format(1x,'Optimal fault plane in element ', &
             I12,': STRESS = ',1PG14.5,' PHI = ',1PG14.5, &
             ' THETA = ',1PG14.6)
            write(luMon,15,err=1000) n,vec(1),vec(2),vec(3)
15               format(1x,'Optimal fault plane normal in element ', &
             I12,': ',3(1PG14.5,1X))
        else
            write(luMon,20,err=1000) R,phi/DEG2RAD,theta/DEG2RAD
            write(luMon,20,err=1000) vec(1),vec(2),vec(3)
20            format(3(1x,1PG14.6))
        endif
    enddo
else
    do n=1,NUMEL
        call eigen (STN(1,n),val,dir)
        call eigensort (val,dir)
        call cartopolar (dir(1,2),dir(2,2),dir(3,2),R,phi,theta)
        call arbitrot (dir,phi,theta,alpha,vec)

        f = (val(1)-val(3))/(TWO*SQRT(ONE+FRMHU*FRMHU))

        call BMULT (vec,3,f)

        XPT(1,n) = REAL(vec(1))
        XPT(2,n) = REAL(vec(2))
        XPT(3,n) = REAL(vec(3))
    enddo
endif
!
return
!
1000    write(stderr,30) luMon
30   format(1x,'coulomb: write error on monitor logical unit ', &
    I3)
call prerr('coulomb')
ierr = EXEC

return
end
!-------------------------------------------------------------------------------
 subroutine arbitrot (X,PHI,THETA,ANGLE,Y)

 USE ALGEBRA
!
! rotates vector X around axis (PHI,THETA) over an angle ANGLE
!
 implicit none
 double precision ZERO,ONE
 parameter (ZERO=0d0, ONE=1d0)
!-pass
 double precision X(3),PHI,THETA,ANGLE,Y(3)
!-locl
 double precision R1(3,3),R2(3,3),RT(3,3)
!
 R1(1,1) =  COS(PHI)
 R1(2,1) = -SIN(PHI)
 R1(3,1) =  ZERO
 R1(1,2) = -R1(2,1)
 R1(2,2) =  R1(1,1)
 R1(3,2) =  ZERO
 R1(1,3) =  ZERO
 R1(2,3) =  ZERO
 R1(3,3) =  ONE
!
 R2(1,1) =  SIN(THETA)
 R2(2,1) =  ZERO
 R2(3,1) = -COS(THETA)
 R2(1,2) =  ZERO
 R2(2,2) =  ONE
 R2(3,2) =  ZERO
 R2(1,3) = -R2(3,1)
 R2(2,3) =  ZERO
 R2(3,3) =  R2(1,1)
 call MATPRD (R2,R1,RT,3,3,3)
 call EQUATE (R1,RT,9)

 R2(1,1) =  ONE
 R2(2,1) =  ZERO
 R2(3,1) =  ZERO
 R2(1,2) =  ZERO
 R2(2,2) =  COS(ANGLE)
 R2(3,2) =  SIN(ANGLE)
 R2(1,3) =  ZERO
 R2(2,3) = -R2(3,2)
 R2(3,3) =  R2(2,2)
 call MATPRD (R2,R1,RT,3,3,3)
 call EQUATE (R1,RT,9)
!
 R2(1,1) =  SIN(THETA)
 R2(2,1) =  ZERO
 R2(3,1) =  COS(THETA)
 R2(1,2) =  ZERO
 R2(2,2) =  ONE
 R2(3,2) =  ZERO
 R2(1,3) = -R2(3,1)
 R2(2,3) =  ZERO
 R2(3,3) =  R2(1,1)
 call MATPRD (R2,R1,RT,3,3,3)
 call EQUATE (R1,RT,9)
!
 R2(1,1) =  COS(PHI)
 R2(2,1) =  SIN(PHI)
 R2(3,1) =  ZERO
 R2(1,2) = -R2(2,1)
 R2(2,2) =  R2(1,1)
 R2(3,2) =  ZERO
 R2(1,3) =  ZERO
 R2(2,3) =  ZERO
 R2(3,3) =  ONE
 call MATPRD (R2,R1,RT,3,3,3)
 call VMPRD (RT,X,Y,3,3)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine eigensort (D,V)
!
! sorts eigenvalues and corresponding vectors in descending order
!
 implicit none
!-pass
 double precision D(3),V(3,3)
!-locl
 integer imin,imax,i,imid
 double precision emin,emax,dl(3),vl(3,3)
!
 imin = 1
 imax = 1
 emin = D(imin)
 emax = D(imax)
 do i=2,3
if (emin.gt.D(i)) then
imin = i
emin = D(imin)
endif
if (emax.lt.D(i)) then
imax = i
emax = D(imax)
endif
 enddo

 if (imin.eq.imax) then
return
 endif
 imid = 6 - imin - imax
!
 dl(1) = D(imax)
 dl(2) = D(imid)
 dl(3) = D(imin)

 do i=1,3
vl(i,1) = V(i,imax)
vl(i,2) = V(i,imid)
vl(i,3) = V(i,imin)
 enddo

 do i=1,3
D(i) = dl(i)
V(i,1) = vl(i,1)
V(i,2) = vl(i,2)
V(i,3) = vl(i,3)
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine eigen (STN,D,V)

use constants, only : dfmin
!
! Computes all eigenvalues and eigenvectors of Cauchy stress tensor STN
! D returns the eigenvalues of A. V is a 3x3 matrix whose columns contain
! on output the normalized eigenvectors of STN.
!
implicit none
integer, parameter :: ZERO=0d0, ONE=1d0
integer, parameter :: NP=3,N=NP,NSTR=6
!-pass
double precision :: STN(NSTR),D(NP),V(NP,NP)
!-locl
character(len=1) :: answer
double precision :: b(NP),z(NP),a(NP,NP)

integer :: ip, iq, nrot, i, j
double precision :: c, s, t, sm, g, h, tau, theta, tresh
!
!     Initialize to the identity matrix and copy STN into a
 a(1,1) = STN(1)
 a(1,2) = STN(4)
 a(1,3) = STN(5)
 a(2,2) = STN(2)
 a(2,3) = STN(6)
 a(3,3) = STN(3)
 a(2,1) = a(1,2)
 a(3,1) = a(1,3)
 a(3,2) = a(2,3)

 do ip=1,N
do iq=1,N
V(ip,iq) = ZERO
enddo
V(ip,ip) = ONE
 enddo
!    Initialize B and D to the diagonal of A
 do ip=1,N
b(ip) = a(ip,ip)
D(ip) = b(ip)
!        Vector Z  will accumulate terms of the form ta(p,q) as in
!        equation (11.1.14)
z(ip) = ZERO
 enddo
 nrot = 0
 do i=1,50
!        Sum off-diagonal elements
sm = ZERO
do ip=1,N-1
do iq=ip+1,N
    sm = sm + ABS(a(ip,iq))
enddo
enddo
!        Case that sum of off-diagonal elements zero branch to end.
if ( sm.le.dfmin ) goto 100
if ( i.lt.4 ) then
!        On the first three sweeps (see notes)
tresh = 0.2d0 * sm / (N*N)
else
tresh = ZERO
endif
do ip=1,N-1
do iq=ip+1,N
    g = 1d3*ABS( a(ip,iq) )
!            After four sweeps, skip the rotation if the off-
!            diagonal element is small.
    if ( (i.gt.4).and.(g.le.dfmin) ) then
    a(ip,iq) = ZERO
    else if ( ABS(a(ip,iq)).gt.tresh ) then
    h = D(iq) - D(ip)
    if ( g.le.dfmin ) then
        t = a(ip,iq) / h
    else
!                Equation (11.1.10)
        theta = 0.5d0 * h / a(ip,iq)
        t = ONE/(ABS(theta)+SQRT(ONE+theta*theta))
        if ( theta.lt.ZERO ) t = -t
    endif
    c        = ONE/SQRT(ONE+t*t)
    s        = t * c
    tau      = s/(ONE+c)
    h        = t * a(ip,iq)
    z(ip)    = z(ip) - h
    z(iq)    = z(iq) + h
    D(ip)    = D(ip) - h
    D(iq)    = D(iq) + h
    a(ip,iq) = ZERO
!            Case of rotations 1 <= j < p
    do j=1,ip-1
        g = a(j,ip)
        h = a(j,iq)
        a(j,ip) = g - s * (h + g * tau)
        a(j,iq) = h + s * (g - h * tau)
    enddo
!            Case of rotations p < j < q
    do j=ip+1,iq-1
        g = a(ip,j)
        h = a(j,iq)
        a(ip,j) = g - s * (h + g * tau)
        a(j,iq) = h + s * (g - h * tau)
    enddo
!            Case of rotations q < j <= N
    do j=iq+1,N
        g = a(ip,j)
        h = a(iq,j)
        a(ip,j) = g - s * (h + g * tau)
        a(iq,j) = h + s * (g - h * tau)
    enddo
    do j=1,n
        g = V(j,ip)
        h = V(j,iq)
        V(j,ip) = g - s * (h + g * tau)
        V(j,iq) = h + s * (g - h * tau)
    enddo
    nrot = nrot + 1
    endif
enddo
enddo
!        Update D with the sum of ta(p,q) and re-initialize Z.
do ip=1,N
b(ip) = b(ip) + z(ip)
D(ip) = b(ip)
z(ip) = ZERO
enddo
 enddo
 write(*,*) 'eigen: 50 iterations should never happen'
read(*,*) answer
!
100    return
 end
