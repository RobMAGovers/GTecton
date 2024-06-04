module ioMod
implicit none
character(len=256) :: infile,polyfile
integer            :: stderr
logical            :: VERBOSE
end module

module epsMod
implicit none
double precision   :: eps
end module

module wrappingMod
implicit none
logical            :: wrapping
end module



program main
!
! Sorts coordinates into an oriented polygon. Uses distance R from first
! point along the polygon. Sorts R and index table NDEX.
! Rob Govers 02-2011 version 1.0
! Wienand Drenth 24-02-2011 version 1.1
! Wienand Drenth 01-11-2011 version 1.2, reverted to Rob's implemention
!       in setdist routine
! Wienand Drenth 02-11-2011 version 1.3, major changes:
!   - correct algorithm to find the segment in which a point lies
!   - bug fix in check on closure of polygon
! Lukas van de Wiel 27-feb-2014 version 1.4
!   - fixed bug in changing error margin argument to number.
!   - added comment in determination of in_xy and in_r
! Lukas van de Wiel 27-jun-2022 version 1.5
!   - Removed common blocks to agree to the Fortran 2018 standard

use ioMod,       only: infile, verbose, polyfile, stderr
use wrappingMod, only: wrapping
use epsMod,      only: eps

implicit none
!-locl
integer MAXN
parameter (MAXN=1e7)
integer iflu,npoly,nx,i,ipc,numc

double precision, allocatable :: x(:,:),r(:)
integer,          allocatable :: ndex(:), np(:)
!-init
allocate(x(2,MAXN))
allocate(r(MAXN))
allocate(ndex(MAXN))
allocate(np(MAXN))

eps      = 1e-4
stderr   = iflu('stderr')
infile   = '<stdin>'
VERBOSE  = .false.
WRAPPING = .false.

!       Gets input file names from command line
call setinput (polyfile,infile)

!       Read polygon data
call rdpoly (x,npoly,MAXN,polyfile)
call testpoly (x,npoly)

!       Read coordinates
nx = npoly
call rdcoord (x,np,nx,MAXN,infile)
ipc = npoly+1
numc = nx-npoly

!       Initialize NDEX
do i=ipc,nx
    ndex(i) = i
enddo
!
if (numc.ge.2) then 
!           Compute R
    call setdist(x,nx,npoly,r)
    if (VERBOSE) then
        write(*,30)
   30           format(1x,'Coordinates and distances before sorting:')
        do i=1,nx
            write(*,40) x(1,i),x(2,i),r(i)
   40               format(1PE20.6,1X,1PE20.6,1X,1PE20.6)
        enddo
    endif
!
!           (Quick)sort index table NDEX based on R.
    call sort (numc,r(ipc),ndex(ipc))
    if (VERBOSE) then
        write(*,50)
   50           format(1x,'Coordinates and distances after 1st sort:')
        do i=ipc,nx
            write(*,40) x(1,ndex(i)),x(2,ndex(i)),r(i)
        enddo
    endif

    if (WRAPPING) then
!               Correct for wrapping and re-sort
!               r(npoly) is the maximum distance (i.e. number of intervals)
        call circular (r(ipc),r(npoly),numc)
        if (VERBOSE) then
            write(*,60)
   60   format(1x,'Coordinates and distances after wrapping:')
            do i=ipc,nx
                write(*,40) x(1,ndex(i)),x(2,ndex(i)),r(i)
            enddo
        endif

        call sort (numc,r(ipc),ndex(ipc))
        if (VERBOSE) then
            write(*,70)
   70   format(1x,'Coordinates and distances after 2nd sort:')
            do i=ipc,nx
                write(*,40) x(1,ndex(i)),x(2,ndex(i)),r(i)
            enddo
        endif
    endif
endif

do i=ipc,nx
    write(*,10) np(ndex(i)),x(1,ndex(i)),x(2,ndex(i))
   10       format(i8,1X,1PG16.6,1X,1PG16.6)
enddo
call exitp(0)

end

!-----------------------------------------------------------------------
subroutine circular (r,rmax,numc)

implicit none
!-pass
integer numc
double precision r(numc),rmax
!-locl
integer imd,i
double precision maxdist,dist

!       wrapped/circular?
maxdist = r(1)
do i=2,numc
    maxdist = MAX(maxdist,ABS(r(i)))
enddo

! check if point is in last segment or not 
! if maxdist > rmax-1 (and maxdist < rmax) than
! a point is in the last segment

if (maxdist.lt.ABS(rmax)-1d0) return

maxdist= ABS(r(2)-r(1))
imd = 1

do i=2,numc-1
    dist = ABS(r(i+1)-r(i))
    if (dist.gt.maxdist) then
        maxdist = dist
        imd = i
    endif
enddo

! adapt the distances by shifting them rmax so the gap is made
! as small as possible
do i=imd+1,numc
    r(i) = r(i) - rmax
enddo

return
end
!-----------------------------------------------------------------------
subroutine setdist (x,nx,npoly,r)

use ioMod,  only: verbose, stderr
use epsMod, only: eps

implicit none
!-pass
integer nx,npoly
double precision x(2,nx),r(nx)
!-locl
logical in_x,in_y
double precision in_xy, in_r
integer i,j
double precision HRT2,xx,yy,rx,ry,dx,dy
!-init
HRT2 = 1d0/SQRT(2d0)
!
!       Initialize R for polygon

do i=1,npoly
    r(i) = FLOAT(i-1)
!            write(*,*) x(1,i), x(2,i), r(i)
enddo

!       Compute R for additional points

do j=npoly+1,nx
    xx = x(1,j)
    yy = x(2,j)

!           Find polygon interval where (xx,yy) occurs
    i = 0
100         i = i + 1

    if (i.ge.npoly) then

        if (VERBOSE) then
            write(stderr,*) 'Where it went wrong'
            write(stderr,*) i-1, j - npoly, npoly
            write(stderr,*) x(1,i-1), x(2,i-1)
            write(stderr,*) x(1,i), x(2,i)
            write(stderr,*) in_xy, in_r, EPS
        endif

        write(stderr,10) xx,yy
   10           format(1x,'polysort: coordinate not in polygon: ',&
        2(1PG20.6))
        write(stderr,*) "Possible cause: too critical search parameter."
                      write(stderr,*)    "Possible solution:"
        write(stderr,*) "Loosen search criteria by adding argument '-e [number].'"
        write(stderr,*) "Increase the number to above the default of 1e-4."
        call exitp(1)
    endif

    if (ABS(xx-x(1,i)).le.EPS .and. &
        ABS(yy-x(2,i)).le.EPS) then
       goto 200
    endif

    if (ABS(xx-x(1,i+1)).le.EPS .and. &
        ABS(yy-x(2,i+1)).le.EPS) then
       goto 200
    endif

    if (.true.) then

    ! when the test point is not one of the segment's ends:
    ! compute in_xy such that the point
    ! x(i) + in_xy(x(i+1) - x(i)) is the point on the line
    ! x(i) -- x(i+1) with shortest distance to (xx,yy) (the test
    ! point),as follows:
    !
    !                        ^   * (xx,yy) 
    !                   in_r |   |
    !  x(i)                  v   |                x(i+1)
    !   *--------------------------------------------*
    !
    ! a  <----------------------->
    ! b  <------------------------------------------->
    !
    !      in_xy = a/b
    !

    in_xy = ((xx - x(1,i))*(x(1,i+1) - x(1,i)) + &
             (yy - x(2,i))*(x(2,i+1) - x(2,i))) / &
      ((x(1,i+1) - x(1,i))*(x(1,i+1) - x(1,i)) + &
       (x(2,i+1) - x(2,i))*(x(2,i+1) - x(2,i)))

    
    if (VERBOSE) then
        write(stderr,*) 'trying to match point to polygon segment',i
        write(stderr,*) 'from',x(1,i),x(2,i),&
                        'to',  x(1,i+1),x(2,i+1)
        write(stderr,*) 'point', xx,yy 
        write(stderr,*) 'giving in_xy', in_xy
    endif


    ! in_r is the square of the distance of the test point (xx,yy) to the
    ! point on the line x(i) -- x(i+1) closest to the test point

    in_r   = (xx - x(1,i) +in_xy*(x(1,i)-x(1,i+1)))*&
             (xx - x(1,i) +in_xy*(x(1,i)-x(1,i+1)))+&
             (yy - x(2,i) +in_xy*(x(2,i)-x(2,i+1)))*&
             (yy - x(2,i) +in_xy*(x(2,i)-x(2,i+1)))

    if (VERBOSE) then
        write(stderr,*) 'giving range', in_r
    endif


    if (SQRT(in_r).le.EPS .and. (in_xy.le.1d0 .and. in_xy.ge.0d0)) then
        ! the point is on the polygon within
        ! numerical precision
        goto 200    
    endif

else  

    in_x = ((xx-x(1,i))*(xx-x(1,i+1)).le.0d0)
    in_y = ((yy-x(2,i))*(yy-x(2,i+1)).le.0d0)

    if (in_x.and.in_y) goto 200

endif


goto 100

200         dx = x(1,i+1)-x(1,i)
    dy = x(2,i+1)-x(2,i)

    in_xy = ((xx - x(1,i))*(x(1,i+1) - x(1,i)) + &
            (yy - x(2,i))*(x(2,i+1) - x(2,i))) / &
      ((x(1,i+1) - x(1,i))*(x(1,i+1) - x(1,i)) + &
       (x(2,i+1) - x(2,i))*(x(2,i+1) - x(2,i)))

    if (ABS(dx).gt.EPS) then
        if (ABS(dy).gt.EPS) then
            rx = HRT2*(xx-x(1,i))/(x(1,i+1)-x(1,i))
            ry = HRT2*(yy-x(2,i))/(x(2,i+1)-x(2,i))
        else
            rx = (xx-x(1,i))/(x(1,i+1)-x(1,i))
            ry = 0d0
        endif
    else
        if (ABS(dy).gt.EPS) then
            rx = 0d0
            ry = (yy-x(2,i))/(x(2,i+1)-x(2,i))
        else
            rx = -1d0
            ry = -1d0
            write(stderr,20) i,i+1
   20               format(1x,'polysort: overlapping polygon ',&
             'coordinates ',I6,' and ',I6)
            call exitp(1)
        endif
    endif
    r(j) = r(i) + SQRT(rx*rx+ry*ry)

!            write(*,*) x(1,i),x(2,i), x(1,i+1),x(2,i+1)
!            write(*,*) i,j, xx,yy,SQRT(rx*rx+ry*ry), in_xy
!            write(*,*) '------- ---------- ---------'
enddo


return
end
!-----------------------------------------------------------------------
subroutine rdcoord (x,np,nx,MAXN,infile)

use ioMod,  only: stderr         

implicit none
!-pass
integer nx,MAXN,np(MAXN)
double precision x(2,MAXN)
character*(*) infile
!-locl
integer luin,nextlu,iflu,ios,ii
double precision xx,yy
external nextlu,iflu
!-init
if (infile.eq.'<stdin>') then
    luin = iflu('stdin')
else
    luin = nextlu(0)
    call openf(luin,infile,'old')
endif

100     read(luin,*,err=1000,end=200,iostat=ios) ii,xx,yy
nx = nx + 1

if (nx.gt.MAXN) then
    write(stderr,10)
   10       format(1x,'polysort: too many coordinates')
    call exitp(1)
endif

np(nx)  = ii
x(1,nx) = xx
x(2,nx) = yy

goto 100
!
200     if (infile.ne.'<stdin>') call closef(luin)

return
1000    write(stderr,1010) ios,nx+1
 1010   format(1x,'polysort: read error ',I5,' coordinate record ',I5/11x,'Expected node number, x- and y-coordinate')
call exitp(1)
end
!-----------------------------------------------------------------------
subroutine testpoly (x,npoly)

! validates whether the points and their sequence indeed form a polygon
! (no lines intersect)
use ioMod, only: stderr

implicit none
!-pass
integer npoly
double precision x(2,npoly)
!-locl
integer i,j,m
double precision x1,y1,x2,y2,x3,y3,x4,y4
logical cross

do i=1,npoly-2
    x1 = x(1,i)
    y1 = x(2,i)
    x2 = x(1,i+1)
    y2 = x(2,i+1)
    m = npoly-1

    if (i.eq.1) m=m-1

    do j=i+2,m
        x3 = x(1,j)
        y3 = x(2,j)
        x4 = x(1,j+1)
        y4 = x(2,j+1)

        if (cross (x1,y1,x2,y2,x3,y3,x4,y4)) then
            write(stderr,10) i,x1,y1,x2,y2,j,x3,y3,x4,y4
   10               format(1x,'polysort: self-intersecting polygon'/1x, &
             'segment ',I5,' with nodes',4(1x,1PG12.2),' cuts'/ &
             ' segment ',I5,' with nodes',4(1x,1PG12.2))
            call exitp(1)
        endif

    enddo
enddo

return

call exitp(1)
end
!-----------------------------------------------------------------------
subroutine rdpoly (x,npoly,MAXN,polyfile)

use ioMod,  only: stderr

implicit none
!-pass
integer npoly,MAXN
double precision x(2,MAXN)
character*(*) polyfile
!-locl
integer lu,nextlu,ios
double precision xx,yy,EPS
external nextlu
!-init
EPS = 2D-4
lu = nextlu(0)
npoly = 0

call openf(lu,polyfile,'old')
100     read(lu,*,err=1000,end=200,iostat=ios) xx,yy
npoly = npoly + 1

if (npoly.gt.MAXN) then
    write(stderr,10)
   10       format(1x,'polysort: too many polygon coordinates')
    call exitp(1)
endif

x(1,npoly) = xx
x(2,npoly) = yy
goto 100

200     call closef(lu)

!       check if first and last coordinates are equal
!        If that is not the case, the first point is added
!        to the list

if ( (ABS(x(1,1)-x(1,npoly)).gt.EPS) .or. &
     (ABS(x(2,1)-x(2,npoly)).gt.EPS) )  then
    npoly = npoly + 1
!           Duplicate first coordinate as last coordinate
    if (npoly.gt.MAXN) then
        write(stderr,10)
        call exitp(1)
    endif
    x(1,npoly) = x(1,1)
    x(2,npoly) = x(2,1)

endif

return

1000    write(stderr,1010) ios,npoly+1
 1010   format(1x,'polysort: read error ',I5,' polygon record ',I4)

call exitp(1)
end
!-----------------------------------------------------------------------
subroutine setinput (polyfile,infile)

use wrappingMod, only: wrapping
use ioMod,       only: stderr, verbose
use epsMod,      only: eps

! Gets input file names from command line

implicit none
!-pass
character*(*) polyfile,infile
!-locl
logical numeric
integer nargs,iargc,i,lnblk
character(len=256) USAGE,arg
double precision chreal
external lnblk,numeric,chreal
!-init

   
call HowToUse('-p polygon_file [-e eps] [-v] [-w] [filename]',USAGE)
nargs = command_argument_count()

if (nargs.lt.2) then
    write(stderr,10) USAGE(1:lnblk(USAGE))
   10       format(256a)
    call exitp(1)
endif

i = 1
100     if (i.gt.nargs) return

call get_command_argument(i, arg)

if (arg.eq.'-p') then
    i = i + 1
    if (i.gt.nargs) then
        write(stderr,10) USAGE(1:lnblk(USAGE))
        call exitp(1)
    endif
    call get_command_argument(i, polyfile)
else if (arg.eq.'-e') then
    i = i + 1
    if (i.gt.nargs) then
        write(stderr,10) USAGE(1:lnblk(USAGE))
        call exitp(1)
    endif
    call get_command_argument(i ,arg)
    if (numeric(arg)) then
        read(arg,*) EPS
!                EPS = chreal(arg)
        if (EPS.le.0.0) then
            write(stderr,20)
   20               format(1x,'polysort: epsilon should be > 0')
            call exitp(1)
        endif
    else
        write(stderr,10) USAGE(1:lnblk(USAGE))
        call exitp(1)
    endif
else if (arg.eq.'-v') then
    VERBOSE = .true.
else if (arg.eq.'-w') then
    WRAPPING = .true.
else
    if (infile.eq.'<stdin>') then
        infile=arg
    else
        write(stderr,10) USAGE(1:lnblk(USAGE))
        call exitp(1)
    endif
endif

i = i + 1
goto 100

end
!-----------------------------------------------------------------------
subroutine sort (n,distance,ndex)
!
! Quick sort algorithm (Numerical Recipes)
!

use ioMod,  only: stderr, verbose
use epsMod, only: eps

!-pass
integer n,ndex(n)
double precision distance(n)
!-locl
integer M,NSTACK
parameter (M=7,NSTACK=50)
character answer
integer i,ir,j,jstack,k,l,istack(NSTACK)
double precision a,temp
integer ib,itemp
jstack=0
l=1
ir=n
1 if(ir-l.lt.M)then
    do j=l+1,ir
        a=distance(j)
        ib=ndex(j)
        do i=j-1,1,-1
            if(distance(i).le.a) then
                goto 2
            endif

            distance(i+1)=distance(i)
            ndex(i+1)=ndex(i)
        enddo
        i=0
2       distance(i+1)=a
        ndex(i+1)=ib
    enddo

    if(jstack.eq.0) then
        return
    endif

    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
else
    k=(l+ir)/2
    temp=distance(k)
    distance(k)=distance(l+1)
    distance(l+1)=temp
    itemp=ndex(k)
    ndex(k)=ndex(l+1)
    ndex(l+1)=itemp

    if(distance(l+1).gt.distance(ir))then
        temp=distance(l+1)
        distance(l+1)=distance(ir)
        distance(ir)=temp
        itemp=ndex(l+1)
        ndex(l+1)=ndex(ir)
        ndex(ir)=itemp
    endif

    if(distance(l).gt.distance(ir))then
        temp=distance(l)
        distance(l)=distance(ir)
        distance(ir)=temp
        itemp=ndex(l)
        ndex(l)=ndex(ir)
        ndex(ir)=itemp
    endif

    if(distance(l+1).gt.distance(l))then
        temp=distance(l+1)
        distance(l+1)=distance(l)
        distance(l)=temp
        itemp=ndex(l+1)
        ndex(l+1)=ndex(l)
        ndex(l)=itemp
    endif

    i=l+1
    j=ir
    a=distance(l)
    ib=ndex(l)
3   continue
    i=i+1

    if(distance(i).lt.a) then
        goto 3
    endif
4   continue

    j=j-1

    if(distance(j).gt.a) then
        goto 4
    endif

    if(j.lt.i) then
        goto 5
    endif

    temp=distance(i)
    distance(i)=distance(j)
    distance(j)=temp
    itemp=ndex(i)
    ndex(i)=ndex(j)
    ndex(j)=itemp
    goto 3
5   distance(l)=distance(j)
    distance(j)=a
    ndex(l)=ndex(j)
    ndex(j)=ib
    jstack=jstack+2
  
    if(jstack.gt.NSTACK) then
        write(stderr,10)
   10       format('polysort PAUSE: NSTACK too small in sort')
        read(*,*) answer
    endif

    if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
    else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
    endif
endif
goto 1
end
!-----------------------------------------------------------------------
logical function cross (x1,y1,x2,y2,x3,y3,x4,y4)

! Coordinates (x1,y1) - (x2,y2) specify line I
! Coordinates (x3,y3) - (x4,y4) specify line II
! Subroutine tests whether intersection of the two lines occurs within
! bracket (x1,y1) - (x2,y2) and bracket (x3,y3) - (x4,y4).

use ioMod,  only: stderr         

implicit none
!-pass
double precision x1,y1,x2,y2,x3,y3,x4,y4
!-locl
double precision, parameter :: eps=1e-4

double precision dxI,dyI,dxII,dyII,aI,bI,aII,bII,xi,yi,da
!-init
cross = .false.

ai  = 0d0
bi  = 0d0
aii = 0d0
bii = 0d0


dxI = (x1-x2)
dyI = (y1-y2)

if (ABS(dxI).le.EPS .and. ABS(dyI).le.EPS) then
    write(stderr,3)
    3       format(1x,'polysort: Line I coordinates overlap')
    call exitp(1)
endif

dxII = (x3-x4)
dyII = (y3-y4)

if (ABS(dxII).le.EPS .and. ABS(dyII).le.EPS) then
    write(stderr,4)
    4       format(1x,'polysort: Line II coordinates overlap')
    call exitp(1)
endif

if (ABS(dxI).le.EPS .and. ABS(dxII).le.EPS) then
    return
endif

if (ABS(dxI).gt.EPS) then
    aI = dyI/dxI
    bI = y1 - aI * x1
endif

if (ABS(dxII).gt.EPS) then
    aII = dyII/dxII
    bII = y3 - aII * x3
endif

if (ABS(dxI).le.EPS) then
    xi = x1
    yi = aII*xi + bII
else if (ABS(dxII).le.EPS) then
    xi = x3
    yi = aI*xi + bI
else
    da = aII - aI
    if (ABS(da).le.EPS) return
    xi = (bI-bII)/da
    yi = (bI*aII-bII*aI)/da
endif

cross = ((x1-xi)*(x2-xi).le.0d0 .and. &
         (y1-yi)*(y2-yi).le.0d0 .and. &
         (x3-xi)*(x4-xi).le.0d0 .and. &
         (y3-yi)*(y4-yi).le.0d0)

return
end
