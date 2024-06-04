module elementSizes
implicit none
double precision :: minelementsize, maxelementsize
end module


program setarea

use elementSizes, only: minelementsize, maxelementsize

implicit none

double precision, parameter :: FACTOR=0.05d0, POW=1.0d0
character(len=80) fna, arg

integer, allocatable :: attr(:)
integer, allocatable :: IEN(:,:)
double precision, allocatable :: DISTS(:)
double precision, allocatable :: X(:)
double precision, allocatable :: Y(:)

integer reftype
integer refmaterial
double precision refarea, area, areamax
double precision reffactor
double precision reftosize
double precision refpow
double precision minrefresult, maxrefresult
double precision dist, absolute_dist, dmax
double precision refpoint(2)
double precision refline(2,2)
logical islineref, ispointref, ismaterialref, isglobalref

double precision :: prd, u
double precision :: xc, xd, xdist, xm, xmin, xmax
double precision :: yc, yd, ydist, ym, ymin, ymax
double precision :: dx1, dy1, dx2, dy2
integer          :: i1, i2, i3
integer          :: narg, iattr, idim, ibound, imark

integer          :: i, im, ip, k, l
integer          :: lu, lue
integer, external :: iflu, nextlu, lnblk ! from tools
integer          :: numel, npoints, numnp

double precision, external :: area0, areaall, area1, area2, area21, area3, area31, area9, areamaterial


integer ERROR

islineref = .false.
ispointref = .false.
ismaterialref = .false.
isglobalref = .false.
refarea = 100.0
refmaterial = -1
reftype = -1
reffactor = -1.0
refpow = 1.0
absolute_dist  = -1.0
minelementsize = -1.0
maxelementsize = -1.0
minrefresult     = -1.0
maxrefresult     = -1.0    

lue = iflu('stderr')
narg = command_argument_count()

if (narg.eq.0) then
    call printusage()
!write(lue,10)
!   10        format(1x,'Usage: setarea basename [options]')
    call exitp(1)
endif

call get_command_argument(1,fna)
! extend so that refinement to a combination of several
! points and lines is possible

if (narg.gt.1) then
    l = 2
    do while (l.le.narg)
        call get_command_argument(l,arg)
        if (arg.eq.'-h') then
            call printusage()
            call exitp(1)
        else if (arg.eq.'-p') then
            call get_command_argument(l+1, arg)
            read (arg,999) refpoint(1)
            call get_command_argument(l+2, arg)
            read (arg,999) refpoint(2)
            ispointref =  .true.
            l = l + 3
        else if (arg.eq.'-l') then
            call get_command_argument(l+1, arg)
            read (arg,999) refline(1,1)
            call get_command_argument(l+2, arg)
            read (arg,999) refline(1,2)
            call get_command_argument(l+3, arg)
            read (arg,999) refline(2,1)
            call get_command_argument(l+4, arg)
            read (arg,999) refline(2,2)
            islineref =  .true.
            l = l + 5
        else if (arg.eq.'-m') then
            call get_command_argument(l+1, arg)
            read (arg,999) refmaterial
            ismaterialref = .true.
            l = l + 2
        else if (arg.eq.'-d') then
            call get_command_argument(l+1, arg)
            read (arg,999) absolute_dist
            l = l + 2
        else if (arg.eq.'-maxsize') then
            call get_command_argument(l+1, arg)
            read (arg,999) maxelementsize
            l = l + 2
        else if (arg.eq.'-minsize') then
            call get_command_argument(l+1, arg)
            read (arg,999) minelementsize
            l = l + 2
        else if (arg.eq.'-minrefresult') then
            call get_command_argument(l+1, arg)
            read (arg,*) minrefresult
            l = l + 2
        else if (arg.eq.'-maxrefresult') then
            call get_command_argument(l+1, arg)
            read (arg,*) maxrefresult
            l = l + 2
        else if (arg.eq.'-f') then
            call get_command_argument(l+1, arg)
            read (arg,999) reffactor
            l = l + 2
        else if (arg.eq.'-r') then
            call get_command_argument(l+1,arg)
            read (arg,*) reftype
            isglobalref = reftype.eq.0
            l = l + 2
        else if (arg.eq.'-a') then
            call get_command_argument(l+1, arg)
            read (arg,999) refarea
            l = l + 2
        else if (arg.eq.'-x') then
            call get_command_argument(l+1, arg)
            read (arg,999) refpow
            l = l + 2
        else
            write(*,*) 'Option:', arg, 'unknown'
            call printusage()
            call exitp(1)
        endif
    enddo

    if (maxrefresult.gt.-1.and. &
        minrefresult.gt.-1.and. &
        minrefresult.gt.maxrefresult) then
        write(lue,*) 'minrefresult must be smaller than maxrefresult.'
        call exitp(1)
    endif

    if (maxelementsize.gt.-1.and. &
        minelementsize.gt.-1.and. &
        minelementsize.gt.maxelementsize) then
        write(lue,*) 'minsize must be smaller than maxsize.'
        call exitp(1)
    endif

    if (reftype.ne.0 .and. &
        reftype.ne.1 .and. &
        reftype.ne.1002 .and. &
        reftype.ne.1021 .and. &
        reftype.ne.1003 .and. &
        reftype.ne.1031 .and. &
        reftype.ne.9) then
          write(lue,*) "Did not recognize refinement type", reftype
          call printusage()
          stop "Leaving seteara..."
    endif

endif

  999 format(G14.0)
lu = nextlu(0)

call openf(lu,fna(1:lnblk(fna))//'.node','old')

read (lu,*,err=1000,end=1000,iostat=ERROR) NUMNP,idim,iattr,ibound
if (NUMNP.le.0) then
    write(lue,*) 'setarea: no nodes in .node file'
    call exitp(1)
endif
if (idim.ne.2) then
    write(lue,*) 'setarea: .node file header should specify dimension=2'
    call exitp(1)
endif
if (iattr.ne.0) then
    write(lue,*) 'setarea is not ready for nodal attributes'
    write(lue,*) '  did you mean to specify node markers in your .poly file?'
    call exitp(1)
endif
if (ibound.ne.1) then
    write(lue,*) 'setarea: expects one node marker in .node file'
    call exitp(1)
endif

ALLOCATE(X(NUMNP),STAT=ERROR)
if (error.ne.0) then
    stop 'setarea: x-coordinate array memory allocation failed'
endif
ALLOCATE(Y(NUMNP),STAT=ERROR)
if (error.ne.0) then
    stop 'setarea: y-coordinate array memory allocation failed'
endif

! find extrema of the domain
xmin = 99999999
xmax = -99999999
ymin = 99999999
ymax = -99999999

do k=1,NUMNP
    read (lu,*,err=1100,end=1100,iostat=ERROR) ip,X(k),Y(k),imark
    if (X(k).gt.xmax) xmax = X(k)
    if (X(k).lt.xmin) xmin = X(k)
    if (Y(k).gt.ymax) ymax = Y(k)
    if (Y(k).lt.ymin) ymin = Y(k)
end do

call closef(lu)
write(lue,20) fna(1:lnblk(fna))//'.node'
   20 format(1x,'EOF ',90a)

!
call openf(lu,fna(1:lnblk(fna))//'.ele','old')
read (lu,*,err=1200,end=1200,iostat=ERROR) NUMEL,npoints,iattr

if (NUMEL.le.0) then
    write(lue,*) 'setarea: no elements in .ele file'
    call exitp(1)
endif
if (npoints.ne.3) then
    write(lue,*) 'setarea: only 3 points per triangle supported'
    call exitp(1)
endif
if (iattr.ne.1) then
    write(lue,*) 'setarea only supports 1 attribute per triangle: material number'
    call exitp(1)
endif

ALLOCATE(IEN(3,NUMEL),STAT=ERROR)
if (error.ne.0) then
    stop 'setarea: element array memory allocation failed'
endif
ALLOCATE(attr(NUMEL),STAT=ERROR)
if (error.ne.0) then
    stop 'setarea: attribute array memory allocation failed'
endif
ALLOCATE(DISTS(NUMEL),STAT=ERROR)
if (error.ne.0) then
    stop 'setarea: distance array memory allocation failed'
endif
attr = 0
!
do k=1,NUMEL
    read (lu,*,err=1300,end=1300,iostat=ERROR) im,(IEN(i,k),i=1,3),attr(k)
    if (attr(k).eq.0) then
        write(lue,*) 'setarea warning: element ',k,' has material number zero'
    endif
    do i=1,3
        IEN(i,k) = IEN(i,k) + 1
    enddo
end do

call closef(lu)
write(lue,20) fna(1:lnblk(fna))//'.ele'

!     default to point refinement in case no refinement (point, line, material) was set
if ((.not.ispointref).and.(.not.islineref) &
    .and.(.not.ismaterialref).and.(.not.isglobalref)) then
    refpoint(1) = 0.5*(xmin+xmax)
    refpoint(2) = 0.5*(ymin+ymax)
    ispointref = .true.
endif

if (reftype.lt.0) then
    reftype = 0
endif

if (reffactor.lt.0.0) then
    reffactor = factor
endif
    
if (ispointref) then
    write(lue,21) refpoint(1), refpoint(2), reftype, reffactor
   21      format(1x,'Performing point refinement towards',1x, &
        '(',G14.6,', ',G14.6,')',/1x,'with refinement type',1x,I4, &
         1x,'and refinement factor',1x,G14.6) 
else if (islineref) then
    write(lue,22) refline(1,1), refline(1,2), refline(2,1), &
      refline(2,2), reftype, reffactor
   22      format(1x,'Performing line refinement towards',1x, &
        '(',G14.6,', ',G14.6,') - (',G14.6,', ',G14.6,')'/1x, &
        'with refinement type',1x,I4, &
        1x,'and refinement factor',1x,F14.6) 
else if (ismaterialref) then
    write(lue,23) refmaterial, reffactor
   23      format(1x,'Performing refinement of elements with',1x, &
       'material number',1x, I4,/1x, &
       'with refinement factor',1x,G14.6)
else if (isglobalref) then
    write(lue,24) reffactor
   24      format(1x,'Performing global refinement of elements with',1x, &
        'with refinement factor',1x,G14.6)
endif

if (minelementsize.gt.-1.0) then
!             write(*,*) 'Only refining elements bigger than', 
!     >         minelementsize
endif

if (maxelementsize.gt.-1.0) then
!            write(*,*) 'Only refining elements smaller than', 
!     >       maxelementsize
endif

    write(*,*) NUMEL
    xmax = MAX(ABS(xmax),ABS(xmin))
    ymax = MAX(ABS(ymax),ABS(ymin))
    areamax = -99999999
    dmax = -999999999

!   compute some maximum distances and areas    
    do k=1,NUMEL
  i1 = IEN(1,k)
  i2 = IEN(2,k)
  i3 = IEN(3,k)
  dx1 = X(i2)-X(i1)
  dy1 = Y(i2)-Y(i1)
  dx2 = X(i3)-X(i1)
  dy2 = Y(i3)-Y(i1)
  prd = 0.5*(dx1*dy2-dy1*dx2)
  if (ABS(prd).lt.1e-20) then
  write(lue,40) k
   40        format(1x,'degenerate triangle ',I5)
  call exitp(1)
  endif
  prd = ABS(prd)
  if (prd.gt.areamax) areamax = prd
!               
  xc = (X(i1)+X(i2)+X(i3))/3d0
  yc = (Y(i1)+Y(i2)+Y(i3))/3d0

  if (islineref) then
      xd = refline(2,1) - refline(1,1)
      yd = refline(2,2) - refline(1,2)

      u = ((xc - refline(1,1))*xd + &
           (yc - refline(1,2))*yd) / &
         (xd*xd + yd*yd)
      if ((u.ge.0.0).AND.(u.le.1.0)) then
          xm = refline(1,1) + u*xd
          ym = refline(1,2) + u*yd
          xdist = xc - xm
          ydist = yc - ym
      else if (u.lt.0.0) then
          xdist = xc - refline(1,1)
          ydist = yc - refline(1,2)
      else if (u.gt.1.0) then
          xdist = xc - refline(2,1)
          ydist = yc - refline(2,2)
      else
          write(0,*) "Setarea could not determine xdist and ydist."
          write(0,*) "Please contact model support."
          stop "Exiting"
      endif
      dist = sqrt(xdist*xdist + ydist*ydist)

  else if (ispointref) then
      xdist = xc - refpoint(1)
      ydist = yc - refpoint(2)
      dist = sqrt(xdist*xdist + ydist*ydist)
  endif
  DISTS(k) = dist

  if (dist.gt.dmax) then
      dmax = dist
  endif

enddo

do k=1,NUMEL
  i1 = IEN(1,k)
  i2 = IEN(2,k)
  i3 = IEN(3,k)
  dx1 = X(i2)-X(i1)
  dy1 = Y(i2)-Y(i1)
  dx2 = X(i3)-X(i1)
  dy2 = Y(i3)-Y(i1)
  prd = 0.5*(ABS(dx1*dy2-dy1*dx2))
  area = prd 
  xc = (X(i1)+X(i2)+X(i3))/3d0
  yc = (Y(i1)+Y(i2)+Y(i3))/3d0

  if (islineref .or. ispointref .or. isglobalref) then

    if ((minelementsize.gt.-1.0.and.area.lt.minelementsize).or. &
       (maxelementsize.gt.-1.0.and.area.gt.maxelementsize))then
      ! element is outside the area range. Do not refine it.
          reftosize = area
    else
      ! call the actual refinement routines, based on type
      if (reftype.eq.0) then
          reftosize = area0 (reffactor,area,DISTS(k),dmax,absolute_dist)
      else if (reftype.eq.1) then
          reftosize = area1(reffactor,area,DISTS(k),dmax,refarea,absolute_dist)
      else if (reftype.eq.1002) then
          reftosize = area2(reffactor,refpow,area,DISTS(k),dmax,absolute_dist)
      else if (reftype.eq.1021) then
!  can be merged with reftype 1002, will be active if the user supplies -a percentage
          reftosize = area21(reffactor,refpow,area,DISTS(k),dmax,refarea,absolute_dist)
      else if (reftype.eq.1003) then
          reftosize = area3(reffactor,refpow,area,DISTS(k),dmax,absolute_dist)
      else if (reftype.eq.1031) then
!  can be merged with reftype 1003, will be active if the user supplies -a percentage
          reftosize = area31(reffactor,refpow,area,DISTS(k),dmax,refarea,absolute_dist)
      else if (reftype.eq.9) then
          reftosize = area9(reffactor,refpow,area,DISTS(k),dmax,absolute_dist)
      else
          write(lue,*) "Did not recognize refinement type", reftype
          call printusage()
          stop "Leaving seteara..."
      endif
    endif
  else if (ismaterialref) then
      reftosize = areamaterial(reffactor, area, attr(k), refmaterial, minelementsize)
  endif

  if (minrefresult.gt.-1.and. &
     reftosize.lt.minrefresult .and. reftosize.gt.0) then
      reftosize=minrefresult
  endif

  if (maxrefresult.gt.-1.and. &
     reftosize.gt.maxrefresult) then
      reftosize=maxrefresult
  endif

  ! filter out accidental wrong material refines
  if(refmaterial.gt.-1) then
      if(refmaterial.ne.attr(k)) then
          reftosize = -1.0d0
      endif
  endif

  write(*,*) k-1, reftosize

enddo
if (reftype.eq.0 .and. absolute_dist.lt.0d0) write(0,*) 'setarea used tapering distance=',dmax

DEALLOCATE(X, STAT=ERROR)
DEALLOCATE(Y, STAT=ERROR)
DEALLOCATE(IEN, STAT=ERROR)
DEALLOCATE(attr, STAT=ERROR)
DEALLOCATE(DISTS, STAT=ERROR)
 
call exitp(0)

1000  write(lue,1001) ERROR
1001  format(1x,'setarea: read error ',I6,' on .node file header')
call exitp(1)
1100  write(lue,1101) ERROR,k
1101  format(1x,'setarea: read error ',I6,' while reading record ',I7,' of .node file')
call exitp(1)
1200  write(lue,1201) ERROR
1201  format(1x,'setarea: read error ',I6,' on .ele file header')
call exitp(1)
1300  write(lue,1301) ERROR,k
1301  format(1x,'setarea: read error ',I6,' while reading record ',I7,' of .ele file')
call exitp(1)

end
!-----------------------------------------------------------------------
double precision function areamaterial(f, area, currentmat, refmat, minimalElementSize)

use elementSizes, only: minelementsize, maxelementsize

implicit none

double precision :: f, area, minimalElementSize
integer          :: refmat,currentmat

if (currentmat.eq.refmat) then 
    areamaterial = area * f
    if (minelementsize.gt.0.0) then
        areamaterial = max (areamaterial,minelementsize)
    endif
    if (maxelementsize.gt.0.0) then
        areamaterial = min (areamaterial,maxelementsize)
    endif
else
    areamaterial = -1.0
endif
return
end
!-----------------------------------------------------------------------
double precision function area0 (f,area,dist,dmax,absolute_dist)

implicit none
!-pass
double precision :: f,area,dist,dmax,absolute_dist
!-locl
double precision, parameter :: ONE=1d0
double precision :: relative_distance,taper

if (absolute_dist.eq.-ONE) then
    relative_distance = dist/dmax
else
    relative_distance = dist/absolute_dist
endif

taper = ONE - EXP(-3d0*relative_distance*relative_distance)
area0 = area * (f + (ONE-f)*taper)

return
end
!-----------------------------------------------------------------------
double precision function area1 (f,area,dist,dmax,refarea,absolute_dist)

implicit none

double precision :: f,area,dist,dmax,refarea,absolute_dist

double precision :: ratio

ratio = refarea * 0.01

if (absolute_dist.eq.-1d0) then
    if (dist/dmax .le. ratio) then
        area1 = area * f
    else
        area1 = -1.0
    endif
else
    if (dist .le. absolute_dist) then
        area1 = area * f
    else
        area1 = -1.0
    endif

endif
return
end
!-----------------------------------------------------------------------
double precision function area2 (f,pow,area,dist,dmax,absolute_dist)

implicit none

double precision :: f,pow,area,dist,dmax,absolute_dist


if (absolute_dist.eq.-1d0) then
    area2 = area * (f + (1d0-f)*(dist/dmax)**pow)
else
    area2 = area * (f + (1d0-f)*(dist/absolute_dist)**pow)
endif
return
end
!-----------------------------------------------------------------------
double precision function area21 (f,pow,area,dist,dmax,refarea,absolute_dist)

implicit none

!-pass
double precision :: f,pow,area,dist,dmax,refarea,absolute_dist
!-local
double precision :: ratio


ratio = refarea * 0.01
if (absolute_dist.eq.-1d0) then
    if (dist/dmax .le. ratio) then
        area21 = area * (f + (1d0-f)*(dist/(ratio * dmax))**pow)
    else
        area21 = -1.0
    endif
else
    if (dist .le. absolute_dist) then    
        area21 = area * (f + (1d0-f)*(dist/absolute_dist)**pow)
    else
        area21 = -1.0
    endif
endif
return
end
!-----------------------------------------------------------------------
double precision function area3 (f,pow,area,dist,dmax,absolute_dist)

implicit none

double precision :: f,pow,area,dist,dmax,absolute_dist


area3 = area * ceiling(((dist/dmax)**pow)/f) * f
return
end
!-----------------------------------------------------------------------
double precision function area31 (f,pow,area,dist,dmax,refarea,absolute_dist)

implicit none

double precision :: f,pow,area,dist,dmax,refarea,absolute_dist

double precision :: ratio


ratio = refarea * 0.01
if (absolute_dist.eq.-1d0) then
    if (dist/dmax .le. ratio) then
        area31 = area * ceiling(((dist/(ratio*dmax))**pow)/f) * f
    else
        area31 = -1.0
    endif
else
    if (dist .le. absolute_dist) then
        area31 = area * ceiling(((dist/absolute_dist)**pow)/f) * f
    else
        area31 = -1.0
    endif
endif
return
end
!-----------------------------------------------------------------------
double precision function area9 (f,pow,area,dist,dmax,absolute_dist)
!area9(reffactor,refpow,area,DISTS(k),dmax,refarea,absolute_dist)
implicit none
double precision :: area, f, pow, dist, dmax, absolute_dist

area9 = pow + f * dist * dist * dist * dist

return
end

!-----------------------------------------------------------------------
subroutine printusage()

implicit none

integer           :: lue
integer, external :: iflu


lue = iflu('stderr')
write(lue,10)
10     format(1x, &
  'Usage: setarea basename [options] [> basename.area]'/1x, &
  'options:'/3x,&
  '-p x y         (point refinement)'/3x,&
  '-l x1 y1 x2 y2 (line refinement)'/3x,&
  '-m mat         (material refinement; elements with'/3x,&
  ' material mat are refined only)'/3x,&
  '-f factor      (< 1 refinement ratio)'/3x,&
  '-r type        (refinement type '/3x,&
  '                0 - tapered element refinement in a domain'/3x,&
  '                1 - uniform refinement'/3x,&
  '                    in part of domain'/3x,&
  '             1002 - refine linearly towards point/line of'/3x,&
  '                    refinement over entire domain'/3x,&
  '             1021 - refine linearly towards point/line of'/3x,&
  '                    refinement over part of domain'/3x,&
  '             1003 - refine linearly towards point/line of'/3x,&
  '                    refinement over entire domain'/3x,&
  '             1031 - refine linearly towards point/line of'/3x,&
  '                    refinement over part of domain'/3x,&
  '                (1002 (1021) and 1003 (1031) differ slightly))'/3x,&
  '             1009 - refinement according to George, 20Jan14'/3x,&
  '-a percentage  (denotes which part of domain should be',&
  ' refined'/3x,&
  '-d absolute distance to which the mesh should be refined.'/3x,&
  '                (100 = entire domain); used only by'/3x,&
  '                refinement types 1, 21 and 31)'/3x,&
  '-minsize min   (minimum element size; specify minimum area'/3x,&
  ' for elements that are to be refined and avoids small'/3x,&
  ' elements when performing successive refinement steps)'/3x,&
  '-maxsize max   (maximum element size; specify maximum area'/3x,&
  ' for elements that are to be refined)'/3x,&
  '-minrefresult min   (estimated minimum element size after'/3x,&
  ' refinement)'/3x,&
  '-maxrefresult max   (estimated maximum element size after'/3x,&
  ' refinement)'/3x,&
  '-x power       (to get sup (> 1) or sub (< 1) linear',&
  ' refinement;'/3x,&
  '                used only by types 2, 21, 3 and 31)'/1x,&
  '(default: point refinement to middle of the domain with'/3x,&
  'refinement type 0, factor 0.05)')
return
end
