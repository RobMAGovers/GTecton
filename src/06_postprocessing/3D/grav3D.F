subroutine gravitytetrahedron (IEN,X,&
                               D,MAT,PRPMAT,STR,TFAULT,&
                               lonmin,lonmax,lonstep,latmin,latmax,latstep,height,&
                               file)
 
USE MODELDEFINITION
use materialsmodule, only: MaterialIDOfElement,matdatactx
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: nsd,NDOF,NSTR,NEN
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif
use iomodule,        only: stderr, verbose
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: coordtype, detdilatation, detgrad, detgrav, &
                           detpot, opgravchange
use errors,          only: exec

! Routine for calculating gravity signals. 
implicit none
double precision, parameter :: REarth=6371000d0

!-pass
integer   :: ien, mat
dimension :: IEN(NEN,*),MAT(*)

 
double precision X,D,PRPMAT,TFAULT,STR
dimension X(NSD,*),D(NDOF,*),TFAULT(NDOF,*),PRPMAT(6,*),STR(NSTR,*)

double precision lonmin,lonmax,lonstep,latmin,latmax,latstep,height
#include "lenchar.i" 
character(len=LENCHAR) file
 



!-locl
integer i,j
integer numnn,numlon,numlat,n,nn,m
double precision lon,lat,rho
double precision xloc,yloc,zloc,gabs,gabs1,gabs0
double precision e,U,dU ! potential
double precision g(3),g1(3),g0(3),dg(3) ! gravity
double precision gEast,gNorth,gUp
double precision gg(3,3),dgg(3,3)
integer lu
logical OpDeform,OpDensChange
logical test
integer, external :: nextlu,ifopen
integer :: ierr
integer :: slipperySequenceNr
!-init
OpDensChange = detdilatation
OpDeform =  OpGravChange 
test = .false.

write(stderr,*) 'coordtype,detpot,detgrav,detgrad'
write(stderr,*) coordtype,detpot,detgrav,detgrad
write(stderr,*) 'verbose',verbose
write(stderr,*) 'effect of dilatation',detdilatation
write(stderr,*) 'calculate gravity change',OpDeform
lu = nextlu(0)
ierr = ifopen(lu,file,'unknown')
if (ierr.ne.0) then
    ierr = EXEC
    return
endif

! set number of elements
if (elmon) then
    numnn = nnel
else
    numnn = NUMEL
endif

if (test) then
    numnn = 1
endif

if (detgrad) then
    write(stderr,*) 'gravity gradients not tested, exit'
    return
endif

if (detgrav) then
    write(stderr,*) 'gravity not tested, exit'
    return
endif

! Calculate number of sample points
numlon = nint(ABS((lonmax-lonmin)/lonstep))+1
numlat = nint(ABS((latmax-latmin)/latstep))+1
write(*,*) 'numlon',numlon,'numlat',numlat
if (numlon.le.0d0.or.numlat.le.0d0) then
    return
endif


! loop on sample points
do i=1,numlon
    do j=1,numlat

! get reference sample point
        lon = lonmin+lonstep*(i-1)
        lat = latmin+latstep*(j-1)
        if (coordtype.eq.'sphe') then
! convert to cartesian coordinates
            call geoc2cart (lon,lat,height,REarth,xloc,yloc,zloc)
 !  write(stderr,*) 'lon',lon,'lat',lat,'height',height,'x y z',xloc,yloc,zloc
        else if (coordtype.eq.'cart') then
   ! do not convert to cartesian coordinates but take lon, lat and heigth as
   ! x,y,z
            xloc=lon
            yloc=lat
            zloc=height
        else
            write(stderr,*) 'coordtype:',coordtype
            write(stderr,*) 'wrong option coordinate type'
            goto 1000
        endif
! initialise
        U = 0d0
        g(1:3) = 0d0
        g1(1:3) = 0d0
        g0(1:3) = 0d0
        gg(1:3,1:3) = 0d0

! loop on all elements
        do nn=1,numnn
!   write(stderr,*) 'element',nn
            if (elmon) then
 ! if monitor
                n = monel(nn)
            else
                n = nn
            endif

 ! get material type
            m = MaterialIDOfElement(n)
 ! get density
            Rho = dble(PRPMAT(5,m))
 !  write(stderr,*) 'density',Rho,'m',m
 !  write(stderr,*) 'all',PRPMAT(:,m)
 ! get gravity contribution from element
   
   if (OpDeform) then
           ! first deformed

     call gravityofelement(IEN,X,D, &
     TFAULT,n,Rho,STR,xloc,yloc,zloc,detpot,detgrav,detgrad,&
     dU,dg,dgg,OpDeform,OpDensChange)

   
     if (detpot) then
       U=U+dU
    !   write(stderr,*) 'potential cumm',n,U
     endif
     if (detgrav) then
       g1=g1+dg
     endif
     if (detgrad) then
       gg=gg+dgg
     endif
  
   ! now undeformed
   !  Deformed = .false.

     call gravityofelement(IEN,X,D, &
     TFAULT,n,Rho,STR,xloc,yloc,zloc,detpot,detgrav,detgrad,&
     dU,dg,dgg,.false.,.false.)


 !    write(stderr,*) 'gravity undeformed',dg
     if (detpot) then
         U=U-dU
     endif
     if (detgrav) then
         g0=g0-dg
 !        write(stderr,*) 'gravity change cummulative',n,g
     endif
     if (detgrad) then
        gg=gg-dgg
     endif

   else
     ! only undeformed
     call gravityofelement(IEN,X,D, &
     TFAULT,n,Rho,STR,xloc,yloc,zloc,detpot,detgrav,detgrad,&
     dU,dg,dgg,OpDeform,.false.)



     if (detpot) then
       U=U+dU
     endif
     if (detgrav) then
       g=g+dg
     endif
     if (detgrad) then
       gg=gg+dgg
     endif
   endif


!   if (elmon) then
! ! write contribution per monitored element
!   if (verbose) then
!    if (detpot) then
!    write(luMon,1) n,lon,lat,height,U
! 1                   format(1x,'potential contribution of element', &
!                 1x,I5,' at (',1PG12.4,',',1PG12.4,',',1PG12.4,') = ', &
!                 1PG11.3)
!    endif
!    if (detgrav) then
!    write(luMon,2) n,lon,lat,height,g
! 2                   format(1x,'gravity contribution of ', &
!                 'element ',I5,' at (',1PG12.4,',',1PG12.4,',',1PG12.4, &
!                 ') = ',1PG11.3)
!    endif
!    if (detgrad) then
!    write(luMon,3) n,lon,lat,height,gg
! 3                   format(1x,'gravity gradient contribution of ', &
!                 'element ',I5,' at (',1PG12.4,',',1PG12.4,',',1PG12.4, &
!                 ') = ',1PG11.3)
!    endif
!   else
!    ! not verbose
!    if (detpot) then
!      write(luMon,*) n,lon,lat,height,U
!    endif
!    if (detgrav) then
!      write(luMon,*) n,lon,lat,height,g
!    endif      
!   endif
!   endif
   
  enddo ! end loop elements

  if (detgrav) then
    if (OpDeform) then
      ! subract gravity from undeformed mesh from gravity of deformed mesh to
      ! get gravity change
      g = g1 - g0
      ! convert gravity change vector to scalar
      gabs1=sqrt(g1(1)*g1(1)+g1(2)*g1(2)+g1(3)*g1(3))
      gabs0=sqrt(g0(1)*g0(1)+g0(2)*g0(2)+g0(3)*g0(3))
      gabs = gabs1-gabs0
    else
      ! convert gravity vector to scalar
      gabs=sqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
    endif
  endif

  ! in case spherical coordinates and gravity or gravity gradient, transform to
  ! NEU coordinates
  if (coordtype.eq.'sphe') then
          if (detgrav) then
             call xyz2enu (g,lon,lat,gEast,gNorth,gUp)
             g(1)=gEast
             g(2)=gNorth
             g(3)=gUp
          endif
          if (detgrad) then
           ! find out how to transform array
           write(stderr,*) 'transformation of gravity gradients not yet implemented'
           return
          endif
  endif

! write result
  if (.not.elmon) then
    if (detpot.and..not.detgrav) then
     write(lu,*) lon,lat,height,gUp
    else if (.not.detpot.and.detgrav) then
     write(lu,*) lon,lat,height,g,gabs
    elseif (detpot.and.detgrav) then
     write(lu,*) lon,lat,height,gUp,g
    endif
  endif
 enddo ! end loop latitude
enddo ! end loop longitude

1000 continue ! weg hier
 return
 end


!-----------------------------------------------------------------------
subroutine gravityofelement (IEN,X,D,TFAULT, &
  N,RHO0,STR,X0,Y0,Z0,detpot,detgrav,detgrad,&
  Potential,Gravity,GravGradient,OpDeform,OpDensChange)

use constants,       only: useDXE
USE MODELDEFINITION
use modeltopology,   only: ndof, nen, nstr, nsd

#ifndef EXPERIMENTAL_ReducedSlipAlloc
use meshdatamodule,  only: meshdatactx
#endif
! use modeldatamodule, only: modeldatactx
! LMF may be deleted from this function
! Routine to calculate gravity signal in (x0,y0,z0) of a single element.
implicit none


!-pass 
 integer N
integer :: ien
 dimension IEN(NEN,*)
 logical detpot,detgrav,detgrad
 logical OpDensChange,OpDeform

double precision :: x0, y0, z0

double precision :: x, d, tfault, str
dimension X(NSD,*),D(NDOF,*),TFAULT(NDOF,*),STR(NSTR,*)

double precision :: Potential,Gravity(3),GravGradient(3,3)
 character(len=4) gravtype


!-locl!
double precision, external :: AJ1
double precision :: xl, dl
 dimension XL(NSD,NEN),DL(NDOF,NEN)
 double precision :: dilatation,RHO0,RHO
! get coordinates of element (XL)
 call LCOORD (X,XL,IEN(1,n))
 ! option to update element with displacements
 if (OpDeform) then
! get displacements (DL)
 call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
!  if (n.eq.39) then
!    write(stderr,*) 'prior fault',n,DL
!  endif
! add faulted node displacements 
 call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN, 180)
! add slippery node displacements

 call ADDSNE (DL,n,NDOF,NEN,useDXE)

! add displacements to coordinates
 call REZONE (XL,DL,1d0)
! if (n.eq.37) then
! write(stderr,*) 'after fault',n,DL
! endif
! if (N.eq.1) then
!   write(stderr,*) 'coord: ',X(1,1),X(2,1),X(3,1)
! endif
 endif

 ! option to include gravity changes due to density changes
 if (OpDensChange) then
! element density changes due to volumetric strain, so update
! get volume change (3 times to small)
! call SHPTET (ZERO,ZERO,ZERO,XL,VOLAFTER,SH,IERR,1)
! meanstrain = dilate(STR(1,n)) dont use function dilate, it returns mean
! strain, not volume change!
 dilatation = AJ1(STR(1,n))
! dilatation2 = 3*meanstrain
! dilatation2 = (VOLAFTER-VOLPRE)/VOLPRE
!  write(*,*) 'dilatation',n,dilatation,dilatation2,dilatation/dilatation2
! update density (rho + drho) (V + dV) = rho V
! drho dV -> 0
! drho = -rho * dV = -rho * divergence(u) = -rho * 1/3 * firstinvariant(strain)
    RHO = RHO0*(1-dilatation)
 else
    RHO = RHO0
 endif

!    Gravity contribution from element based on method of Werner and Scheeres 1996
 call WernerScheeres(XL,X0,Y0,Z0,RHO,detpot,detgrav,detgrad,Potential,Gravity,GravGradient)
! 

 return
 end

!-----------------------------------------------------------------------
  subroutine WernerScheeres(X,X0,Y0,Z0,Rho,detpot,detgrav,detgrad,Potential,Gravity,GravGradient)
!    Gravity contribution from element based on method of Werner and Scheeres x

! functions to calculate potential, gravity or gravity gradient at (X0,Y0,Z0)
! due to a tetrahedron with corner coordinates X and mass density rho.

! Ref: Werner & Scheeres, EXTERIOR GRAVITATION OF A POLYHEDRON DERIVED AND
! COMPARED WITH HARMONIC AND MASCON GRAVITATION REPRESENTATIONS 
! OF ASTEROID 4769 CASTALIA, Celestial Mechanics and Dynamical Astronomy, 1996
!
! dbt broerse May 2016

 USE MODELDEFINITION
 USE CONSTANTS, only: grav_const
 USE ALGEBRA
use modeltopology, only: NSD,NEN,NDOF,NSTR
use iomodule,           only: stderr

 implicit none
!-pass
 double precision X
 dimension X(NSD,NEN)
 double precision :: X0,Y0,Z0
 double precision :: Rho
 logical detpot,detgrav,detgrad

!-locl
  logical test
  integer NEdges,NFaces,NEdgesFace,NCornersFace
  parameter (NEdges=6,NFaces=4,NEdgesFace=3,NCornersFace=3)
  integer i,j
  integer side(2),localnodes(2),localedges(2)
  integer lien(3)
  double precision nf(NSD,NFaces),nef(NSD,NEdgesFace,NFaces)
  double precision rf(NSD,NCornersFace),rfl(NSD)
  double precision re(NSD,2),Ee(NSD,NSD),Le
  double precision Ff(NSD,NSD),omf
  double precision rel(2),sidelength,sidevec(3)
  double precision prodEEre(3),prodFFrf(3),crossprodrf(3)
  double precision Potential,Gravity(3),GravGradient(3,3)
  double precision edgecontrp,facecontrp
  double precision edgecontrg(3),facecontrg(3)
  double precision edgecontrgg(3,3),facecontrgg(3,3)

  test = .false.
!  write(stderr,*) 'grav const',grav_const
!  write(stderr,*)
!  write(stderr,*) 'density',RHO
 if (NEN.ne.4) then
 ! since code has been written for tetrahedrons NEN has to be four, since there are four corners
   write(stderr,1) NEN
 1        format(1x,'error with element type, code only written for tetrahedra, current element corners ',I5)
 endif
 
 ! potential is determined using equation (10)
 ! gravity (attraction) is determined using equation (15)
 ! gravity gradient is determined using equation (16)
 
 ! all equations have the same elements in common:
 ! (a)  re : a vector from the field point to any point on the edge
 ! (b)  Ee : edge dyad (section 2.1.8)
 ! (c)  Le : the dimensionless per-edge factor, equation (7)
 ! (d)  rf : a vector from the field point to some fixed point in the face plane
 ! (e)  Ff : face dyad (section 2.1.8)
 ! (f)  omf : the dimensionless perface factor, equation (27)
 !
 ! with variables
 ! sigma = density
 ! G = gravitional constant
 ! nf = outward-pointing face normal vector 
 ! nef = edge-normal vector (depending on element face)
 !
 ! furthermore
 ! rel = length of re
 ! el = length of edge
 ! rfl = length of rf

!-init
!  if (gravtype.eq.'pote') then 
    edgecontrp=0d0
    facecontrp=0d0
    Potential=0d0
!  elseif (gravtype.eq.'grav') then
    do i=1,3
      edgecontrg(i)=0d0
      facecontrg(i)=0d0
      Gravity(i)=0d0
    enddo
!  elseif (gravtype.eq.'grad') then
    do i=1,3
     do j=1,3
       edgecontrgg(i,j)=0d0
       facecontrgg(i,j)=0d0
       GravGradient(i,j)=0d0
     enddo
    enddo
!  endif
 ! 
 if (test) then 
   X(1,1)=0d0
   X(2,1)=0d0
   X(3,1)=0d0
   X(1,2)=1d4
   X(2,2)=0d0
   X(3,2)=0d0
   X(1,3)=0d0
   X(2,3)=1d4
   X(3,3)=0d0
   X(1,4)=0d0
   X(2,4)=0d0
   X(3,4)=1d4
 endif

! determine per element, the normal vectors  
 call getnormals(X,nf,nef)
 
! do i=1,NFaces
!   write(stderr,*) 'nf',i,nf(:,i)
!   do j=1,NEdgesFace
!     write(stderr,*) 'nef',nef(:,j,i)
!   enddo
! enddo

! start with contribution from tetrahedron edges 
! loop on element edges
 do i=1,NEdges
  ! get local side number for current edge
  call EDGES(i,side,localnodes,localedges)
!  write(stderr,*) 'side nrs',side
  ! Get local nodal point numbers of the side (lien)
  call SIDENP(side(1),.true.,lien)
  
  ! calculate vector re from field point to edge
  ! first point on edge 
!  write(stderr,*) 'x',x(1,lien(localnodes(1))),'x0',X0
  re(1,1) = x(1,lien(localnodes(1)))-X0
  re(2,1) = x(2,lien(localnodes(1)))-Y0
  re(3,1) = x(3,lien(localnodes(1)))-Z0
  ! second point on edge 
  re(1,2) = x(1,lien(localnodes(2)))-X0
  re(2,2) = x(2,lien(localnodes(2)))-Y0
  re(3,2) = x(3,lien(localnodes(2)))-Z0
!  write(stderr,*) 're',re(1,1),re(2,1),re(3,1),re(1,2),re(2,2),re(3,2)
!  write(stderr,*) 'x1',x(:,lien(localnodes(1))),'x2',x(:,lien(localnodes(2)))
!  write(stderr,*) 'nf', nf(:,side(1)),'nef',nef(:,localedges(1),side(1))
!  write(stderr,*) 'nf', nf(:,side(2)),'nef',nef(:,localedges(2),side(2))
!  write(stderr,*) 'sides',side(1),side(2)

  ! calculate Ee edge diad 
  call edgediad(nf(:,side(1)),nef(:,localedges(1),side(1)),nf(:,side(2)),nef(:,localedges(2),side(2)),EE)
!  write(stderr,*) 'edge diad',Ee
  ! calculate Le the dimensionless per-edge factor
  
  ! first calculate absolute lengths of re
  do j=1,2
    rel(j)=sqrt(re(1,j)*re(1,j)+re(2,j)*re(2,j)+re(3,j)*re(3,j))
!    write(stderr,*) 'r vector length',rel(j)
  enddo
  ! than length of edge
  do j=1,3
    sidevec(j)=x(j,lien(localnodes(1)))-x(j,lien(localnodes(2)))
!    write(stderr,*) 'side vector',sidevec(j)
  enddo
  sidelength=sqrt(sidevec(1)*sidevec(1)+sidevec(2)*sidevec(2)+sidevec(3)*sidevec(3))
!  write(stderr,*) 'side length',sidelength
  ! now per-edge factor
  Le = dlog((rel(1)+rel(2)+sidelength)/(rel(1)+rel(2)-sidelength))
  ! check for infinity
  if (rel(1)+rel(2)-sidelength==0) then
          Le=0
  endif
  if (detpot) then
   ! potential
   call VMPRD(Ee,re(:,1),prodEEre,3,3)
   edgecontrp=edgecontrp+dot_product(re(:,1),prodEEre)*Le
!   write(stderr,*) 'edge contribution',edgecontrp
  endif
  if (detgrav) then
    ! gravitational attraction
    call VMPRD(Ee,re(:,1),prodEEre,3,3)
    edgecontrg=edgecontrg+prodEEre*Le
!    write(stderr,*) 'edge contribution',edgecontrg
  endif
  if (detgrad) then
    ! gravity gradient
    edgecontrgg=edgecontrgg+Ee*Le
  endif
 
 enddo ! end loop on edges
 
 ! proceed with face contribution
 ! loop on element edges
 do i=1,NFaces

  ! determine vector rf a vector from the field point to
  ! some fixed point in the face plane (i.e. any point of the face)
!  write(stderr,*) 'x',x(1,lien(1)),'x0',X0
!  write(stderr,*) 'point',X0,Y0,Z0
  call SIDENP(i,.true.,lien)
  do j=1,NCornersFace
   rf(1,j) = x(1,lien(j))-X0
   rf(2,j) = x(2,lien(j))-Y0
   rf(3,j) = x(3,lien(j))-Z0
   rfl(j) = sqrt(rf(1,j)*rf(1,j)+rf(2,j)*rf(2,j)+rf(3,j)*rf(3,j))
!   write(stderr,*) 'x',x(:,lien(j)),lien(j)
!   write(stderr,*) 'length rf',rfl(j)
  enddo

  ! face dyad (outer product nf with itself)
  call outerproduct(nf(:,i),nf(:,i),Ff) 
!  write(stderr,*) 'Ff',Ff
  ! omf : the dimensionless perface factor
  call CROSSP (rf(:,2),rf(:,3),crossprodrf)
!  write(stderr,*) 'rf1',rf(:,1)
!  write(stderr,*) 'rf2',rf(:,2)
!  write(stderr,*) 'rf3',rf(:,3)
!  write(stderr,*) 'cross',crossprodrf
  omf = 2*atan2(dot_product(rf(:,1),crossprodrf),(rfl(1)*rfl(2)*rfl(3)+rfl(1)*dot_product(rf(:,2),rf(:,3))+rfl(2)*dot_product(rf(:,1),rf(:,3))+rfl(3)*dot_product(rf(:,1),rf(:,2))))
!  write(stderr,*) 'omf',omf
  if (detpot) then
    ! potential
      call VMPRD(Ff,rf(:,1),prodFFrf,3,3)
      facecontrp=facecontrp+dot_product(rf(:,1),prodFFrf)*omf
  endif
  if (detgrav) then
    ! gravitational attraction
      call VMPRD(Ff,rf(:,1),prodFFrf,3,3)
      facecontrg=facecontrg+prodFFrf*omf
 !     write(stderr,*) 'face contribution',facecontrg
  endif
  if (detgrad) then
      ! gravity gradient
      facecontrgg=facecontrgg+Ff*omf
  endif
 enddo

 ! now sum all contributions
 if (detpot) then
    Potential=0.5d0*grav_const*Rho*(edgecontrp-facecontrp)
 endif
 if (detgrav) then

    Gravity=grav_const*Rho*(facecontrg-edgecontrg)
!    write(stderr,*) 'gravity ',Gravity
!    write(stderr,*) 'grav_const',grav_const
!    write(stderr,*) 'Rho',Rho
!    write(stderr,*) 'face',facecontrg
!    write(stderr,*) 'edge',edgecontrg
!    write(stderr,*)
 endif
 if (detpot) then
    GravGradient=grav_const*Rho*(edgecontrgg-facecontrgg)
 endif
! 
!   if (test) then
!    do i=1,4
!      write(stderr,*) 'face normal',i,nf(1,i),nf(2,i),nf(3,i)
!      do j=1,3
!        write(stderr,*) 'edge normal',i,j,nef(1,j,i),nef(2,j,i),nef(3,j,i)
!      enddo
!    enddo
!   endif


 

 return
 end
!-----------------------------------------------------------------------
  subroutine getnormals(xl,nf,nef)
! computes face normals and edge per face normals for given element
  implicit none
  integer NSD,NEN
  parameter (NSD=3,NEN=4)
!-pass
  double precision nf(NSD,NEN),nef(NSD,3,NEN)
  double precision xl(NSD,NEN),xn(NSD),xne(NSD)
!-locl
  integer lien(3),k,l
  double precision area
 ! loop on each face
 do k=1,4
 ! Get local nodal point numbers of the side (lien)
  call SIDENP(k,.true.,lien)
!  calculate face normal vector nf
  call FCGEOM (xl,lien,area,xn)
  ! put in vector
  nf(1,k)=xn(1)
  nf(2,k)=xn(2)
  nf(3,k)=xn(3)
  
  ! calculate edge normal vector (per face) nef, each face has 3 edges
  do l=1,3
    call edgenormal(xl,xn,lien,l,xne)
    nef(1,l,k)=xne(1)
    nef(2,l,k)=xne(2)
    nef(3,l,k)=xne(3)
  enddo
 enddo
 return
 end
 
!-----------------------------------------------------------------------
 subroutine edgenormal(xl,xn,lien,edge,xne)
 ! computes the normal vector of an edge per face. 
 ! this vector is perpendicular to the edge and lies
 ! in the face plane, and points outwards
 
 ! to do so, take the cross product between face normal vector and side vector
 ! xn: face normal vector
 ! xne: edge normal vector, in plane of face, outward pointing
 ! edge: number of edge
 ! lien: face nodes
 ! xl: element coordinates
 implicit none
 integer NSD,NEN
 parameter (NSD=3,NEN=4)
!-pass
 integer edge
 integer lien(3)
 double precision xl(NSD,NEN),xn(NSD),xne(NSD)

!-locl
 integer node(2)
 integer i
 double precision X(NSD,2),e(NSD)
 double precision area

 ! select edge nodes 
 call edgenp(LIEN,edge,node)

 ! select node coordinates
 do i=1,2
   X(1,i) = XL(1,node(i))
   X(2,i) = XL(2,node(i))
   X(3,i) = XL(3,node(i))
 enddo

 ! edge vector 
 do i=1,3
   e(i)= X(i,2)-X(i,1)
 enddo 
 
 ! cross product edge vector with face normal vector 
 ! this should give edge perpendicular vector, in face plane,
 ! and outward pointing.
 ! outward pointingness is assured by taking edge vector point 2 - point 1
 
 xne(1) = e(2)*xn(3)-e(3)*xn(2)
 xne(2) = e(3)*xn(1)-e(1)*xn(3)  
 xne(3) = e(1)*xn(2)-e(2)*xn(1)
 
 ! normalize 
 area=sqrt(xne(1)*xne(1)+xne(2)*xne(2)+xne(3)*xne(3))
 ! note: area is not the area of the triangular face
 xne(1)=xne(1)/area
 xne(2)=xne(2)/area
 xne(3)=xne(3)/area
 
 return
 end
 
!-----------------------------------------------------------------------
 subroutine edgenp(LIEN,edge,node)

use iomodule,           only: stderr

 ! return nodes for a given side edge
 implicit none
!-pass
 integer LIEN(3)
 integer edge
 integer node(2)


 if (edge.eq.1) then
   node(1)=lien(1)
   node(2)=lien(2)
 else if (edge.eq.2) then
   node(1)=lien(2)
   node(2)=lien(3) 
 else if (edge.eq.3) then
   node(1)=lien(3)
   node(2)=lien(1) 
 else
      write(stderr,1) edge
 1        format(1x,'edgenp: edge number outside range [1,3]:',I5)
      return 
 endif
 
 
 return
 end
 
!-----------------------------------------------------------------------
     subroutine geoc2cart (Lon,Lat,height,RE,X,Y,Z)
! convert geocentric coordinates cartesian 
! geocentric coordinates are to be given in degrees
     USE CONSTANTS, only: deg2rad
     implicit none
! !-pass
     double precision X,Y,Z
     double precision Lon,Lat,RE,height,r
!- comm

!-locl
     integer i
     double precision xylen,latr,lonr,rcoslat
! !-init
      r = RE + height
      latr = lat*deg2rad
      lonr = lon*deg2rad
      z = r*sin(latr)
      rcoslat = r*cos(latr)
      x = rcoslat*cos(lonr);
      y = rcoslat*sin(lonr);
!      write(*,*) sqrt(x*x+y*y+z*z)
!      write(*,*) 'x,y,z',x,y,z
!      write(*,*) r,RE,height
      return
      end
      
!-----------------------------------------------------------------------
  subroutine EDGES(edge,side,localnodes,localedges)

use iomodule,           only: stderr

  ! gives neighbouring sides and local nodes for each edge
  implicit none
!-pass
  integer edge
  integer side(2)
  integer localnodes(2)
  integer localedges(2)

! edge numbering follows from side numbering and node numbering
! edge | nodes | side/local edge/local node
! 1    : 1-2   |   1 / 1 / 1-2
! 2    : 2-4   |   1 / 2 / 2-3
! 3    : 4-1   |   1 / 3 / 3-1
! 4    : 1-3   |   2 / 1 / 1-2
! 5    : 3-2   |   2 / 2 / 2-3
! 6    : 4-3   |   3 / 2 / 2-3
! 

! localnodes gives local nodes that belong to edge in side(1)
! localedges gives local edge in side(1) and side(2)

  if (edge.eq.1) then
  ! 1-2
    side(1)=1
    localnodes(1)=1
    localnodes(2)=2
    localedges(1)=1 ! local edge of side(1)
    side(2)=2
    localedges(2)=3 ! local edge of side(2) (reverse point order: 2-1)
  else if (edge.eq.2) then
  ! 2-4
    side(1)=1
    localnodes(1)=2
    localnodes(2)=3
    localedges(1)=2
    side(2)=4
    localedges(2)=3
  else if (edge.eq.3) then
  ! 4-1
    side(1)=1
    localnodes(1)=3
    localnodes(2)=1
    localedges(1)=3
    side(2)=3
    localedges(2)=1
  else if (edge.eq.4) then
  ! 1-3
    side(1)=2
    localnodes(1)=1
    localnodes(2)=2
    localedges(1)=1
    side(2)=3
    localedges(2)=3
  else if (edge.eq.5) then
  ! 3-2
    side(1)=2
    localnodes(1)=2
    localnodes(2)=3
    localedges(1)=2
    side(2)=4
    localedges(2)=1
  else if (edge.eq.6) then
  ! 4-3
    side(1)=3
    localnodes(1)=2
    localnodes(2)=3
    localedges(1)=2
    side(2)=4
    localedges(2)=2
  else
      write(stderr,1) edge
 1        format(1x,'edges: edge number outside range [1,6]:',I5)
      return
    endif

return
end

!-----------------------------------------------------------------------
!  function outerprod(a,b)
  subroutine outerproduct(a,b,outerprod)
  ! outer product for vectors of length 3
  implicit none
  !-pass 
  double precision, dimension(3), intent(in) ::  a
  double precision, dimension(3), intent(in) :: b
  double precision, dimension(3,3), intent(inout) ::  outerprod
  
  !-locl
  integer i,j
  do j=1,3
   do i=1,3
!    write(stderr,*) 'i',i,'j',j
!    write(stderr,*) 'a',a(i),'b',b(j)
    outerprod(i,j)=a(i)*b(j)
   enddo
  enddo 

  return
  end

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine addmat(a,b,sumab)
  ! add two 3x3 matrices
  implicit none
  !-pass 
  double precision, dimension(3,3), intent(in) ::  a
  double precision, dimension(3,3), intent(in) ::  b
  double precision, dimension(3,3), intent(inout) ::  sumab
  
  !-locl
  integer i,j
!  write(stderr,*) 'a',a
!  write(stderr,*) 'b',b
  do j=1,3
   do i=1,3
    sumab(i,j) = 0d0
!write(stderr,*) sumab(i,j)
!    write(stderr,*) i,j
!    write(stderr,*) a(i,j)
!    write(stderr,*) b(i,j)
    sumab(i,j)=a(i,j)+b(i,j)
!    write(stderr,*) sumab(i,j)
   enddo
  enddo 

  return
  end
  
!-----------------------------------------------------------------------
  subroutine edgediad(nfa,nefa,nfb,nefb,EE)
  ! determine EE, the edge diad
  implicit none
!-pass
  double precision nfa(3),nefa(3),nfb(3),nefb(3)
  double precision, intent(inout), dimension(3,3) :: EE

!-locl
  double precision parta(3,3),partb(3,3)
!  double precision, external :: addmat,outerprod
!  write(stderr,*) 'in edgediad function'
!  write(stderr,*) 'nfa',nfa
!  write(stderr,*) 'nefa',nefa
  call outerproduct(nfa,nefa,partA)
  call outerproduct(nfb,nefb,partB)
  !partA = outerprod(nfa,nefa)
  !partB = outerprod(nfb,nefb)
!  write(stderr,*) 'partA and B',partA,partB
  ! add the two parts
!   write(stderr,*) 'edgediad'
  call  addmat(partA,partB,EE)
  !write(stderr,*) edgediad
 ! write(stderr,*) EE
  return
  end
!-------------------------------------------------------------------------------
      subroutine xyz2enu (X,Lon,Lat,E,N,U)
      ! transform unit vectors in global cartesian to east, north, up directions

      implicit none
! -pass
      double precision, intent(in) :: Lon,Lat,X
      double precision, intent(out) :: E,N,U
      dimension X(3)
! -locl
      real sinlon,sinlat,coslon,coslat

      ! equations, see:
      ! http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
        sinlon=sin(Lon)
        sinlat=sin(Lat)
        coslon=cos(Lon)
        coslat=cos(Lat)
        E = -sinlon*X(1) + coslon*X(2)
        N = -coslon*sinlat*X(1) - sinlon*sinlat*X(2) + coslat*X(3)
        U = coslon*coslat*X(1)  + sinlon*coslat*X(2)  + sinlat*X(3)
       return
      end
!-------------------------------------------------------------------------------
subroutine gravityOfDomain(point, gravityVec, Potential, GravGradient)

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use materialsmodule, only: matdatactx
use plotControl,     only: detpot,detgrav,detgrad

implicit none


! pass
double precision :: point(3)   ! compute gravity here
double precision :: gravityVec(3)  ! x, y, z
double precision :: Potential
double precision :: GravGradient(3,3)


!local
double precision :: PotentialOfOneElem
double precision :: GravGradientOfOneElem(3,3)
double precision :: gravityVecOfOneElem(3)


double precision :: rho
integer          :: iElem

! computes the gravity of the complete mesh at a certain point.
! by calling WernerScheeres for every element


gravityVec = 0d0
Potential = 0d0
GravGradient = 0d0

do iElem=1,meshdatactx%neglobal

    rho = matdatactx%PRPMAT(5,meshdatactx%mat(iElem))

    gravityVecOfOneElem = 0d0

    ! .false.,.true.,.false.   =  not compute potential, do compute gravity, do not compute gravity gradient

    call gravityofelement(meshdatactx%IEN, &
                          meshdatactx%X, &
                          modeldatactx%D, &
                          modeldatactx%TFAULT, &
                          iElem, &
                          Rho,modeldatactx%STR, &
                          point(1),point(2),point(3), &
                          detpot,detgrav,detgrad, &
                          PotentialOfOneElem,gravityVecOfOneElem,GravGradientOfOneElem, &
                          .false.,.false.)

    gravityVec   = gravityVec   + gravityVecOfOneElem
    Potential    = Potential    + PotentialOfOneElem
    GravGradient = GravGradient + GravGradientOfOneElem
enddo

end subroutine

