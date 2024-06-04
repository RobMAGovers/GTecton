subroutine SYDATA ()
 
#ifdef SPARSE
!               Define structs for sparse solution of mechanical pdes.
#else
!               Computes stiffness matrix topology and bandwidth
#endif

USE MODELDATAMODULE, only: modeldatactx, &
                           idiag
USE MODELDEFINITION, only: MODE, iword, NEQlocal, NEQ, &
                           NEQext, MTOT, NUMNP, &
                           NSIZEA
USE MODELTOPOLOGY,   only: ndof, nen
USE ALGEBRA,         only: clear
use debugmodule,     only: allocateError
#ifdef SPARSE
USE AOMODULE,        only: nequations, nvertices, &
                           equationsmask, equationsmaskx, verticesmask
#endif
use iomodule, only: outfil, luout
implicit none


integer ERROR
integer, ALLOCATABLE :: d_nz(:)
integer, ALLOCATABLE :: o_nz(:)
integer :: IPNNZ
integer :: meanbn
 
#ifdef SPARSE
  if (MODE.eq.1 .or. MODE.eq.2 .or. MODE.eq.6) then

    ALLOCATE(d_nz(NEQlocal), STAT=ERROR)
    call AllocateError("d_nz",ERROR)

    ALLOCATE(o_nz(NEQlocal), STAT=ERROR)
    call AllocateError("o_nz",ERROR)

    d_nz = 0
    o_nz = 0

    ! everything in local per processor numbering
    call SetNonZeroPattern(d_nz, o_nz,NDOF,NEN)
    call InitializeMatrix(d_nz,o_nz,NEN,NDOF,0)

    DEALLOCATE(d_nz, STAT=ERROR)
    DEALLOCATE(o_nz, STAT=ERROR)

    if(allocated(equationsmask)) then
        deallocate(equationsmask)
    endif

    if(allocated(equationsmaskx)) then
        deallocate(equationsmaskx)
    endif

!    call CLEAR (A(IPNNZ),NEQ+1,"IPNNZ)")
  endif

#else 

  NSIZEA = 0
  meanbn = 0

    call DIAG(IDIAG,NEQ,NSIZEA)
 
  meanbn = NSIZEA/NEQ
  if (OUTFIL(1)) WRITE(luout(1),1) NEQ,NSIZEA,meanbn
1   format(/// &
  ' E Q U A T I O N   S Y S T E M   D A T A               ',  //5X, &
  ' NUMBER OF EQUATIONS  . . . . . . . . . . . . . (NEQ) =',I12//5X, &
  ' NUMBER OF TERMS IN STIFFNESS . . . . . . . . . (NA ) =',I12//5X, &
  ' MEAN HALF BANDWIDTH. . . . . . . . . . . . . . (MB ) =',I12)
 

#endif /* SPARSE */

#ifdef SPARSE

  ! nequations is the GLOBAL number of equations. So B is available for every equation throughout the domain
  ALLOCATE(modeldatactx%B(nequations), STAT=ERROR)
  call AllocateError("B",ERROR)
  modeldatactx%B = 0d0
  
  ALLOCATE(modeldatactx%BRES(nequations), STAT=ERROR)
  call AllocateError("BRES",ERROR)
  modeldatactx%BRES = 0d0
  
  ALLOCATE(modeldatactx%BTOT(nequations), STAT=ERROR)
  call AllocateError("BTOT",ERROR)
  modeldatactx%BTOT = 0d0
  
  NEQext = nequations

  if (allocated(verticesmask)) then
    deallocate(verticesmask)
  endif

#else

  ALLOCATE(modeldatactx%B(NEQlocal), STAT=ERROR)
  call AllocateError("B",ERROR)
  
  ALLOCATE(modeldatactx%BRES(NEQlocal), STAT=ERROR)
  call AllocateError("BRES",ERROR)
  
  ALLOCATE(modeldatactx%BTOT(NEQlocal), STAT=ERROR)
  call AllocateError("BTOT",ERROR)
  
  modeldatactx%BTOT = 0d0
  modeldatactx%BRES = 0d0
  modeldatactx%B = 0d0 
  
  NEQext = NEQlocal
  
#endif /* SPARSE */


return
end subroutine
!-------------------------------------------------------------------------------
#ifndef SPARSE
subroutine DIAG (IDIAG,NEQ,NA)
!
!    Changes column heights into diagonal addresses
!
implicit none

integer :: IDIAG,NEQ,NA
integer :: i


dimension IDIAG(NEQ)
!-init
NA = 1
IDIAG(1) = 1
if (NEQ.eq.1) then
    return
endif

do i=2,NEQ
    IDIAG(i) = IDIAG(i) + IDIAG(i-1) + 1
enddo

NA = IDIAG(NEQ)

return
end
#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
 subroutine SETAIJ (I,J,IA,JA,NEQ,MAXJA)
!
! builds IA and JA, assuming that the matrix topology is symmetric
!
!    IA   - integer one-dimensional array containing pointers to delimit
!           rows in JA and A.  Dimension = NEQ+1
!
!    JA   - integer one-dimensional array containing the column indices
!           corresponding to the elements of (storage array) A.
!           Dimension = number of nonzero entries in the full matrix M. In
!        case of a symmetric M, only upper triangle storage is indexed.
!
!  GENERAL STORAGE SYSTEM OF SPARSE MATRICES (Yale Sparse Matrix Format or
!  Compressed Row Storage)
!
!    The nonzero entries of the matrix M are stored row-by-row in the
!    array A. To identify the individual nonzero entries in each row,
!    we need to know in which column each entry lies. These column
!    indices are stored in the array JA i.e., if  A(k) = M(i,j),  then
!    JA(k) = j. To identify the individual rows, we need to know where
!    each row starts. These row pointers are stored in the array IA.
!    i.e., if M(i,j) is the first nonzero entry (stored) in the i-th row
!    and  A(k) = M(i,j),  then  IA(i) = k.  moreover, IA(n+1) points to
!    the first location following the last element in the last row.
!    thus, the number of entries in the i-th row is  IA(i+1) - IA(i),
!    the nonzero entries in the i-th row are stored consecutively in
!
!            A(IA(i)),  A(IA(i)+1),  ..., A(IA(i+1)-1),
!
!    and the corresponding column indices are stored consecutively in
!
!            JA(IA(i)), JA(IA(i)+1), ..., JA(IA(i+1)-1).
!
!    when the coefficient matrix is symmetric, only the nonzero entries
!    in the upper triangle need be stored.  for example, the matrix
!
!             ( 1  0  2  3  0 )
!             ( 0  4  0  0  0 )
!         M = ( 2  0  5  6  0 )
!             ( 3  0  6  7  8 )
!             ( 0  0  0  8  9 )
!
!    could be stored as
!
!            | 1  2  3  4  5  6  7  8  9 10 11 12 13
!         ---+--------------------------------------
!         IA | 1  4  5  8 12 14
!         JA | 1  3  4  2  1  3  4  1  3  4  5  4  5
!          A | 1  2  3  4  2  5  6  3  6  7  8  8  9
!
!    or (symmetrically) as
!
!            | 1  2  3  4  5  6  7  8  9
!         ---+--------------------------
!         IA | 1  4  5  7  9 10
!         JA | 1  3  4  2  3  4  4  5  5
!          A | 1  2  3  4  5  6  7  8  9          .
!
 
use debugmodule, only: xit
use iomodule,    only: stderr

 implicit none
!-pass
 integer I,J,IA,JA,NEQ,MAXJA
 dimension IA(*),JA(MAXJA)

!-locl
 integer k,m
!
 if (I.ge.J) return
!
 if (IA(I+1).gt.IA(I)) then
k = IA(I)
do while (k.le.IA(I+1)-1)
if (JA(k).eq.J) return
k = k + 1
enddo
if (J.lt.JA(IA(I))) then
if (IA(NEQ+1).ge.MAXJA) goto 100
do k=IA(NEQ+1),IA(I),-1
    if (k.gt.0) JA(k+1) = JA(k)
enddo
JA(IA(I)) = J
else if (J.gt.IA(I+1)-1) then
if (IA(NEQ+1).ge.MAXJA) goto 100
do k=IA(NEQ+1),IA(I+1),-1
    if (k.gt.0) JA(k+1) = JA(k)
enddo
JA(IA(I+1)) = J
else
m = IA(I)
do while (m.lt.IA(I+1)-1.and.(J-JA(m))*(J-JA(m+1)).gt.0)
    m = m + 1
enddo
if (IA(NEQ+1).ge.MAXJA) goto 100
do k=IA(NEQ+1),m+1,-1
    if (k.gt.0) JA(k+1) = JA(k)
enddo
JA(m+1) = J
endif
 else
if (IA(NEQ+1).ge.MAXJA) goto 100
do k=IA(NEQ+1),IA(I),-1
if (k.gt.0) JA(k+1) = JA(k)
enddo
JA(IA(I)) = J
 endif
 do k=I+1,NEQ+1
IA(k) = IA(k) + 1
 enddo
!
 return
!
100    write(stderr,1) MAXJA
 1   format(/1x,'SETAIJ: need more temporary memory than ',I9)
 call xit(1," ")
 end
#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
subroutine SETNZ_NEIGHBOURS(IA,JA,myJA,ID,IDX,doidx,NDOF, &
              nextnum,eqnno,numnb,num,node,NEQ,MAXJA)

implicit none
!-pass
integer :: NDOF,nextnum,eqnno,numnb,num,node
integer :: NEQ, MAXJA
logical :: doidx
integer :: IA(NEQ+1),JA(MAXJA),myJA(numnb),ID(NDOF,*),IDX(NDOF,*)
!-locl
integer :: m,k
integer :: nxtnb

IA(eqnno) = nextnum

do m=1,NDOF
    eqnno = 0
    if (doidx) then
        eqnno = IDX(m,node)
    endif
    if (eqnno.gt.0) then
       JA(nextnum) = eqnno
       nextnum = nextnum + 1
    endif
    eqnno = ID(m,node)
    if (eqnno.gt.0) then
        JA(nextnum) = eqnno
        nextnum = nextnum + 1
    endif
enddo

do k=1,numnb
    nxtnb = myJA(k)
    do m=1,NDOF
        eqnno = 0
        if (doidx) then
            eqnno=IDX(m,nxtnb)
        endif
        if (eqnno.gt.0) then
            JA(nextnum) = eqnno
            nextnum = nextnum + 1
        endif
        eqnno=ID(m,nxtnb)
        if (eqnno.gt.0) then
            JA(nextnum) = eqnno
            nextnum = nextnum + 1
        endif
    enddo
enddo

return
end
#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
subroutine SetNonZeroNeighbours(d_nz, o_nz, NEQlocal,&
                                node, i, next, doslip, NDOF, NEN)

! knowing the nonzero pattern and passing this on to PETSc
! gives a significant (1 to 2 orders of magnitude) improvement in performance. Hence 
! we make the effort to compute it.

! it is called once for every point of this partition, for every degree of freedom.
! it is called      for vertex i                 

   USE MESHDATAMODULE
   USE MODELDATAMODULE
   USE AOMODULE
   USE MODELCTX

implicit none

   integer :: NEQlocal
   integer NDOF, NEN, node, i , next
   logical doslip
   integer d_nz(NEQlocal)
   integer o_nz(NEQlocal)
   integer jj, k, l, nb, eqnno


   next = next + 1

if(next.gt.size(d_nz,1)) then
    return
endif


   d_nz(next) = 0
   o_nz(next) = 0

   ! first own node
   do jj=1,NDOF
  eqnno = 0

if (doslip) then
    eqnno = equationsmaskx(jj,node)
    if (eqnno.ne.0) then

!       write(*,*) "accessing ", next, "of ", NEQlocal


           d_nz(next) = d_nz(next) + 1

 !                  write(*,*) 'a slip) d_nz', next, ' increased to', d_nz(next)
 
!                  if (d_nz(next).gt.NEQglobal) then
 !                      d_nz(next) = NEQglobal
 !                  endif

      endif
  endif

  eqnno = equationsmask(jj,node)

  if (eqnno.ne.0) then

       d_nz(next) = d_nz(next) + 1

 !                  if (d_nz(next).gt.NEQglobal) then
 !                      d_nz(next) = NEQglobal
 !                  endif

!               write(*,*) 'a noslip) d_nz', next, ' increased to', d_nz(next)
  endif
   enddo

   ! do neighbors

!write(*,*) "size of itot: ", size(meshdatactx%itot,1)
!write(*,*) "size of AdjM: ", size(meshdatactx%AdjM,1), size(meshdatactx%AdjM,2)

do k=1,meshdatactx%itot(i)
  nb = meshdatactx%AdjM(k,i)
  ! a neighbour in own range of nodal points
  if (nb.le.meshdatactx%Nvlocal) then
      do l=1,NDOF
          eqnno = 0

          if (doslip) then
              eqnno = equationsmaskx(l,vertices(nb))

              if (eqnno.ne.0) then
                  d_nz(next) = d_nz(next) + 1
!                          if (d_nz(next).gt.NEQglobal) then
!                              d_nz(next) = NEQglobal
!                          endif

!                          write(*,*) 'b slip) d_nz', next, ' increased to', d_nz(next)
              endif

          endif

          eqnno = equationsmask(l,vertices(nb))
          if (eqnno.ne.0) then
              d_nz(next) = d_nz(next) + 1
!                      if (d_nz(next).gt.NEQglobal) then
!                           d_nz(next) = NEQglobal
!                      endif

!                      write(*,*) 'b noslip) d_nz', next, ' increased to', d_nz(next)
          endif
      enddo

  ! a neighbour in off diagonal block
  else
      do l=1,NDOF
          eqnno = 0
          if (doslip) then

              eqnno = equationsmaskx(l,vertices(nb))
              if (eqnno.ne.0) then
                  o_nz(next) = o_nz(next) + 1
!                          write(*,*) 'c slip) o_nz', next, ' increased to', o_nz(next)
              endif
          endif

          eqnno = equationsmask(l,vertices(nb))
          if (eqnno.ne.0) then
              o_nz(next) = o_nz(next) + 1
          endif
      enddo
  endif
   enddo




   return
   end
#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
   subroutine SetNonZeroPattern(d_nz, o_nz,NDOF,NEN)

! knowing the nonzero pattern and passing this on to PETSc
! gives a significant improvement in performance. Hence 
! we make the effort to compute it.


   USE MESHDATAMODULE
   USE MODELDEFINITION
   USE MODELDATAMODULE
   USE AOMODULE
   USE MODELCTX
use iomodule

implicit none

   integer NDOF,NEN, meanb, meanbn
   integer i, j,jj,k,l,next, node, eqnno, m
   integer d_nz(NEQlocal)
   integer o_nz(NEQlocal)
   logical doslip
! variables for correcting for linked nodes
   logical indexFound

   integer masterEQ,slaveEQ

!      do i=1,modeldatactx%NEQlocal

  next = 0
   doslip = (NUMSLP.gt.0)

!      write(*,*) 'SetNonZeroPattern says: equationsmask', equationsmask

do i=1,meshdatactx%Nvlocal
    node = vertices(i)
    do j=1,NDOF
 
       eqnno = 0

       if (doslip) then
          eqnno =  equationsmaskx(j,node)

          if (eqnno.gt.0) then
              ! 'next' is being increased in this sub
              call SetNonZeroNeighbours(d_nz, o_nz, NEQlocal,&
                   node, i, next, doslip, NDOF, NEN)
          endif

      endif

      ! check if the node is a linked slave node.
      ! if so, no equation is added here.

      m = INLINK(j, node)

      if (m.gt.0) then
!                  write(*,*) 'hi, a linked node', node, m, dof 
      endif

!              write(*,*) '******', i, j, m, '**************'

      eqnno = equationsmask(j,node)

      if (eqnno.gt.0 .and. m.lt.0) then  
          ! m lt 0 -> not linked node
          call SetNonZeroNeighbours(d_nz, o_nz, NEQlocal,&
            node, i, next, doslip, NDOF, NEN)
      endif


  enddo
!          write(*,*) 'after point ', i, 'd_nz',d_nz
   enddo

   ! modify for linked nodes. A master node gets extra nonzeros from it slave nodes.
   ! Currently, all the nodes have there number of own equations determined.
   ! However, the problem is that it is not known how many nonzero a d.o.f.
   ! outside this partition has.
   ! sit in corner and cry and just set to 6 for 2D (three nodes with two equations each, 
   ! or 12 for 3D (four nodes with three d.o.f.'s each). Just make sure we do not make them 
   ! larger than the total number of equations.
   
!      write(*,*) 'hoi, NEQlocal: ', NEQlocal


   do i=1,NLINK
!          write(*,*) 'tecdat says: dof',modeldatactx%LINK(1,i), &
!                     'of node ',  modeldatactx%LINK(2,i), &
!                     'linked to ', modeldatactx%LINK(3,i)

  ! link has the global vertex numbers
  ! we have to translate this to local, to
  ! find the matching equation number.

!----------------------------------- first attempt
!          if (meshdatactx%v2p(modeldatactx%LINK(3,i)).eq.getrank()) then
!              write(*,*) 'yay, we have a master node in this partition' 

!              if (meshdatactx%v2p(modeldatactx%LINK(2,i)).eq.getrank()) then

!                  masterEQ = equationsmask(modeldatactx%LINK(1,i),modeldatactx%LINK(3,i))

!                  ! for the slave idea, the number is not required.
!!                  write(*,*) 'yay, we have the slave node also in this partition, increase d_nz'
!                  if (NDOF.eq.2) then
!                      d_nz(masterEQ) = d_nz(masterEQ) + 6
!                  else if (NDOF.eq.3) then
!                      d_nz(masterEQ) = d_nz(masterEQ) + 12
!                  else
!                      write(*,*) 'NDOF should be 2 or 3. This should not happen. Contact model support'
!                  endif

!                  if (d_nz(masterEQ).gt.NEQlocal) then
!                      d_nz(masterEQ) = NEQlocal
!                  endif

!              else
!                  if (getrank().gt.1) then
!                  ! slave node is in another partition
!                      if (NDOF.eq.2) then
!                          o_nz(masterEQ) = o_nz(masterEQ) + 6
!                      else if (NDOF.eq.3) then
!                          o_nz(masterEQ) = o_nz(masterEQ) + 12
!                      else
!                          write(*,*) 'NDOF should be 2 or 3. This should not happen. Contact model support'
!                      endif

!                      if (o_nz(masterEQ).gt.NEQglobal-NEQlocal) then
!                          o_nz(masterEQ) = NEQglobal-NEQlocal
!                      endif
!                  endif
!              endif
!          endif
!------------------------------------

!   do j=1,NEQlocal
!       just to be sure... really... that is bad
!       d_nz(j) = d_nz(j) + 6
!       o_nz(j) = o_nz(j) + 6
!   enddo
   enddo

   do i=1,NEQlocal
  if (d_nz(i).gt.NEQlocal) then
      d_nz(i) = NEQlocal  
  endif
  if (getrank().gt.1) then
      if (o_nz(i).gt.NEQlocal) then
          o_nz(i) = NEQlocal
      endif
  endif
   enddo


!      write(*,*) 'after accomodating linked nodes: d_nz: ', d_nz

   NSIZEA = 0

   do i=1,NEQlocal
  NSIZEA = NSIZEA + d_nz(i)
  NSIZEA = NSIZEA + o_nz(i)
   enddo

   if (OUTFIL(1)) then
  meanbn = NSIZEA
!
!   bandwidth / 2
!   <------->
!  (4 1 1 0 1 0 0 0 )
!  (1 4 1 1 0 1 0 0 )
!  (1 1 4 1 1 0 1 0 )
!  (0 1 1 4 1 1 0 1 )
!  (1 0 1 1 4 1 1 0 )
!  (0 1 0 1 1 4 1 1 )
!  (0 0 1 0 1 1 4 1 )
!  (0 0 0 1 0 1 1 4 )

  meanbn = NINT( DBLE(meanbn+1)/DBLE(2*NEQlocal) )
  write(luout(1),2) NEQlocal,NSIZEA,meanbn
 2       format(///1X, &
     'E Q U A T I O N   S Y S T E M   D A T A',//6X, &
     'NUMBER OF EQUATIONS  . . . . . . . . . . . . . (NEQ) =', &
     I12//6X, &
     'NUMBER OF TERMS IN STIFFNESS . . . . . . . . . (NA ) =', &
     I12//6X, &
     'MEAN HALF BANDWIDTH. . . . . . . . . . . . . . (MB ) =', &
     I12)
   endif
!

!      write(*,*) 'SetNonZeroPattern says: d_nz: ', d_nz


   return
   end
#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
   subroutine SETNZPWD(IA,JA,IEN,ID,IDX,NDOF,NEN,MAXJA)
   USE MODELDEFINITION
   USE ALGEBRA
use debugmodule, only: iecho, xit
use iomodule
!use petscksp
!
! this defines the CRS (compressed row storage) for the matrix based on a partition graph computed by MeTis.

implicit none

!-pass
   integer NDOF, nen
integer :: MAXJA
integer :: IA, JA, ID, IDX, IEN
   dimension IA(*),JA(MAXJA),ID(NDOF,*),IDX(NDOF,*),IEN(NEN,*)


#include "petsc/finclude/petscsys.h" 
! for PetscMPIInt and PetscErrorCode

!-locl
   logical there
   integer,external :: NINT
   PetscMPIInt irank
   PetscErrorCode ierr
   integer i,j,jj,k,m
   integer, allocatable :: nxadj(:),nadjncy(:),elmts(:)
   integer nextnum, num, numnb, nextnb, etype, numflag, eqnno,ERROR
   logical doidx
integer :: ios, meanbn
!-init

write(*,*) "setnzpwd ntype: ", NTYPE



   if (MAXJA.lt.NEQ) then
  goto 1100
   endif

   if (iecho.ne.0) then
 write(stdout,1,advance='no')
 1    format(1x,'Establishing stiffness matrix topology ... ')
   endif
!

   doidx = .false.
   ALLOCATE (nxadj(NUMNP+1), STAT=ERROR)
   if (ERROR.ne.0) then
  write(*,*) 'Could not allocate nxadj. Error ', ERROR
  call xit(1," ")
   endif

   if (ERROR.ne.0) then
  goto 1000
   endif

   numflag = 1  !fortran style numbering
   if (NTYPE.eq.1) then !triangle
  etype = 1 !triangles

  ALLOCATE (nadjncy(10*NUMNP), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate nadjncy. Error ', ERROR
      call xit(1," ")
  endif

  ALLOCATE (elmts(3*NUMEL), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate elmts. Error ', ERROR
      call xit(1," ")
  endif

  do i=1,NUMEL
      elmts((i-1)*3 + 1) = IEN(1,i)
      elmts((i-1)*3 + 2) = IEN(2,i)
      elmts((i-1)*3 + 3) = IEN(3,i)
  enddo

   else if (NTYPE.eq.6) then !tetrahedra
  etype = 2 !tetrahedra
  ALLOCATE (nadjncy(20*NUMNP), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate nadjncy. Error ', ERROR
      call xit(1," ")
  endif

  ALLOCATE (elmts(4*NUMEL), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate elmts. Error ', ERROR
      call xit(1," ")
  endif

  do i=1,NUMEL
      elmts((i-1)*4 + 1) = IEN(1,i)
      elmts((i-1)*4 + 2) = IEN(2,i)
      elmts((i-1)*4 + 3) = IEN(3,i)
      elmts((i-1)*4 + 4) = IEN(4,i)
  enddo

   else if (NTYPE.eq.5) then !hexahedra
  etype = 3 !hexahedra
  ALLOCATE (nadjncy(10*NUMNP), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate nadjncy. Error ', ERROR
      call xit(1," ")
  endif

  ALLOCATE (elmts(8*NUMEL), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate elmts. Error ', ERROR
      call xit(1," ")
  endif

  do i=1,NUMEL
      elmts((i-1)*8 + 1) = IEN(1,i)
      elmts((i-1)*8 + 2) = IEN(2,i)
      elmts((i-1)*8 + 3) = IEN(3,i)
      elmts((i-1)*8 + 4) = IEN(4,i)
      elmts((i-1)*8 + 5) = IEN(5,i)
      elmts((i-1)*8 + 6) = IEN(6,i)
      elmts((i-1)*8 + 7) = IEN(7,i)
      elmts((i-1)*8 + 8) = IEN(8,i)
  enddo

   else
  write(stderr,31) NTYPE
31     format(1x,'SETNZP: NTYPE = ',I1,' not supported')
  call xit(1," ")
   endif

   if (ERROR.ne.0) goto 1000

   call METIS_MeshToNodal(NUMEL, NUMNP, elmts, etype, numflag, &
     nxadj, nadjncy)

!
!     now reform the METIS graph data into one suited for gtecton.
!     This takes care of boundary conditions and the fact that there are
!     multiple equations per node

   call CLEAR(IA,NEQ+1,"IA")

   doidx = (NUMSLP.ne.0)
   nextnum=1

   do i=1,NUMNP
  num=nxadj(i)!+1
  numnb = nxadj(i+1) - nxadj(i)

  do j=1,NDOF
      eqnno =0
      if (doidx) then
          eqnno = IDX(j,i)
      endif
      if (eqnno.gt.0) then
          call SETNZ_NEIGHBOURS(IA,JA,nadjncy(num:num+numnb-1), &
             ID,IDX,doidx,NDOF, &
             nextnum,eqnno,numnb,num,i,NEQ,MAXJA)
      endif
      eqnno = ID(j,i)
      if (eqnno.gt.0) then
          call SETNZ_NEIGHBOURS(IA,JA,nadjncy(num:num+numnb-1), &
             ID,IDX,doidx,NDOF, &
             nextnum,eqnno,numnb,num,i,NEQ,MAXJA)
      endif
  enddo

   enddo

   IA(NEQ+1) = nextnum

950   NSIZEA = IA(NEQ+1)-1

   if (OUTFIL(1)) then
  meanbn = NSIZEA
!
!   bandwidth / 2
!   <------->
!  (4 1 1 0 1 0 0 0 )
!  (1 4 1 1 0 1 0 0 )
!  (1 1 4 1 1 0 1 0 )
!  (0 1 1 4 1 1 0 1 )
!  (1 0 1 1 4 1 1 0 )
!  (0 1 0 1 1 4 1 1 )
!  (0 0 1 0 1 1 4 1 )
!  (0 0 0 1 0 1 1 4 )
  meanbn = NINT( DBLE(meanbn+1)/DBLE(2*NEQ) )
  write(luout(1),2) NEQ,NSIZEA,meanbn
 2       format(///1X, &
     'E Q U A T I O N   S Y S T E M   D A T A',//6X, &
     'NUMBER OF EQUATIONS  . . . . . . . . . . . . . (NEQ) =', &
     I12//6X, &
     'NUMBER OF TERMS IN STIFFNESS . . . . . . . . . (NA ) =', &
     I12//6X, &
     'MEAN HALF BANDWIDTH. . . . . . . . . . . . . . (MB ) =', &
     I12)
   endif
!

   return
!
1000  write(stderr,4) ERROR
14   format(/1x,'SETNZP: error in allocating arrays ',I9)
   call xit(1," ")

1100  write(stderr,4) MAXJA
 4   format(/1x,'SETNZP: need more temporary memory than ',I9)
   call xit(1," ")

1200  write(stderr,5) ios
 5   format(/1x,'SETNZP: read error ',I4,' on "MSPARSE.DAT"')
   call xit(1," ")

1300  write(stderr,6)
 6   format(/1x,'SETNZP: premature EOF error on "MSPARSE.DAT"')
   call xit(1," ")

1500  write(stderr,8) ios
 8   format(/1x,'SETNZP: write error ',I4,' on "SPARSE.DAT"')
   call xit(1," ")

   end
#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
subroutine SETNZP (IA,JA,IEN,ID,LINK,LM,LMX, &
                   NDOF,NEN,MAXJA)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: inlink, &
                           modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: inlink, &
                           modeldatactx
#endif

USE MODELDEFINITION
USE ALGEBRA
use debugmodule,     only: iecho, xit
use iomodule
!
! first defines the complete non-zero topology of the matrix to infer
! the number of non-zeroes on each row from it. Returns the result in IA
! JA is an array for intermediate storage.

implicit none

!-pass
integer :: MAXJA
integer NDOF,NEN!,NUMEL
integer :: IA, JA, LM, LMX, IEN, ID, LINK
dimension IA(*),JA(MAXJA),LM(NDOF,NEN,*),LMX(NDOF,NEN,*), &
   IEN(NEN,*),ID(NDOF,*),LINK(3,*)

!#include "petsc/finclude/petscsysdef.h" 
! for PetscMPIInt and PetscErrorCode
!-locl
integer, parameter :: PETSC_COMM_WORLD = 0
logical there
integer, external :: NINT
PetscMPIInt irank
PetscErrorCode ierr
integer myIA(NEQ+1),myJA(MAXJA),cntr(NEQ),i,ii,j,jj,k,m
integer :: ios
integer :: kode
integer :: lu
integer, external :: nextlu
integer :: meanbn
integer :: mumeq, numeq
integer :: node
!-init
if (MAXJA.lt.NEQ) then
    goto 1100
endif
!
if (iecho.ne.0) then
    write(stdout,1,advance='no')
endif
 1   format(1x,'Establishing stiffness matrix topology ... ')

!
call MPI_Comm_rank (PETSC_COMM_WORLD,irank,ierr)

inquire (file='MTOPO.DAT',exist=there)

if (there) then
    lu = nextlu(0)
    call fbopen(lu,'MTOPO.DAT','old')
    read(lu,err=1200,end=1300,iostat=ios)(IA(i),i=1,NEQ)
    call closef(lu)
    NSIZEA = 0
    do i=1,NEQ
        NSIZEA = NSIZEA + IA(i)
    enddo
    goto 950
endif

call CLEAR (JA,MIN(NEQ*NEQ,MAXJA),"JA")

do k=1,NEQ+1
    IA(k) = k
enddo

do k=1,NEQ
    JA(k) = k
enddo
!
do k=1,NUMEL
    do j=1,NEN
        node = IEN(j,k)
        do i=1,NDOF
!lwi linktest
            m = INLINK (i,node)
!            write(*,*) ''

            if (m.gt.0) then
                numeq = ID(i,LINK(3,m))
            else
                numeq = LM(i,j,k)
            endif

            if (numeq.eq.0) then
                goto 300
            endif

            do jj=1,NEN
                kode = IEN(jj,k)
                do ii=1,NDOF

!lwi linktest
                    m = INLINK (ii,kode)
                    if (m.gt.0) then
                        mumeq = ID(ii,LINK(3,m))
                    else
                        mumeq = LM(ii,jj,k)
                    endif

                    if (mumeq.gt.0) then
                        call SETAIJ(numeq,mumeq,IA, JA,NEQ,MAXJA)
                    endif

                    if (NUMSN.eq.0) then
                        cycle
                    endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
                    if(elementSlipperyID(k,90).gt.0) then
                        mumeq = IABS(LMX(ii,jj,elementSlipperyID(k,90)))
                    else
                        mumeq=0 
                    endif
#else
                    mumeq = IABS(LMX(ii,jj,k))
#endif                    
! original here mumeq = IABS(LMX(ii,jj,k))

                    if (mumeq.gt.0) then
                        call SETAIJ(numeq,mumeq,IA, &
                                    JA,NEQ,MAXJA)
                    endif

                enddo
            enddo

300         if (NUMSN.eq.0) then 
                cycle
            endif
#ifdef EXPERIMENTAL_ReducedSlipAlloc
            if(elementSlipperyID(k,91).gt.0) then
                numeq = IABS(LMX(i,j,elementSlipperyID(k,91)))
            else
                numeq=0
            endif
#else
            numeq = IABS(LMX(i,j,k))
#endif
            if (numeq.eq.0) then
                cycle
            endif
            do jj=1,NEN
                kode = IEN(jj,k)
                do ii=1,NDOF

!lwi linktest
                    m = INLINK (ii,kode)
                    if (m.gt.0) then
                        mumeq = ID(ii,LINK(3,m))
                    else
                        mumeq = LM(ii,jj,k)
                    endif

                    if (mumeq.gt.0) then
                        call SETAIJ(numeq,mumeq,IA, &
                                    JA,NEQ,MAXJA)
                    endif
#ifdef EXPERIMENTAL_ReducedSlipAlloc
                    if(elementSlipperyID(k,92).gt.0) then
                        mumeq = IABS(LMX(ii,jj,elementSlipperyID(k,92)))
                        if (mumeq.gt.0) then  ! is this check still necessary?
                            call SETAIJ(numeq,mumeq,IA, &
                                        JA,NEQ,MAXJA)
                        endif
                    endif
#else
                    mumeq = IABS(LMX(ii,jj,k))
                    if (mumeq.gt.0) then
                        call SETAIJ(numeq,mumeq,IA, &
                                    JA,NEQ,MAXJA)
                    endif
#endif
                enddo
            enddo
        enddo
    enddo
enddo

!
NSIZEA = 2*IA(NEQ+1)-NEQ-2
! WIENAND here write non symmetric crs for later use
!          write(stdout,*) "IA(NEQ): ", IA(NEQ),
!     >     " and IA(NEQ+1): ", IA(NEQ+1)
do i=1,NEQ+1
    myIA(i) = 0
enddo

do i=1,MAXJA
    myJA(i) = 0
enddo

do i=1,NEQ
    cntr(i) = 0
enddo

do i=1,NEQ
    myIA(i)=myIA(i)+IA(i)
    do j=IA(i)+1,IA(i+1)-1
        jj = JA(j)
        do k=jj+1,NEQ+1
            myIA(k)=myIA(k) + 1
        enddo
    enddo
enddo

myIA(NEQ+1)=myIA(NEQ+1)+IA(NEQ+1)

do i=1,NEQ
    j=IA(i+1)
    jj=j-IA(i)

    do k=1,jj
        myJA(myIA(i+1)-k) = JA(j-k)
    enddo
enddo

do i=1,NEQ
    do j=myIA(i)+1,myIA(i+1)-1
        if(myJA(j).gt.i) then
            myJA(myIA(myJA(j)) + cntr(myJA(j))) = i;
            cntr(myJA(j)) = cntr(myJA(j)) + 1
        endif
    enddo
enddo

if (irank.eq.0) then
    lu =nextlu(0)
    inquire (file='MSPARSE.DAT',exist=there)
    if (there) then
        call fbopen(lu, 'MSPARSE.DAT', 'old')
!              close(unit=lu,status='delete')
    endif
    call fbopen(lu, 'MSPARSE.DAT', 'new')
!          write(stdout,*) "writing sparsity file"
!          write(stdout,*) "IA(NEQ): ", IA(NEQ),
!     >     " and IA(NEQ+1): ", IA(NEQ+1)
    write(lu,err=1500,iostat=ios) (myIA(i),i=1,NEQ+1)
    write(lu,err=1500,iostat=ios) (myJA(i),i=1,myIA(NEQ+1))
    call closef(lu)
endif

123 format(I5)
! WIENAND end write
do i=1,NEQ
    JA(IA(i)) = 0
enddo

do i=1,NEQ-1
    do jj=IA(i)+1,IA(i+1)-1
        j = JA(jj)
        JA(IA(j)) = JA(IA(j)) + 1
    enddo
enddo

do i=1,NEQ
    IA(i) = JA(IA(i)) + IA(i+1)-IA(i)
enddo

call CLEAR (JA,MAXJA,"JA")

if (irank.eq.0) then
    lu = nextlu(0)
    call fbopen(lu,'MTOPO.DAT','new')
    write(lu,err=1400,iostat=ios) (IA(i),i=1,NEQ)
    call closef(lu)
endif
!
950 if (OUTFIL(1)) then
    meanbn = 0
    do i=1,NEQ
        meanbn = meanbn + IA(i)
    enddo
!
!   bandwidth / 2
!   <------->
!  (4 1 1 0 1 0 0 0 )
!  (1 4 1 1 0 1 0 0 )
!  (1 1 4 1 1 0 1 0 )
!  (0 1 1 4 1 1 0 1 )
!  (1 0 1 1 4 1 1 0 )
!  (0 1 0 1 1 4 1 1 )
!  (0 0 1 0 1 1 4 1 )
!  (0 0 0 1 0 1 1 4 )
   meanbn = NINT( DBLE(meanbn+1)/DBLE(2*NEQ) )
   write(luout(1),2) NEQ,NSIZEA,meanbn
 2       format(///1X, &
     'E Q U A T I O N   S Y S T E M   D A T A',//6X, &
     'NUMBER OF EQUATIONS  . . . . . . . . . . . . . (NEQ) =', &
     I12//6X, &
     'NUMBER OF TERMS IN STIFFNESS . . . . . . . . . (NA ) =', &
     I12//6X, &
     'MEAN HALF BANDWIDTH. . . . . . . . . . . . . . (MB ) =', &
     I12)
endif

call MPI_Barrier(PETSC_COMM_WORLD,ierr)

return
!
1100    write(stderr,4) MAXJA
 4   format(/1x,'SETNZP: need more temporary memory than ',I9)
 call xit(1," ")

1200    write(stderr,5) ios
 5   format(/1x,'SETNZP: read error ',I4,' on "MTOPO.DAT"')
 call xit(1," ")

1300    write(stderr,6)
 6   format(/1x,'SETNZP: premature EOF error on "MTOPO.DAT"')
 call xit(1," ")

1400    write(stderr,7) ios
 7   format(/1x,'SETNZP: write error ',I4,' on "MTOPO.DAT"')
 call xit(1," ")

1500    write(stderr,8) ios
 8   format(/1x,'SETNZP: write error ',I4,' on "SPARSE.DAT"')
 call xit(1," ")

 end
#endif
!-------------------------------------------------------------------------------

#ifdef SPARSE
subroutine InitSolve (IFLAG)

USE MODELDEFINITION
use modelctx,    only: getrank, &
                       getsize
use debugmodule, only: iecho
use petscksp
use spetscmodule, only: M_A, M_init, M_sles, &
                        T_A, T_init, T_sles
!
! initiates PETSc solver context, sets preconditioner and Krylov subspace method
! IFLAG=0: mechanical
! IFLAG=1: thermal
!
!
!
!
! possible solver types:
!#define KSPRICHARDSON "richardson"
!#define KSPCHEBYSHEV  "chebyshev"
!#define KSPCG         "cg"
!#define KSPGROPPCG    "groppcg"
!#define KSPPIPECG     "pipecg"
!#define   KSPCGNE       "cgne"
!#define   KSPNASH       "nash"
!#define   KSPSTCG       "stcg"
!#define   KSPGLTR       "gltr"
!#define KSPGMRES      "gmres"
!#define   KSPFGMRES     "fgmres"
!#define   KSPLGMRES     "lgmres"
!#define   KSPDGMRES     "dgmres"
!#define   KSPPGMRES     "pgmres"
!#define KSPTCQMR      "tcqmr"
!#define KSPBCGS       "bcgs"
!#define   KSPIBCGS      "ibcgs"
!#define   KSPFBCGS      "fbcgs"
!#define   KSPFBCGSR     "fbcgsr"
!#define   KSPBCGSL      "bcgsl"
!#define KSPCGS        "cgs"
!#define KSPTFQMR      "tfqmr"
!#define KSPCR         "cr"
!#define KSPPIPECR     "pipecr"
!#define KSPLSQR       "lsqr"
!#define KSPPREONLY    "preonly"
!#define KSPQCG        "qcg"
!#define KSPBICG       "bicg"
!#define KSPMINRES     "minres"
!#define KSPSYMMLQ     "symmlq"
!#define KSPLCD        "lcd"
!#define KSPPYTHON     "python"
!#define KSPGCR        "gcr"
!#define KSPSPECEST    "specest"



implicit none

#include "petsc/finclude/petscsys.h" 
! for PETSC_TRUE and PETSC_COMM_WORLD

#include "petsc/finclude/petscmat.h"
! for SAME_NONZERO_PATTERN

#include "petsc/finclude/petscpc.h"
! for type PC

#include "petsc/finclude/petscksp.h"
! for type KSP


!



!-locl
PC   MyPc
KSP  MyKsp
MatType mattype
KSPType ksptype
PCType pctype
PetscErrorCode ierr
PetscReal tol

integer :: iflag

data tol/1.0D-07/
!
if (IFLAG.eq.0) then  ! iflag 0: mechanical
                    ! iflag 1: thermal

    if (M_Init) then

        if (iecho.eq.5) then
            write(*,*) 'Solver already exists, just set initial guess'
        endif

        ! the check for M_init ensures that the solver M_SLES already exists.
        ! before m_init, it is not there yet.
        ! Because BRES is used for (temporary) storage:

        if (IRESDU.eq.0) then
            call KSPSetInitialGuessNonzero (M_SLES, PETSC_TRUE,ierr)
        endif

    else

        if (iecho.eq.5) then
            write(*,*) 'Solver does not exists, create'
        endif

        ! make the solver context
        call KSPCreate (PETSC_COMM_WORLD,M_SLES,ierr)

        ! PCSOR is the only preconditioner that works with MPIAIJ matrix type :\
        ! TODO: experiment with other matrix types that may be compatible with better preconditioners.

        call KSPGetPC (M_SLES,MyPc,ierr)

        if (getsize().eq.1) then
            ! PCLU is effectively a direct solve,
            ! using a LU-computation and elimation.
            ! Even for runs on many processors, PCLU on 
            ! one single processor could still be faster.
            ! so only use multiple processor runs if memory is an issue.
            ! call PCSetType (MyPc,PCLU,ierr)
            call PCSetType (MyPc,PCLU,ierr)
!            call PCSetType (MyPc,PCNONE,ierr)

        else
            ! PCSOR is a very bad preconditioner.
            ! better is to overwrite it using the HYPRE preconditioner,
            ! by adding to petsc_options.txt, the line:
            ! -pc_type hypre -pc_hypre_type pilut
            ! or one of the other hypre types:
            ! parasails, boomeramg, euclid
             call PCSetType (MyPc,PCSOR,ierr)
!             call PCSetType (MyPc,PCNONE,ierr)
            ! call PCSetType (MyPc,PCILU,ierr)
        endif

!       make ifdef MUMPS + ifndef PETSC232
#ifdef MUMPS
#ifndef PETSC232
!    write(*,*) 'setting solver type to mumps'
!      call PCFactorSetMatSolverPackage(MyPC, MAT_SOLVER_MUMPS, ierr)
#endif
#endif

!write(*,*) "setting tolerances to", tol, &
!PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER

        call KSPSetTolerances (M_SLES,&       ! solver
          tol, &                              ! relative tol
          PETSC_DEFAULT_REAL, &
          PETSC_DEFAULT_REAL, &
          99999999, &                         ! maximum iterations
          ierr)

        if (ierr.ne.0) then
            write(*,*) "ERROR: Failed to set PETSc tolerances"
            stop "Leaving"
        endif

        call KSPSetType(M_SLES,KSPCG,ierr)
!        call KSPSetType(M_SLES,KSPRICHARDSON, ierr)   <- less accurate. do not use

        call KSPSetFromOptions (M_SLES,ierr)

        M_Init = .true.

        call MatGetType(M_A, mattype,ierr)
        call KSPGetType (M_SLES,ksptype,ierr)
        call PCGetType (MyPc,pctype,ierr)

    endif

else if (IFLAG.eq.1) then

    if (T_Init) then

        if (iecho.eq.5) then
            write(*,*) 'Thermal solver already exists, just set initial guess'
        endif

        ! the check for M_init ensures that the solver M_SLES already exists.
        ! before m_init, it is not there yet.
        ! Because BRES is used for (temporary) storage:

        if (IRESDU.eq.0) then
            call KSPSetInitialGuessNonzero (T_SLES, PETSC_TRUE,ierr)
        endif

    else

        if (iecho.eq.5) then
            write(*,*) 'Thermal solver does not exists, create'
        endif

        ! make the solver context
        call KSPCreate (PETSC_COMM_WORLD,T_SLES,ierr)

        ! PCSOR is the only preconditioner that works with MPIAIJ matrix type :\
        ! TODO: experiment with other matrix types that may be compatible with better preconditioners.

        call KSPGetPC (T_SLES,MyPc,ierr)

        if (getsize().eq.1) then
            ! PCLU is effectively a direct solve,
            ! using a LU-computation and elimation.
            ! Even for runs on many processors, PCLU on 
            ! one single processor could still be faster.
            ! so only use multiple processor runs if memory is an issue.
            ! call PCSetType (MyPc,PCLU,ierr)

            call PCSetType (MyPc,PCLU,ierr)
        else
            ! PCSOR is a very bad preconditioner.
            ! better is to overwrite it using the HYPRE preconditioner,
            ! by adding to petsc_options.txt, the line:
            ! -pc_type hypre -pc_hypre_type pilut
            ! or one of the other hypre types:
            ! parasails, boomeramg, euclid
            ! call PCSetType (MyPc,PCSOR,ierr)
            call PCSetType (MyPc,PCNONE,ierr)   
            ! call PCSetType (MyPc,PCILU,ierr)   ! not available in parallel
        endif

!       make ifdef MUMPS + ifndef PETSC232
#ifdef MUMPS
#ifndef PETSC232
!    write(*,*) 'setting solver type to mumps'
!      call PCFactorSetMatSolverPackage(MyPC, MAT_SOLVER_MUMPS, ierr)
#endif
#endif

!write(*,*) "setting tolerances to", tol, &
!PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER


        call KSPSetTolerances (T_SLES,&       ! solver
          tol, &                              ! relative tol
          PETSC_DEFAULT_REAL, &
          PETSC_DEFAULT_REAL, &
          99999999, &                         ! maximum iterations
          ierr)

        if (ierr.ne.0) then
            write(*,*) "ERROR: Failed to set PETSc tolerances"
            stop "Leaving"
        endif

!        call KSPSetType(T_SLES,KSPGMRES,ierr)
!        call KSPSetType(T_SLES,KSPBCGS,ierr)
!        call KSPSetType(T_SLES,KSPCG,ierr)

        if (IADVEC.eq.0) then
            call KSPSetType(T_SLES,KSPCG,ierr)
        else
            ! GMRES makes an enormous difference for advection models!
            call KSPSetType(T_SLES,KSPGMRES,ierr)
        endif

        call KSPSetFromOptions (T_SLES,ierr)

        T_Init = .true.

        call MatGetType(T_A, mattype,ierr)
        call KSPGetType (T_SLES,ksptype,ierr)
        call PCGetType (MyPc,pctype,ierr)

    endif


endif


!
!     Effective upon successive calls:
if (IFLAG.eq.0) then
    call KSPSetOperators (M_SLES,M_A,M_A,ierr)
else if (IFLAG.eq.1) then
    call KSPSetOperators (T_SLES,T_A,T_A,ierr)
endif 

! show the solver stuff
!      write(*,*) 'arg to kspview: ', PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_STDOUT_WORLD
!      call KSPView(M_SLES,PETSC_VIEWER_STDOUT_WORLD)

!      if (iecho.eq.5 .and. getrank().eq.0 .and m_init) then
!               write(*,*) '******** Solver properties: ***********'
!               call kspview(M_SLES, PETSC_VIEWER_STDOUT_SELF, ierr)
!          endif
!      endif



!
   return
   end subroutine

#endif
!SPARSE

!-------------------------------------------------------------------------------
#ifdef SPARSE
   subroutine PetscInit ()

USE MODELCTX,    only: initializemodelctx
use filemodule,  only: WORKPATH
use debugmodule, only: iecho, &
                       usedPetscOptionsTxt, xit
use spetscmodule, only: PetSc_Init
use petscksp

implicit none

! Initialized PETSc

#include "petsc/finclude/petscsys.h"   
! for the PETSC_NULL_CHARACTER

!#include "petsc/finclude/petscksp.h"

!#include "petsc/finclude/petscsys.h"
!#include "mpif.h"
!#include "mpif-mpi-io.h"


!-locl
!    Not accounted for possibilty that integer differes from PetscInt
!    PetscMPIInt getsize!nprocs,irank
PetscInt mloc
PetscErrorCode ierr
logical :: file_exists

usedPetscOptionsTxt = .false.

if (.not.PetSc_Init) then
  
inquire(file=trim(WORKPATH)//'/'//"petsc_options.txt", exist=file_exists)       

    if (file_exists) then
        if (iecho.eq.2) then
            write(*,*) 'Initializing using options from petsc_options.txt'
        endif
        usedPetscOptionsTxt = .true.
        call PetscInitialize (trim(WORKPATH)//'/'//"petsc_options.txt",ierr)
    else
        if (iecho.eq.2) then
            write(*,*) 'petsc_options.txt not found, initializing from default'
        endif
!        write(*,*) "Calling petsc initialize without file", PETSC_NULL_CHARACTER
        call PetscInitialize (PETSC_NULL_CHARACTER,ierr)
!        write(*,*) "Finished petsc initialize without file", ierr
    endif

    if (ierr.ne.0) then
        write(*,*) 'PetscInitialize failed with error: ', ierr
        call xit(1," ")
    endif

!    write(*,*) "Calling ini modelctx"
    call initializemodelctx()
!    write(*,*) "Finished ini modelctx"

    PetSc_Init = .true.
endif

return
end

#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
   subroutine GetMemoryAllocationWD(IA,JA,Istart,Iend,D_NZ,O_NZ, &
         NUMEQ,MAXJA,IFLAG)
! IFLAG = 0: mechanical
! IFLAG = 1: thermal

use debugmodule, only: xit

implicit none


!-pass
   integer IFLAG, NUMEQ, MAXJA, D_NZ, O_NZ
   PetscInt  Istart,Iend
   integer IA(NUMEQ+1), JA(MAXJA)
!-locl
   logical there
   integer tempdmax, maxnnz
   integer i,j
!
   if(IFLAG.eq.0) then
      D_NZ=0
      O_NZ=0
      maxnnz=0
      do i=Istart + 1,Iend
          tempdmax=0
          do j=IA(i),IA(i+1)-1
              if ((JA(j).le.Istart).or. &
                 (JA(j).gt.Iend)) then
                   tempdmax = tempdmax + 1;
              endif
          enddo
          O_NZ=max(O_NZ,tempdmax)
          D_NZ=max(D_NZ,IA(i+1) - IA(i) - tempdmax)
      enddo
!          endif
   else
      D_NZ=0
      O_NZ=0
      maxnnz=0
      do i=Istart + 1,Iend
          tempdmax=0
          do j=IA(i), IA(i+1)-1
              if ((JA(j).le.Istart).or. &
               (JA(j).gt.Iend))then
                  tempdmax = tempdmax + 1;
              endif
          enddo
          O_NZ=max(O_NZ,tempdmax)
          D_NZ=max(D_NZ,IA(i+1)-IA(i) - tempdmax)
      enddo
   endif
   return
   end
#endif

!-------------------------------------------------------------------------------
#ifdef SPARSE
   subroutine GetMemoryAllocation(Istart,Iend,D_NZ,O_NZ, &
         NEQ,MAXJA,IFLAG)
! IFLAG = 0: mechanical
! IFLAG = 1: thermal
use debugmodule, only: iecho, xit
use iomodule

implicit none


!-pass
   PetscInt  Istart,Iend
   integer IFLAG, NEQ, MAXJA, D_NZ, O_NZ
!-locl
integer :: ios
   logical there
   integer tempdmax, maxnnz
   integer i,j
integer :: RPOINTERS(NEQ+1),CMARKERS(MAXJA)
integer :: lu
integer, external :: nextlu
!
   if(IFLAG.eq.0) then
  inquire(file='MSPARSE.DAT', exist=there)
  if(there) then

      if (iecho.ne.0) then
          write(stdout,*) 'reading MSPARSE.DAT'
      endif

      lu = nextlu(0)

      call fbopen(lu,'MSPARSE.DAT','old')
      read(lu,err=1200,end=1300,iostat=ios) &
          (RPOINTERS(i),i=1,NEQ+1)
      read(lu,err=1200,end=1300,iostat=ios) &
          (CMARKERS(i),i=1,RPOINTERS(NEQ+1)-1)
      call closef(lu)

      D_NZ=0
      O_NZ=0
      maxnnz=0

      do i=Istart + 1,Iend
          tempdmax=0
          do j=RPOINTERS(i),RPOINTERS(i+1)-1
              if ((CMARKERS(j).le.Istart).or. &
                 (CMARKERS(j).gt.Iend)) then
                   tempdmax = tempdmax + 1;
              endif
          enddo
          O_NZ=max(O_NZ,tempdmax)
          D_NZ=max(D_NZ,RPOINTERS(i+1)-RPOINTERS(i) - tempdmax)
      enddo

  endif
   else
  inquire(file='TSPARSE.DAT', exist=there)
  if(there) then

      if (iecho.ne.0) then
          write(stdout,*) 'reading TSPARSE.DAT'
      endif

      lu = nextlu(0)
      call fbopen(lu,'TSPARSE.DAT','old')
      read(lu,err=1201,end=1301,iostat=ios) &
          (RPOINTERS(i),i=1,NEQ+1)
      read(lu,err=1201,end=1301,iostat=ios) &
          (CMARKERS(i),i=1,RPOINTERS(NEQ+1)-1)
      call closef(lu)

      D_NZ=0
      O_NZ=0
      maxnnz=0

      do i=Istart + 1,Iend
          tempdmax=0
          do j=RPOINTERS(i),RPOINTERS(i+1)-1
              if ((CMARKERS(j).le.Istart).or. &
               (CMARKERS(j).gt.Iend))then
                  tempdmax = tempdmax + 1;
              endif
          enddo
          O_NZ=max(O_NZ,tempdmax)
          D_NZ=max(D_NZ,RPOINTERS(i+1)-RPOINTERS(i) - tempdmax)
      enddo

  endif
   endif
   return

1200  write(stderr,5) ios
 5   format(/1x,'GetMemoryAllocation: read error ',I4, ' on MSPARSE.DAT')
   call xit(1," ")

1300  write(stderr,6)
 6   format(/1x,'GetMemoryAllocation: premature EOF error on MSPARSE.DAT')
   call xit(1," ")

1201  write(stderr,7) ios
 7   format(/1x,'GetMemoryAllocation: read error ',I4, ' on TSPARSE.DAT')
   call xit(1," ")

1301  write(stderr,8)
 8   format(/1x,'GetMemoryAllocation: premature EOF error on TSPARSE.DAT')
   call xit(1," ")

   end
#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
   subroutine InitializeMatrix(d_nz, o_nz,NEN,NDOF,IFLAG)
!      USE MODELDATAMODULE
   USE MODELDEFINITION
!      USE MESHDATAMODULE
   USE AOMODULE
USE MODELCTX,            only: getsize, &
                               getrank, file_outputf, parout
use spetscmodule,     only: M_A, M_X, M_B, M_Xseq, M_Bseq, M_init, M_sles, &
                            T_A, T_X, T_B, T_Xseq, T_Bseq, T_init, T_sles, &
                            PetSc_Init
use petscksp


#ifdef PLOTNONZEROES
   use debugmodule, only: debug, iecho, parallelLog, &
                          matrixEntryNonZero
#else
   use debugmodule, only: debug, iecho, parallelLog
#endif

implicit none


#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscksp.h"

!



integer NEN,NDOF, i,j,k,ns,nb,jj,kk, mmm,nnn
integer IFLAG
!      PetscMPIInt getsize
   PetscErrorCode ierr
   PetscInt d_nz(*)
   PetscInt o_nz(*)
   PetscInt row(3)
   PetscInt col(3)

if (IFLAG.eq.0) then
    ! initialize mechanical matrix
    goto 100
endif

if (IFLAG.eq.1) then
    ! initialize thermal matrix
    goto 200
endif


 100  if (debug) then
!  if (parout) then
!      write(FILE_outputf,*) 'start mechanical matrix assembly'
!  endif
   endif



call MatCreate(PETSC_COMM_WORLD, M_A, ierr)
if (ierr.ne.0) then
    write(*,*) "The mechanical matrix could not be initialized."
    write(*,*) "MatCreate returned error code", ierr
    stop "Leaving GTecton."
endif

!if (getsize().eq.1) then
     call MatSetType(M_A, MATAIJ, ierr)
!else
!    call MatSetType(M_A, MATSBAIJ , ierr)
!endif

if (ierr.ne.0) then
    write(*,*) "The type of the mechanical matrix could not be set."
    write(*,*) "MatSetType returned error code", ierr
    stop "Leaving GTecton."
endif


call MatSetSizes(M_A, NEQlocal,NEQlocal, NEQglobal, NEQglobal, ierr)
if (ierr.ne.0) then
    write(*,*) "The size of the mechanical matrix could not be set."
    write(*,*) "MatSetSizes returned error code", ierr
    stop "Leaving GTecton."
endif



#ifdef PLOTNONZEROES
allocate(matrixEntryNonZero(NEQglobal,NEQglobal))
do iEq=1,NEQglobal
    do jEq=1,NEQglobal
        matrixEntryNonZero(iEq, jEq) = .false.
    enddo
enddo
#endif

! if (getrank().eq.0) then
!    write(*,*) 'degrees of freedom: ', NEQglobal
! endif

if (iecho.eq.5) then
!    write(*,*) 'mechanical matrix size set to ', NEQglobal
!    write(*,*) 'of which rank ',getrank(),' has ', NEQlocal
endif

!      call MatSetSizes(M_A, PETSC_DECIDE, PETSC_DECIDE, NEQglobal, NEQglobal, ierr)

#ifdef MUMPS
!      write(*,*) 'MUMPS defined: matsettype'
!      write(*,*) 'M_A',M_A
!      write(*,*) 'MATAIJMUMPS',MATAIJMUMPS
!      call MatSetType(M_A, MATAIJMUMPS, ierr)
#else
!      write(*,*) 'no mumps, no type'
!      call MatSetType(M_A, MATAIJMUMPS, ierr)
#endif

!      call MatSetType(M_A, MATMPIMAIJ, ierr)

   call MatSetFromOptions(M_A, ierr)

   if (getsize().gt.1) then

  ! temporary fix to give room for the parallel nodes
  ! this is ugly, but it gets shit working
  if (NLINK.gt.0) then
      do i=1,NEQlocal
!                  d_nz(i) = d_nz(i) + 2*NLINK
!                  if (d_nz(i) > NEQlocal) then
!                      d_nz(i) = NEQlocal
!                  endif
!                  o_nz(i) = o_nz(i) + 2*NLINK
!                  if (o_nz(i) > NEQglobal - NEQlocal) then
!                      o_nz(i) = NEQglobal - NEQlocal
!                  endif
      enddo
  endif


  call MatMPIAIJSetPreallocation(M_A, 0, d_nz, 0, o_nz, ierr)

   else
!          if (iecho.eq.5) then

  ! temporary fix to give room for the parallel nodes
  ! this is ugly, but it gets shit working
  if (NLINK.gt.0) then
      do i=1,NEQlocal
                  d_nz(i) = d_nz(i) + 2*NLINK
                  if (d_nz(i) > NEQlocal) then
                      d_nz(i) = NEQlocal
                  endif
                  o_nz(i) = o_nz(i) + 2*NLINK
                  if (o_nz(i) > NEQlocal) then
                      o_nz(i) = NEQlocal
                  endif
      enddo
  endif
!               write(*,*) 'd_nz', d_nz(1:7)
!               do i=1,7
!                   d_nz(i) = 7
!               enddo
!          endif
  call MatSeqAIJSetPreallocation(M_A, 0, d_nz, ierr)
   endif

   call MatSetOption (M_A,MAT_SYMMETRIC,PETSC_TRUE,ierr)
if (ierr.ne.0) then
    write(*,*) "Setting M_A option symmetry returned an error:", ierr
    stop "Leaving Gtecton"
endif

   call MatSetOption (M_A,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr)
if (ierr.ne.0) then
    write(*,*) "Setting M_A option symmetry eternal returned an error:", ierr
    stop "Leaving Gtecton"
endif


   ! this prevents errors when our preallocation was insufficient.
   ! It should not be necessary, but for linked nodes, it still is.
   call MatSetOption(M_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
if (ierr.ne.0) then
    write(*,*) "Setting M_A option no new non zeroes returned an error:", ierr
    stop "Leaving Gtecton"
endif



   ! make sure the run does not crash when a new nonzero has been set.
   ! This occurs with slippery nodes, and should be handled nicely. TODO
!      call MatSetOption (M_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,pierr)


if (iecho.eq.8) then
    call ParallelLog("","Creating mechanical load vector")
endif


call VecCreate(PETSC_COMM_WORLD, M_B, ierr)
if (ierr.ne.0) then
    write(*,*) "Veccreate M_B returned an error:", ierr
    stop "Leaving Gtecton"
endif

call VecSetSizes(M_B, NEQlocal, NEQglobal, ierr)
if (ierr.ne.0) then
    write(*,*) "VecSetSize M_B returned an error:", ierr
    stop "Leaving Gtecton"
endif

call VecSetFromOptions(M_B,ierr)
if (ierr.ne.0) then
    write(*,*) "VecSetFromOptions M_B returned an error:", ierr
    stop "Leaving Gtecton"
endif

!call vecView(M_B,PETSC_VIEWER_STDOUT_WORLD , iErr)


call VecDuplicate(M_B, M_X, ierr)
if (ierr.ne.0) then
    write(*,*) "Vecduplicate M_B -> M_X returned an error:", ierr
    stop "Leaving Gtecton"
endif

call createlocaltoglobalmapping(M_A,M_B,0)

if (iecho.eq.8) then
    call ParallelLog("","Create mechanical load vector")
    write(*,*) "M_B has pointer ", M_B
endif


   goto 300

 200  if (debug) then
   if (parout) write(FILE_outputf,*) 'start thermal matrix assembly'
   endif


   call MatCreate(PETSC_COMM_WORLD, T_A, ierr)
   call MatSetSizes(T_A, NTEQlocal, &
                    NTEQlocal, &
                    NTEQglobal, &
                    NTEQglobal, &
                    ierr)

if (iecho.eq.5) then
    write(*,*) 'temperature matrix size set to ', NTEQglobal
    write(*,*) 'of which rank ',getrank(),' has ', NTEQlocal
endif


!if (getsize().eq.1) then
    call MatSetType(T_A, MATAIJ, ierr)
!else
!    call MatSetType(T_A, MATSBAIJ , ierr)
!endif



#ifdef MUMPS
!      call MatSetType(T_A, MATAIJMUMPS, ierr)
#endif

   call MatSetFromOptions(T_A, ierr)



   if (getsize().gt.1) then
  call MatMPIAIJSetPreallocation(T_A, 0, &
           d_nz, 0, o_nz, ierr)
   else
  call MatSeqAIJSetPreallocation(T_A,0, &
           d_nz,ierr)
   endif

! turn off new zero allocation error
!   call MatSetOption(T_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)

   if (IADVEC.eq.0) then
  call MatSetOption (T_A,MAT_SYMMETRIC,PETSC_TRUE,ierr)
  call MatSetOption (T_A,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr)
   else
  call MatSetOption (T_A,MAT_STRUCTURALLY_SYMMETRIC, &
                     PETSC_TRUE,ierr)
   endif

   call VecCreate(PETSC_COMM_WORLD, T_B, ierr)
!      write(*,*) 'VecSetSizes with: ', NTEQlocal
   call VecSetSizes(T_B, NTEQlocal, NTEQglobal, ierr)
   call VecSetFromOptions(T_B,ierr)
   call VecDuplicate(T_B, T_X, ierr)

   call createlocaltoglobalmapping(T_A,T_B,1)

 300  if (debug) then
  if (parout) then
      write(FILE_outputf,*) 'finished matrix assembly'
  endif
   endif

   return
   end
#endif


!-------------------------------------------------------------------------------
#ifdef SPARSE
 subroutine UserInitializeWD (IA,JA,NEQ,IFLAG,MAXJA)

use debugmodule, only: iecho, parallelLog
use spetscmodule,     only: M_A, M_X, M_B, M_Xseq, M_Bseq, M_init, M_sles, &
                            T_A, T_X, T_B, T_Xseq, T_Bseq, T_init, T_sles, &
                            PetSc_Init
use petscksp
use modeldefinition,  only: ITMODE,IADVEC,IDIFT
!
! Initializes PETSc, matrix and array structs
! IFLAG = 0: mechanical
! IFLAG = 1: thermal

implicit none


!



#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscsys.h"
!-pass
integer :: MAXJA
 integer IA(*),JA(MAXJA),NEQ,IFLAG
!-locl
!    Not accounted for possibilty that integer differes from PetscInt
 PetscMPIInt nprocs,irank
 PetscErrorCode ierr
   PetscInt  Istart,Iend
   integer i, max_dnnz, max_onnz, tempdmax, maxnnz
   logical there
!
   call MPI_Comm_rank (PETSC_COMM_WORLD,irank,ierr)
   call MPI_Comm_size (PETSC_COMM_WORLD,nprocs,ierr)

!#ifdef PETSC232
   if (IFLAG.eq.0) then
  call MatCreate(PETSC_COMM_WORLD, M_A, ierr)
  call MatSetSizes(M_A, PETSC_DECIDE, PETSC_DECIDE, NEQ, NEQ, ierr)
#ifdef MUMPS
!          call MatSetType(M_A, MATAIJMUMPS, ierr)
#endif
  call MatSetFromOptions(M_A, ierr)
!     WD: preallocate memory for the matrix
  if (nprocs.gt.1) then
      call MatGetOwnershipRange(M_A, Istart, Iend, ierr)
      call GetMemoryAllocationWD(IA,JA,Istart,Iend,max_dnnz, &
       max_onnz, NEQ,MAXJA,IFLAG)

      call MatMPIAIJSetPreallocation(M_A, max_dnnz, &
          PETSC_NULL_INTEGER,max_onnz, PETSC_NULL_INTEGER, ierr)

  else
      do i=1,NEQ
          IA(i) = IA(i+1)-IA(i)
      enddo

!              lu =nextlu(0)
!              open(lu,FILE="nonzerosRowNewMech.dat",STATUS='unknown')
!              do i=1,NEQ
!                  write(lu,*) "NZZ[", i,"]: ", IA(i)
!              enddo
!              close(lu)
      call MatSeqAIJSetPreallocation(M_A,0, &
          IA,ierr)
  endif
  call MatSetOption (M_A,MAT_SYMMETRIC,PETSC_TRUE,ierr)
  call MatSetOption (M_A,MAT_SYMMETRY_ETERNAL, &
                     PETSC_TRUE,ierr)
!     WD: preallocate memory for the vectors
  if (nprocs.gt.1) then
      call VecCreateMPI(PETSC_COMM_WORLD, Iend-Istart, &
          PETSC_DETERMINE, M_B, ierr)
  else
      call VecCreateSeq(PETSC_COMM_SELF,NEQ,M_B,ierr)
  endif

if (iecho.eq.8) then
    call parallelLog("UserInitializeWD","Made M_B")
    write(*,*) "M_B has pointer: ", M_B
endif

  call VecDuplicate (M_B,M_X,ierr)
  call VecCreateSeq(PETSC_COMM_SELF, NEQ, M_Bseq, ierr)
  call VecDuplicate(M_Bseq,M_Xseq,ierr)
   else
!           TODO WD: find a smart way to initialize the thermal matrix
  call MatCreate(PETSC_COMM_WORLD, T_A, ierr)
  call MatSetSizes(T_A, PETSC_DECIDE, PETSC_DECIDE, &
     NEQ, NEQ, ierr)
#ifdef MUMPS
!          call MatSetType(T_A, MATAIJMUMPS, ierr)
#endif
  if (nprocs.gt.1) then
      call MatGetOwnershipRange(T_A, Istart, Iend, ierr)
!              initialize mpi matrix
      call GetMemoryAllocationWD(IA,JA,Istart,Iend,max_dnnz, &
        max_onnz,NEQ,MAXJA,IFLAG)
      call MatMPIAIJSetPreallocation(T_A, max_dnnz, &
          PETSC_NULL_INTEGER,max_onnz, PETSC_NULL_INTEGER, ierr)
  else
      do i=1,NEQ
          IA(i) = IA(i+1)-IA(i)
      enddo
!              lu =nextlu(0)
!              open(lu,FILE="nonzerosRowNewTemp.dat",STATUS='unknown')
!              do i=1,NEQ
!                  write(lu,*) "NZZ[", i,"]: ", IA(i)
!              enddo
!              close(lu)
      call MatSeqAIJSetPreallocation(T_A,0, &
          IA,ierr)
!        call MatCreateSeqAIJ (PETSC_COMM_SELF,NEQ,NEQ,
!     >         0,NNZ,T_A,ierr)
  endif
  if (IADVEC.eq.0) then
      call MatSetOption (T_A,MAT_SYMMETRIC,PETSC_TRUE,ierr)
      call MatSetOption (T_A,MAT_SYMMETRY_ETERNAL, &
                         PETSC_TRUE,ierr)
  else
      call MatSetOption (T_A,MAT_STRUCTURALLY_SYMMETRIC, &
                         PETSC_TRUE,ierr)
  endif
!           TODO WD: find something to allocate storage smartly
  if (nprocs.gt.1) then
      call VecCreateMPI(PETSC_COMM_WORLD, Iend-Istart, &
          PETSC_DETERMINE, T_B, ierr)
  else
      call VecCreateSeq(PETSC_COMM_SELF,NEQ,T_B,ierr)
  endif
!        call VecCreateSeqWithArray (PETSC_COMM_SELF,NEQ,
!     >         PETSC_NULL_SCALAR,T_B,ierr)
  call VecDuplicate (T_B,T_X,ierr)
  call VecCreateSeq(PETSC_COMM_SELF,NEQ,T_Bseq, ierr)
  call VecDuplicate(T_Bseq,T_Xseq,ierr)
   endif

 return
 end
#endif
!SPARSE
!-------------------------------------------------------------------------------
#ifdef SPARSE
 subroutine UserInitialize (NNZ,NEQ,IFLAG,MAXJA)
!
! Initializes PETSc, matrix and array structs
! IFLAG = 0: mechanical
! IFLAG = 1: thermal

use spetscmodule, only: M_A, M_B, M_X, M_Bseq, M_Xseq, &
                        T_A, T_B, T_X, T_Bseq, T_Xseq
use petscksp
use modeldefinition, only: ITMODE,IADVEC,IDIFT

implicit none


!



!#include "petsc/finclude/petscsysdef.h"
!#include "petsc/finclude/petscmatdef.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscksp.h"

#include "petsc/finclude/petscmat.h"  
! for MAT_STRUCTURALLY_SYMMETRIC

#include "petsc/finclude/petscsys.h" 
! for PETSC_COMM_WORLD and PETSC_COMM_SELF

!-pass
 integer NNZ(*),NEQ,IFLAG
!-locl
!    Not accounted for possibilty that integer differes from PetscInt
 PetscMPIInt nprocs,irank
 PetscErrorCode ierr
   PetscInt  Istart,Iend
   integer i, max_dnnz, max_onnz, tempdmax, maxnnz
   logical there
integer :: maxja
   integer RPOINTERS(NEQ+1),CMARKERS(MAXJA)
!
   call MPI_Comm_rank (PETSC_COMM_WORLD,irank,ierr)
   call MPI_Comm_size (PETSC_COMM_WORLD,nprocs,ierr)
!#ifdef PETSC232
   if (IFLAG.eq.0) then
  call MatCreate(PETSC_COMM_WORLD, M_A, ierr)
  call MatSetSizes(M_A, PETSC_DECIDE, PETSC_DECIDE, &
         NEQ, NEQ, ierr)
#ifdef MUMPS
!          call MatSetType(M_A, MATAIJMUMPS, ierr)
#endif
  call MatSetFromOptions(M_A, ierr)
!     WD: preallocate memory for the matrix
  if (nprocs.gt.1) then
      call MatGetOwnershipRange(M_A, Istart, Iend, ierr)
      call GetMemoryAllocation(Istart,Iend,max_dnnz,max_onnz, &
          NEQ,MAXJA,IFLAG)
      call MatMPIAIJSetPreallocation(M_A, max_dnnz, &
          PETSC_NULL_INTEGER,max_onnz, PETSC_NULL_INTEGER, ierr)
  else
      call MatSeqAIJSetPreallocation(M_A,0, &
          NNZ,ierr)
  endif
  call MatSetOption (M_A,MAT_SYMMETRIC,PETSC_TRUE,ierr)
  call MatSetOption (M_A,MAT_SYMMETRY_ETERNAL, &
                     PETSC_TRUE,ierr)
!     WD: preallocate memory for the vectors
  if (nprocs.gt.1) then
      call VecCreateMPI(PETSC_COMM_WORLD, Iend-Istart, &
          PETSC_DETERMINE, M_B, ierr)
  else
      call VecCreateSeq(PETSC_COMM_SELF,NEQ,M_B,ierr)
  endif
  call VecDuplicate (M_B,M_X,ierr)
  call VecCreateSeq(PETSC_COMM_SELF, NEQ, M_Bseq, ierr)
  call VecDuplicate(M_Bseq,M_Xseq,ierr)
   else
!           TODO WD: find a smart way to initialize the thermal matrix
  call MatCreate(PETSC_COMM_WORLD, T_A, ierr)
  call MatSetSizes(T_A, PETSC_DECIDE, PETSC_DECIDE, &
     NEQ, NEQ, ierr)
#ifdef MUMPS
!          call MatSetType(T_A, MATAIJMUMPS, ierr)
#endif
  if (nprocs.gt.1) then
      call MatGetOwnershipRange(T_A, Istart, Iend, ierr)
!              initialize mpi matrix
      call GetMemoryAllocation(Istart,Iend,max_dnnz,max_onnz, &
          NEQ,MAXJA,IFLAG)
      call MatMPIAIJSetPreallocation(T_A, max_dnnz, &
          PETSC_NULL_INTEGER,max_onnz, PETSC_NULL_INTEGER, ierr)
  else
      call MatSeqAIJSetPreallocation(T_A,0, &
          NNZ,ierr)
!        call MatCreateSeqAIJ (PETSC_COMM_SELF,NEQ,NEQ,
!     >         0,NNZ,T_A,ierr)
  endif
  if (IADVEC.eq.0) then
      call MatSetOption (T_A,MAT_SYMMETRIC,PETSC_TRUE,ierr)
      call MatSetOption (T_A,MAT_SYMMETRY_ETERNAL, &
                         PETSC_TRUE,ierr)
  else
      call MatSetOption (T_A,MAT_STRUCTURALLY_SYMMETRIC, &
                         PETSC_TRUE,ierr)
  endif
!           TODO WD: find something to allocate storage smartly
  if (nprocs.gt.1) then
      call VecCreateMPI(PETSC_COMM_WORLD, Iend-Istart, &
          PETSC_DETERMINE, T_B, ierr)
  else
      call VecCreateSeq(PETSC_COMM_SELF,NEQ,T_B,ierr)
  endif
!        call VecCreateSeqWithArray (PETSC_COMM_SELF,NEQ,
!     >         PETSC_NULL_SCALAR,T_B,ierr)
  call VecDuplicate (T_B,T_X,ierr)
  call VecCreateSeq(PETSC_COMM_SELF,NEQ,T_Bseq, ierr)
  call VecDuplicate(T_Bseq,T_Xseq,ierr)
   endif

!
 return
 end
#endif
!SPARSE
!-------------------------------------------------------------------------------
!#ifdef SPARSE
! block data NITPETSC

! PetSc_Init = .false.
!M_Init=.false.
!T_Init=.false.

! end
!#endif
