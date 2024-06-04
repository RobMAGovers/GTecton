subroutine tempdat ()

USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELDEFINITION
USE MODELTOPOLOGY
USE MODELCTX
USE ALGEBRA
use debugmodule, only: iecho, parallelLog
use iomodule          
#ifdef SPARSE
USE AOMODULE,   only: equationsmaskt, &
                       equationsmaskxt, &
                       nequationst
#endif

#ifdef SPARSE
! Define structs for sparse solution of thermal pdes.
#else
! Compute bandwidth of stiffness matrix and allocate space in the main array
#endif
!
implicit none

!-locl
integer              :: meanbn
integer              :: ERROR
integer, ALLOCATABLE :: d_nz(:)
integer, ALLOCATABLE :: o_nz(:)

!
if (iecho.eq.8) then
     call ParallelLog("tempdat", "welcome")
endif

if (MODE.lt.4) then
    return
endif
!
#ifdef SPARSE
!

if (iecho.eq.8) then
     call ParallelLog("tempdat", "initializing thermal matrix")
endif


ALLOCATE(d_nz(NTEQlocal), STAT=ERROR)
ALLOCATE(o_nz(NTEQlocal), STAT=ERROR)
call SetNonZeroPatternThermal(d_nz, o_nz, NEN)

if(allocated(equationsmaskt)) then
    deallocate(equationsmaskt)
endif

if(allocated(equationsmaskxt)) then
    deallocate(equationsmaskxt)
endif



if (iecho.eq.8) then
     call ParallelLog("tempdat", "nonzero pattern determined")
endif

!        write(FILE_outputf, *) (d_nz(i),i=1,NTEQlocal)
!        write(FILE_outputf, *) (o_nz(i),i=1,NTEQlocal)
call InitializeMatrix(d_nz,o_nz,NEN,NDOF,1)
DEALLOCATE(d_nz, STAT=ERROR)
DEALLOCATE(o_nz, STAT=ERROR)

if (iecho.eq.8) then
     call ParallelLog("tempdat", "initialized thermal matrix")
endif


if (iecho.eq.8) then
     call ParallelLog("tempdat", "clearing thingy")
endif

ALLOCATE(modeldatactx%AA(NTEQlocal), STAT=ERROR)

#else
! not sparse 
NAA = 0
!        establish diagonal addresses in ITDIAG
call DIAG (itdiag,NTEQ,NAA)

 meanbn = NAA/NTEQ
 if (OUTFIL(1)) WRITE(luout(1),1) NTEQ,NAA,meanbn
 1   format(///1x, &
 'T H E R M A L   E Q U A T I O N   S Y S T E M   D A T A',//5X, &
 'NUMBER OF EQUATIONS  . . . . . . . . . . . .  (NTEQ) =',I12// &
 5X, &
 'NUMBER OF TERMS IN STIFFNESS . . . . . . . . . (NAA) =',I12// &
 5X, &
 'MEAN HALF BANDWIDTH. . . . . . . . . . . . . . . . . =',I12)
!
!    allocate space for storage of stiffness matrix
#endif
!    allocate space for load vector
#ifdef SPARSE
   ALLOCATE(modeldatactx%BT(nequationst), STAT=ERROR)
   NTEQext = nequationst
#else
   ALLOCATE(modeldatactx%BT(NTEQlocal), STAT=ERROR)
   NTEQext = NTEQlocal
#endif
!
 return
 end
!-------------------------------------------------------------------------------
#ifdef SPARSE
subroutine SetNonZeroNeighborsThermal(d_nz, o_nz, node, &
                                      i, next, doidtx)
USE MESHDATAMODULE
USE MODELDATAMODULE
USE AOMODULE
USE MODELCTX
use iomodule          

integer NEN, node, i , next
logical doidtx
integer d_nz(*)
integer o_nz(*)
integer k, nb, eqnno

next = next + 1
d_nz(next) = 0
o_nz(next) = 0
! first own node
eqnno = 0

if (doidtx) then
    eqnno = equationsmaskxt(node)
endif

if (eqnno.gt.0) then
    d_nz(next) = d_nz(next) + 1
endif

eqnno = equationsmaskt(node)

if (eqnno.gt.0) then
    d_nz(next) = d_nz(next) + 1
endif

   ! neighbors
do k=1,meshdatactx%itot(i)
    nb = meshdatactx%AdjM(k,i)

    ! neighbour in own range of nodal points
    if (nb.le.meshdatactx%Nvlocal) then

        eqnno = 0
        if (doidtx) then
            eqnno = equationsmaskxt(vertices(nb))
        endif
        if (eqnno.gt.0) then
            d_nz(next) = d_nz(next) + 1
        endif
        eqnno = equationsmaskt(vertices(nb))
        if (eqnno.gt.0) then
            d_nz(next) = d_nz(next) + 1
        endif
        ! neighbour in off diagonal block

    else
        eqnno = 0
        if (doidtx) then
            eqnno = equationsmaskxt(vertices(nb))
        endif
        if (eqnno.gt.0) then
            o_nz(next) = o_nz(next) + 1
        endif
        eqnno = equationsmaskt(vertices(nb))
        if (eqnno.gt.0) then
            o_nz(next) = o_nz(next) + 1
        endif
    endif
enddo

return
end
#endif

!-------------------------------------------------------------------------------
#ifdef SPARSE
subroutine SetNonZeroPatternThermal(d_nz, o_nz, NEN)

USE MESHDATAMODULE,  only: meshdatactx
USE MODELDEFINITION, only: numsn, numfn, idift, nAA, NTEQlocal
USE AOMODULE,       only: vertices, equationsmaskxt, equationsmaskt
USE MODELCTX,        only: getrank
use debugmodule,     only: iecho, parallelLog
use iomodule          

implicit none


integer NEN, meanbn
integer i, next, node, eqnno
integer d_nz(*)
integer o_nz(*)
logical doidtx

next = 0
doidtx =  (NUMSN+NUMFN.gt.0 .and.IDIFT.eq.1)

do i=1,meshdatactx%Nvlocal

    node = vertices(i)
    eqnno = 0

    if (doidtx) then
        eqnno = equationsmaskxt(node)
    endif

    if (eqnno.gt.0) then

        if (iecho.eq.8) then
             call ParallelLog("SetNonZeroPatternThermal", "set slippery neighbor nonzeroes")
        endif

        call SetNonZeroNeighborsThermal(d_nz, o_nz, node, i, &
           next, doidtx)
    endif

    eqnno = equationsmaskt(node)

    if (eqnno.gt.0) then

        if (iecho.eq.8) then
             call ParallelLog("SetNonZeroPatternThermal", "set normal neighbor nonzeroes")
        endif

        call SetNonZeroNeighborsThermal(d_nz, o_nz, node, i, &
           next, doidtx)

        if (iecho.eq.8) then
             call ParallelLog("SetNonZeroPatternThermal", &
                              "finished setting normal neighbor nonzeroes")
        endif


    endif
enddo


NAA = 0

do i=1,NTEQlocal
    NAA = NAA + d_nz(i)
    NAA = NAA + o_nz(i)
enddo

if (OUTFIL(1)) then
    meanbn = NAA
    meanbn = NINT( DBLE(meanbn+1)/DBLE(2*NTEQlocal) )
    write(luout(1),2) NTEQlocal,NAA,meanbn
 2       format(///1x, &
     'T H E R M A L   E Q U A T I O N   S Y S T E M   D A T A', &
     //5X, &
     'NUMBER OF EQUATIONS  . . . . . . . . . . . .  (NTEQ) =', &
     I12//5X, &
     'NUMBER OF TERMS IN STIFFNESS . . . . . . . . . (NAA) =', &
     I12//5X, &
     'MEAN HALF BANDWIDTH. . . . . . . . . . . . . . . . . =', &
     I12)
endif



return
end subroutine
#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE

subroutine SETNZTWD (IA,JA,IEN,IDT,IDTX, NEN,MAXJA)

USE MODELDEFINITION
USE ALGEBRA
use debugmodule, only: iecho, xit
use iomodule          

implicit none

!
!     defines the CRS for the thermal stiffness matrix based on graph from MeTiS
!
!-pass
integer NEN,MAXJA

integer :: IA(*),JA(MAXJA),IDT(*),IDTX(*),IEN(NEN,*)

#include "petsc/finclude/petscsys.h"


!-locl
integer              :: node,meanbn
integer, external    :: NINT
integer              :: i,k
integer, allocatable :: nxadj(:),nadjncy(:),elmts(:)
integer              :: nextnum, num, numnb, etype, numflag,eqnno,ERROR
integer              :: nxtnb
logical              :: doidtx

if (MAXJA.lt.NTEQ) then
    goto 1100
endif

if (iecho.ne.0) then
    write(stdout,1)
 1   format(1x,'Establishing thermal stiffness matrix topology ... ')
endif

doidtx = .false.
ALLOCATE (nxadj(NUMNP+1), STAT=ERROR)

if (ERROR.ne.0) then
    goto 1000
endif

numflag = 1  !fortran style numbering

if (NTYPE.eq.1) then !triangle
    etype = 1 !triangles

    ALLOCATE (nadjncy(10*NUMNP), STAT=ERROR)
    ALLOCATE (elmts(3*NUMEL), STAT=ERROR)

    do i=1,NUMEL
        elmts((i-1)*3 + 1) = IEN(1,i)
        elmts((i-1)*3 + 2) = IEN(2,i)
        elmts((i-1)*3 + 3) = IEN(3,i)
    enddo

else if (NTYPE.eq.6) then !tetrahedra

    etype = 2 !tetrahedra

    ALLOCATE (nadjncy(20*NUMNP), STAT=ERROR)
    ALLOCATE (elmts(4*NUMEL), STAT=ERROR)

    do i=1,NUMEL
        elmts((i-1)*4 + 1) = IEN(1,i)
        elmts((i-1)*4 + 2) = IEN(2,i)
        elmts((i-1)*4 + 3) = IEN(3,i)
        elmts((i-1)*4 + 4) = IEN(4,i)
    enddo

else if (NTYPE.eq.5) then !hexahedra

    etype = 3 !hexahedra

    ALLOCATE (nadjncy(10*NUMNP), STAT=ERROR)
    ALLOCATE (elmts(8*NUMEL), STAT=ERROR)

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

if (ERROR.ne.0) then
    goto 1000
endif

call METIS_MeshToNodal(NUMEL, NUMNP, elmts, etype, numflag, &
                       nxadj, nadjncy)


!     now reform the METIS graph data into one suited for gtecton
!     this takes care of boundary conditions and the fact there are
!     multiple equations per node

call CLEAR(IA,NEQ+1,"IA")

doidtx =  (NUMSN+NUMFN.gt.0.and.IDIFT.eq.1)
nextnum=1

do i=1,NUMNP
    num=nxadj(i)
    numnb = nxadj(i+1) - nxadj(i)
    eqnno = 0
    if (doidtx) then
        eqnno = IDTX(i)
    endif
    if (eqnno.gt.0) then
        IA(eqnno) = nextnum
        JA(nextnum) = eqnno
        nextnum = nextnum + 1

        do k=1,numnb
            nxtnb = nadjncy(num + k - 1)
            eqnno = 0

            if (doidtx) then
                eqnno = IDTX(nxtnb)
            endif

            if (eqnno.gt.0) then
                JA(nextnum) = eqnno
                nextnum = nextnum + 1
            endif

            eqnno = IDT(nxtnb)

            if (eqnno.gt.0) then
                JA(nextnum) = eqnno
                nextnum = nextnum + 1
            endif

        enddo
    endif
    eqnno = IDT(i)
    if (eqnno.gt.0) then

        IA(eqnno) = nextnum
        JA(nextnum) = eqnno
        nextnum = nextnum + 1

        do k=1,numnb
            nxtnb = nadjncy(num + k - 1)! + 1
            eqnno = 0

            if (doidtx) then
                eqnno = IDTX(nxtnb)
            endif

            if (eqnno.gt.0) then
                JA(nextnum) = eqnno
                nextnum = nextnum + 1
            endif

            eqnno = IDT(nxtnb)

            if (eqnno.gt.0) then
                JA(nextnum) = eqnno
                nextnum = nextnum + 1
            endif

        enddo
    endif
enddo

IA(NTEQ+1) = nextnum

NAA = IA(NTEQ+1)-1

if (OUTFIL(1)) then
    meanbn = NAA
    meanbn = NINT( DBLE(meanbn+1)/DBLE(2*NTEQ) )
    write(luout(1),2) NTEQ,NAA,meanbn
 2       format(///1x, &
     'T H E R M A L   E Q U A T I O N   S Y S T E M   D A T A', &
     //5X, &
     'NUMBER OF EQUATIONS  . . . . . . . . . . . .  (NTEQ) =', &
     I12//5X, &
     'NUMBER OF TERMS IN STIFFNESS . . . . . . . . . (NAA) =', &
     I12//5X, &
     'MEAN HALF BANDWIDTH. . . . . . . . . . . . . . . . . =', &
     I12)
endif

return


1000  write(stderr,4) ERROR
14   format(/1x,'SETNZT: error in allocating arrays ',I9)
call xit(1," ")

1100  write(stderr,4) MAXJA
 4   format(/1x,'SETNZT: need more temporary memory than ',I9)
call xit(1," ")
call xit(1," ")

end subroutine

#endif
!-------------------------------------------------------------------------------
#ifdef SPARSE
subroutine SETNZT (IA,JA,IEN,IDT,LMT,LMTX,NEN,MAXJA)

USE MODELDEFINITION
USE ALGEBRA
use debugmodule, only: iecho, xit
use iomodule          
use AOmodule,    only: PETSC_COMM_WORLD

!
! first defines the complete non-zero topology of the thermal stiffness
! matrix to infer the number of non-zeroes on each row from it. Returns the
! result in IA. JA is an array for intermediate storage.
!

implicit none

!-pass
integer NEN,MAXJA

integer     :: IA(*),JA(MAXJA),LMT(NEN,*),LMTX(NEN,*),IEN(NEN,*),IDT(*)

! dimension IA(*),JA(MAXJA),LMT(NEN,*),LMTX(NEN,*), &
!   IEN(NEN,*),IDT(*)


!-locl
logical          :: nodift,there
integer          :: k,j,node,numeq,jj,kode,mumeq,i,meanbn,lu, ios
integer,external :: NINT,nextlu
PetscMPIInt      :: irank
PetscErrorCode   :: ierr
integer          :: myIA(NTEQ+1),myJA(MAXJA),cntr(NTEQ)

!-init
if (MAXJA.lt.NTEQ) then
    goto 1100
endif

nodift = (NUMSLP+NUMFN.eq.0.or.IDIFT.eq.0)
!
if (iecho.ne.0) then
    write(stdout,1)
 1   format(1x,'Establishing thermal stiffness matrix topology ... ')
endif

call MPI_Comm_rank (PETSC_COMM_WORLD,irank,ierr)

inquire (file='TTOPO.DAT',exist=there)
if (there) then
    lu = nextlu(0)

    call fbopen(lu,'TTOPO.DAT','old')
    read(lu,err=1200,end=1300,iostat=ios)(IA(i),i=1,NTEQ)
    call closef(lu)

    NAA = 0
    do i=1,NTEQ
        NAA = NAA + IA(i)
    enddo
    goto 950
endif

call CLEAR (JA,MIN(NTEQ*NTEQ,MAXJA),"JA")

do k=1,NTEQ+1
    IA(k) = k
enddo

do k=1,NTEQ
    JA(k) = k
enddo

do k=1,NUMEL
    do j=1,NEN

        node = IEN(j,k)
        numeq = LMT(j,k)

        if (numeq.eq.0) then
             cycle
        endif

        do jj=1,NEN
            kode = IEN(jj,k)
            mumeq = LMT(jj,k)
            if (mumeq.gt.0) then
                call SETAIJ(numeq,mumeq,IA, JA,NTEQ,MAXJA)
            endif
            if (nodift) then
                cycle
            endif

            mumeq = IABS(LMTX(jj,k))
            if (mumeq.gt.0) then
                call SETAIJ(numeq,mumeq,IA, JA,NTEQ,MAXJA)
            endif
        enddo

        if (nodift) then
            cycle
        endif

        numeq = IABS(LMT(j,k))
        if (numeq.eq.0) then
            cycle
        endif

        do jj=1,NEN
            kode = IEN(jj,k)
            mumeq = LMT(jj,k)
            if (mumeq.gt.0) then
                call SETAIJ(numeq,mumeq,IA, JA,NTEQ,MAXJA)
            endif
            mumeq = IABS(LMTX(jj,k))
            if (mumeq.gt.0) then
                call SETAIJ(numeq,mumeq,IA, JA,NTEQ,MAXJA)
            endif
        enddo
    enddo
enddo

!
!    determine total storage requirements, taking into account that
!    SETAIJ determined upper triangle indices only
! WIENAND here write non symmetric crs for later use
!          write(stdout,*) "IA(NTEQ): ", IA(NTEQ),
!     >     " and IA(NTEQ+1): ", IA(NTEQ+1)

do i=1,NTEQ+1
    myIA(i) = 0
enddo

do i=1,MAXJA
    myJA(i) = 0
enddo

do i=1,NTEQ
    cntr(i) = 0
enddo

do i=1,NTEQ
    myIA(i)=myIA(i)+IA(i)
    do j=IA(i)+1,IA(i+1)-1
        jj = JA(j)
        do k=jj+1,NTEQ+1
            myIA(k)=myIA(k) + 1
        enddo
    enddo
enddo

myIA(NTEQ+1)=myIA(NTEQ+1)+IA(NTEQ+1)

do i=1,NTEQ
    j=IA(i+1)
    jj=j-IA(i)

    do k=1,jj
        myJA(myIA(i+1)-k) = JA(j-k)
    enddo
enddo

do i=1,NTEQ
    do j=myIA(i)+1,myIA(i+1)-1
        if(myJA(j).gt.i) then
            myJA(myIA(myJA(j)) + cntr(myJA(j))) = i;
            cntr(myJA(j)) = cntr(myJA(j)) + 1
        endif
    enddo
enddo

if (irank.eq.0) then
    lu =nextlu(0)
    inquire (file='TSPARSE.DAT',exist=there)
    if (there) then
        call fbopen(lu, 'TSPARSE.DAT', 'old')
!              close(unit=lu,status='delete')
    endif

    call fbopen(lu, 'TSPARSE.DAT', 'new')
!          write(stdout,*) "writing sparsity file"
!          write(stdout,*) "IA(NTEQ): ", IA(NTEQ),
!     >     " and IA(NTEQ+1): ", IA(NTEQ+1)
    write(lu,err=1500,iostat=ios) (myIA(i),i=1,NTEQ+1)
    write(lu,err=1500,iostat=ios) (myJA(i),i=1,myIA(NTEQ+1))

    call closef(lu)

endif

  123 format(I5)
! WIENAND end write

do i=1,NTEQ
    JA(IA(i)) = 0
enddo

do i=1,NTEQ-1
    do jj=IA(i)+1,IA(i+1)-1
        j = JA(jj)
        JA(IA(j)) = JA(IA(j)) + 1 ! account for contribution M(j,i)
    enddo
enddo

do i=1,NTEQ
    IA(i) = JA(IA(i)) + IA(i+1)-IA(i)
enddo

call CLEAR (JA,MAXJA,"JA")
NAA = 2*IA(NTEQ+1)-NTEQ-2
if (irank.eq.0) then
    lu = nextlu(0)
    call fbopen(lu,'TTOPO.DAT','new')
    write(lu,err=1400,iostat=ios) (IA(i),i=1,NTEQ)
    call closef(lu)
endif
!
950    if (OUTFIL(1)) then
    meanbn = 0

    do i=1,NTEQ
        meanbn = meanbn + IA(i)
    enddo

    meanbn = NINT( DBLE(meanbn+1)/DBLE(2*NTEQ) )
    write(luout(1),2) NTEQ,NAA,meanbn
 2       format(///1x, &
     'T H E R M A L   E Q U A T I O N   S Y S T E M   D A T A', &
     //5X, &
     'NUMBER OF EQUATIONS  . . . . . . . . . . . .  (NTEQ) =', &
     I12//5X, &
     'NUMBER OF TERMS IN STIFFNESS . . . . . . . . . (NAA) =', &
     I12//5X, &
     'MEAN HALF BANDWIDTH. . . . . . . . . . . . . . . . . =', &
     I12)
endif


call MPI_Barrier(PETSC_COMM_WORLD,ierr)


return


1100    write(stderr,4) MAXJA
 4   format(/1x,'SETNZT: need more temporary memory than ',I9)
 call xit(1," ")

1200    write(stderr,5) ios
 5   format(/1x,'SETNZT: read error ',I4,' on "TTOPO.DAT"')
 call xit(1," ")

1300    write(stderr,6)
 6   format(/1x,'SETNZT: premature EOF error on "TTOPO.DAT"')
 call xit(1," ")

1400    write(stderr,7) ios
 7   format(/1x,'SETNZT: write error ',I4,' on "TTOPO.DAT"')
 call xit(1," ")

1500  write(stderr,8) ios
 8   format(/1x,'SETNZT: write error ',I4,' on "TSPARSE.DAT"')
 call xit(1," ")

end subroutine
#endif
