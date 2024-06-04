#ifdef EXPERIMENTAL_FRICTION

module FrictionModule

  type FrictionEntry
    integer :: faultNumber
    integer :: point(3)
    double precision :: staticFriction
    double precision :: Winkler
    double precision :: sideArea
    double precision :: NodalFrictionContribution(3)
    double precision :: NodalStrengthContribution(3)
  end type FrictionEntry

  type(FrictionEntry), allocatable :: friction(:)
  
contains


subroutine frictionTraction (fault_number,DIFORC,nFrictionEntriesGlobal,NTYPE,NDOF,&
                             writeToFile, WriteToFileName)

use modeldatamodule, only: modeldatactx


! Computes tractions on fault segments from differential forces on slippery nodes.
! Finds the solution to the overdetermined system of equations using LAPACK 
! routine DGELS, which is based on least squares minimization of 
! || B - A*X || using QR or LQ factorization.
! RG 2018-7-20

implicit none

! pass
integer fault_number,nFrictionEntriesGlobal,NTYPE,NDOF
double precision DIFORC(:,:)

logical :: writeToFile
character(len=*) :: WriteToFileName

! local
double precision ZERO,TWO,THREE
parameter (ZERO=0d0, TWO=2d0, THREE=3d0)
logical, parameter :: test = .FALSE.,condition = .TRUE.
integer, allocatable :: EQN(:,:),nodes(:),GlobalNode(:)
double precision, allocatable :: A(:,:),B(:,:),WORK(:)
integer numeqn,ientry,jentry,node,NSN,NSEGM,LWORK,inode,ieq,idof,info,max_node,&
 min_node,numnodes,eqnentry,eqncount
double precision area,AverageArea,AverageB(3)

integer, external :: nextlu
integer           :: fileID, iError


!***** determine number of segments NSEGM in fault fault_number
NSEGM=0
do ientry=1,nFrictionEntriesGlobal
  if (friction(ientry)%faultNumber .eq. fault_number) NSEGM = NSEGM + 1
enddo
if (test) print*,'NSEGM=',NSEGM

!***** determine number of slippery nodes NSN in this fault 
ALLOCATE ( nodes(3*NSEGM) )
nodes = 0
numnodes = 0
do ientry=1,nFrictionEntriesGlobal
  if (friction(ientry)%faultNumber .eq. fault_number) then
    numnodes = numnodes + 1
    nodes(numnodes) = friction(ientry)%point(1)
    numnodes = numnodes + 1
    nodes(numnodes) = friction(ientry)%point(2)
    if (NTYPE.eq.6) then
      numnodes = numnodes + 1
      nodes(numnodes) = friction(ientry)%point(3)
    endif
  endif
enddo
max_node = maxval(nodes)
do ientry=numnodes+1,3*NSEGM
  nodes(ientry) = max_node
enddo
min_node = minval(nodes)-1
NSN = 0
do while (min_node < max_node)
  NSN = NSN + 1
  min_node = minval(nodes, mask=nodes>min_node)
enddo
DEALLOCATE (nodes)

if (test) print*,'NSN=',NSN
ALLOCATE ( EQN(3,NSN) ); EQN = 0
ALLOCATE ( GlobalNode(NSN) ); GlobalNode = 0

!***** loop FrictionEntry array to set local equation numbers
numeqn = 0
eqnentry = 0
do ientry=1,nFrictionEntriesGlobal
  if (friction(ientry)%faultNumber .eq. fault_number) then
    if (test) write(*,81) ientry,fault_number
81  format('friction entry',I3,' is for fault ',I3)
    if (numeqn.eq.0) then
      eqnentry = eqnentry + 1
      numeqn = numeqn + 1
      EQN(1,eqnentry) = numeqn
      GlobalNode(numeqn) = friction(ientry)%point(1)
      if (test) write(*,82) numeqn,GlobalNode(numeqn),1,eqnentry,EQN(1,eqnentry)
82    format('GlobalNode(',I3,')=',I3,' EQN(',I1,',',I3,')=',I3)
      numeqn = numeqn + 1
      EQN(2,eqnentry) = numeqn
      GlobalNode(numeqn) = friction(ientry)%point(2)
      if (test) write(*,82) numeqn,GlobalNode(numeqn),2,eqnentry,EQN(2,eqnentry)
      if (NTYPE.eq.6) then
        numeqn = numeqn + 1
        EQN(3,eqnentry) = numeqn
        GlobalNode(numeqn) = friction(ientry)%point(3)
        if (test) write(*,82) numeqn,GlobalNode(numeqn),3,eqnentry,EQN(3,eqnentry)
      endif
    else
      node = friction(ientry)%point(1)
      eqnentry = eqnentry + 1
      eqncount = 0
      !**** search previous records for the same global node number
      do jentry=1,ientry-1
        if (friction(jentry)%faultNumber .eq. fault_number) then
          eqncount = eqncount + 1
          if (node.eq.friction(jentry)%point(1)) then
            EQN(1,eqnentry) = EQN(1,eqncount)
            if (test) write(*,83) 1,eqnentry,1,eqncount,EQN(1,eqncount)
83          format('set EQN(',I1,',',I3,')=EQN(',I1,',',I3,')=',I3)
            goto 100
          endif
        endif
      enddo
      eqncount = 0
      do jentry=1,ientry-1
        if (friction(jentry)%faultNumber .eq. fault_number) then
          eqncount = eqncount + 1
          if (node.eq.friction(jentry)%point(2)) then
            EQN(1,eqnentry) = EQN(2,eqncount)
            if (test) write(*,83) 1,eqnentry,2,eqncount,EQN(2,eqncount)
            goto 100
          endif
        endif
      enddo
      if (NTYPE.eq.6) then
        eqncount = 0
        do jentry=1,ientry-1
          if (friction(jentry)%faultNumber .eq. fault_number) then
            eqncount = eqncount + 1
            if (node.eq.friction(jentry)%point(3)) then
              EQN(1,eqnentry) = EQN(3,eqncount)
              if (test) write(*,83) 1,eqnentry,3,eqncount,EQN(3,eqncount)
              goto 100
            endif
          endif
        enddo
      endif
      numeqn = numeqn + 1
      EQN(1,eqnentry) = numeqn
      GlobalNode(numeqn) = node
      if (test) write(*,82) numeqn,GlobalNode(numeqn),1,eqnentry,EQN(1,eqnentry)

100   node = friction(ientry)%point(2)
      eqncount = 0
      do jentry=1,ientry-1
        if (friction(jentry)%faultNumber .eq. fault_number) then
          eqncount = eqncount + 1
          if (node.eq.friction(jentry)%point(1)) then
            EQN(2,eqnentry) = EQN(1,eqncount)
            if (test) write(*,83) 2,eqnentry,1,eqncount,EQN(1,eqncount)
            goto 200
          endif
        endif
      enddo
      eqncount = 0
      do jentry=1,ientry-1
        if (friction(jentry)%faultNumber .eq. fault_number) then
          eqncount = eqncount + 1
          if (node.eq.friction(jentry)%point(2)) then
            EQN(2,eqnentry) = EQN(2,eqncount)
            if (test) write(*,83) 2,eqnentry,2,eqncount,EQN(2,eqncount)
            goto 200
          endif
        endif
      enddo
      if (NTYPE.eq.6) then
        eqncount = 0
        do jentry=1,ientry-1
          if (friction(jentry)%faultNumber .eq. fault_number) then
            eqncount = eqncount + 1
            if (node.eq.friction(jentry)%point(3)) then
              EQN(2,eqnentry) = EQN(3,eqncount)
              if (test) write(*,83) 2,eqnentry,3,eqncount,EQN(3,eqncount)
              goto 200
            endif
          endif
        enddo
      endif
      numeqn = numeqn + 1
      EQN(2,eqnentry) = numeqn
      GlobalNode(numeqn) = node
      if (test) write(*,82) numeqn,GlobalNode(numeqn),2,eqnentry,EQN(2,eqnentry)

200   if (NTYPE.ne.6) goto 300
      node = friction(ientry)%point(3)
      eqncount = 0
      do jentry=1,ientry-1
        if (friction(jentry)%faultNumber .eq. fault_number) then
          eqncount = eqncount + 1
          if (node.eq.friction(jentry)%point(1)) then
            EQN(3,eqnentry) = EQN(1,eqncount)
            if (test) write(*,83) 3,eqnentry,1,eqncount,EQN(1,eqncount)
            goto 300
          endif
        endif
      enddo
      eqncount = 0
      do jentry=1,ientry-1
        if (friction(jentry)%faultNumber .eq. fault_number) then
          eqncount = eqncount + 1
          if (node.eq.friction(jentry)%point(2)) then
            EQN(3,eqnentry) = EQN(2,eqncount)
            if (test) write(*,83) 3,eqnentry,2,eqncount,EQN(2,eqncount)
            goto 300
          endif
        endif
      enddo
      eqncount = 0
      do jentry=1,ientry-1
        if (friction(jentry)%faultNumber .eq. fault_number) then
          eqncount = eqncount + 1
          if (node.eq.friction(jentry)%point(3)) then
            EQN(3,eqnentry) = EQN(3,eqncount)
            if (test) write(*,83) 3,eqnentry,3,eqncount,EQN(3,eqncount)
            goto 300
          endif
        endif
      enddo
      numeqn = numeqn + 1
      EQN(3,eqnentry) = numeqn
      GlobalNode(numeqn) = node
      if (test) write(*,82) numeqn,GlobalNode(numeqn),3,eqnentry,EQN(3,eqnentry)

300   continue
    endif
  endif
enddo
if (numeqn.ne.NSN) then
  write(0,*) 'traction: numeqn=',numeqn,'NSN=',NSN
  stop 'fatal'
endif
if (test) print*,'GlobalNode: ',GlobalNode
if (test) print*,'EQN: ',EQN

! start setting up matrix system
! A is the rectangular matrix as defined in DGELS (NSEGM<NSN)
! B is the RHS vector with NDOF RHSs as defined in DGELS
! WORK is a work array as defined in DGELS
! LWORK is 2*NSEGM as defined in DGELS
! NDOF is equal to the number of RHS vectors for which a solution is sought

ALLOCATE( A(NSN,NSEGM) ); A = ZERO
ALLOCATE( B(NSN,NDOF) ) ; B = ZERO
LWORK = 2*NSEGM
ALLOCATE( WORK(LWORK) ) 

if (test) print*,'setting up A'
if (condition) then
    AverageArea = 0d0
    eqncount = 0
    do ientry=1,nFrictionEntriesGlobal
      if (friction(ientry)%faultNumber .eq. fault_number) then
        eqncount = eqncount + 1
        area = friction(ientry)%sideArea
        if (area.lt.1D-6) then
            write(0,*) 'frictionTraction: area of entry ',ientry,' (close to) zero'
            stop 'fatal'
        endif
      endif
      AverageArea = AverageArea + area
    enddo
    if (eqncount.lt.1) then
      write(0,*) 'frictionTraction: no segments in the fault with label ',fault_number
      stop 'fatal'
    endif
    AverageArea = AverageArea/DBLE(eqncount)
endif
eqncount = 0
do ientry=1,nFrictionEntriesGlobal
  if (friction(ientry)%faultNumber .eq. fault_number) then
    eqncount = eqncount + 1
    area = friction(ientry)%sideArea
    if (area.lt.1D-6) then
        write(0,*) 'frictionTraction: area of entry ',ientry,' (close to) zero'
        stop 'fatal'
    endif
    area = area / AverageArea
    if (NTYPE.eq.6) then
      do inode=1,3
        ieq = EQN(inode,eqncount)
        A(ieq,eqncount) = A(ieq,eqncount) + area/THREE
        if (test) write(*,9) ientry,ieq,eqncount
9       format('   Add area/3 of fault entry ',I1,' to A(',I1,',',I1,')')
      enddo
    else
      do inode=1,2
        ieq = EQN(inode,eqncount)
        A(ieq,eqncount) = A(ieq,eqncount) + area/TWO
        if (test) write(*,85) ientry,ieq,eqncount
85      format('   Add area/2 of fault entry ',I1,' to A(',I1,',',I1,')')
      enddo
    endif
  endif
enddo


if (test) then
    print*,'setting up B'
endif

if (condition) then
  AverageB = 0d0

  do ieq=1,NSN
    node = GlobalNode(ieq)
    do idof=1,NDOF
      AverageB(idof) = AverageB(idof) + ABS(DIFORC(idof,node))
    enddo
  enddo

  do idof=1,NDOF

    AverageB(idof) = AverageB(idof) / DBLE(NSN)
    if (AverageB(idof).lt.1D-9) then
        AverageB(idof) = 1d0
    endif
  enddo

endif

do ieq=1,NSN
  node = GlobalNode(ieq)
  do idof=1,NDOF
    if (test) then
        write(*,13) ieq,idof,idof,node
    endif
13  format('   setting B(',I1,',',I1,') to DIFORC(',I1,',',I3,')')
    B(ieq,idof) = DIFORC(idof,node) / AverageB(idof)
  enddo
enddo

DEALLOCATE( EQN )
DEALLOCATE( GlobalNode )
!---------------- solve -------------
call DGELS ('N',NSN,NSEGM,NDOF,A,NSN,B,NSN,WORK,LWORK,info)
INFO=-1
!-------------------------------------
if (info.lt.0) then
  write(0,*) 'traction: dgels argument ',-info,' has illegal value'
  stop 'fatal'
else if (info.gt.0) then
  write(0,*) 'traction: dgels the ',info,'-th diagonal element of the triangular&
   & factor of A is zero'
  stop 'fatal'
endif
DEALLOCATE( A )
DEALLOCATE( WORK )

if (condition) then
  jentry = 0
  do ientry=1,nFrictionEntriesGlobal
    if (friction(ientry)%faultNumber .eq. fault_number) then
      jentry = jentry + 1
      do idof=1,ndof
! inverse preconditioning
!        B(jentry,idof) = B(jentry,idof)*AverageArea/AverageB(idof)
        B(jentry,idof) = B(jentry,idof) / AverageArea * AverageB(idof)

      enddo
    endif
  enddo
endif


write(*,*) " "

fileID=nextlu()
open(unit=fileID, file=WriteToFileName, iostat=iError)
if (iError.ne.0) then
    write(*,*) "Could not open file ", trim(WriteToFileName)
    write(*,*) "to write fictions to. Error value ", iError
    stop "Leaving GTecton..."
endif

!write(*,*) "frictra has ntype: ", ntype, fileID
!write(*,*) "frictra has file stuff", fileID, writeToFile, trim(WriteToFileName)

! print results
jentry = 0
do ientry=1,nFrictionEntriesGlobal
    if (friction(ientry)%faultNumber .eq. fault_number) then
        jentry = jentry + 1
        if (NTYPE.eq.6) then
            ! NTYPE 6 -> tetrahedrons
            if(writeToFile) then
                write(fileID,*) friction(ientry)%faultNumber, &
                                friction(ientry)%point(1:3), &
                                (B(jentry,idof),idof=1,NDOF)
            else
                print*,friction(ientry)%faultNumber,&
                       friction(ientry)%point(1:3),&
                       (B(jentry,idof),idof=1,NDOF)
            endif
        else
            ! likely triangles.
            if(writeToFile) then
                write(fileID,*) friction(ientry)%faultNumber, &
                                friction(ientry)%point(1:2), &
                                (B(jentry,idof),idof=1,NDOF)
            else
                print*,friction(ientry)%faultNumber,&
                       friction(ientry)%point(1:2),&
                       (B(jentry,idof),idof=1,NDOF)
            endif
        endif
    endif
enddo
DEALLOCATE( B )

close(fileID)

return
end subroutine

end module FrictionModule

#endif
