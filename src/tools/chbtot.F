 program main
 USE ALGEBRA, only: CLEAR
!
! modifies BTOT for new boundary conditions for restart
!
 implicit double precision (a-h,o-z)
 parameter (NUMNP=5624, NDOF=2, MAX=NUMNP*NDOF)
 dimension BTOT1(MAX),BTOT2(MAX)
 dimension IBOND1(NDOF,NUMNP),IBOND2(NDOF,NUMNP)
 logical there,echo
!-init
 data echo /.true./
 call CLEAR(BTOT1,MAX,"BTOT1")
 call CLEAR(BTOT2,MAX,"BTOT2")
 call CLEAR(IBOND1,MAX,"IBOND1")
 call CLEAR(IBOND2,MAX,"IBOND2")
!
!    define bcs for phase1 model and determine # of eqn's
 NEQ1 = MAX
 do i=1,37
IBOND1(1,i) = 1
NEQ1 = NEQ1 - 1
 enddo
 IBOND1(2,37) = 1
 NEQ1 = NEQ1 - 1
 do i=5588,5624
    IBOND1(1,i) = 1
    NEQ1 = NEQ1 - 1
 enddo
!
!    define bcs for phase2 model and determine # of eqn's
 NEQ2 = MAX
 IBOND2(2,37) = 1
 NEQ2 = NEQ2 - 1
 do i=5588,5624
    IBOND2(1,i) = 1
    NEQ2 = NEQ2 - 1
 enddo
!
 inquire(file='BTOT.DAT',exist=there)
if (.not.there) then
write(0,1)
 1        format(1x,'"BTOT.DAT" not found')
call exitp(1)
 endif
lutmp = nextlu(0)
call fbopen (lutmp,'BTOT.DAT','old')
read(lutmp,err=100,end=200,iostat=ios) (BTOT1(i),i=1,NEQ1)
call closef(lutmp)
!
 n1 = 0
 n2 = 0
 do n=1,NUMNP
do i=1,NDOF
if (IBOND1(i,n).eq.0) n1 = n1 + 1
if (IBOND2(i,n).eq.0) n2 = n2 + 1
if (IBOND1(i,n).eq.0 .and. IBOND2(i,n).eq.0) &
       BTOT2(n2) = BTOT1(n1)
if (echo) then
    if (IBOND1(i,n).eq.0) then
    if (IBOND2(i,n).ne.0) write(0,2) i,n,n1
 2            format(1x,'Zeroed BTOT1(',I1,',',I4,') (neq=', &
           I5,')')
    else
    if (IBOND2(i,n).eq.0) write(0,3) i,n,n2
 3            format(1x,'Zeroed BTOT2(',I1,',',I4,') (neq=', &
           I5,')')
    endif
endif
enddo
 enddo
!
call fbopen (lutmp,'BTOT2.DAT','unknown')
write(lutmp,err=300,iostat=ios) (BTOT2(i),i=1,NEQ2)
call closef(lutmp)
!
 call exitp(0)
!
100     write(stderr,101) ios
  101   format(1x,'BTOT1: read error ',I5)
 call exitp(1)
200     write(stderr,201) ios
  201   format(1x,'BTOT1: EOF error ',I5)
 call exitp(1)
300     write(stderr,301) ios
  301   format(1x,'BTOT2: write error ',I5)
 call exitp(1)
 end
!-----------------------------------------------------------------------
 subroutine CLEAR(A,N)
 double precision A(N)
 do i=1,N
A(i) = 0d0
 enddo
 return
 end
!-----------------------------------------------------------------------
 subroutine CLEAR(IA,N)
 dimension IA(N)
 do i=1,N
IA(i) = 0
 enddo
 return
 end
