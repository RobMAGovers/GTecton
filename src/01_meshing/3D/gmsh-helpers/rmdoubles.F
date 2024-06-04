program main

! removes double element/node entries from input

implicit none
integer MAX
parameter (MAX=10000000)
!-locl
real                 :: xn,yn,zn
integer              :: N,i,element,node
integer, allocatable :: LIST(:,:)
real,    allocatable :: V(:,:)

allocate(LIST(2,MAX))
allocate(MAX(3,MAX))

N = 0
100    read(*,10,end=200,err=1000) element,node,xn,yn,zn

10    format(2I8,3E14.0)
if (N.gt.0) then
    do i=1,N
        if (LIST(1,i).eq.element .and. &
            LIST(2,i).eq.node) then
			goto 150
		endif
    enddo
endif

N = N + 1

if (N.gt.MAX) then
	stop 'overflow'
endif

LIST(1,N) = element
LIST(2,N) = node

V(1,N) = xn
V(2,N) = yn
V(3,N) = zn

150    goto 100

200    do i=1,N
    write(*,20) LIST(1,i),LIST(2,i),V(1,i),V(2,i),V(3,i)
   20        format(2I8,3(1PG14.6))
enddo




stop "Completed" 

1000    stop "read error"
end program
        
