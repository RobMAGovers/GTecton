    program main
*
c removes double element/node entries from input
*
    implicit none

    integer MAX
    parameter (MAX=1000000)
c-locl
    real V,xn,yn,zn
    integer LIST,N,i,element,node
    dimension LIST(2,MAX),V(3,MAX)
*
    N = 0
100    read(*,10,end=200,err=1000) element,node,xn,yn,zn
   10    format(2I5,3E14.0)

    ! check whether element/node combo is already in list.
    if (N.gt.0) then
        do i=1,N
            ! if so, skip it
            if (LIST(1,i).eq.element.and.LIST(2,i).eq.node) then
                goto 150
            endif
        enddo
    endif

    N = N + 1

    if (N.gt.MAX) then
        stop 'Too many node/element combinations'
    endif

    ! if not, add it to the list.
    LIST(1,N) = element
    LIST(2,N) = node
    V(1,N) = xn
    V(2,N) = yn
    V(3,N) = zn

150    goto 100

    
200    do i=1,N
        write(*,20) LIST(1,i),LIST(2,i),V(1,i),V(2,i),V(3,i)
   20        format(2I5,3(1PG14.6))
    enddo

    stop
1000    stop 'read error'
    end
        
