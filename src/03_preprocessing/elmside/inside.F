        program main
*
        real px(3),py(3),x,y
        integer i
        real dx1,dy1,dx2,dy2,dx3,dy3,O,w1,w2,w3
        data px(1),py(1) /0.0, 0.0/
        data px(2),py(2) /1.0, 0.0/
        data px(3),py(3) /0.0, 1.0/
*
100     do i=1,3
            write(*,*) px(i),py(i)
        enddo
        write(*,1)
   1    format(1x,'x,y > ',$)
        read(*,*) x,y
*
        call PNPOLY (x,y,px,py,3,inout)
        if (inout.eq.-1) goto 100
*
        dx1 = px(2) - px(1)
        dy1 = py(2) - py(1)
        dx2 = px(3) - px(2)
        dy2 = py(3) - py(2)
        O = dx1*dy2 - dy1*dx2
*
        dx1 = px(2) - x
        dy1 = py(2) - y
        dx2 = px(3) - px(2)
        dy2 = py(3) - py(2)
        w1 = dx1*dy2 - dy1*dx2
*
        dx1 = x - px(1)
        dy1 = y - py(1)
        dx2 = px(3) - x
        dy2 = py(3) - y
        w2 = dx1*dy2 - dy1*dx2
*
        dx1 = px(2) - px(1)
        dy1 = py(2) - py(1)
        dx2 = x - px(2)
        dy2 = y - py(2)
        w3 = dx1*dy2 - dy1*dx2
*
        print*,'O,w1,w2,w3=',O,w1,w2,w3
        print*,'L1,L2,L3=',w1/O,w2/O,w3/O
        goto 100
*
        end
c-------------------------------------------------------------------------------
        subroutine PNPOLY(PX,PY,XX,YY,N,INOUT)
*
C Determines whether a point is inside a polygon 
c The vertices may be listed clockwise or anticlockwise.
c The first may optionally be repeated; if so N may 
c optionally be increased by 1. 
c The input polygon may be a compound polygon consisting 
c of several separate subpolygons. If so, the first vertex 
c of each subpolygon must be repeated, and when calculating 
c N, these first vertices must be counted twice. 
c INOUT is the only parameter whose value is changed. 
c written by Randolph Franklin, university of ottawa, 7/70. 
c Method: a vertical line is drawn thru the point in question. If it 
c crosses the polygon an odd number of times, then the 
c point is inside of the polygon.
*
        IMPLICIT NONE
c-pass
        INTEGER N       ! Number of vertices of polygon
        INTEGER INOUT   ! Return value: -1 outside, 0 = on edge, 1 = inside
        REAL XX(N)      ! polygon x-coordinates
        REAL YY(N)      ! polygon y-coordinates
        REAL PX,PY      ! point coordinate
c-local
        INTEGER MAXDIM
        PARAMETER (MAXDIM=200)
        LOGICAL mx,my,nx,ny
        INTEGER i,j
        REAL x(MAXDIM),y(MAXDIM)
*
        if (N.gt.MAXDIM) then
            write(*,1)
1           format(1x,'PNPOLY overflow')
            return
        endif
        do i=1,N
            x(i) = XX(i) - PX
            y(i) = YY(i) - PY
        enddo
        INOUT = -1 
        do 2 i=1,N
            j = 1 + MOD(i,N)
            MX = x(i).GE.0.0
            NX = x(j).GE.0.0 
            MY = y(i).GE.0.0
            NY = y(j).GE.0.0 
            if (.not.((my.or.ny).and.(mx.or.nx)).or.(mx.and.nx)) goto 2
            if (.not.(my.and.ny.and.(mx.or.nx).and..not.(mx.and.nx)))
     >       goto 3 
            INOUT=-INOUT 
            goto 2
3           if((y(i)*x(j)-x(i)*y(j))/(x(j)-x(i))) 2,4,5
4           INOUT = 0
            return
5           INOUT = -INOUT
2       continue
        return
        end 
