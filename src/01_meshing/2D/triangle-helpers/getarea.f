!     getarea returns some information regarding the mesh
!     1. the area of the largest element
!     2. the area for the smallest element
!     3. the minimum and maximum distance between nodal points
!     4. the minimal and maximum height of an element
      program getarea


      implicit double precision (a-h,o-z)
      parameter (MAXEL=999999, MAXNP=MAXEL)
      parameter (FACTOR=0.05d0, POW=1.0d0)

      double precision, parameter :: dfmax = 1.79d+308
      double precision, parameter :: eps = 2.22d-16


      character*80 fna
      dimension attr(3)
      integer, allocatable :: IEN(:,:)
      double precision, allocatable :: X(:)
      double precision, allocatable :: Y(:)

!      dimension IEN(3,MAXEL),X(MAXNP),Y(MAXNP)
      integer ERROR

*
      lue = iflu('stderr')
      narg = iargc()
      if (narg.eq.0) then
          write(lue,10)
   10          format(1x,'Usage: setarea basename')
          call exitp(1)
      endif
      call getarg(1,fna)
*
      lu = nextlu(0)
      call openf(lu,fna(1:lnblnk(fna))//'.node','old')
      read (lu,*) NUMNP,idim,iattr,ibound
      if (NUMNP.gt.MAXNP) stop 'overflow on NUMNP'
      ALLOCATE(X(NUMNP),STAT=ERROR)
      ALLOCATE(Y(NUMNP),STAT=ERROR)

      xmin = dfmax
      xmax = -dfmax
      ymin = dfmax
      ymax = -dfmax
      do k=1,NUMNP
          read (lu,*) ip,X(k),Y(k),imark
          if (X(k).gt.xmax) xmax = X(k)
          if (X(k).lt.xmin) xmin = X(k)
          if (Y(k).gt.ymax) ymax = Y(k)
          if (Y(k).lt.ymin) ymin = Y(k)
      end do
      call closef(lu)
      write(lue,20) fna(1:lnblnk(fna))//'.node'
   20      format(1x,'EOF ',90a)
*
      call openf(lu,fna(1:lnblnk(fna))//'.ele','old')
      read (lu,*) NUMEL,npoints,iattr
      if (NUMEL.gt.MAXEL) stop 'overflow on NUMEL'
      if (npoints.ne.3) then
          write(lue,30)
   30          format('3 points per triangle supported')
          call exitp(1)
      endif
      ALLOCATE(IEN(3,NUMEL),STAT=ERROR)

*
      attr(1) = 0.
      attr(2) = 0.
      attr(3) = 0.
      do k=1,NUMEL
          if (iattr.gt.0) then
             read (lu,*) im,(IEN(i,k),i=1,3),(attr(i),i=1,iattr)
          else
             read (lu,*) im,(IEN(i,k),i=1,3)
          endif
          do i=1,3
            IEN(i,k) = IEN(i,k) + 1
          enddo
      end do
      call closef(lu)
      write(lue,20) fna(1:lnblnk(fna))//'.ele'
*
      xmax = MAX(ABS(xmax),ABS(xmin))
      ymax = MAX(ABS(ymax),ABS(ymin))
      areamax = -dfmax
      areamin = dfmax
      distmax = -dfmax
      distmin = dfmax
      heightmin = dfmax
      heightmax = -dfmax
!      distxmax = -dfmax
!      distxmin = dfmax
!      distymax = -dfmax
!      distymin = dfmax
      numelmax = 0
      numelmin = 0
      numelhmin = 0
      numelhmax = 0
      do k=1,NUMEL
          i1 = IEN(1,k)
          i2 = IEN(2,k)
          i3 = IEN(3,k)
          dx1 = X(i2)-X(i1)
          dy1 = Y(i2)-Y(i1)
          dx2 = X(i3)-X(i1)
          dy2 = Y(i3)-Y(i1)
          dx3 = X(i3)-X(i2)
          dy3 = Y(i3)-Y(i2)
          prd = 0.5*(dx1*dy2-dy1*dx2)
          if (ABS(prd).lt.eps) then
            write(lue,40) k
   40            format(1x,'degenerate triangle ',I5)
            call exitp(1)
          endif
          prd = ABS(prd)
          if (prd.gt.areamax) then
            areamax = prd
            numelmax = k
        endif
          if (prd.lt.areamin) then
            areamin = prd
            numelmin = k
        endif
!       next compute min and max dist between nodal points
        dx12 = SQRT(dx1*dx1 + dy1*dy1)
        dx13 = SQRT(dx2*dx2 + dy2*dy2)
        dx23 = SQRT(dx3*dx3 + dy3*dy3)
        distmax = max(max(dx12, max(dx13,dx23)),distmax)
        distmin = min(min(dx12, min(dx13,dx23)),distmin)
!        write(*,*) 'For element ', k, ' with area ', prd
!        write(*,*) ' and distmax ', max(dx12, max(dx13,dx23)),
!     .    ' and distmin ', min(dx12, min(dx13,dx23))
!        write(*,*) ' minimal height ', 
!     .    2.0 * prd/max(dx12, max(dx13,dx23))
        
        heighttmp = 2.0 * prd/max(dx12, max(dx13,dx23))
        if (heighttmp.lt.heightmin) then
            heightmin = heighttmp
            numelhmin = k
        endif
        heighttmp = 2.0 * prd/min(dx12, min(dx13,dx23))
        if (heighttmp.gt.heightmax) then
            heightmax = heighttmp
            numelhmax = k
        endif

      enddo
      write(*,*) 'Maximum element area: ', areamax
      write(*,*) ' for element ', numelmax, ' with vertices: '
      write(*,*) IEN(1,numelmax),': ', X(IEN(1,numelmax)), 
     . Y(IEN(1,numelmax))
      write(*,*) IEN(2,numelmax),': ', X(IEN(2,numelmax)),
     . Y(IEN(2,numelmax))
      write(*,*) IEN(3,numelmax),': ', X(IEN(3,numelmax)),
     . Y(IEN(3,numelmax))
      write(*,*) 'Minimum element area: ', areamin
      write(*,*) ' for element ', numelmin, ' with vertices: '
      write(*,*) IEN(1,numelmin),': ', X(IEN(1,numelmin)),
     . Y(IEN(1,numelmin))
      write(*,*) IEN(2,numelmin),': ', X(IEN(2,numelmin)),
     . Y(IEN(2,numelmin))
      write(*,*) IEN(3,numelmin),': ', X(IEN(3,numelmin)),
     . Y(IEN(3,numelmin))
      write(*,*) 'Maximum distance between nodal points: ', distmax
      write(*,*) 'Minimum distance between nodal points: ', distmin
      write(*,*) 'Maximum height of an element         : ', heightmax
      write(*,*) ' for element ', numelhmax
      write(*,*) 'Minimum height of an element         : ', heightmin
      write(*,*) ' for element ', numelhmin
!      write(*,*) 'Max x-distance between nodal points: ', distxmax
!      write(*,*) 'Min x-distance between nodal points: ', distxmin
!      write(*,*) 'Max y-distance between nodal points: ', distymax
!      write(*,*) 'Min y-distance between nodal points: ', distymin
c      write(*,*) (areamax + 3.0*areamin) *0.25
c      areamax2 = -dfmax
c      areamin2 = dfmax
c      do k=1,NUMEL
c          i1 = IEN(1,k)
c          i2 = IEN(2,k)
c          i3 = IEN(3,k)
c          dx1 = X(i2)-X(i1)
c          dy1 = Y(i2)-Y(i1)
c          dx2 = X(i3)-X(i1)
c          dy2 = Y(i3)-Y(i1)
c          prd = ABS(dx1*dy2-dy1*dx2)
c          area = (areamax + prd*3d0) / 4d0
c          if (area.gt.areamax2) areamax2 = area
c          if (area.lt.areamin2) areamin2 = area
c      enddo
c      write(*,*) 'Max area (2): ', areamax2
c      write(*,*) 'Min area (2): ', areamin2
*
      DEALLOCATE(X, STAT=ERROR)
      DEALLOCATE(Y, STAT=ERROR)
      DEALLOCATE(IEN, STAT=ERROR)

      end

c--------------------------------------------------------------------------------    
      double precision function EnsureNotZero(x1,x2,x3,eps,xmax)


!      double precision x1,x2,x3,eps,xmax, EnsureNotZero
      double precision x1,x2,x3,eps,xmax

      if (x1.lt.eps) then
          if (x2.lt.eps) then
              if (x3.lt.eps) then
                  EnsureNotZero = xmax
              else
                  EnsureNotZero = x3
              endif
          else
              if (x3.lt.eps) then
                  EnsureNotZero = x2
              else
                  EnsureNotZero = min(x2,x3)
              endif
          endif
      else
          if ((x2.lt.eps).AND.(x3.lt.eps)) then
              EnsureNotZero = x1
          else if (x2.lt.eps) then
              EnsureNotZero = min(x1,x3)
          else if (x3.lt.eps) then
              EnsureNotZero = min(x1,x2)
          else
              EnsureNotZero = min(x1, min(x2,x3))
          endif
      endif
      return 
      end




