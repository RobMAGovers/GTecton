    program main

! computes relative velocity 

    implicit none
    double precision one,gr,RE
    parameter (one=1d0,gr=3.141592653589793116d0/1.8d2,RE=6375d3)
!-locl
    double precision lonp,latp,wp,lonx,latx,w(3),x(3),v(3),dotpr
    double precision vl(3),azimuth,vel
    external dotpr

    call shell (lonp,latp,wp,lonx,latx)

    lonp = lonp * gr
    latp = latp * gr
    wp = wp * gr
    lonx = lonx * gr
    latx = latx * gr

    w(1) = SIGN(one,wp)*COS(latp)*COS(lonp)
    w(2) = SIGN(one,wp)*COS(latp)*SIN(lonp)
    w(3) = SIGN(one,wp)*SIN(latp)
    x(1) = COS(latx)*COS(lonx)
    x(2) = COS(latx)*SIN(lonx)
    x(3) = SIN(latx)
    call crossp (w,x,v)

    call gtol(lonx,latx,v,vl)
    azimuth = ATAN2(vl(2),vl(1))/gr
    if (azimuth.lt.0d0) azimuth = azimuth + 360d0
    vel = ABS(wp)*RE*SQRT(dotpr(v,v))

    write(*,10) vel*1e-3,azimuth
   10    format(1x,'Velocity = ',F7.3,' mm/yr, azimuth = ',F8.3,' deg')

    end
!-----------------------------------------------------------------------
    subroutine shell (lonp,latp,wp,lonx,latx)

    implicit none
!-pass
    double precision lonp,latp,wp,lonx,latx
!-locl
    logical numeric
    character(len=80) usage,arg(7),dummy
    integer iflu,stderr,narg,iargc,lnblk,i
    double precision chreal
    external iflu,numeric,chreal,lnblk,adjustl

    call HowToUse ('-p lon lat rate -x lon lat',usage)
    stderr = iflu('stderr')

    narg = command_argument_count()
    if (narg.lt.7) then
        write(stderr,'(80a)') usage(1:lnblk(usage))
        call exitp(1)
    endif
    do i=1,7
        call get_command_argument(i,arg(i))
        call adjustl(arg(i))
    enddo

    if (arg(1).ne.'-p') then
        call strswp (arg(1),arg(5))
        call strswp (arg(2),arg(6))
        call strswp (arg(3),arg(7))
        call strswp (arg(3),arg(4))
        call strswp (arg(2),arg(3))
        call strswp (arg(1),arg(2))
    endif

    if (arg(1)(1:2).ne.'-p') then
        write(stderr,'(80a)') usage(1:lnblk(usage))
        call exitp(1)
    endif
    if (numeric(arg(2))) then
        lonp = chreal(arg(2))
    else
        write(stderr,'(80a)') usage(1:lnblk(usage))
        call exitp(1)
    endif
    if (numeric(arg(3))) then
        latp = chreal(arg(3))
    else
        write(stderr,'(80a)') usage(1:lnblk(usage))
        call exitp(1)
    endif
    if (numeric(arg(4))) then
        wp = chreal(arg(4))
    else
        write(stderr,'(80a)') usage(1:lnblk(usage))
        call exitp(1)
    endif
    if (arg(5)(1:2).ne.'-x') then
        write(stderr,'(80a)') usage(1:lnblk(usage))
        call exitp(1)
    endif
    if (numeric(arg(6))) then
        lonx = chreal(arg(6))
    else
        write(stderr,'(80a)') usage(1:lnblk(usage))
        call exitp(1)
    endif
    if (numeric(arg(7))) then
        latx = chreal(arg(7))
    else
        write(stderr,'(80a)') usage(1:lnblk(usage))
        call exitp(1)
    endif

    write(*,10) lonp,latp,wp
   10    format(1x,'Rotation pole: lon = ',F8.3,' deg, lat = ',F7.3, &
          ' deg, rate = ',F7.3,' deg/My')
    write(*,20) lonx,latx
   20    format(1x,'Coordinate: lon = ',F8.3,' deg, lat = ',F7.3)

    return
    end
!-----------------------------------------------------------------------
    subroutine strswp(a,b)
    character*(*) a,b
    character(len=80) c
    c = a
    a = b
    b = c
    return
    end
!-----------------------------------------------------------------------
        subroutine crossp(a,b,c)
        implicit none
        double precision a(3),b(3),c(3)
        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
        return
        end
!-----------------------------------------------------------------------
        double precision function dotpr(a,b)
        implicit none
        double precision a(3),b(3)
        dotpr=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
        return
        end
!-----------------------------------------------------------------------
        subroutine gtol(lon,lat,u,v)

! Transforms from global frame (x1: 0N,0E; x2: 0N,90E; x3: 90N) to local frame
! (x1: north; x2: east; x3: down).

        implicit none
!-pass
    double precision lon,lat,u(3),v(3)
!-locl
    integer i,j
        double precision t(3,3)

!    Set up rotation matrix
        t(1,1)= -SIN(lat)*COS(lon)
        t(1,2)= -SIN(lat)*SIN(lon)
        t(1,3)=  COS(lat)
        t(2,1)= -SIN(lon)
        t(2,2)=  COS(lon)
        t(2,3)=  0d0
        t(3,1)= -COS(lat)*COS(lon)
        t(3,2)= -COS(lat)*SIN(lon)
        t(3,3)= -SIN(lat)

        do i=1,3
           v(i)=0.d0
           do j=1,3
              v(i) = v(i) + t(i,j)*u(j)
           end do
        end do

        return
        end
