    program main

! tensor transformation

    implicit none
    integer i,j
    double precision T(3,3),R(3,3),dr,lon,lat,phi,lambda
!-init
    dr = DATAN(1D0)/45d0
    do i=1,3
        do j=1,3
        T(i,j) = 0d0
        enddo
    enddo

    T(3,2) = 1d0
    T(2,3) = -T(3,2)
    lon = 50.0
    lat = 20.0

    phi = lon*dr
    lambda = (90d0-lat)*dr ! co-latitude

    R(1,1) = COS(lambda)*COS(phi)
    R(1,2) = -COS(lambda)*SIN(phi)
    R(1,3) = SIN(lambda)
    R(2,1) = SIN(phi)
    R(2,2) = COS(phi)
    R(2,3) = 0d0
    R(3,1) = -SIN(lambda)*COS(phi)
    R(3,2) = SIN(lambda)*SIN(phi)
    R(3,3) = COS(lambda)

    call dyade (R,T)

    write(*,1) lon,lat
    1    format(1x,'@ lon,lat=',F7.2,1x,F6.2,':')
    write(*,2) T(1,1),T(1,2),T(1,3)
    2    format(3(1X,F9.5))
    write(*,2) T(2,1),T(2,2),T(2,3)
    write(*,2) T(3,1),T(3,2),T(3,3)

    end
!-------------------------------------------------------------------------------
    subroutine dyade (R,T)

! transforms tensor T by rotation matrix R

    implicit none
!-pass
    double precision R(3,3),T(3,3)
!-locl
    double precision TL(3,3)
    integer i,j,k,l

    do i=1,3
        do j=1,3
        TL(i,j) = 0d0
        do k=1,3
            do l=1,3
            TL(i,j) = TL(i,j) + R(i,k)*R(j,l)*T(k,l)
            enddo
        enddo
        enddo
    enddo
    do i=1,3
        do j=1,3
        T(i,j) = TL(i,j)
        enddo
    enddo

    return
    end
