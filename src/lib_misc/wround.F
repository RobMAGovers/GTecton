subroutine wround (S,x,maxdig,peqn)


! writes x in readible format to string S. Max. # digits is MAXDIG

implicit none
double precision, parameter :: TEN=10.0D0
double precision, parameter :: dfmin = 2.23d-308


!-pass
integer maxdig
character*(*) S
logical peqn

double precision  :: x
!-locl
double precision  :: xx, xp, xxp
integer           :: nfields, ndig, ne
character(len=80) :: f
character(len=1)  :: q
!-external
integer, external :: isigdig,LOGRND
!-init
q = CHAR(39)

xx = ABS(x)

if (xx.le.dfmin) then
    write(S,1)
1        format('0')
else
!   XX e [0.1, 10000> ?
    if (xx.ge.0.1 .and. xx.lt.1D4) then
!       write in "I"- or "F"-format
        nfields = MAX0(INT(LOG10(xx)+0.01)+1,1)
!       determine number of significant digits in xx
        ndig = isigdig (xx,maxdig,dfmin)

           if (ndig.gt.0) then
               nfields = nfields+ndig+1
        endif

        if (x.lt.0.) then
            nfields=nfields+1
        endif

        if (ndig.eq.0) then
            write(f,2,err=1000) nfields
    2            format('(I',I1,')')
            write(S,fmt=f,err=1000) INT(x)
        else
            write(f,4,err=1000) nfields,ndig
    4            format('(F',I1,'.',I1,')')
            write(S,fmt=f,err=1000) x
        endif
    else
!        write "E"-format
        xp = LOGRND(xx,dfmin)
        xx = xx / TEN**xp
        if (.not.peqn) then
            ndig = isigdig (xx,maxdig,dfmin)
            nfields = ndig+4

            if (x.lt.0.) then
                nfields = nfields + 1
            endif

            xxp = ABS(xp)
            ne = 1

            if (xxp.ge.10.)  then
                ne = ne + 1
            endif

            if (xxp.ge.100.) then
                ne = ne + 1
            endif

            nfields= nfields + ne
            if (nfields.lt.10) then
                write(f,6,err=1000) nfields,ndig,ne
    6                format('(1PE',I1,'.',I1,'E',I1,')')
            else
                write(f,8,err=1000) nfields,ndig,ne
    8                format('(1PE',I2,'.',I1,'E',I1,')')
            endif
            write(S,fmt=f,err=1000) x
        else
            nfields = MAX0(INT(LOG10(xx))+1,1)
!            determine number of significant digits in xx
            ndig = isigdig (xx,maxdig,dfmin)
            nfields = nfields+ndig+1
            if (x.lt.0.) nfields=nfields+1
            xxp = ABS(xp)
            ne = 1
            if (xxp.ge.10.)  ne = ne + 1
            if (xxp.ge.100.) ne = ne + 1
            if (xp.lt.0.)    ne = ne + 1
            if (ne.lt.10) then
                write(f,16,err=1000) q,q,nfields,ndig,q,q,ne,q,q
   16            format('(',a1,'@',a1,',F',i1,'.',i1,',',a1, &
             ' ~10 sup ',a1,',I',i1,',',a1,'@',a1,')')
                    else if (ne.lt.100) then
                write(f,18,err=1000) q,q,nfields,ndig,q,q,ne,q,q
   18            format('(',a1,'@',a1,',F',i1,'.',i1,',',a1,&
             ' ~10 sup ',a1,',I',i2,',',a1,'@',a1,')')
            else
                write(f,20,err=1000) q,q,nfields,ndig,q,q,ne,q,q
   20            format('(',a1,'@',a1,',F',i1,'.',i1,',',a1,&
             ' ~10 sup ',a1,',I',i3,',',a1,'@',a1,')')
            endif
            write(S,fmt=f,err=1000) x/TEN**xp,INT(xp)
        endif
        endif
    endif
    return
1000    call prerr('wround')
    call exitp(1)
    end
!-------------------------------------------------------------------------------
    integer function isigdig (x,maxdig,dfmin)

! determines the number of significant digits in x

    implicit none
    double precision, parameter :: TEN=10d0
    double precision x,dfmin
!-pass
    integer maxdig
!-local
double precision :: y,yy, Ypower
    integer :: i

    y = ABS(x)
    i = 0
100    if (i.le.maxdig) then
        Ypower = TEN**DBLE(i)
        yy = DBLE(INT(y*Ypower))/Ypower
        if (ABS(y-yy).gt.dfmin) then
        i = i + 1
        goto 100
        endif
    else
        i = maxdig
    endif
    isigdig = i
    return
    end
!-------------------------------------------------------------------------------
    integer function LOGRND(x,dfmin)

    implicit none

    double precision :: xx, x, dfmin
    integer          :: ix

    xx = LOG10(ABS(x))
    ix = INT(xx)
    if (ABS(xx-DBLE(ix)).lt.dfmin) then
        LOGRND = ix
    else
        if (xx.lt.0.0) then
            LOGRND = ix - 1
        else
            LOGRND = ix
        endif
    endif
    return
    end
