integer function indx (a1,a2)

implicit none

! returns an integer value indicating the starting position within
! character string a1 of substring a2. If a2 occurs more than once in
! a1, the starting position of the first occurence is returned.
! Backslashed characters will not be tested and will therefore not be
! recognized!


!-pass
character    :: a1*(*),a2*(*)
!local
character(1) :: backslash
logical      :: backslashed
integer      :: inx, i1, i2, l1, l2


inx       = 0
backslash = CHAR(92)
l1        = LEN(a1)
l2        = LEN(a2)
if (l1.ge.l2) then
    i1 = 1
    i2 = 1
    backslashed = .false.
10  if (i1.le.l1) then
           if (backslashed) then
!            do not compare
            backslashed = .false.
            i2 = 1
        else
            if (a2(i2:i2).eq.backslash) then
!            do not interpret backslashes in search string
                if (a1(i1:i1).eq.a2(i2:i2)) then
                    i2 = i2 + 1
                    if (i2.gt.l2) then
                        inx = i1 - i2 + 2
                        goto 20
                    endif
                    backslashed = .false.
                else
                    i2 = 1
                    backslashed = .false.
                endif
            else
                if (a1(i1:i1).eq.backslash) then
                    backslashed = .true.
                else if (a1(i1:i1).eq.a2(i2:i2)) then
                    i2 = i2 + 1
                    if (i2.gt.l2) then
                        inx = i1 - i2 + 2
                        goto 20
                    endif
                    backslashed = .false.
                else
                    if (i2.gt.1) then
                        i1 = i1 - 1
                    endif
                    i2 = 1
                    backslashed = .false.
                endif
            endif
        endif
        i1 = i1 + 1
        goto 10
    endif
20        continue
endif
indx = inx
return
end function
