integer function lnblk(string)

implicit none

!     returns position of last non-blank in "string"

!-pass
character string*(*)
!-local
character(len=1)  :: space, tab
integer           :: i, lue
integer, external :: iflu

space = CHAR(32)
tab   = CHAR(9)

i = LEN(string)

if (i.gt.0) then

100 if (i.gt.0) then
        if (string(i:i).eq.space .or. &
            string(i:i).eq.tab) then

            i=i-1
            goto 100
        endif
    else
        i=0
    endif

else
    lue = iflu('stderr')
    write(lue,1)
    1     format(1x,'lnblk: string length error')
    i = -1
endif

lnblk=i
return
end function
