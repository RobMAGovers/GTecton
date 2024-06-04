character(len=255) function dirname (string)

implicit none

!-pass
character(len=255) :: string
!-locaL
character(len=1)   :: slash
integer, external  :: lnblk, iflu
integer            :: i, l, m
integer                  :: lue

data slash /'/'/
lue=iflu('stderr')

l = lnblk(string)
if (l.gt.1) then
    i = l-1

10  if (i.ge.1) then
        if (string(i:i).eq.slash) goto 20
        i = i - 1
        goto 10
    endif

20  if (i.eq.1) then
        dirname = slash
    else if (i.gt.1) then
        m = LEN(dirname)
        if (i-1.le.m) then
            dirname = string(1:i-1)
        else
            write(lue,1)
1            format(1x,'dirname: insufficient length error')
            call exitp(1)
        endif
    else
        dirname = '.'
    endif

else

   if (string(1:1).eq.slash) then
       dirname = slash
   else
       dirname = '.'
   endif

endif
return
end function
