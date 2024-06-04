logical function AsciF (fileName)

implicit none

!-pass
character         :: fileName*(*)
!-locl
logical           :: there,isbinary
integer, external :: lnblk, iflu
integer           :: l, m, lue
!-init

lue = iflu('stderr')

l = LEN(fileName)
if (l.le.0) then
    write(lue,2)
2        format(1x,'AsciF: file argument error')
    call exitp(1)
endif

m=lnblk(fileName)
if (m.le.0) then
    write(lue,4)
4        format(1x,'AsciF: no file argument error')
    call exitp(1)
endif

inquire (file=fileName,exist=there)
if (.not.there) then
    write(lue,6) CHAR(34),fileName(1:m),CHAR(34)
6        format(1x,'AsciF: non-existing file: ',80A)
    call exitp(1)
endif

AsciF = (.not.isbinary(fileName,m))
return
end function
