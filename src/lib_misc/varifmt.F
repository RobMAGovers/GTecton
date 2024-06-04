character*(*) function varifmt(i)

implicit none
!-pass
integer i
!-locl
character(len=10) :: cnumber
integer           :: n, m, ierr, lue
integer, external :: fnblnk, lnblk, iflu
!-init
lue=iflu('stderr')
varifmt = ' '

n = LEN(varifmt)

if (n.lt.13) then
    write(lue,1)
1        format(1x,'varifmt: length error')
    return
endif

write(cnumber,2,err=20,iostat=ierr) i
2    format(I10)
n = fnblnk(cnumber)
m = lnblk(cnumber)

if (n.ge.m .or. m.le.0) then
    return
endif

varifmt = '(A' // cnumber(n:m) // ')'

10    return

20    write(lue,3) ierr
   3    format(1x,'varifmt: internal write error ',i3) 
call prerr('varifmt')
goto 10
end function
