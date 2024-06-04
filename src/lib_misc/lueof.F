integer function lueof (lu)

! system-independent routine for setting file pointer to EndOfFile

implicit none

!-pass
integer            :: lu
!-local
logical open
character(len=1)   :: a
character(len=11)  :: fm
character(len=255) :: fname
integer            :: ios, lue
integer, external  :: iflu, lnblk

lue = iflu('stderr')

inquire(unit=lu,opened=open,err=1000,iostat=ios)

if (.not.open) then
    write(lue,1)
1        format(1x,'lueof: unit=',I3,' not open')
    lueof = 1
    return
endif

inquire(unit=lu,form=fm,err=2000,iostat=ios)

if (fm.eq.'FORMATTED' .or. fm.eq.'formatted') then
100        read(lu,'(a)',end=200,err=2000,iostat=ios) a
    goto 100
200        continue
else
    write(lue,2)
2        format(1x,'lueof: cannot skip to EOF of unformatted file')
    lueof = 1
endif

lueof = 0
return

1000    write(lue,3) lu,ios
3    format(1x,'lueof: unit=',I3,' I/O error ',I3)
lueof = 1
return
2000    inquire(unit=lu,name=fname,err=1000,iostat=ios)
write(lue,4) ios,'"',fname(1:lnblk(fname)),'"'
4    format(1x,'lueof: I/O error ',I3,' on file ',255a)
lueof = 1
return
end
