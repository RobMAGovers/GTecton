    integer function intchr(cnumber)
    character*(*) cnumber
    external lnblk
    lue=iflu('stderr')
    read(cnumber,1,err=10,end=20,iostat=ierr) intchr
    1    format(I10)
    return
10    write(lue,2) ierr,cnumber(1:lnblk(cnumber))
    2    format('intchr: internal read error ',i3,' in ',80a)
    call exitp(1)
20    write(lue,3) ierr,cnumber(1:lnblk(cnumber))
    3    format('intchr: internal read EOF error ',i3,' in ',80a)
    call exitp(1)
    end
