program main
! This program aims to convert a grid as described by triangle into a 
! grid as described by gtecton. For this it read the the domain.[node/ele/poly]
! files and it creates tecin.dat.[nps/ele]
!

implicit none

character(len=80) Nfile,Efile
logical :: writeSerialMode

!double precision :: attr(3)
integer           :: attr(3)
integer           :: m(3)
integer, external :: NINT,iflu,nextlu
double precision  :: x, y
integer           :: i, iattr, ibound, idim
integer           :: im, imark, ios, ip, np, ne, k
integer           :: luo, lui, lue
integer           :: npoints

writeSerialMode = .false.

call shell (Nfile,Efile,writeSerialMode)

lue = iflu('stderr')
lui = nextlu(0)

call openf(lui,Nfile,'old')
read (lui,*,err=990,end=1000,iostat=ios) np,idim,iattr,ibound

if (iattr.gt.0) then
    write(lue,1)
endif

1    format('node attributes not carried along')

if (ibound.eq.0) then
    write(lue,2)
2        format('external boundary nodes should be labeled <> 0')
    call exitp(1)
endif

!********* Do the nodes

luo = nextlu(lui+1)
call openf(luo,'tecin.dat.nps','unknown')

do k=1,np
    read (lui,*,err=990,end=1000,iostat=ios) ip,x,y,imark

    if (writeSerialMode) then
        write (luo,33) k,imark,x,y
    else
        write (luo,3) k,imark,x,y
    endif
enddo

!3      format (I12,I12,2E14.6)
3      format (I12,I12,2E25.17)  ! to match the output precision of triangle.
!    3  format(I12,I12,2E14)
33     format (I5,I5,2(1PG14.6)) ! serial old sk00l

write(luo,'(a)') 'end nps'
call closef(luo)
call closef(lui)
write(*,*) 'NUMNP = ',np

!******** Do the elements

call openf(lui,Efile,'old')
read (lui,*,err=990,end=1000,iostat=ios) ne,npoints,iattr

if (npoints.ne.3) then
    write(lue,4)
4        format('3 points per triangle supported')
    call exitp(1)
endif

if (iattr.gt.3) then
    write(lue,5)
5        format('a maximum of 3 attributes per triangle is supported')
endif
  
call openf(luo,'tecin.dat.elm','unknown')

attr(1) = 1
attr(2) = 0
attr(3) = 0

do k=1,ne

   if (iattr.gt.0) then
       read (lui,*,err=990,end=1000,iostat=ios) im,(m(i),i=1,3), (attr(i),i=1,iattr)
   else
       read (lui,*,err=990,end=1000,iostat=ios) im,(m(i),i=1,3)
   endif

   do i=1,3
         m(i) = m(i) + 1
         if (m(i).gt.np) then
            write(lue,7)
7           format('node numbering not zero-based')
            call exitp(1)
         endif
   enddo

   if (writeSerialMode) then 
       write (luo,66) k, attr(1), (m(i),i= 1,3),m(3)
   else
       write (luo,6) k, attr(1), (m(i),i= 1,3),m(3)
   endif
end do

6       format (6I12)
66       format (6I5)


write(luo,'(a)') 'end elm'

call closef(luo)
call closef(lui)

write(*,*) 'NUMEL = ',ne

call exitp(0)

990    write(lue,*) 'read error ',ios
call exitp(1)
1000    write(lue,*) 'EOF error ',ios
call exitp(1)
end
!-----------------------------------------------------------------------
subroutine shell (Nfile,Efile,writeSerialMode)
!
implicit none
!-pass
character*(*) Efile,Nfile
logical :: writeSerialMode
!-locl
character(len=80) usage,arg
integer iflu,stderr,narg,lnblk
external iflu,lnblk
!
!        call HowToUse ('-n TriangleNodeFile -e TriangleElmFile',usage)
stderr = iflu('stderr')
!
Nfile = ''
Efile = ''

narg = command_argument_count()

if (narg.ne.4 .and. narg.ne.5) then
write(stderr,'(80a)') usage(1:lnblk(usage))
call exitp(1)
endif

call get_command_argument(1,arg)

if (arg.eq.'-n') then
call get_command_argument(2,Nfile)
else if (arg.eq.'-e') then
call get_command_argument(2,Efile)
else
write(stderr,'(80a)') usage(1:lnblk(usage))
call exitp(1)
endif

call get_command_argument(3,arg)

if (arg.eq.'-n') then
call get_command_argument(4,Nfile)
else if (arg.eq.'-e') then
call get_command_argument(4,Efile)
else
write(stderr,'(80a)') usage(1:lnblk(usage))
call exitp(1)
endif

if (narg.eq.5) then
call get_command_argument(5,arg)
if (arg.eq.'serial') then
    writeSerialMode = .true.
endif
endif

if (lnblk(Nfile).eq.0 .or. lnblk(Efile).eq.0) then
write(stderr,'(80a)') usage(1:lnblk(usage))
call exitp(1)
endif

return
end
