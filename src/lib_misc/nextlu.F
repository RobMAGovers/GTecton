integer function nextlu (lu0)

! Gives back an unused logical unit number
! Search starts at lu0

implicit none

!-pass
integer :: lu0, lumax
integer :: numtry
!-locl
logical :: isOpen
integer :: stderr,stdin,stdout,lu

integer, external :: iflu

stdin  = iflu('stdin')
stdout = iflu('stdout')
stderr = iflu('stderr')

lumax  = 1000
numtry = 0

lu = MOD(abs(lu0),lumax)


! ***** keep iterating through numbers until we found a 
! ***** positive one that has not yet been opened,
! ***** or until we have reched a 100, indicating that
! ***** there are too many open files.

100 if (lu.le.0 .or. lu.ge.lumax) then
    lu = 1
endif

if (lu.ne.stdin  .and. &
    lu.ne.stdout .and. &
    lu.ne.stderr) then

    inquire (unit=lu,opened=isOpen)

    if (.not.isOpen) then
        ! Yes, got free one! Let's get out of here
        goto 300
    endif
else
!    write(*,*) 'same as a standard io'
endif

lu = lu + 1
numtry = numtry + 1


if (numtry.gt.lumax) then
    write(stderr,1) lumax
1        format(1x,'nextlu: more than ',I3,' open files')
    call exitp(1)
endif

! bugger, file handler was alreay occupied. Try the next one.
goto 100

! whiiiiiiii :-)
300 nextlu = lu


return
end
