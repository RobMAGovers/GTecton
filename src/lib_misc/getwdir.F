! called only by Plt3d/xeclib3D.F

subroutine getwdir(path)

implicit none

!    overlay routine for call to C-routine

!-pass
character*(*) :: path
!-local
integer       :: n

n = LEN(path)

if (n.le.0) then
    return
endif

call cwd(path)
return
end subroutine
