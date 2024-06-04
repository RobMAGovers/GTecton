integer function iprinfo(cmd,ncmd)

use documentation, only: plt3d_manual

implicit none

!     Print command info

!-pass
integer ncmd
character cmd*(*)
dimension cmd(ncmd)


!-locl
#include "lenchar.i"

if (ncmd.le.1) then
    !       single "?": give list of commands
    call plt3d_manual('list')
else
    !       information on a specific command requested
    call plt3d_manual(trim(cmd(2)))
endif

return
end
