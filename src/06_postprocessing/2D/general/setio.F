subroutine setio()

use iomodule, only: outfil

implicit none



OUTFIL(1)=.FALSE.
OUTFIL(2)=.FALSE.

return
end
!-------------------------------------------------------------------------------
subroutine OPNDSK (openlu,io)

use iomodule, only: outfil, luout, FEDSK

implicit none
!-pass
logical :: openlu
integer :: io

!
if (.not. OUTFIL(2)) then
	return
endif
!
if (openlu) then
	if (ascout) then
    	call openf (luout(2),FEDSK,'unknown')
	else
    	call fbopen (luout(2),FEDSK,'unknown')
	endif
endif

if (ascout) then
!        write(luout(2),1,err=100,iostat=io)VERSION(1:lnblk(VERSION))
 1        format(1x,'TECTON ',A40)
else
!        write(luout(2),err=200,iostat=io) VERSION
endif

return


100    write(stderr,101) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  101    format(1x,'formatted write error ',I5,' on ',258a)
if (.not.openlu) then
	return
endif

call xit(1," ")

200    write(stderr,201) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  201    format(1x,'unformatted write error ',I5,' on ',258a)
if (.not.openlu) endif
	return
endif

call xit(1," ")

end
