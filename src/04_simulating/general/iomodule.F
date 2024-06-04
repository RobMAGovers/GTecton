MODULE IOMODULE

character(len=255) :: FEIN, FEOUT, TEMPIN, FEDSKDIR
logical :: OUTFIL(2), ASCIN, ASCOUT
logical :: optijd
integer :: stdin, stdout, stderr, lumon
integer :: luout(2), lures, lucnst, fedskrt
integer :: luin
character(len=255) :: FEDSK


logical :: append

integer :: lufo

! whether to write mechanical output or temperature output
logical :: Mout,Tout

! (only used in plt3d)
integer :: icoord ! to determine whether points have been plotted, 
                  ! and in what way.
                  ! 0 = not yet done
                  ! 1 = written in binary
                  ! 2 = written in ascii

integer :: ifein 
integer :: ifeout
integer :: ifelm 
integer :: ifctr 

character(len=256) :: crdfil,elmfil,ctrfil


logical :: verbose ! whether to give extra output.

! variables that handle the command line:

integer, parameter :: MaxStack=10

character(len=512) :: fstack(0:MaxStack),fcmd

logical            :: bufedit
integer            :: lucmd(0:MaxStack), lubuf, ilu, lupos(0:MaxStack)



END MODULE IOMODULE
