subroutine triostat (ierr)

! g77 iostat runtime error messages

implicit none
integer ierr,lue,iflu
character(len=35) list
dimension list(100:131)
external iflu
data list/ &
          'error in format',&
          'illegal unit number',&
          'formatted io not allowed',&
          'unformatted io not allowed',&
          'direct io not allowed',&
          'sequential io not allowed',&
          'can not backspace file',&
          'null file name',&
          'can not stat file',&
          'unit not connected',&
          'off end of record',&
          'truncation failed in endfile',&
          'incomprehensible list input',&
          'out of free space',&
          'unit not connected',&
          'read unexpected character',&
          'bad logical input field',&
          'bad variable type',&
          'bad namelist name',&
          'variable not in namelist',&
          'no end record',&
          'variable count incorrect',&
          'subscript for scalar variable',&
          'invalid array section',&
          'substring out of bounds',&
          'subscript out of bounds',&
          'can not read file',&
          'can not write file',&
          '"new" file exists',&
          'can not append to file',&
          'non-positive record number',&
          'I/O started while already doing I/O'/

if (ierr.ge.100.and.ierr.le.131) then
    lue = iflu('stderr')
    write(lue,'(1x,a35)') list(ierr)
endif
return
end
