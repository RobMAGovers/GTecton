integer function iflu (arg)

implicit none

character*(*) :: arg
integer       :: stderr,stdin,stdout

stderr = 0
stdin  = 5
stdout = 6

if (arg.eq.'stderr' .or. arg.eq.'STDERR') then
    iflu = stderr
else if (arg.eq.'stdin' .or. arg.eq.'STDIN') then
    iflu = stdin
else if (arg.eq.'stdout' .or. arg.eq.'STDOUT') then
    iflu = stdout
else
    write(stderr,1)
1        format(1x,'iflu: argument error')
    iflu = -1
endif

return
end

