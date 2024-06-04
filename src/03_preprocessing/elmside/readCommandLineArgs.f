subroutine readcommandlineargs()

    use globals, only: opn, &
                       printcoords, &
                       sideonly, &
                       printweights, &
                       invertsign, &
                       tecindatpartfnps, &
                       tecindatpartfelm, &
                       loopfile, &
                       vertexlist, &
                       ndimensions, &
                       endpoints, &
                       algorithm, &
                       iecho

    implicit none

    character (len=128) :: arg

    logical, external :: stringToLogical
    integer, external :: stringToInteger

    integer :: nargs, iarg

    ! set defaults
    opn = .false.
    printcoords = .false.
    sideonly = .false.
    printweights = .false.
    invertsign = .false.

    tecindatpartfnps = ''
    tecindatpartfelm = ''
    vertexlist = ''
    loopfile = ''

    algorithm = 1
    iecho = 0

    ! no default for the vertexlist

    nargs = command_argument_count()





    ! fortran parses command line arguments as
    ! elmside a=something b=anotherthing
    ! while c, and the unix standard is
    ! elmside -a something -b anotherthing.

    ! the issue is that in Fortran, 
    ! -a something
    ! is seen as two command line arguments, so in parsing,
    ! we have do two steps in one

    iarg = 1

    do while (iarg.le.nargs)

        call get_command_argument(iarg , arg)       ! tecin.dat.partf.nps

!           write(*,*) "argument: ", trim(argQualifier), " with ", trim(argContent)
            ! strings (filenames)
        if (trim(arg).eq."-n") then
            iarg = iarg + 1
            call get_command_argument(iarg , arg)

            tecindatpartfnps = trim(arg)

        else if (trim(arg).eq."-e") then
            iarg = iarg + 1
            call get_command_argument(iarg , arg)
            tecindatpartfelm = trim(arg)
        else if (trim(arg).eq."-l") then
            iarg = iarg + 1
            call get_command_argument(iarg , arg)
            vertexlist = trim(arg)

        else if (trim(arg).eq."-L") then
            iarg = iarg + 1
            call get_command_argument(iarg , arg)
            loopfile = trim(arg)

        ! integers
        else if (trim(arg).eq."-d") then
            iarg = iarg + 1
            call get_command_argument(iarg , arg)
            ndimensions = stringToInteger(trim(arg))

        else if (trim(arg).eq."-a") then
            iarg = iarg + 1
            call get_command_argument(iarg , arg)
            algorithm = stringToInteger(trim(arg))

        else if (trim(arg).eq."-v") then
            iecho = 1

       ! logicals/booleans
        else if (trim(arg).eq."-w") then
            printweights = .true.
        else if (trim(arg).eq."-c") then
            printcoords = .true.
        else if (trim(arg).eq."-o") then
            opn = .true.
        else if (trim(arg).eq."-s") then
            sideonly = .true.
        else if (trim(arg).eq."-i") then
            invertsign = .true.
        else if (trim(arg).eq."-p") then
            endpoints = .true.
        else
            write(*,*) "Command line argument ", trim(arg), \
                       " not recognized"
        endif

        iarg = iarg + 1

    enddo


    if (iecho.gt.0) then
        write(*,*) "continuing with command line arguments:"
        write(*,*) "partitioned nodal points file:  ", tecindatpartfnps
        write(*,*) "partitioned element file:       ", tecindatpartfelm
        if (algorithm.eq.3) then
            write(*,*) "loopfile                        ", loopfile
        endif
        write(*,*) "vertexlist                      ", vertexlist
        write(*,*) "number of dimensions            ", ndimensions
        write(*,*) "(verbosity) echo                ", iecho
        write(*,*) "add weights to output:          ", printweights
        write(*,*) "add coordinates to output:      ", printcoords
        write(*,*) "output in opn format:           ", opn
        write(*,*) "give only the side numbers:     ", sideonly
        write(*,*) "invert the sign of the weights: ", invertsign
        write(*,*) "write the endpoints:            ", endpoints
    endif

end subroutine

logical function stringToLogical(string)

implicit none

character(len=128) :: string

if ((string(1:1).eq.'1').or.&
    (string(1:4).eq.'true').or.&
    (string(1:1).eq.'t').or.&
    (string(1:1).eq.'y').or.&
    (string(1:3).eq.'yes')) then
    stringToLogical=.true.
else
    stringToLogical=.false.
endif

end function


integer function stringToInteger(string)

implicit none

character(len=1) :: string

!write(*,*) " converting ", string(1:1)

read(string(1:1),*) stringToInteger

!write(*,*) " to integer ", stringToInteger

end function



subroutine validatecommandlineargs()

    use globals, only: ndimensions, &
                       tecindatpartfnps, &
                       tecindatpartfelm, &
                       vertexlist, &
                       algorithm, loopfile

    implicit none

    if ((ndimensions.ne.2) .and. (ndimensions.ne.3)) then
        write(*,*) 'dimensions should be 2 or 3; -d 2 or -d 3; exiting'
        stop
    endif

    if (tecindatpartfnps.eq.'') then
        write(*,*) 'Nodal point file not supplied; -n <filename>; exiting'
        stop
    endif

    if (tecindatpartfelm.eq.'') then
        write(*,*) 'Element file not supplied; -e <filename>; exiting'
        stop
    endif

    if (vertexlist.eq.'') then
        write(*,*) 'Vertex list file not supplied; -l <filename>; exiting'
        stop
    endif

    if (algorithm.eq.3 .and. loopfile.eq.'') then
        write(*,*) 'loop file not supplied; -L <filename>; exiting'
        stop
    endif

end subroutine

