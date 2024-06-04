!***********************************************************
!         Check whether the domain is one
!         piece. If not, partition tools may
!         run into trouble.
!***********************************************************

module globalVars
    implicit none

! to store file data
    integer                            :: ndimensions
    integer                            :: nvertices   ! number of match points

    integer                            :: nelems
    integer, allocatable, dimension(:,:) :: vertices
    integer                            :: d1,d2    ! dummies

! filenames passed in arguments
    character (len=128)                :: connectivityFileName  ! optional

! level of output
    integer                            :: verbose

! subs


end module

!***********************************************************
!           Main program  
!***********************************************************

program domainUnity

    use globalVars, only: verbose, &
                          connectivityFileName, &
                          vertices, &
                          d1,d2

    implicit none

    integer :: nelems, iElem
    integer :: nPoints, iPoint, updatePoint

    integer, allocatable :: reach(:)
    logical :: reachChanged
    logical :: oneChunk
    integer :: iUpdate, iTest, nUpdates1, nUpdates2


! commmand line argument variables are all in the globalVars module, hence no arguments
! (line numbered to possible skip or insert testing stuff)

    call readcommandlineargs()
    call validatecommandlineargs()

!***** read tecin.dat.elm ***
    if (verbose.eq.1) then
        write(0,*) "Reading connectivity file",trim(connectivityFileName)
    endif

    call countlines(trim(connectivityFileName),nelems, 32)
    nelems = nelems - 1 ! remove the 'end' line

    if (nelems.eq.-1) then
        write(0,*) "Partitioned element file empty: ", trim(connectivityFileName)
    endif


    if (verbose.eq.1) then
        write(0,*) "nlines in tecin elem ", nelems
    endif

    allocate(vertices(4,nelems))

    open(unit = 233, file=trim(connectivityFileName))
    do iElem=1,nelems
        read (233,*) d1, d2, vertices(:,iElem)
    enddo
    close(233)        

    if (verbose.eq.1) then
        write(0,*) 'Read ', nelems, ' elements'
    endif


!* create array and fill it. Start at nodal point 1

    ! find highest noe number
    nPoints = maxval(vertices)


    allocate(reach(nPoints))
    reach = 0
    reach(1) = 1

    ! the reach array will get a nonzero for every node that can can be reached
    ! from point 1, by stepping through the elements
    ! At some point, during this iteration, this array will no longer
    ! change. It will be either full of nonzeros then (the domain is 
    ! a single connected block) or there will still be zeros
    ! (the domain is split into multiple disconnected subdomains)

!    write(*,*) "nmeuk", nPoints, nelems

    nUpdates1=0
    nUpdates2=0

    reachChanged = .true.
    do while (reachChanged)
!        write(*,*) "A new round!", nUpdates1, nUpdates2
        nUpdates1=0
        nUpdates2=0
        reachChanged = .false.
        do iPoint=1,nPoints
            
            if (reach(iPoint).eq.2) then
                  	! this points is connected to point 1,
                ! but all connecting points have already
                ! been registered. Do nothing

            else if (reach(iPoint).eq.1) then
                ! this points is connected to point 1,
                ! see whether more points are connected.
                ! it is on the edge of current chunk.
                ! find new neighbors
                do iTest=1,4
                    do iElem=1,nelems
                        if (vertices(iTest,iElem).eq.iPoint) then
                            ! found an element that contains this point
                            do iUpdate = 1,4
                                updatePoint = vertices(iUpdate,iElem)
                                if (reach(updatePoint).eq.0) then
                                    ! new point to the chunk
                                    reach(updatePoint) = 1
                                    nUpdates1 = nUpdates1 + 1
                                    reachChanged = .true.
                                else
                                    ! the point was already included
                                endif
                            enddo

                        endif
                    enddo
                enddo
                reach(iPoint) = 2
                nUpdates2 = nUpdates2 + 1
                reachChanged = .true.
            else
                  	! this points is not yet connected to the domain

            endif

        enddo
    enddo

    ! the whole chunk attached to point 1 has been marked.
    ! see if there are any empty points

    oneChunk=.true.

    do iPoint=1,nPoints
        if (reach(iPoint).eq.0) then
            write(*,*) "Point", iPoint, "is not connected to point 1"
            oneChunk=.false.
        endif
    enddo

    if (oneChunk) then
        write(*,*) "all Points forms a single connected chunk"
    else
           write(*,*) "domain consists of multiple chunks."
    endif

end program


!***** subroutines for reading the command line arguments

subroutine readcommandlineargs()

    use globalVars, only: verbose,&
                          connectivityFileName

    implicit none

    character (len=128) :: arg

    logical, external :: stringToLogical
    integer, external :: stringToInteger

    integer :: nargs, iarg

    ! set defaults
    verbose = 0

    connectivityFileName = ''

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

        call get_command_argument(iarg , arg)      ! tecin.dat.partf.nps

        if (trim(arg).eq."-e") then    
            iarg = iarg + 1
            call get_command_argument(iarg , arg)
            connectivityFileName = trim(arg)

        ! integers
        else if (trim(arg).eq."-v") then
            verbose = 1
        else
            write(*,*) "Command line argument ", trim(arg), \
                       " not recognized"
        endif

        iarg = iarg + 1

    enddo


    if (verbose.gt.0) then
        write(*,*) "continuing with command line arguments:"
        write(*,*) "element file:                   ", connectivityFileName
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

    use globalVars, only: connectivityFileName

    implicit none

    if (connectivityFileName.eq.'') then
        write(*,*) 'Element file not supplied; -e <filename>; exiting'
        stop
    endif

end subroutine


subroutine countlines(filename, nlines, filehandle)
    use globalVars, only: verbose

    implicit none

    character(len=*)   :: filename
    integer            :: nlines, filehandle, ierr

    nlines = 0
    open(unit=filehandle, file=trim(filename), iostat=ierr)

    if (verbose.eq.1) then
        write(0,*) "file ", trim(filename), " openend with status ", ierr
    endif

    if (ierr.ne.0) then
        write(0,*) "elmside wanted to count lines in ", trim(filename)
        write(0,*) "but could not open it, error", ierr
    endif


    do while (.true.)
        read (filehandle,*, end=10)
        nlines = nlines + 1
    enddo
10  close(filehandle)

end subroutine
