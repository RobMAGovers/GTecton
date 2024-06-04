subroutine readInput

use globals

implicit none


integer :: i,j


!****** read node list *****

    if (iecho.eq.1) then
        write(0,*) "Reading vertex list ",trim(vertexlist)
    endif

    call countlines(vertexlist,nvertices, 22)
    allocate(vertices(nvertices))

    if (iecho.eq.1) then
        write(0,*) "nlines in node list ", nvertices
    endif

    if (nvertices.eq.-1) then
        write(0,*) "nodelist file empty: ", trim(vertexlist)
    endif

    open(unit = 23, file=trim(vertexlist))
    do i=1,nvertices
        read (23,'(I12)') vertices(i)
    enddo
    close(23)

    if (nvertices.eq.0) then
        write(0,*) 'vertexlist empty. Please supply vertices to match.'
        stop
    endif

    if (nvertices.eq.1) then
        write(0,*) 'vertexlist contains only a single point. This does not define a side.'
        stop
    endif

    if (iecho.eq.1) then
        write(0,*) 'Read ', nvertices, ' vertices'
    endif


!***** read tecin.dat.part.elm ***
    if (iecho.eq.1) then
        write(0,*) "Reading partitioned element file ",trim(tecindatpartfelm)
    endif

    call countlines(trim(tecindatpartfelm),nelems, 32)
    nelems = nelems - 1 ! remove the 'end' line

    if (nelems.eq.-1) then
        write(0,*) "Partitioned element file empty: ", trim(tecindatpartfelm)
    endif


    if (iecho.eq.1) then
        write(0,*) "nlines in tecin elem ", nelems
    endif

    allocate(v(4,nelems))

    open(unit = 233, file=trim(tecindatpartfelm))
    do i=1,nelems
        read (233,*) d1, d2, d3, v(1,i), v(2,i), v(3,i), v(4,i)
    enddo
    close(233)

    if (iecho.eq.1) then
        write(0,*) 'Read ', nelems, ' elements'
    endif

!***** read tecin.dat.part.nps ***
    if (iecho.eq.1) then
        write(0,*) "Reading partitioned vertex file ",trim(tecindatpartfnps)
    endif

    call countlines(trim(tecindatpartfnps),nallvertices, 42)
    nallvertices = nallvertices - 1 ! -1 to remove the 'end' line

    if (nallvertices.eq.-1) then
        write(0,*) "Partitioned nodal point file empty: ", trim(tecindatpartfnps)
    endif


    if (iecho.eq.1) then
         write(0,*) "nlines in tecin nodes ", nallvertices
    endif

    d1 = maxval(vertices)
    if (d1 .gt. nallvertices) then
        write(0,*) "vertex selection list goes up to index", d1
        write(0,*) "while mesh has only ",nallvertices, "vertices"
        write(0,*) "Exiting..."
        stop
    endif

    if (ndimensions.eq.2) then
        allocate(coords(2,nallvertices))
        allocate(neighbourIDs(nallvertices, maxneighbours2D))
    else
        allocate(neighbourIDs(nallvertices, maxneighbours3D))
        allocate(coords(3,nallvertices))
    endif
    allocate(nneighbours(nallvertices))

    if (algorithm.eq.3) then
        ! ***** loopfile ***
        if (iecho.eq.1) then
            write(0,*) "Reading loop file ",trim(loopfile)
        endif
        call countlines(trim(loopfile),numloop, 50)
        if (numloop.eq.0) then
            write(0,*) "loop file empty: ", trim(loopfile)
        endif
        if (iecho.eq.1) then
             write(0,*) "loop consists of ", numloop," vertices"
        endif
        allocate(loopcoords(2,numloop))
        open(unit = 50, file=trim(loopfile))
        do i=1,numloop
            read (50,*) loopcoords(1,i),loopcoords(2,i)
        enddo
        close(50)
    endif


!***** determine extrema of the domain ****

    xmin =  1e10
    ymin =  1e10
    zmin =  1e10
    xmax = -1e10
    ymax = -1e10
    zmax = -1e10


    open(unit = 43, file=trim(tecindatpartfnps))
    if (ndimensions.eq.2) then
        do i=1,nallvertices
               ! dummies are partition, index, and marker
            read (43,*) d1, d2, d3, coords(1,i), coords(2,i),   nneighbours(i), (neighbourIDs(i,j),j=1,nneighbours(i))

            if (coords(1,i).lt.xmin) then
                xmin = coords(1,i)
            endif
            if (coords(2,i).lt.ymin) then
                ymin = coords(2,i)
            endif

            if (coords(1,i).gt.xmax) then
                xmax = coords(1,i)
            endif
            if (coords(2,i).gt.ymax) then
                ymax = coords(2,i)
            endif

        enddo
    else
        do i=1,nallvertices
            ! dummies are partition, index, and marker
            read (43,*) d1, d2, d3, coords(1,i), coords(2,i), coords(3,i),  nneighbours(i), (neighbourIDs(i,j),j=1,nneighbours(i))

            if (coords(1,i).lt.xmin) then
                xmin = coords(1,i)
            endif
            if (coords(2,i).lt.ymin) then
                ymin = coords(2,i)
            endif
            if (coords(3,i).lt.zmin) then
                zmin = coords(3,i)
            endif

            if (coords(1,i).gt.xmax) then
                xmax = coords(1,i)
            endif
            if (coords(2,i).gt.ymax) then
                ymax = coords(2,i)
            endif
            if (coords(3,i).gt.zmax) then
                zmax = coords(3,i)
            endif

        enddo
    endif
    close(43)

    if (iecho.eq.1) then
        write(0,*) 'Read ', nallvertices, ' vertices'
        write(0,*) 'Finished reading input data'
    endif


end subroutine

subroutine countlines(filename, nlines, filehandle)
    use globals, only: iecho
    implicit none
    character(len=*)   :: filename
    integer            :: nlines, filehandle, ierr

    nlines = 0
    open(unit=filehandle, file=trim(filename), iostat=ierr)

    if (iecho.eq.1) then
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

