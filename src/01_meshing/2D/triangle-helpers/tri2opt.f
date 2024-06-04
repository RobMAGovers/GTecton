	program triangle2opt

c rewrite output of "triangle" in TOAST format for "optimize"
c 2/1997 Paul Meijer
c Modified 3/97 by Govers
c convention: external boundary points marked with an 1
c
	character*80 fna
	dimension m(3),attr(3)

	lue = iflu('stderr')
	read (*,'(a)') fna

c nodes
	open (10,file= fna(1:lnblnk(fna))//'.node')
	rewind 10

	read (10,*) np,idim,iattr,ibound
	if (iattr.gt.0) write(lue,1)
    1	format('node attributes not carried along')
	if (ibound.eq.0) then
	    write(lue,2)
    2	    format('external boundary nodes should be labeled with 1')
	    call exitp(1)
	endif
	
	write (*,'(a)') 'MeshData 4.0'
	write (*,'(/,"NodeList ",i4,2x,i1)') np,2

	do k= 1,np
	    read (10,*) ip,x,y,imark
	    if (imark.eq.1) then
		write (*,3) 2,'B',x,y,imark
    3		format (i1,a1,2x,1PG14.6,2x,1PG14.6,2x,I3)
	    else
		write (*,3) 2,'N',x,y,imark
	    endif
	end do
	close (10)

c elements
	open (10,file= fna(1:lnblnk(fna))//'.ele')
	rewind 10

	read (10,*) ne,npoints,iattr
	if (npoints.ne.3) then
	    write(lue,4)
    4	    format('3 points per triangle supported')
	    call exitp(1)
	endif
	if (iattr.gt.3) then
	    write(lue,5)
    5	    format('a maximum of 3 attributes per triangle is supported')
	endif

	write (*,'(/,"ElementList ",i4,2x,i1)') ne,2

	attr(1) = 0.
	attr(2) = 0.
	attr(3) = 0.
	do k= 1,ne
	   if (iattr.gt.0) then
	       read (10,*) im,(m(i),i= 1,3),(attr(i),i=1,iattr)
	   else
	       read (10,*) im,(m(i),i= 1,3)
	   endif
	   do i=1,3
	     m(i) = m(i) + 1
             if (m(i).gt.ne) then
		write(lue,7)
    7		format('node numbering not zero-based')
		call exitp(1)
	     endif
	   enddo
	   write (*,6) 'a',(NINT(attr(i)),i=1,3),(m(i),i= 1,3)
    6	   format (a1,3(2x,I5),3(2x,i4))
	end do

	write (*,'(/,"PointSource  0")')
	write (*,'("0  0")')
	write (*,'(/,"TimeSpec  0  0  0  0")')

	end
