subroutine wGMTfile (pdata,ndata,itype,lu,ascii,lue,ierr)

! Routine for writing a "GMT" data file to file lu

implicit none

!-pass
integer ndata,itype,lu,lue,ierr
logical ascii
double precision :: pdata
dimension pdata(*)

!-init
ierr = 0

!    Perform scaling on small and large numbers
!    call wGMTreal (pdata,ndata,itype)

if (itype.eq.4) then
    call wPSXYMULTI (pdata,ndata,lu,ierr)
else if (itype.eq.7) then
    call wPSCONTOUR (pdata,ndata,lu,ascii,lue,ierr)
else 
    write(lue,1) itype
1        format(1x,'wPfile: type ',I4,' data not implemented')
    ierr = 1
endif

return
end

!-------------------------------------------------------------------------------

subroutine wPSXYMULTI (pdata,ndata,lu,ierr)

! Routine for writing gmt data for PSXY multisegment file 

implicit none
!-pass
integer ndata,lu,ierr
double precision :: pdata
dimension pdata(3,*)
!-locl
logical ready
integer :: j,k,n,m
!-init
ready  = .false.

if (ndata.le.0) return

n = 1

100    m = n + 3

if (m.eq.ndata) then
    ready = .true.
endif

if (m.gt.ndata) goto 400


write(lu,4,err=400)
4    format('> ',1x)

write(lu,6,err=400) ((pdata(j,k),j=1,2),k=n,m)
6    format(E15.6,1X,E15.6)

n = m + 1

if (.not.ready) goto 100

return

400    call prerr('wPSXYMULTI')

ierr = 1

return
end

    
!-------------------------------------------------------------------------------
subroutine wPSCONTOUR (pdata,ndata,lu,ascii,lue,ierr)

! Routine for writing data for gmt command pscontour 

implicit none
!-pass
integer ndata,lu,lue,ierr
logical ascii
double precision :: pdata
dimension pdata(3,*)
!-local
integer :: j,k


if (ndata.le.0) then
    write(*,*) 'wPSCONTOUR has no data to write. Leaving'
    return
endif

!    Check if data complete
if (MOD(ndata,3).ne.0) then
    write(lue,2) ndata
2        format(1x,'wGMTCONTOUR: ndata not multiple of 3; ndata = ',I5)
    ierr = 1
    return
endif

!    Count number of triangles in data

if (lu.eq.0) then
    write(*,*) 'wPSCONTOUR has file handle 0. Cannot write data'
endif

if (ascii) then
    write(lu,6,err=200) ((REAL(pdata(j,k)),j=1,3),k=1,ndata)
6        format(E15.6,1X,E15.6,1X,E15.6,1X)
else
   write(lu,err=200) ((REAL(pdata(j,k)),j=1,3),k=1,ndata)
endif


return

200    call prerr('wPSCONTOUR')
ierr = 1
return
end
!-------------------------------------------------------------------------------
subroutine wGMTreal (pdata,ndata,itype)

! scales data within limits

implicit none
!-pass

integer ndata,itype
double precision :: pdata
dimension pdata(*)

! todo huh?

return
end
