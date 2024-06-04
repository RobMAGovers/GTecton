   MODULE REGRIDDATA
  integer :: numnewelements
  integer, allocatable :: newelements(:)
  integer, allocatable :: newtriangles(:,:)
  double precision, allocatable :: newstress(:,:)
  double precision, allocatable :: newstrain(:,:)
  double precision, allocatable :: newslpdxe(:,:,:)
  integer, allocatable :: newminnum(:)

character(len=255) MKGRID, NEWELMS, NEWSTR, NEWSTN, NEWDXE, NEWNUM



  contains

  subroutine writeregriddata(NSTR,NEN,NDOF)

#include "regrid.i"
  integer i,NSTR,NEN,NDOF
  integer lu,nextlu,lustn,lustr,ludxe,lunum
  external nextlu

#ifdef EXPERIMENTAL_ReducedSlipAlloc
stop "Regridding in reduced memory mode not yet supported"
#else


  NEWELMS = 'newelms'
  NEWSTN = 'newstress'!tmpnm('newstress')
  NEWSTR = 'newstrain'!tmpnm('newstrain')
  NEWDXE = 'newdxe'!tmpnm('newdxe')
  NEWNUM = 'newnum'!tmpnm('newnum')

  lu = nextlu(0)
  call openf(lu,NEWELMS,'unknown')
  lustn = nextlu(lu+1)
  call openf(lustn,NEWSTN,'unknown')
  lustr = nextlu(lustn+1)
  call openf(lustr,NEWSTR,'unknown')
  if (allocated(newslpdxe)) then
      ludxe = nextlu(lustr+1)
      call openf(ludxe,NEWDXE,'unknown')
      lunum = nextlu(ludxe+1)
      call openf(lunum,NEWNUM,'unknown')
  endif

  do i=1,numnewelements
      write(lu,65) newelements(i), &
       (newtriangles(j,i),j=1,3)
65         format(4I8)
      write(lustn,66) newelements(i), &
        (newstress(j,i),j=1,NSTR)
      write(lustr,66) newelements(i), &
        (newstrain(j,i),j=1,NSTR)
66         format(I8,3E26.18)
      if (allocated(newslpdxe)) then
          write(ludxe,67) newelements(i), &
          ((newslpdxe(k,j,i),j=1,NEN),k=1,NDOF)
67             format(I8,8E26.18)
          write(lunum,68) newelements(i),newminnum(i)
68             format(2I8)
      endif
  enddo

  close(lu)
  close(lustn)
  close(lustr)
  if (allocated(newslpdxe)) then
      close(ludxe)
      close(lunum)
  endif

#endif

  end subroutine
   END MODULE REGRIDDATA
