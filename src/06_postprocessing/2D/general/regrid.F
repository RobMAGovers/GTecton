   subroutine regrid (IEN,X,D,MAT,NSLIP,MATNO,NP,CON,MAXCON,ierr)
   USE REGRIDDATA
 USE MODELDEFINITION, only: NUMNP




!
!----------------------------------------------------------------------
! creates Delauney tesselation of elements in a domain with material number
! MATNO. New elements are written to ./newelms
!
! SNGL implicit real (a-h,o-z)
   implicit double precision (a-h,o-z)
   parameter(third=1.d0/3)



#ifdef EXPERIMENTAL_ReducedSlipAlloc
stop "Regridding not yet supported in reduced memory mode. Leaving"
#else

!-pass
   integer MAXCON,MATNO,ierr
   dimension X(NSD,*),D(NDOF,*),IEN(NEN,*),MAT(*),NSLIP(5,*)
   integer NP(NUMNP) ! list of unique node numbers that need be tesselated
   integer CON(2,MAXCON) ! list of region boundary connections
!-locl
!#include "lenchar.i"
   logical hit1,hit2,hit3,quad
   integer me,mnp,mcn,i,m,j,k,num1,num2,num3,lu,nextlu!,csystem,
!     >   rien(3)
!     some storage and counters for materials other than MATNO
   double precision HOLEPOINT(NSD,20)
   integer          HOLEPASSED(2,20), ii,jj,kk

!      character(len=255 command
   external nextlu!,csystem
!      character(len=)tmpnm*255
!      integer lnblk
   integer ERROR
   TYPE(triangulateio) :: tri_in, tri_out
   character(len=40) triswitch, MKGRID
!-init
!      ierr = OK
   MKGRID = 'mkgrid'!tmpnm('mkgrid')
!      NEWELMS = 'newelms'!tmpnm('newelms')
!
   me = 0
   mnp = 0
   mcn = 0
   do i=1,20
  do j=1,NSD
      HOLEPOINT(j,i) = .0d0
  enddo
  HOLEPASSED(1,i)=0
  HOLEPASSED(2,i)=0
   enddo
   ii = 1
   do i=1,NUMAT
    if (i.ne.MATNO) then
        HOLEPASSED(1,ii) = i
        ii = ii+1
    endif
   enddo

   do i=1,NUMEL
  m = MAT(i)
  if (m.eq.MATNO) then
      me = me + 1
      if (IEN(3,i).ne.IEN(4,i)) then
          write(stderr,10)
10             format(1x,'regrid cannot work on quads')
          ierr = EXEC
          return
      endif
!       add element nodes to list if not already in there
      do 100 j=1,3
          do k=1,mnp
              if (NP(k).eq.IEN(j,i)) goto 100
          enddo
          mnp = mnp + 1
          if (mnp.gt.NUMNP) then
              write(stderr,20) NUMNP
20                 format(1x,'regrid: overflow on NUMNP=',I8)
              ierr = EXEC
              return
          endif
          NP(mnp) = IEN(j,i)
100           continue
!
!       Boundary test
      hit1 = .false.
      hit2 = .false.
      hit3 = .false.
      do 200 j=1,NUMEL
          if (i.eq.j) goto 200
!      if (MAT(j).eq.m) goto 200
          if (IEN(3,j).ne.IEN(4,j)) then
              write(stderr,10)
              ierr = EXEC
              return
          endif
          num1 = 0
          num2 = 0
          num3 = 0
          do jj=1,3
              if (IEN(1,i).eq.IEN(jj,j)) then
                  num1 = num1 + 1
                  num3 = num3 + 1
              endif
              if (IEN(2,i).eq.IEN(jj,j)) then
                  num1 = num1 + 1
                  num2 = num2 + 1
              endif
              if (IEN(3,i).eq.IEN(jj,j)) then
                  num2 = num2 + 1
                  num3 = num3 + 1
              endif
          enddo
          if (num1.eq.2.and.MAT(j).ne.m) then
!          material boundary between nodes 1 and 2
              mcn = mcn + 1
              if (mcn.gt.MAXCON) then
                  write(stderr,30) MAXCON
30                     format(1x,'regrid: overflow on MAXCON=',I8)
                  ierr = EXEC
                  return
              endif
              CON(1,mcn) = IEN(1,i)
              CON(2,mcn) = IEN(2,i)
          endif
          if (num2.eq.2.and.MAT(j).ne.m) then
!          material boundary between nodes 2 and 3
              mcn = mcn + 1
              if (mcn.gt.MAXCON) then
                  write(stderr,30) MAXCON
                  ierr = EXEC
                  return
              endif
              CON(1,mcn) = IEN(2,i)
              CON(2,mcn) = IEN(3,i)
          endif
          if (num3.eq.2.and.MAT(j).ne.m) then
!          material boundary between nodes 3 and 1
              mcn = mcn + 1
              if (mcn.gt.MAXCON) then
                  write(stderr,30) MAXCON
                  ierr = EXEC
                  return
              endif
              CON(1,mcn) = IEN(3,i)
              CON(2,mcn) = IEN(1,i)
          endif
          hit1 = (hit1.or.(num1.eq.2))
          hit2 = (hit2.or.(num2.eq.2))
          hit3 = (hit3.or.(num3.eq.2))
200           continue
      if (.not.hit1) then
!    side 1 is external
          mcn = mcn + 1
          if (mcn.gt.MAXCON) then
              write(stderr,30) MAXCON
              ierr = EXEC
              return
          endif
          CON(1,mcn) = IEN(1,i)
          CON(2,mcn) = IEN(2,i)
      endif
      if (.not.hit2) then
!    side 2 is external
          mcn = mcn + 1
          if (mcn.gt.MAXCON) then
              write(stderr,30) MAXCON
              ierr = EXEC
              return
          endif
          CON(1,mcn) = IEN(2,i)
          CON(2,mcn) = IEN(3,i)
      endif
      if (.not.hit3) then
!    side 3 is external
          mcn = mcn + 1
          if (mcn.gt.MAXCON) then
              write(stderr,30) MAXCON
              ierr = EXEC
              return
          endif
          CON(1,mcn) = IEN(3,i)
          CON(2,mcn) = IEN(1,i)
      endif
  else
      do ii=1,NUMAT-1
          if (HOLEPASSED(1,ii).eq.m.and. &
           HOLEPASSED(2,ii).eq.0) then
              HOLEPASSED(1,ii) = m
              HOLEPASSED(2,ii) = i
              HOLEPOINT(1,ii) = third*( & !0.3333*( &
              X(1,IEN(1,i))+D(1,IEN(1,i)) + &
              X(1,IEN(2,i))+D(1,IEN(2,i)) + &
              X(1,IEN(3,i))+D(1,IEN(3,i)))
              HOLEPOINT(2,ii) = third*( & !0.3333*( &
              X(2,IEN(1,i))+D(2,IEN(1,i)) + &
              X(2,IEN(2,i))+D(2,IEN(2,i)) + &
              X(2,IEN(3,i))+D(2,IEN(3,i)))
              if (iecho.eq.1) then
                  write(stdout,*) 'Found other material elt:'
                  write(stdout,*) 'Element: ', i, &
                   ' with material: ', m
                  write(stdout,*) 'Nodal point numbers: ', &
                   IEN(1,i), IEN(2,i), IEN(3,i)
                  write(stdout,*) 'Nodal point coordinates: '
                  write(stdout,*) X(1,IEN(1,i))+D(1,IEN(1,i)), &
              X(2,IEN(1,i))+D(2,IEN(1,i))
                  write(stdout,*) X(1,IEN(2,i))+D(1,IEN(2,i)), &
              X(2,IEN(2,i))+D(2,IEN(2,i))
                  write(stdout,*) X(1,IEN(3,i))+D(1,IEN(3,i)), &
              X(2,IEN(3,i))+D(2,IEN(3,i))
              endif
          endif
      enddo
  endif
   enddo
!
!     BEGIN new triangle interface
   triswitch = 'pQ'
   call init_triangle_input(triswitch(1:lnblk(triswitch)), tri_in, &
  mnp, 0,0,0,0,mcn, NUMAT-1, 0)
!      call init_triangle_output(triswitch(1:lnblk(triswitch)),
!     .  tri_in,tri_out)
   call init_triangle_output_regrid(triswitch(1:lnblk(triswitch)), &
tri_in,tri_out,me)
   do i=1,tri_in%numberofpoints
  xup = X(1,NP(i))+D(1,NP(i))
  yup = X(2,NP(i))+D(2,NP(i))
  tri_in%pointlist(2*i-1) = xup
  tri_in%pointlist(2*i) = yup
   enddo
   do i=1,tri_in%numberofsegments
  j = 1
  do while (j.le.mnp .and. NP(j).ne.CON(1,i))
      j = j + 1
  enddo
  k = 1
  do while (k.le.mnp .and. NP(k).ne.CON(2,i))
      k = k + 1
  enddo
  tri_in%segmentlist(2*i-1) = j!-1
  tri_in%segmentlist(2*i) = k!-1
   enddo
   do i=1,tri_in%numberofholes
  tri_in%holelist(2*i-1) = HOLEPOINT(1,i)
  tri_in%holelist(2*i)   = HOLEPOINT(2,i)
   enddo
!      write(stdout,*) 'nu potentieel wegschrijven'
!     write the .poly to disk for inspection / debugging purposes
   if (iecho.eq.1) then
  write(stdout,*) 'writing triangle file ', &
   MKGRID(1:lnblk(MKGRID)),'.poly to disk'
!          lu = nextlu(0)
!          call openf (lu,MKGRID(1:lnblk(MKGRID))//'.poly','unknown')
!          write(lu,*) tri_in%numberofpoints,2,0,0
!          do i=1,tri_in%numberofpoints
!              write(lu,*) i-1,tri_in%pointlist(2*i-1),
!     .           tri_in%pointlist(2*i)
!          enddo
!          write(lu,*) tri_in%numberofsegments,0
!          do i=1,tri_in%numberofsegments
!              write(lu,*) i-1,tri_in%segmentlist(2*i-1)-1,
!     .         tri_in%segmentlist(2*i)-1
!          enddo
!   ensure that areas for other materials are considered as holes
!          write(lu,*) tri_in%numberofholes
!          do i=1,tri_in%numberofholes
!              write(lu,*) i-1, tri_in%holelist(2*i-1),
!     .               tri_in%holelist(2*i)
!          enddo
!          call closef(lu)
  lu = nextlu(0)
  call openf (lu,MKGRID(1:lnblk(MKGRID))//'.poly','unknown')
  call dump_triangle_info(tri_in,lu)
  call closef(lu)
   endif
   call triangulate(tri_in, tri_out, triswitch(1:lnblk(triswitch)))
!     ENDENDEND new triangle interface
!      lu = nextlu(0)
!    call openf (lu,'/tmp/mkgrid.poly','unknown')
!      call openf (lu,MKGRID(1:lnblk(MKGRID))//'.poly','unknown')
!      write(lu,*) mnp,2,0,0
!      do i=1,mnp
!          xup = X(1,NP(i))+D(1,NP(i))
!          yup = X(2,NP(i))+D(2,NP(i))
!          write(lu,*) i-1,xup,yup
!      enddo
!      write(lu,*) mcn,0
!      do i=1,mcn
!          j = 1
!          do while (j.le.mnp .and. NP(j).ne.CON(1,i))
!              j = j + 1
!          enddo
!          k = 1
!          do while (k.le.mnp .and. NP(k).ne.CON(2,i))
!              k = k + 1
!          enddo
!          write(lu,*) i-1,j-1,k-1
!      enddo
!   ensure that areas for other materials are considered as holes
!      write(lu,*) NUMAT - 1
!      do i=1,NUMAT-1
!          write(lu,*) i-1, HOLEPOINT(1,i), HOLEPOINT(2,i)
!      enddo
!      call closef(lu)
!
!    command = 'cd /tmp; triangle -pQ mkgrid.poly'
!      command = 'triangle -pQ '
!     .      //MKGRID(1:lnblk(MKGRID))//'.poly'
!      command = 'cd /tmp; triangle -pQ '
!     .      //MKGRID(1:lnblk(MKGRID))//'.poly'
!      ierr = csystem(command)
!      if (ierr.ne.0) then
!          write(stderr,40) MKGRID(1:lnblk(MKGRID)), ierr
!   40     format(1x,'regrid: "triangle -pQ ',A,'.poly" exit ',I5)
!          ierr = EXEC
!          return
!      endif
!      luo = nextlu(lu+1)
!      call openf(luo,NEWELMS,'unknown')
!      call openf(lu,MKGRID(1:lnblk(MKGRID))//'.1.ele','old')
!    call openf(luo,'/tmp/newelms','unknown')
!    call openf(lu,'/tmp/mkgrid.1.ele','old')
!      read(lu,*,err=300,end=400,iostat=ios) nne,i,j
!      if (i.ne.3) then
!          write(stderr,50) MKGRID(1:lnblk(MKGRID))
!   50        format(1x,'regrid: no triangles in "/tmp/mkgrid.1.ele" ')
!   50     format(1x,'regrid: no triangles in "',A,'.1.ele" ')
!          ierr =  EXEC
!          return
!      endif
!      if (nne.ne.me) write(stderr,60) me,nne,MKGRID(1:lnblk(MKGRID)),
!     >   MKGRID(1:lnblk(MKGRID))
!   60   format(1x,'regrid: original # of elements = ',I5/8x,
!     >   'new # of elements = ',I5/8x,
!     >   'Poly file: ',A,'.poly',/8x,
!     >   'Ele  file: ',A,'.1.ele')
!     after this the program simply contiues, pretending as if nothing has been
!     done
   n = 0
   call NXTNDX(MAT,MATNO,NUMEL,n)
   ALLOCATE(newelements(tri_out%numberoftriangles),STAT=ERROR)
   ALLOCATE(newtriangles(3,tri_out%numberoftriangles),STAT=ERROR)
   numnewelements = tri_out%numberoftriangles
   do i=1,tri_out%numberoftriangles!nne
!          RIEN(1) = tri_out%trianglelist(3*i - 2)
!          RIEN(2) = tri_out%trianglelist(3*i - 1)
!          RIEN(3) = tri_out%trianglelist(3*i)
!          rien(1) = NP(rien(1))
!          rien(2) = NP(rien(2))
!          rien(3) = NP(rien(3))
!    read(lu,*,err=300,end=400,iostat=ios) j,RIEN(1),RIEN(2), &
!              RIEN(3)
!        rien(1) = NP(rien(1)+1)
!        rien(2) = NP(rien(2)+1)
!        rien(3) = NP(rien(3)+1)
  newelements(i) = n
  newtriangles(1, i) = NP(tri_out%trianglelist(3*i - 2))!rien(1)
  newtriangles(2, i) = NP(tri_out%trianglelist(3*i - 1))!rien(2)
  newtriangles(3, i) = NP(tri_out%trianglelist(3*i))!rien(3)
!          write(luo,65) n,(rien(j),j=1,3)
!   65     format(4I8)
  if (i.lt.numnewelements) call NXTNDX(MAT,MATNO,NUMEL,n)
   enddo
!      close(lu)
!      close(lu,status='delete')
!      close(luo)
!    command = 'rm -f /tmp/mkgrid.*'
!      command = 'rm -f '//MKGRID(1:lnblk(MKGRID))//'.*'
   call dealloc(tri_in)
   call dealloc(tri_out)

!      ierr = csystem(command)
!      if (ierr.ne.0) then
!          write(stderr,70) ierr, command
!   70     format(1x,'regrid: "'A'" exit ',I5)
!          ierr = EXEC
!          return
!      endif
!
   return

#endif
!
!300    write(stderr,301) ios, MKGRID(1:lnblk(MKGRID))
! 301    format(1x,'regrid: read error ',I4,' on "',A,'.1.ele"')
!     ierr = EXEC
!    return
!400    write(stderr,401) MKGRID(1:lnblk(MKGRID))
! 401    format(1x,'regrid: premature EOF error on "',A,'.1.ele"')
!     ierr = EXEC
!    return
 end subroutine
!-----------------------------------------------------------------------------
 subroutine NXTNDX(MAT,MATNO,NUMEL,N)
 implicit none
 integer MATNO,MAT(*),N,NUMEL
 N = N + 1
 do while (MAT(N).ne.MATNO)
N = N + 1
if (N.gt.NUMEL) stop 'fail'
 enddo
 return
 end
!-----------------------------------------------------------------------------
   subroutine mapelq (IEN,X,D,STN,STR,DXE,ierr)
!
! maps element based quantities (stress, strain, slippery node displacements)
! to newly tesselated elements
!
   USE REGRIDDATA
! SNGL    implicit real (a-h,o-z)
   implicit double precision (a-h,o-z)
   parameter(third=1.d0/3)



!
!-pass
   integer ierr
dimension X(NSD,*),D(NDOF,*),IEN(NEN,*),STN(NSTR,*),STR(NSTR,*), &
   DXE(NDOF,NEN,*)
!-locl
!#include "lenchar.i"
   double precision xl(2,4),rc(2),c(2),sh(3,4),xs(4),det,dist
!      integer lu,nextlu,lustn,lustr,ludxe,rien(3)
!      external nextlu
!      character(len=)tmpnm*255
   integer ERROR, cnt
!-init
   ierr = OK
!      NEWSTN = 'newstress'!tmpnm('newstress')
!      NEWSTR = 'newstrain'!tmpnm('newstrain')
!      NEWDXE = 'newdxe'!tmpnm('newdxe')
!      NEWNUM = 'newnum'!tmpnm('newnum')
!
!lu = nextlu(0)
!call openf(lu,NEWELMS,'old')
!lustn = nextlu(lu+1)
!call openf(lustn,NEWSTN,'unknown')
!lustr = nextlu(lustn+1)
!call openf(lustr,NEWSTR,'unknown')
!if (NUMSLP.gt.0) then
!    ludxe = nextlu(lustr+1)
!    call openf(ludxe,NEWDXE,'unknown')
!    lunum = nextlu(ludxe+1)
!    call openf(lunum,NEWNUM,'unknown')
!endif
   ALLOCATE(newstress(NSTR,numnewelements),STAT=ERROR)
   ALLOCATE(newstrain(NSTR,numnewelements),STAT=ERROR)
   if (NUMSLP.gt.0) then
  ALLOCATE(newslpdxe(NDOF,NEN,numnewelements),STAT=ERROR)
  ALLOCATE(newminnum(numnewelements),STAT=ERROR)
   endif
   cnt = 0
   do cnt = 1,numnewelements
!read(lu,'(4I8)',end=200,err=300,iostat=ios) n,rien(1),
!     >     rien(2),rien(3)
!      write(*,*) "old ", n, rien(1), rien(2), rien(3)
!      write(*,*) "new ", newelements(cnt), newtriangles(1, cnt),
!     .      newtriangles(2, cnt), newtriangles(3, cnt)

  do j=1,3
      k = newtriangles(j,cnt)!rien(j)
      do i=1,NSD
          xl(i,j) = X(i,k) + D(i,k)
      enddo
  enddo
  xl(1,4) = xl(1,3)
  xl(2,4) = xl(2,3)
!      call SHAP20 (0d0,-1d0/3d0,xl,det,sh,xs,.false.,ierr)
!      if (ierr.ne.0) return
  rc(1) = third * (xl(1,1) + xl(1,2) + xl(1,3))
!       rc(1) =  sh(3,1)*xl(1,1)+sh(3,2)*xl(1,2)+sh(3,3)*xl(1,3)
!     >   +sh(3,4)*xl(1,4)
  rc(2) = third * (xl(2,1) + xl(2,2) + xl(2,3))
!        rc(2) =  sh(3,1)*xl(2,1)+sh(3,2)*xl(2,2)+sh(3,3)*xl(2,3)
!     >   +sh(3,4)*xl(2,4)
  dist = FLMAX()
  do m=1,NUMEL
      do j=1,3
          k = IEN(j,m)
          do i=1,NSD
              xl(i,j) = X(i,k) + D(i,k)
          enddo
      enddo
      xl(1,4) = xl(1,3)
      xl(2,4) = xl(2,3)
!         call SHAP20 (0d0,-1d0/3d0,xl,det,sh,xs,.false.,ierr)
!          if (ierr.ne.0) return
      c(1) = third * (xl(1,1) + xl(1,2) + xl(1,3))
      c(2) = third * (xl(2,1) + xl(2,2) + xl(2,3))

!          c(1) = sh(3,1)*xl(1,1)+sh(3,2)*xl(1,2)+sh(3,3)*xl(1,3)
!     >          +sh(3,4)*xl(1,4)
!          c(2) = sh(3,1)*xl(2,1)+sh(3,2)*xl(2,2)+sh(3,3)*xl(2,3)
!     >          +sh(3,4)*xl(2,4)
      dd = (rc(1)-c(1))**2 + (rc(2)-c(2))**2
      if (dd.lt.dist) then
          dist = dd
          minno = m
      endif
  enddo
!    write(lustn,10) n,(STN(i,minno),i=1,NSTR)
!   10    format(I8,3E26.18)
  do i=1,NSTR
      newstress(i, cnt) = STN(i,minno)
  enddo
!    write(lustr,10) n,(STR(i,minno),i=1,NSTR)
  do i=1,NSTR
      newstrain(i, cnt) = STR(i,minno)
  enddo
  if (NUMSLP.gt.0) then
!        write(ludxe,20) n,((DXE(i,j,minno),j=1,NEN),i=1,NDOF)
!   20        format(I8,8E26.18)
    do j=1,NEN
        do i=1,NDOF
            newslpdxe(i,j,cnt) = DXE(i,j,minno)
        enddo
    enddo
    newminnum(cnt) = minno
!        write(lunum,'(2I8)') n,minno
  endif
   enddo!goto 100
!
!200    close(lu)
!    close(lustn)
!    close(lustr)
!    if (NUMSLP.gt.0) then
!        close(ludxe)
!        close(lunum)
!    endif
   return
!
!300    write(stderr,301) ios, NEWELMS(1:lnblk(NEWELMS))
! 301    format(1x,'regrid: read error ',I4,' on "',A,'"')
!     ierr = EXEC
!    return
 end
!------------------------------------------------------------------------------
 subroutine renew (IEN,STN,STR,DXE,NSLIP,NSELSD,ierr)
!
! stores regridding quantities back into their new array locations
!
   USE REGRIDDATA
 implicit double precision (a-h,o-z)



!-pass
 integer IEN(NEN,*),NSLIP(5,*),NSELSD(2,*),ierr
 double precision STN(NSTR,*),STR(NSTR,*),DXE(NDOF,NEN,*)
!-locl
!    integer luelm,nextlu,lustn,lustr,ios,n,rien(3),m,
   integer lien(2)
!    double precision sxx,sxy,syy,dxel(2,4)
   integer ERROR, cnt

#ifdef EXPERIMENTAL_ReducedSlipAlloc
write(*,*) "Regridding not yet supported with reduced memory alloc"
#else

!
!    ierr = OK
!    luelm = nextlu(0)
!    call openf(luelm,NEWELMS,'old')
!    lustn = nextlu(luelm+1)
!    call openf(lustn,NEWSTN,'old')
!    lustr = nextlu(lustn+1)
!    call openf(lustr,NEWSTR,'old')
!    if (NUMSLP.gt.0) then
!        ludxe = nextlu(lustr+1)
!        call openf(ludxe,NEWDXE,'old')
!        lunum = nextlu(ludxe+1)
!        call openf(lunum,NEWNUM,'old')
!    endif
!
   cnt = 0
   do cnt=1,numnewelements!cnt = cnt + 1
   n = newelements(cnt)
!read(luelm,'(4I8)',end=200,err=300,iostat=ios) n,rien(1),
!     >     rien(2),rien(3)
!    read(lustn,5,end=400,err=500,iostat=ios)
!     >     m,sxx,syy,sxy
!    5    format(I8,3E26.18)
!      write(*,*) "old ",  m,sxx,syy,sxy
!      write(*,*) "new ", newelements(cnt), newstress(1, cnt),
!     .  newstress(2,cnt),newstress(3, cnt)
!    if (m.ne.n) then
!        write(stderr,10)
!   10        format(1x,'regrid synchronization error')
!        ierr = EXEC
!        return
!    endif
 STN(1,n) = newstress(1,cnt)!sxx
 STN(2,n) = newstress(2,cnt)!syy
 STN(3,n) = newstress(3,cnt)!sxy
!    read(lustr,5,end=600,err=700,iostat=ios)
!     >     m,sxx,syy,sxy
!    if (m.ne.n) then
!        write(stderr,10)
!        call exitp(1)
!    endif
!      write(*,*) "old strain ",  m,sxx,syy,sxy
!      write(*,*) "new strain", newelements(cnt), newstrain(1, cnt),
!     .  newstrain(2,cnt),newstrain(3, cnt)
 STR(1,n) = newstrain(1,cnt)!sxx
 STR(2,n) = newstrain(2,cnt)!syy
 STR(3,n) = newstrain(3,cnt)!sxy
 if (NUMSLP.gt.0) then
!        read(ludxe,'(I8,8E26.18)',end=800,err=900,iostat=ios)
!     >         m,((dxel(i,j),i=1,NDOF),j=1,NEN)
!        if (m.ne.n) then
!            write(stderr,10)
!            call exitp(1)
!        endif
do j=1,NEN
    DXE(1,j,n) = newslpdxe(1,j,cnt)!dxel(1,j)
    DXE(2,j,n) = newslpdxe(2,j,cnt)!dxel(2,j)
enddo
!        read(lunum,'(2I8)',end=1000,err=1100,iostat=ios)
!     >         m,mold
!        if (m.ne.n) then
!            write(stderr,10)
!            call exitp(1)
!        endif
    mold = newminnum(cnt)
do i=1,NUMSLP
    if (NSLIP(1,i).eq.mold) NSLIP(1,i) = n
enddo
if (NSLSKEW.gt.0) then
    do i=1,NSLSKEW
    if (NSELSD(1,i).eq.mold) then
        NSELSD(1,i) = n
    call SIDENP(NSELSD(2,i),.false.,lien)
    m = kside (rien,IEN(lien(1),mold), &
           IEN(lien(2),mold),NEN,.false.)
    if (m.eq.0) then
        write(stderr,20) n,i
20                       format(1x,'regrid: cannot establish side ', &
                     'number of new element ',I8,', i.e. fault', &
                     ' parallel element entry ',I5)
        ierr = EXEC
        return
    endif
    NSELSD(2,i) = m
    endif
    enddo
endif
 endif
 IEN(1,n) = newtriangles(1,cnt)!rien(1)
 IEN(2,n) = newtriangles(2,cnt)!rien(2)
 IEN(3,n) = newtriangles(3,cnt)!rien(3)
 IEN(4,n) = newtriangles(3,cnt)!rien(3)
   enddo!goto 100

   if (iecho.eq.1) then
  write(stdout,*) 'writing new post-regrid data to disk'
  call writeregriddata(NSTR,NEN,NDOF)
   endif
   DEALLOCATE(newelements, STAT=ERROR)
   DEALLOCATE(newstress, STAT=ERROR)
   DEALLOCATE(newstrain, STAT=ERROR)
   if (NUMSLP.gt.0) then
  DEALLOCATE(newslpdxe, STAT=ERROR)
  DEALLOCATE(newminnum, STAT=ERROR)
   endif
!
!200    close(luelm,status='delete')
!    close(lustn,status='delete')
!    close(lustr,status='delete')
!    if (NUMSLP.gt.0) then
!        close(ludxe,status='delete')
!        close(lunum,status='delete')
!    endif
 return

#endif

!
!300    write(stderr,301) ios, NEWELMS(1:lnblk(NEWELMS))
! 301    format(1x,'regrid: read error ',I4,' on "',A,'"')
!     ierr = EXEC
!    return
!400    write(stderr,401) NEWSTN(1:lnblk(NEWSTN))
! 401    format(1x,'regrid: premature EOF on "',A,'"')
!     ierr = EXEC
!    return
!500    write(stderr,501) ios, NEWSTN(1:lnblk(NEWSTN))
! 501    format(1x,'regrid: read error ',I4,' on "',A,'"')
!     ierr = EXEC
!    return
!600    write(stderr,601)  NEWSTR(1:lnblk(NEWSTR))
! 601    format(1x,'regrid: premature EOF on "',A,'"')
!     ierr = EXEC
!    return
!700    write(stderr,701) ios, NEWSTR(1:lnblk(NEWSTR))
! 701    format(1x,'regrid: read error ',I4,' on "',A,'"')
!     ierr = EXEC
!    return
!800    write(stderr,801) NEWDXE(1:lnblk(NEWDXE))
! 801    format(1x,'regrid: premature EOF on "',A,'"')
!     ierr = EXEC
!    return
!900    write(stderr,901) ios, NEWDXE(1:lnblk(NEWDXE))
! 901    format(1x,'regrid: read error ',I4,' on "',A,'"')
!     ierr = EXEC
!    return
!1000    write(stderr,1001) NEWNUM(1:lnblk(NEWNUM))
! 1001    format(1x,'regrid: premature EOF on "',A,'"')
!     ierr = EXEC
!    return
!1100    write(stderr,1101) ios, NEWNUM(1:lnblk(NEWNUM))
! 1101    format(1x,'regrid: read error ',I4,' on "',A,'"')
!     ierr = EXEC
!    return
 end
!--------------------------------------------------------------------------
subroutine writeelm(IEN,MAT)
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)



!-pass
dimension IEN(NEN,*),MAT(*)
!-locl
 character(len=12)fmt
integer i,n,lu,nextlu
 external nextlu
!-init
lu = nextlu(50)
call openf(lu,'tecin.dat.elm_reg','unknown')

 fmt = '(I12,I5,9I8)'


!
do n=1,NUMEL
   write(lu,fmt=fmt,err=100,iostat=ios) n,MAT(n), &
      (IEN(i,n),i=1,NEN)
enddo
write(lu,10)
  10    format('end element definitions')
 close(lu)
!
 return
100    write(stderr,101) ios
 101    format(1x,'regrid: write error ',I5,' on "tecin.dat.elm_reg"')
 ierr = EXEC
 return
end
!--------------------------------------------------------------------------
subroutine writeslp(NSLIP,DIFORC,NSELSD,IWINX,WINX)
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)



!-pass
dimension NSLIP(5,*),DIFORC(NDOF,*),IWINX(NDOF,*),WINX(NDOF,*), &
   NSELSD(2,*)
!-locl
 character(len=20)fmt
integer i,n,lu,nextlu
 external nextlu
!-init
lu = nextlu(50)
call openf(lu,'tecin.dat.sns_reg','unknown')

fmt = '(I12,I12,3I5,3G14.8)'


!
do n=1,NUMSLP
   write(lu,fmt=fmt,err=100,iostat=ios) (NSLIP(i,n),i=1,5), &
      (DIFORC(i,n),i=1,NDOF)
enddo
write(lu,10)
  10    format('end slippery nodes')
!
 if (NSLSKEW.ne.0) then

fmt = '(I12,I5)'

do n=1,NSLSKEW
    write(lu,fmt=fmt,err=100,iostat=ios) (NSELSD(i,n),i=1,2)
enddo

write(lu,20)
20        format('end fault parallel elements')
 endif
!
 if (NWINKX.ne.0) then
fmt = '(2I12,3I5)  '

do n=1,NUMNP
if (IWINX(1,n).ne.0 .or. IWINX(2,n).ne.0) &
       write(lu,fmt=fmt,err=100,iostat=ios) (IWINX(i,n),i=1, &
       NDOF)
enddo
 endif
 write(lu,30)
30    format('end differential Winkler codes')
 if (NWINKX.ne.0) then
fmt = '(2I12,3G14.8)'

do n=1,NUMNP
    if (IWINX(1,n).ne.0 .or. IWINX(2,n).ne.0) then
        write(lu,fmt=fmt,err=100,iostat=ios) (WINX(i,n),i=1, NDOF)
    endif
enddo
 endif

 write(lu,40)
40    format('end differential Winkler magnitudes')
!
 close(lu)
!
 return
100    write(stderr,101) ios
 101    format(1x,'regrid: write error ',I5,' on "tecin.dat.sns_reg"')
 ierr = EXEC
 return
end
!------------------------------------------------------------------------------
 integer function kside (IEN,N,M,NEN,QUAD)
!
!    sets element side number defined by global nodes N and M
!
!    returns 0 (fail) or int (element side)
!-pass
 logical QUAD
 integer IEN(*),N,M,NEN
!-locl
 integer i
!
 i = 1
100    if (QUAD) then
j = MOD(i,NEN)+1
 else
j = MOD(i,NEN-1)+1
 endif
 if ((IEN(i).eq.N .and. IEN(j).eq.M).or. &
   (IEN(i).eq.M .and. IEN(j).eq.N)) goto 200
 if (QUAD) then
i = MOD(i+1,NEN)
 else
i = MOD(i+1,NEN-1)
 endif
 if (i.ne.0) goto 100
!
200    kside = i
 return
 end
