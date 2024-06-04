subroutine setmon (cmd,ncmd,IEN,NUMEL,NUMNP,ierr)

  USE ALGEBRA
  use documentation,   only: plt3d_manual
  use modeltopology,   only: nen
  use iomodule,        only: icoord, ifelm, stderr
  use monitorVars,     only: elmon, npmon, maxmon, monel, monnp, nnel, nnp, &
                             alloc_node_monitor_data,alloc_elem_monitor_data
  use plotControl,     only: npbased
  use errors,          only: ok, exec, syntax
  
  implicit none
  !-pass
  integer :: ncmd
  character cmd*(*)
  integer NUMEL,NUMNP,ierr
  integer ien
  dimension cmd(*),IEN(NEN,*)
  
  !-locl
  logical, external :: numeric
  double precision, external :: chreal
  integer :: i, j, k, l
  integer :: nprev

  !-init
  ierr = OK
1 format(1x,'Usage: mo[nitor] n|e numbers')

  if (cmd(2)(1:1).eq.'n') then

      if (ncmd.lt.3) then
          call plt3d_manual('mo')
          ierr = SYNTAX
      else
          call alloc_node_monitor_data
          call thenum (cmd,3,ncmd,monnp,nnp,1,NUMNP,NUMNP,ierr)
          npmon = (ierr.eq.OK)
      endif
      npmon = (ierr.eq.OK)

  else if (cmd(2)(1:1).eq.'e') then

      if (ncmd.lt.3) then
          call plt3d_manual('mo')
          ierr = SYNTAX
      else
          call alloc_elem_monitor_data
          call thenum (cmd,3,ncmd,monel,nnel,1,NUMEL,NUMEL,ierr)
      endif
      elmon = (ierr.eq.OK)
  endif

  return

end subroutine
!-------------------------------------------------------------------------------
subroutine thenum (cmd,iicmd,ncmd,iarray,n,m,max,num,ierr)
 
  ! parses cmd(iicmd) to cmd(ncmd) and puts resulting monitor numbers
  ! into IARRAY(M,n), IARRAY(M,n+1), etc... . Range checking on monitor
  ! numbers is performed: NUM is NUMEL or NUMNP.

  use iomodule,        only: stderr
  use errors,          only: exec, syntax

  implicit none
  !-pass
  integer :: iicmd,ncmd,n,m,max,num,ierr
  character(len=*) :: cmd(*)
  integer :: iarray(M,max)

  !-locl
  logical :: numeric
  logical, parameter :: thenum_debug=.false.
  integer,external :: lnblk
  integer :: i,j,k,l,icmd,ibase,istep,iend
  double precision, external :: chreal

  if (thenum_debug) then
      write(0,*) 'entered thenum.  iicmd,ncmd,n,m,max,num=',iicmd,ncmd,n,m,max,num
  endif

  icmd = iicmd
  do while (icmd.le.ncmd) 

    l = lnblk(cmd(icmd))                ! find string index L of end of string
    i = index(cmd(icmd)(1:l),',')       ! find string index of a comma
    if (thenum_debug) then
        write(0,"(80a)") cmd(icmd)(1:l)
        write(0,*) 'iindex=',i,' l=',l
    endif
    if (i.eq.1 .or. i.eq.l) then
        write(stderr,1) 
1       format(1x,'monitor incrementing error')
        ierr = EXEC
        return
    endif

    if (i.le.0) then
        ! syntax: the string does not contain a comma

        if (numeric(cmd(icmd)(1:l))) then

            if (thenum_debug) then
                write(0,*) 'the string contains a number'
            endif

            if (n+1.le.max) then

                if (thenum_debug) then
                    write(0,*) 'there is space for storing the number in IARRAY'
                endif

                k = INT(chreal(cmd(icmd)(1:l)))   ! convert the string to an integer
                if (thenum_debug) then
                    write(0,*) 'the number k=',k
                endif

                if (k.ge.1 .and. k.le.num) then
                    if (thenum_debug) then
                        write(0,*) 'the number fits in the expected range 1-',num
                    endif
                    n = n + 1

                    iarray(M,n) = k
                    if (thenum_debug) then
                        write(0,*) 'set iarray(',m,',',n,') to ',k
                    endif
                else

                    write(stderr,2) num
2                   format(1x,'number out of range [1,',I7,']')
                    ierr = EXEC
                    return

                endif

            else

                write(stderr,3) max
3               format(1x,'overflow error: max=',I5)
                ierr = EXEC
                return
            endif
        else
            write(stderr,4)
4           format(1x,'syntax error')
            ierr = SYNTAX
            return
        endif

    else

        if (thenum_debug) then
            write(0,*) 'the string contains one or more comma''s'
        endif

        ! i is the index of the first comma from the left
        if (numeric(cmd(icmd)(1:i-1))) then

            if (thenum_debug) then
                write(0,*) 'the string before the comma contains a number'
            endif

            ibase = INT(chreal(cmd(icmd)(1:i-1)))       ! convert the first substring to an integer
            if (thenum_debug) then
                write(0,*) 'the number before the first comma=',ibase
            endif

        else

            write(stderr,4)
            ierr = SYNTAX 
            return

        endif

        j = i + index(cmd(icmd)(i+1:l),',')     ! find the index of a possible second comma
        if (thenum_debug) then
            write(0,*) 'the index of the second comma j=',j
        endif
        if (j.eq.i+1) then
            ! the string reads "number,,"
            write(stderr,4)
            ierr = SYNTAX
            return
        endif

        if (j.gt.i+1) then

            ! syntax: ibase,istep,iend
            if (thenum_debug) then
                write(0,*) 'the string contains two or more comma''s'
            endif

            if (numeric(cmd(icmd)(i+1:j-1))) then

                if (thenum_debug) then
                    write(0,*) 'the string before the comma contains a number'
                endif

                istep = INT(chreal(cmd(icmd)(i+1:j-1))) ! convert the second string to an integer
                if (thenum_debug) then
                    write(0,*) 'the number between the 1st and 2nd comma=',istep
                endif

            else

                write(stderr,4)
                ierr = SYNTAX 
                return

            endif

            if (numeric(cmd(icmd)(j+1:l))) then

                if (thenum_debug) then
                    write(0,*) 'the string after the 2nd comma contains a number'
                endif

                iend = INT(chreal(cmd(icmd)(j+1:l))) ! convert the 3rd string to an integer

            else

                write(stderr,4)
                ierr = SYNTAX 
                return

            endif

        else

            ! syntax: ibase,iend
            if (thenum_debug) then
                write(0,*) 'the string contains a single comma'
            endif

            if (numeric(cmd(icmd)(i+1:l))) then

                if (thenum_debug) then
                    write(0,*) 'the string after the comma contains a number'
                endif

                iend  = INT(chreal(cmd(icmd)(i+1:l)))
                if (thenum_debug) then
                    write(0,*) 'the number after the comma=',iend
                endif

                if (iend.gt.ibase) then
                    istep=1
                    if (thenum_debug) then
                        write(0,*) 'iend>ibase: set istep=',istep
                    endif
                else if (iend.lt.ibase) then
                    istep=-1
                    if (thenum_debug) then
                        write(0,*) 'iend<ibase: set istep=',istep
                    endif
                else
                    write(stderr,1)
                    ierr = SYNTAX
                    return
                endif

            else

                write(stderr,4)
                ierr = SYNTAX 
                return
            endif

        endif

        if (MOD(IABS(iend-ibase),IABS(istep)).ne.0) then
            write(stderr,5)
5           format(1x,'WARNING: increment rounded')
        endif

        ! generate the numbers
        do i=ibase,iend,istep
            if (n+1.gt.max) then
                write(stderr,3) max
                ierr = EXEC
                return
            endif
            if (i.lt.1 .or. i.gt.num) then
                write(stderr,2) num
                ierr = EXEC
                return
            endif
            n = n + 1
            iarray(M,n) = i

            if (thenum_debug) then
                write(0,*) 'set iarray(',M,',',n,') to ',i
            endif
        enddo

    endif

    icmd = icmd + 1

  enddo

  return

end subroutine thenum
!--------------------------------------------------------------------------------
! RG 2023-7-28 subroutine thenum (cmd,iicmd,ncmd,iarray,n,m,max,num,ierr)
! RG 2023-7-28
! RG 2023-7-28!    parses cmd(iicmd) to cmd(ncmd) and puts resulting monitor numbers
! RG 2023-7-28!    into IARRAY(M,n), IARRAY(M,n+1), etc... . Range checking on monitor
! RG 2023-7-28!    numbers is performed: NUM is NUMEL or NUMNP.
! RG 2023-7-28
! RG 2023-7-28use iomodule,        only: stderr
! RG 2023-7-28use errors,          only: exec, syntax
! RG 2023-7-28
! RG 2023-7-28implicit none
! RG 2023-7-28
! RG 2023-7-28!-pass
! RG 2023-7-28integer :: iicmd, ncmd, n, m, max, num, ierr
! RG 2023-7-28integer :: iarray(m,max)
! RG 2023-7-28 character(len=*) cmd
! RG 2023-7-28 dimension cmd(*)
! RG 2023-7-28
! RG 2023-7-28
! RG 2023-7-28
! RG 2023-7-28
! RG 2023-7-28!-locl
! RG 2023-7-28logical numeric
! RG 2023-7-28integer, external :: lnblk
! RG 2023-7-28real, external :: chreal
! RG 2023-7-28integer :: icmd
! RG 2023-7-28
! RG 2023-7-28integer :: i, l, k, j
! RG 2023-7-28integer :: ibase, istep, iend
! RG 2023-7-28
! RG 2023-7-28icmd = iicmd
! RG 2023-7-28
! RG 2023-7-28100    if (icmd.gt.ncmd) then
! RG 2023-7-28    return
! RG 2023-7-28endif
! RG 2023-7-28
! RG 2023-7-28    l = lnblk(cmd(icmd))
! RG 2023-7-28    i = index(cmd(icmd)(1:l),',')
! RG 2023-7-28    !    write(*,"(80a)") cmd(icmd)(1:l)
! RG 2023-7-28    !    write(*,*) 'iindex=',i,' l=',l
! RG 2023-7-28
! RG 2023-7-28    if (i.eq.1 .or. i.eq.l) then
! RG 2023-7-28        write(stderr,1)
! RG 2023-7-28 1        format(1x,'monitor incrementing error')
! RG 2023-7-28        ierr = EXEC
! RG 2023-7-28        return
! RG 2023-7-28    endif
! RG 2023-7-28
! RG 2023-7-28    if (i.le.0) then
! RG 2023-7-28        ! syntax contains no comma.
! RG 2023-7-28!        write(*,*) "no comma", cmd(icmd)(1:l), numeric(cmd(icmd)(1:l))
! RG 2023-7-28        if (numeric(cmd(icmd)(1:l))) then
! RG 2023-7-28            ! argument is a single number
! RG 2023-7-28!            write(*,*) "n+1, max", n+1, max
! RG 2023-7-28            if (n+1.le.max) then
! RG 2023-7-28!                write(*,*) "fuckadoo1", cmd(icmd)(1:l)
! RG 2023-7-28                read(cmd(icmd)(1:l),*) k
! RG 2023-7-28
! RG 2023-7-28!                write(*,*) "fuckadoo2", chreal(cmd(icmd)(1:l))
! RG 2023-7-28!                write(*,*) "fuckadoo3", INT(chreal(cmd(icmd)(1:l)))
! RG 2023-7-28!                k = INT(chreal(cmd(icmd)(1:l)))
! RG 2023-7-28!                write(*,*) "monitoring node: ", k
! RG 2023-7-28                if (k.ge.1 .and. k.le.num) then
! RG 2023-7-28                    n = n + 1
! RG 2023-7-28                    iarray(M,n) = k
! RG 2023-7-28    !    write(*,*) 'iarray(',n,')=',k
! RG 2023-7-28                else
! RG 2023-7-28!                    write(stderr,2) num
! RG 2023-7-28 2                   format(1x,'Setmon; number out of range [1,',I12, ']')
! RG 2023-7-28                    write(stderr,*) "Setmon wants to monitor node", k
! RG 2023-7-28                    if (k.lt.1) then
! RG 2023-7-28                        write(stderr,*) "node number must be positive"
! RG 2023-7-28                    else
! RG 2023-7-28                        write(stderr,*) "but the mesh has only ", num, "nodes"
! RG 2023-7-28                    endif
! RG 2023-7-28
! RG 2023-7-28                    ierr = EXEC
! RG 2023-7-28                    return
! RG 2023-7-28                endif
! RG 2023-7-28            else
! RG 2023-7-28                write(stderr,3) max
! RG 2023-7-28 3            format(1x,'Setmon; overflow error: max=',I12)
! RG 2023-7-28                ierr = EXEC
! RG 2023-7-28                return
! RG 2023-7-28            endif
! RG 2023-7-28        else
! RG 2023-7-28            write(stderr,4)
! RG 2023-7-28 4        format(1x,'Setmon; syntax error')
! RG 2023-7-28            ierr = SYNTAX
! RG 2023-7-28            return
! RG 2023-7-28        endif
! RG 2023-7-28    else
! RG 2023-7-28        if (numeric(cmd(icmd)(1:i-1))) then
! RG 2023-7-28            ibase = INT(chreal(cmd(icmd)(1:i-1)))
! RG 2023-7-28        else
! RG 2023-7-28            write(stderr,4)
! RG 2023-7-28            ierr = SYNTAX
! RG 2023-7-28            return
! RG 2023-7-28        endif
! RG 2023-7-28        j = i + index(cmd(icmd)(i+1:l),',')
! RG 2023-7-28
! RG 2023-7-28        !    write(*,*) 'jindex=',j
! RG 2023-7-28        if (j.eq.i+1) then
! RG 2023-7-28            write(stderr,4)
! RG 2023-7-28            ierr = SYNTAX
! RG 2023-7-28            return
! RG 2023-7-28        endif
! RG 2023-7-28
! RG 2023-7-28        if (j.gt.i+1) then
! RG 2023-7-28        ! syntax: ibase,istep,iend
! RG 2023-7-28            if (numeric(cmd(icmd)(i+1:j-1))) then
! RG 2023-7-28                istep = INT(chreal(cmd(icmd)(i+1:j-1)))
! RG 2023-7-28            else
! RG 2023-7-28                write(stderr,4)
! RG 2023-7-28                ierr = SYNTAX
! RG 2023-7-28                return
! RG 2023-7-28            endif
! RG 2023-7-28
! RG 2023-7-28            if (numeric(cmd(icmd)(j+1:l))) then
! RG 2023-7-28                iend = INT(chreal(cmd(icmd)(j+1:l)))
! RG 2023-7-28            else
! RG 2023-7-28                write(stderr,4)
! RG 2023-7-28                ierr = SYNTAX
! RG 2023-7-28                return
! RG 2023-7-28            endif
! RG 2023-7-28        else
! RG 2023-7-28            ! syntax: ibase,iend
! RG 2023-7-28            if (numeric(cmd(icmd)(i+1:l))) then
! RG 2023-7-28                iend  = INT(chreal(cmd(icmd)(i+1:l)))
! RG 2023-7-28                if (iend.gt.ibase) then
! RG 2023-7-28                    istep=1
! RG 2023-7-28                else if (iend.lt.ibase) then
! RG 2023-7-28                    istep=-1
! RG 2023-7-28                else
! RG 2023-7-28                    write(stderr,1)
! RG 2023-7-28                    ierr = SYNTAX
! RG 2023-7-28                    return
! RG 2023-7-28                endif
! RG 2023-7-28            else
! RG 2023-7-28                write(stderr,4)
! RG 2023-7-28                ierr = SYNTAX
! RG 2023-7-28                return
! RG 2023-7-28            endif
! RG 2023-7-28        endif
! RG 2023-7-28
! RG 2023-7-28        if (MOD(IABS(iend-ibase),IABS(istep)).ne.0) then
! RG 2023-7-28            write(stderr,5)
! RG 2023-7-28 5        format(1x,'WARNING: increment rounded')
! RG 2023-7-28        endif
! RG 2023-7-28
! RG 2023-7-28        do i=ibase,iend,istep
! RG 2023-7-28            if (n+1.gt.max) then
! RG 2023-7-28                write(stderr,3) max
! RG 2023-7-28                ierr = EXEC
! RG 2023-7-28                return
! RG 2023-7-28            endif
! RG 2023-7-28            if (i.lt.1 .or. i.gt.num) then
! RG 2023-7-28                write(stderr,2) num
! RG 2023-7-28                ierr = EXEC
! RG 2023-7-28                return
! RG 2023-7-28            endif
! RG 2023-7-28            n = n + 1
! RG 2023-7-28            iarray(M,n) = i
! RG 2023-7-28        enddo
! RG 2023-7-28
! RG 2023-7-28    endif
! RG 2023-7-28
! RG 2023-7-28    icmd = icmd + 1
! RG 2023-7-28    goto 100
! RG 2023-7-28
! RG 2023-7-28!endif
! RG 2023-7-28
! RG 2023-7-28end subroutine
