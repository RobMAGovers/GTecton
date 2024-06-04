subroutine itst1 (A,N,istep)

use iomodule,        only: stderr

implicit none

integer N,A, istep
dimension A(N)

integer    :: ii,jj


write(stderr,'(''A('',I3,'')'')') N
jj = MAX0(istep,1)

do ii=1,N
 if (MOD(ii,jj).eq.0) then
write(stderr,'(I4,I12)') ii,A(ii)
 endif
enddo

return
end
!-------------------------------------------------------------------------------
subroutine itst2 (A,N,M,istep)

use iomodule,        only: stderr

implicit none

integer n,m,A, istep
dimension A(N,M)

integer :: ii,jj,nn

write(stderr,'(''A('',I3,'','',I3,'')'')') N,M
jj = MAX0(istep,1)

do ii=1,M
 if (MOD(ii,jj).eq.0) then
write(stderr,'(I4,8I12)') ii,(A(nn,ii),nn=1,N)
 endif
enddo

return
end
!-------------------------------------------------------------------------------
subroutine rdalsy ()

use symbols

implicit none



call rdsymb (IBcFix,SYFIX,MFIX,MAXSYM)
call rdsymb (IBcSld,SYSLID,MSLID,MAXSYM)
call rdsymb (IBcWnk,SYWINK,MWINK,MAXSYM)
call rdsymb (IBcVel,SYVELO,MVELO,MAXSYM)
call rdsymb (IBcStn,SYSTRS,MSTRS,MAXSYM)
call rdsymb (IBcFor,SYFORC,MFORC,MAXSYM)
call rdsymb (IBcSh,SYSHEAR,MSHEAR,MAXSYM)
call rdsymb (IBcOpn,SYOPEN,MOPEN,MAXSYM)
call rdsymb (IBcShr,SYSHR,MSHR,MAXSYM)
return
end
!-------------------------------------------------------------------------------
subroutine rdsymb (numsy,SYMB,NSYM,MAXSYM)

USE ALGEBRA
use debugmodule,     only: iecho
use iomodule,        only: stderr

implicit none
!-pass
integer numsy,NSYM,MAXSYM
double precision :: symB(2,MAXSYM)

!-locl
character(len=2) charnu
character(len=255) fname
character(len=78) PU
integer luloc
integer, external :: lnblk,fnblnk, nextlu
integer :: i, j, l, n
integer :: ios, itype

NSYM = 0
call get_environment_variable('PU',PU)
l = lnblk(PU)

if (l.le.0) then
 write(stderr,1)
1        format(1x,'rdsymb: P-util directory $PU undefined')
 call exitp(1)
endif

write(charnu,2,err=1000,iostat=ios) numsy
2    format(I2)
fname = PU(1:l) // '/symbols/' // &
  charnu(fnblnk(charnu):lnblk(charnu))
luloc = nextlu(0)
call openf(luloc,fname,'old')
call clear(SYMB,2*MAXSYM,"symb")
l = 0
100    read(luloc,*,err=2000,end=200,iostat=ios) itype

if (itype.eq.-998.or. &
 itype.eq.-999.or. &
 itype.eq.-997) then
 goto 100
endif

backspace(luloc,err=2000,iostat=ios)
read(luloc,*,err=2000,end=2000,iostat=ios) itype,n

if (itype.ne.2) then
 write(stderr,3)
3        format(1x,'rdsym: only data type 2 implemented')
 call exitp(1)
endif

l = l + 1

SYMB(1,l) = DBLE(itype)
SYMB(2,l) = DBLE(n)

read(luloc,*,err=2000,end=2000,iostat=ios) &
  ((SYMB(i,j+l),i=1,2),j=1,n)

l = l + n

if (l.gt.MAXSYM) then
 write(stderr,4) numsy
4        format(1x,'rdsym: overflow for symbol ',I2)
 call exitp(1)
endif

goto 100

200    call closef(luloc)

if (iecho.ge.4) then
 write(stderr,5) l,numsy
5       format(1x,'rdsym: read ',I2,' entries for symbol number ', &
 I2)
endif

NSYM = l

if (numsy.eq.8 .or. &
 numsy.eq.9 .or. &
 numsy.eq.13.or. &
 numsy.eq.37.or. &
 numsy.eq.47.or. &
 numsy.eq.48.or. &
 numsy.eq.49) then 
 call ROTSYM(SYMB, NSYM,90d0)
endif
!
return

1000    write(stderr,1001) ios,numsy
1001    format(1x,'rdsymb: internal read error ',I3,' while reading ', &
  I7)

call prerr('rdsymb')
call exitp(1)

2000    write(stderr,2001) ios
2001    format(1x,'rdsymb: input error ',I3)

call prerr('rdsymb')
call exitp(1)

end
!-------------------------------------------------------------------------------
logical function RDTOO (IPRINT,ISTEP,NPRT,LMODE)

! determines if output is on file
implicit none

!-pass
integer ISTEP,NPRT
logical LMODE
integer :: iprint, iprt
dimension IPRINT(*)
!-init

if (.not. LMODE) then
 RDTOO = .false.
 goto 300
endif

if (ISTEP.eq.0) then
 RDTOO = .true.
 goto 300
endif

RDTOO = .false.

if (NPRT.eq.0) goto 300

if (NPRT.lt.0) then
 if (MOD(ISTEP,IABS(NPRT)).eq.0) RDTOO = .true.
 goto 300
endif

!    check if ISTEP occurs in IPRINT
IPRT = 1

100    if (IPRT.le.NPRT) then
 if (IPRINT(IPRT).eq.ISTEP) goto 200
 IPRT = IPRT + 1
 goto 100
endif

goto 300
200    RDTOO = .true.
300    return
end
!-------------------------------------------------------------------------------
subroutine readij (IJSUR,NUMEL,NSUR,IJFILE,ierr)

use iomodule,        only: stdout
use errors,          only: exec, ok


implicit none
!-pass
integer NUMEL,NSUR,ierr
double precision :: ijSur
dimension IJSUR(NSUR,*)
character IJFILE*(*)


!-locl
integer luij,i,j
integer, external :: ifbopen, nextlu
!-init
luij  = nextlu(0)
ierr  = OK

!     WIENAND begin
  write(stdout, *) 'inREADIJ: ', IJFILE
!     WIENAND end

ierr = ifbopen (luij,IJFILE,'old')

if (ierr.ne.0) then
 ierr = EXEC
 return
endif

read(luij,err=100) &
  ((IJSUR(i,j),i=1,NSUR),j=1,NUMEL)

call closef(luij)

50    return

100    call prerr('readij')
ierr = EXEC
goto 50
end
!-------------------------------------------------------------------------------
subroutine shell(prgname,ierr)

use debugmodule,     only: iecho
use iomodule,        only: stderr, fcmd, fstack, stdin
use errors,          only: ok, prgerr, syntax

implicit none
!-pass
integer ierr
character prgname*(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) arg
integer, external :: lnblk
integer :: iarg, narg
ierr = OK

narg  = command_argument_count()
fcmd  = ' '
1    format(1x,80A)

call get_command_argument(0,prgname)

if (narg.gt.0) then
    do iarg=1,narg
        call get_command_argument(iarg,arg)

        if (arg.eq.'echo') then
            iecho = 1
        elseif (arg.eq.'noecho') then
            iecho = 0
        else
            if (lnblk(fcmd).le.0) then
                if (lnblk(arg).gt.LEN(fcmd)) then
                    write(stderr,2)
2                   format(1x,'shell: LEN(fcmd) insufficient')
                    ierr = PRGERR
                    return
                endif
                fcmd = arg(1:lnblk(arg))
            else
                write(stderr,1) 'Usage: "', &
                    prgname(1:lnblk(prgname)), &
                  ' [no[echo]] [commandfile]'

                ierr = SYNTAX
                fcmd = ' '
                goto 20
            endif
        endif
    enddo
endif

20    fstack(0) = fcmd

if (lnblk(fcmd).le.0) fstack(0) = 'stdin'

return
end
!-------------------------------------------------------------------------------
subroutine tst1 (A,N,istep)

use iomodule, only: stderr

implicit none

integer          :: N, istep, ii, jj
double precision :: A(N)


write(stderr,'(''A('',I3,'')'')') N

jj = MAX0(istep,1)

do ii=1,N
 if (MOD(ii,jj).eq.0) then
write(stderr,'(I4,1PG12.4)') ii,A(ii)
 endif
enddo

return
end
!-------------------------------------------------------------------------------
subroutine tst2 (A,N,M,istep)

use iomodule, only: stderr

implicit none

integer          :: N, M, istep, ii, jj, nn
double precision :: A(N,M)


write(stderr,'(''A('',I3,'','',I3,'')'')') N,M

jj = MAX0(istep,1)

do ii=1,M
    if (MOD(ii,jj).eq.0) then
        write(stderr,'(I4,4(1PG12.4))') ii,(A(nn,ii),nn=1,N)
    endif
enddo

return
end

!-------------------------------------------------------------------------------
subroutine tst2Integers (A,N,M,istep)

use iomodule, only: stderr

implicit none

integer          :: N, M, istep, ii, jj, nn
integer          :: A(N,M)


write(stderr,'(''A('',I3,'','',I3,'')'')') N,M

jj = MAX0(istep,1)

do ii=1,M
    if (MOD(ii,jj).eq.0) then
        write(stderr,'(I4,4(1PG12.4))') ii,(A(nn,ii),nn=1,N)
    endif
enddo

return
end


!-------------------------------------------------------------------------------
subroutine tstinp (ITST)

USE MESHDATAMODULE
use modeldatamodule, only: modeldatactx
use materialsmodule, only: grav
use debugmodule,     only: iecho
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldefinition, only: mode, itmode, numsn, numel, numnp, &
                           nElemsWithSlipperyNodes
#else
use modeldefinition, only: mode, itmode, numsn, numel, numnp
#endif
use modeltopology,   only: nen, ndof, nsd, nstr
use iomodule,        only: stderr

! subroutine to write input data
!    ITST =  0: write FE-input file data
!    ITST <> 0: write FE-output file data

implicit none



!-local
integer :: itst, increm

!-init
increm = 1

if (iecho.ne.2) then
    return
endif

if (ITST.eq.0) then
 write(stderr,'(''FE-input file contents'')')

 write(stderr,'(''nodal coordinates:'')')
 call tst2(meshdatactx%X,NSD,NUMNP,increm)

 write(stderr,'(/''boundary condition codes:'')')
 call tst2integers(modeldatactx%ibond,NDOF,NUMNP,increm)

 write(stderr,'(/''boundary conditions:'')')
 call tst2(modeldatactx%bond,NDOF,NUMNP,increm)

 write(stderr,'(/''gravity:'')')
 call tst1(grav,NDOF,1)

 write(stderr,'(/''element definitions:'')')
 call itst2(meshdatactx%ien,NEN,NUMEL,increm)

 if (MODE.ge.3) then
if (ITMODE.eq.0) then
    write(stderr,'(/''initial temperatures:'')')
    call tst1(modeldatactx%T,NUMNP,increm)
endif
write(stderr,'(/''thermal boundary condition codes:'')')

call itst1(modeldatactx%ITBC,NUMNP,increm)
write(stderr,'(/''thermal boundary conditions:'')')

call tst1(modeldatactx%TBC,NUMNP,increm)
 endif
else
 if (MODE.eq.2 .or. MODE.eq.6) then
write(stderr,'(''Mechanical FE-output files content'')')
write(stderr,'(''displacements:'')')
call tst2 (modeldatactx%D,NDOF,NUMNP,increm)
write(stderr,'(/''differential displacements:'')')

if (NUMSN.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
    call tst2 (modeldatactx%DXE,NDOF,NEN*nElemsWithSlipperyNodes,1)
#else
    call tst2 (modeldatactx%DXE,NDOF,NEN*NUMEL,1)
#endif
endif

write(stderr,'(/''stresses:'')')
call tst2 (modeldatactx%stn,NSTR,NUMEL,increm)
 endif

 if (MODE.eq.5 .or. MODE.eq.6) then
write(stderr,'(''Thermal FE-output files contents'')')
write(stderr,'(''temperatures:'')')
call tst1 (modeldatactx%T,NUMNP,increm)
write(stderr,'(/''heat fluxes:'')')
call tst2 (modeldatactx%Thermalflux, &
           NSD,NUMEL,increm)
 endif
endif
return
end


!-------------------------------------------------------------------------------

#ifdef EXPERIMENTAL_NORATESTORAGE
subroutine femout(D,STN,STR,DXE,PRPMAT,PRPLAS,T,TX, &
    FLUX,SURF,IDSLE,DELD,DELDXE,IREQ,mechrd,temprd, &
    ierr)
#else
subroutine femout(D,STN,STR,DXE,PRPMAT,PRPLAS,RATE,T,TX, &
    FLUX,SURF,IDSLE,DELD,DELDXE,IREQ,mechrd,temprd, &
    ierr)
#endif

use constants,       only: dfmax
USE ALGEBRA
USE MODELDEFINITION
uSE MESHDATAMODULE
USE MODELDATAMODULE
use formatsmodule,   only: FMT_fedsk2D_displacement, &  
                           FMT_fedsk2D_velocities, &
                           FMT_fedsk2D_diff_displacement, &
                           FMT_fedsk3D_displacement, &
                           FMT_fedsk3D_velocities, &
                           FMT_fedsk_stressStrainRate_nstr3, &
                           FMT_fedsk_stressStrainRate_nstr4, &
                           FMT_fedsk_stressStrainRate_nstr5, &
                           FMT_fedsk_stressStrainRate_nstr6, &
                           FMT_fedsk3D_diff_displacement
use debugmodule,     only: iecho
use modeltopology,   only: NSD,NDOF,NSTR,NEN

use extrema,         only: stressPositiveSmallest,       stressPositiveBiggest, &
                           stressNegativeSmallest,       stressNegativeBiggest, &
                           strainPositiveSmallest,       strainPositiveBiggest, &
                           strainNegativeSmallest,       strainNegativeBiggest, &
                           strainRatePositiveSmallest,   strainratePositiveBiggest, &
                           strainRateNegativeSmallest,   strainRateNegativeBiggest, &
                           displacementPositiveSmallest, displacementPositiveBiggest, &
                           displacementNegativeSmallest, displacementNegativeBiggest
use materialsmodule, only: LMAT
use iomodule
use errors,          only: exec, ok


!
! routine to read FEDSK files
!

implicit none
!-pass
logical mechrd,temprd
integer IREQ

integer          :: idsle
double precision :: prpmat, t, prplas, surf, tx, deld, deldxe
#ifndef EXPERIMENTAL_NORATESTORAGE
double precision :: rate
#endif

double precision :: d(NDOF,NUMNP)
double precision :: DXE(NDOF,NEN,NUMEL)
double precision :: FLUX(NSD,NUMEL)
double precision :: STN(NSTR,NUMEL)
double precision :: STR(NSTR,NUMEL)

dimension PRPMAT(6,*), &
  T(*),PRPLAS(9,*),SURF(NSD,*), DELD(NDOF,*),&
    IDSLE(*),TX(*), DELDXE(NDOF,NEN,*)

#ifndef EXPERIMENTAL_NORATESTORAGE
dimension   RATE(NSTR,*)
#endif

!-locl
#ifdef EXPERIMENTAL_ReducedSlipAlloc
character(len=255) :: dataRecord
#endif
character(len=255) :: frmt
character(len=255) :: datatype ! to pass the datatype on the a possible error
double precision   :: effectiveSum
integer            :: iVertex, iDOF

integer            :: i, iElem, j, jj, k, m, n
integer            :: io
integer            :: ierr

!-init
ierr  = OK
!
stressPositiveSmallest = dfmax
stressNegativeSmallest = -dfmax
stressPositiveBiggest = 0d0
stressNegativeBiggest = 0d0

strainPositiveSmallest = dfmax
strainNegativeSmallest = -dfmax
strainPositiveBiggest = 0d0 
strainNegativeBiggest = 0d0

strainRatePositiveSmallest = dfmax
strainRateNegativeSmallest = -dfmax
strainRatePositiveBiggest = 0d0
strainRateNegativeBiggest = 0d0

displacementPositiveSmallest = dfmax
displacementNegativeSmallest = -dfmax
displacementPositiveBiggest = 0d0 
displacementNegativeBiggest = 0d0


if ((.not.temprd).and.(.not.mechrd)) then
    return
endif

!*****************************************************
!    read temperatures data block
!*****************************************************

if (iecho.eq.2) then
     write(stderr,*) 'Reading output data'
endif

if (temprd) then
    if (IREQ.gt.0 .and. IADVEC.ne.0) then
        call clear(D,NDOF*NUMNP,"D")
        if (ascin) then
            write(frmt,80) NDOF
80          format('(I12,1X,',I1,'E15.0)')
            do n=1,NUMNP
                read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
                 j,(D(i,n),i=1,NDOF)
                if (j.ne.n) then
                    write(stderr,85) n
                    ierr = EXEC
                    return
                endif
            enddo
        else  ! read as binary
            read(luin,err=1000,end=2000,iostat=io) &
               ((D(i,n),i=1,NDOF),n=1,NUMNP)
        endif
        if (iecho.eq.2) then
            write(stderr,90)
        endif
    endif
    call clear (T,NUMNP,"T")
    if (ascin) then
        frmt='(I12,1X,E15.0)'
        do n=1,NUMNP
            read(luin,frmt,err=1000,end=2000,iostat=io) &
           j,T(n)
            if (j.ne.n) then
                write(stderr,10) n
10              format(1x,'femout [a]: synchronization error in', &
                1x,'temperatures data block for element ',I12)
                ierr = EXEC
                return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) (T(n),n=1,NUMNP)
    endif

    if (iecho.eq.2) then 
        write(stderr,15)
15        format(1x,'temperatures have been read')
    endif

    if (NUMFN+NUMSN.gt.0.and.IDIFT.eq.1) then

        call clear (TX,NUMNP,"TX")
        frmt='(I12,1X,F11.0)'
        if (ascin) then
            read(luin,frmt,err=1000,end=2000,iostat=io) &
               (j,TX(j),m=1,NUMNP)
        else
            read(luin,err=1000,end=2000,iostat=io) (TX(j),j=1, &
               NUMNP)
        endif

        if (iecho.eq.2) then 
            write(stderr,30)
30            format(1x,'differential temperatures have been read')
        endif
    endif

!************ read the thermal fluxes

    call clear (FLUX,NSD*NUMEL,"flux")

    if (ascin) then
        write(frmt,35) NSD
35      format('(I12,1X,',I1,'E15.0)')
        do n=1,NUMEL
            read(luin,fmt=frmt,err=1000,end=2000,iostat=io) j,(FLUX(i,n),i=1,NSD)
            if (iecho.eq.22) then
                write(*,*) 'read fluxes: ', n, (FLUX(1:NSD,n))
            endif
            if (j.ne.n) then
                write(stderr,40) n
40                   format(1x,'femout [b]: synchronization error in', &
                 1x,'heat flux data block for element ',I12)
                 ierr = EXEC
                 return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) &
             ((FLUX(i,n),i=1,NSD),n=1,NUMEL)
    endif

    if (iecho.eq.2) then
        write(stderr,45)
45        format(1x,'heat fluxes have been read')
    endif

endif

!*****************************************************
!  Read mechanical output
!*****************************************************

if (mechrd) then

    !************ PRPMAT **********************
    if (IREQ.gt.0) then
        if (LMAT) then
            m = NUMEL
        else
            m = NUMAT
        endif

        if (IDIG(ICVIS,4,1).ne.0 .and. MODE.eq.6) then
            if (ascin) then
                do n=1,m
                    read(luin,fmt=50,err=1000,end=2000,iostat=io) &
                     j,PRPMAT(3,n),PRPMAT(4,n)
50                  format(I12,1X,E15.0,1X,F4.0)
                    if (j.ne.n) then
                        write(stderr,55) n
55                      format(1x,'femout [c]: synchronization error', &
                         1x,'in viscosity data block for material ', &
                         I12)
                        ierr = EXEC
                        return
                    endif
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                (PRPMAT(3,n),PRPMAT(4,n),n=1,m)
            endif
            if (iecho.eq.2) write(stderr,60)
60          format(1x,'viscosities have been read')
        endif

        !************ PRPLAS **********************
        if (IDIG(ICVIS,4,2).ne.0) then
            if (ascin) then
                do n=1,m
                    read(luin,fmt=65,err=1000,end=2000,iostat=io) &
                     j,PRPLAS(3,n),PRPLAS(4,n)
65                  format(I12,1X,E15.0,1X,E15.0)
                    if (j.ne.n) then
                        write(stderr,70) n
70                      format(1x,'femout [d]: synchronization error', &
                         1x,'in plastic data block for material ', &
                         I12)
                        ierr = EXEC
                        return
                    endif
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                (PRPLAS(3,n),PRPLAS(4,n),n=1,m)
            endif
            if (iecho.eq.2) write(stderr,75)
75          format(1x,'plastic parameters have been read')
        endif
    endif


    !****************** displacement ********************
    call clear_matrix(modeldatactx%DELD, NDOF, NUMNP)
    datatype = "displacements"
    call clear_matrix(modeldatactx%D, NDOF, NUMNP)
    if (ascin) then
        if (NDOF.eq.2) then
            frmt = FMT_fedsk2D_displacement
        else
            frmt = FMT_fedsk3D_displacement
        endif
        do n=1,NUMNP
            read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
              j,(modeldatactx%D(i,n),i=1,NDOF)
            if (j.ne.n) then
                write(stderr,85) n
85              format(1x,'femout [e]: synchronization error in', &
                 1x,'displacements data block for node ',I12)
                ierr = EXEC
                return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) &
        ((modeldatactx%D(i,n),i=1,NDOF),n=1,NUMNP)
    endif
    if (iecho.eq.2) write(stderr,90)
90  format(1x,'displacements have been read')

    ! determine displacement extrema
    do iVertex=1,NUMNP
        effectiveSum = 0d0
        do iDOF=1,NDOF
            effectiveSum = effectiveSum + modeldatactx%D(iDOF,iVertex)**2
        enddo
        effectiveSum = sqrt(effectiveSum)
        if (effectiveSum.ge.0d0) then
            if (effectiveSum.gt.displacementPositiveBiggest) then
                displacementPositiveBiggest = effectiveSum
            endif
            if (effectiveSum.lt.displacementPositiveSmallest) then
                displacementPositiveSmallest = effectiveSum
            endif
        else if (effectiveSum.le.0d0) then
            if (effectiveSum.lt.displacementNegativeBiggest) then
                displacementNegativeBiggest = effectiveSum
            endif
            if (effectiveSum.gt.displacementNegativeSmallest) then
                displacementNegativeSmallest = effectiveSum
            endif
        endif
    enddo

#ifdef EXPERIMENTAL_FRICTION

    !******************** Differential forces  *******************
    if (numslpglobal.gt.0) then
        datatype = "diff forces"
        call clear_matrix(modeldatactx%Diforc, NDOF, NUMNP)
        if (ascin) then
            if (NDOF.eq.2) then
                frmt = FMT_fedsk2D_displacement
            else
                frmt = FMT_fedsk3D_displacement
            endif
        
            do n=1,NUMNP
                read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%Diforc(i,n),i=1,NDOF)
                if (j.ne.n) then
                    write(stderr,86) n
86                  format(1x,'femout [e]: synchronization error in', &
                     1x,'differential forces data block for node ',I12)
                    ierr = EXEC
                    return
                endif
            enddo
        else
            read(luin,err=1000,end=2000,iostat=io) &
             ((modeldatactx%Diforc(i,n),i=1,NDOF),n=1,NUMNP)
        endif
    endif
    if (iecho.eq.2) write(stderr,91)
91  format(1x,'differential forces have been read')

#endif

    if (IREQ.gt.0) then
        !******************** velocity  *******************
        datatype = "velocities"
        if (ascin) then
            do n=1,NUMNP
                if (NDOF.eq.2) then
                    read(luin,fmt=FMT_fedsk2D_velocities,err=1000,end=2000,iostat=io) &
                     j,(modeldatactx%DELD(i,n),i=1,NDOF)
                else
                    read(luin,fmt=FMT_fedsk3D_velocities,err=1000,end=2000,iostat=io) &
                     j,(modeldatactx%DELD(i,n),i=1,NDOF)
                endif
                if (j.ne.n) then
                    write(stderr,95) n
95                  format(1x,'femout [e]: synchronization error', &
                     'in velocities data block for node ',I12)
                    ierr = EXEC
                    return
                endif
            enddo
        else
            read(luin,err=1000,end=2000,iostat=io) &
             ((modeldatactx%DELD(i,n),i=1,NDOF),n=1,NUMNP)
        endif
        if (iecho.eq.2) write(stderr,96)
96      format(1x,'velocities have been read')
    endif 

    !******************** DXE (differential displacement)  *******************
    datatype = "differential displacements"
    if (NUMSN.gt.0) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc
        call CLEAR (DXE,NDOF*NEN*nElemsWithSlipperyNodes)
        if (ascin) then
            if (NDOF.eq.2) then
                frmt = FMT_fedsk2D_diff_displacement
            else if (NDOF.eq.3) then
                frmt = FMT_fedsk3D_diff_displacement
            else
                write(stderr,*) "diff disp says: NDOF not 2 or 3; contact model support"
            endif
    
            if (iecho.eq.2) then
                write(stderr,*) 'Going to read diff displ with format: ', frmt
            endif
    
            do n=1,nElemsWithSlipperyNodes
                do k=1,NEN  ! quick and dirty fix
                    read(luin,"(a)",err=1000,end=2000, iostat=io) dataRecord
    
                    read(dataRecord,fmt=frmt,err=1000,end=2000, iostat=io) j,jj,(DXE(i,k,n),i=1,NDOF)
    
                enddo
            enddo
        else
            ! read as binaries
            read(luin,err=1000,end=2000,iostat=io) &
            (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
        endif
#else
        call CLEAR (DXE,NDOF*NEN*NUMEL,"DXE")
        if (ascin) then
            if (NDOF.eq.2) then
                frmt = FMT_fedsk2D_diff_displacement
            else if (NDOF.eq.3) then
                frmt = FMT_fedsk3D_diff_displacement
            endif
            do n=1,NUMEL
                do k=1,NEN
                    read(luin,fmt=frmt,err=1000,end=2000, iostat=io) j,jj,(DXE(i,k,n),i=1,NDOF)
                    if (j.ne.n) then
                        write(stderr,100) n
  100                   format(1x,'femout [f]: synchronization',1x, &
                        'error in slippery nodes data block', &
                        ' for node ',I12)
                        ierr = EXEC
                        return
                    endif
                enddo
            enddo
        else
            ! read as binaries
            read(luin,err=1000,end=2000,iostat=io) &
            (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
        endif
#endif

        if (iecho.eq.2) write(stderr,105)
  105   format(1x,'slippery displacements have been read')

!******************** DELDXE (differential velocities)  *******************
        datatype = "differential velocities"

        if (IREQ.gt.0) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc
            call CLEAR (modeldatactx%DELDXE,NDOF*NEN*nElemsWithSlipperyNodes,"modeldatactx%DELDXE")
#else
            call CLEAR (modeldatactx%DELDXE,NDOF*NEN*NUMEL,"modeldatactx%DELDXE")
#endif
            modeldatactx%DELDXE = 0d0

#ifdef EXPERIMENTAL_ReducedSlipAlloc
            if (ascin) then 
                do n=1,nElemsWithSlipperyNodes
                    do k=1,NEN

                        read(luin,*) dataRecord

!                        write(*,*) "read record     ***"// dataRecord(1:50)//"***"


                        read(dataRecord,fmt=frmt,err=1000,end=2000, iostat=io) &
                         j,jj, (modeldatactx%DELDXE(i,k,n),i=1,NDOF)

!                        write(*,*) "read diff. velo ", j, jj, modeldatactx%DELDXE(1:NDOF,k,n)

                    enddo
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                (((modeldatactx%DELDXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,nElemsWithSlipperyNodes)
            endif
#else
            if (ascin) then
                do n=1,NUMEL
                    do k=1,NEN
                        read(luin,fmt=frmt,err=1000,end=2000, &
                         iostat=io) j,jj, &
                         (modeldatactx%DELDXE(i,k,n),i=1,NDOF)
                        if (j.ne.n) then
                            write(stderr,107) n
  107                       format(1x,'femout [g]: synchronization', &
                             ' error in slippery velocities', &
                             ' for node ',I12)
                            ierr = EXEC
                            return
                        endif
                    enddo
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                (((modeldatactx%DELDXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
            endif
#endif
            if (iecho.eq.2) write(stderr,108)
  108       format(1x,'differential velocities have been read')
        endif

    endif   ! NUMSN>0


    !********************** stress ************
    datatype = "stress"
    call clear (STN,NSTR*NUMEL,"STN")
    if (ascin) then
        write(frmt,130) NSTR
  130   format('(I12,1X,',I1,'E15.0)')
        do n=1,NUMEL
            if (nstr.eq.3) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STN(i,n),i=1,NSTR)
            else if (nstr.eq.4) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STN(i,n),i=1,NSTR)
            else if (nstr.eq.5) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STN(i,n),i=1,NSTR)
            else if (nstr.eq.6) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STN(i,n),i=1,NSTR)
            endif
    
            if (j.ne.n) then
                write(stderr,135) n
  135           format(1x,'femout [j]: synchronization error in', &
                1x,'stresses data block for element ',I12)
                ierr = EXEC
                return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) &
        ((modeldatactx%STN(i,n),i=1,NSTR),n=1,NUMEL)
    endif
    if (iecho.eq.2) write(stderr,140)
140 format(1x,'stresses have been read')

    ! determine stress extrema
    do iElem=1,NUMEL
        effectiveSum = 0d0
        do iDOF=1,NDOF
            effectiveSum = effectiveSum + modeldatactx%STN(iDOF,iElem)
        enddo
        if (effectiveSum.ge.0d0) then
            if (effectiveSum.gt.stressPositiveBiggest) stressPositiveBiggest = effectiveSum
            if (effectiveSum.lt.stressPositiveSmallest) stressPositiveSmallest = effectiveSum
        else if (effectiveSum.le.0d0) then
            if (effectiveSum.lt.stressNegativeBiggest) stressNegativeBiggest = effectiveSum
            if (effectiveSum.gt.stressNegativeSmallest) stressNegativeSmallest = effectiveSum
        endif
    enddo


    !****************** strains ****************************
    datatype = "strain"
    call clear (STR,NSTR*NUMEL,"STR")
    if (ascin) then
        do n=1,NUMEL
            ! format adapts to long or short data format
            if (nstr.eq.3) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                j,(modeldatactx%STR(i,n),i=1,NSTR)
            else if (nstr.eq.4) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                j,(modeldatactx%STR(i,n),i=1,NSTR)
            else if (nstr.eq.5) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                j,(modeldatactx%STR(i,n),i=1,NSTR)
            else if (nstr.eq.6) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                j,(modeldatactx%STR(i,n),i=1,NSTR)
            endif
            if (j.ne.n) then
                write(stderr,142) n
  142           format(1x,'femout [k]: synchronization error in', &
                 1x,'strains data block for element ',I12)
                ierr = EXEC
                return
            endif
        enddo
     else
        read(luin,err=1000,end=2000,iostat=io) &
         ((modeldatactx%STR(i,n),i=1,NSTR),n=1,NUMEL)
     endif
     if (iecho.eq.2) write(stderr,143)
143  format(1x,'strains have been read')

     ! Scaling off-diagonal strain components (engineering strain to infinitesimal)
     if (NTYPE.eq.1) then
         ! pln
         do n=1,NUMEL
             modeldatactx%STR(3,n) = modeldatactx%STR(3,n)*5D-1
         enddo
     elseif (NTYPE.eq.2) then
         ! opn
         do n=1,NUMEL
             modeldatactx%STR(3,n) = modeldatactx%STR(3,n)*5D-1
             modeldatactx%STR(4,n) = modeldatactx%STR(4,n)*5D-1
             modeldatactx%STR(5,n) = modeldatactx%STR(5,n)*5D-1
         enddo
     elseif (NTYPE.eq.3) then
         ! axi
         do n=1,NUMEL
             modeldatactx%STR(4,n) = modeldatactx%STR(4,n)*5D-1
         enddo
     elseif (NTYPE.eq.6) then
         ! f3d
         do n=1,NUMEL
             modeldatactx%STR(4,n) = modeldatactx%STR(4,n)*5D-1
             modeldatactx%STR(5,n) = modeldatactx%STR(5,n)*5D-1
             modeldatactx%STR(6,n) = modeldatactx%STR(6,n)*5D-1
         enddo
     else
        write(stderr,*) 'femout: unknown NTYPE=',NTYPE
        call exitp(1)
     endif

     do iElem=1,NUMEL
         effectiveSum = 0d0
         do iDOF=1,NDOF
             effectiveSum = effectiveSum + modeldatactx%STR(iDOF,iElem)
         enddo
         if (effectiveSum.ge.0d0) then
             if (effectiveSum.gt.strainPositiveBiggest) then
                 strainPositiveBiggest = effectiveSum
             endif
             if (effectiveSum.lt.strainPositiveSmallest) then
                 strainPositiveSmallest = effectiveSum
             endif
         else if (effectiveSum.le.0d0) then
             if (effectiveSum.lt.strainNegativeBiggest) then
                 strainNegativeBiggest = effectiveSum
             endif
             if (effectiveSum.gt.strainNegativeSmallest) then
                 strainNegativeSmallest = effectiveSum
             endif
         endif
     enddo

#ifndef EXPERIMENTAL_NORATESTORAGE

     !****************** strain rates ****************************
     datatype = "strain rates"
     if (IREQ.gt.0) then
        if (ascin) then
            write(frmt,145) NSTR
  145       format('(I12,1X,',I1,'E15.0)')
            do n=1,NUMEL
                if (nstr.eq.3) then
                    read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                        j,(modeldatactx%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.4) then
                    read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                        j,(modeldatactx%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.5) then
                    read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                        j,(modeldatactx%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.6) then
                    read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                        j,(modeldatactx%RATE(i,n),i=1,NSTR)
                endif
                if (j.ne.n) then
                    write(stderr,150) n
  150               format(1x,'femout [l]: synchronization error', &
                     1x,'in rates data block for element ',I12)
                    ierr = EXEC
                    return
                endif
            enddo
        else
            read(luin,err=1000,end=2000,iostat=io) &
             ((modeldatactx%RATE(i,n),i=1,NSTR),n=1,NUMEL)
        endif


        if (iecho.eq.2) write(stderr,155)
  155   format(1x,'total strain rates have been read')

        ! Scaling off-diagonal strain rate components 
        if (NTYPE.eq.1) then
            ! pln
            do n=1,NUMEL
                modeldatactx%RATE(3,n) = modeldatactx%RATE(3,n)*5D-1
            enddo
        elseif (NTYPE.eq.2) then
            ! opn
            do n=1,NUMEL
                modeldatactx%RATE(3,n) = modeldatactx%RATE(3,n)*5D-1
                modeldatactx%RATE(4,n) = modeldatactx%RATE(4,n)*5D-1
                modeldatactx%RATE(5,n) = modeldatactx%RATE(5,n)*5D-1
            enddo
        elseif (NTYPE.eq.3) then
            ! axi
            do n=1,NUMEL
                modeldatactx%RATE(4,n) = modeldatactx%RATE(4,n)*5D-1
            enddo
        elseif (NTYPE.eq.6) then
            ! f3d
            do n=1,NUMEL
                modeldatactx%RATE(4,n) = modeldatactx%RATE(4,n)*5D-1
                modeldatactx%RATE(5,n) = modeldatactx%RATE(5,n)*5D-1
                modeldatactx%RATE(6,n) = modeldatactx%RATE(6,n)*5D-1
            enddo
        endif
    
        ! determine strain rate extrema
        do iElem=1,NUMEL
            effectiveSum = 0d0
            do iDOF=1,NDOF
                effectiveSum = effectiveSum + modeldatactx%RATE(iDOF,iElem)
            enddo
            if (effectiveSum.ge.0d0) then
                if (effectiveSum.gt.strainratePositiveBiggest) then
                    strainratePositiveBiggest = effectiveSum
                endif
                if (effectiveSum.lt.strainratePositiveSmallest) then
                    strainratePositiveSmallest = effectiveSum
                endif
            else if (effectiveSum.le.0d0) then
                if (effectiveSum.lt.strainrateNegativeBiggest) then
                    strainrateNegativeBiggest = effectiveSum
                endif
                if (effectiveSum.gt.strainrateNegativeSmallest) then
                    strainrateNegativeSmallest = effectiveSum
                endif
            endif
        enddo
    endif

#endif

    !****************** surface uplifts ****************************

    datatype = "surface uplifts"
    if (NSURF.gt.0.and.NSED.gt.0) then
        if (ascin) then
            write(frmt,160) NDOF
  160       format('(I12,1X,',I1,'E15.0)')
            do n=1,NSURF
                read(luin,fmt=frmt,err=1000,end=2000, &
                iostat=io) j,(SURF(i,n),i=1,NSD)
                if (j.ne.n) then
                    write(stderr,165) n
  165               format(1x,'femout [m]: synchronization err', &
                   1x,'in block for surface node ',I12)
                    ierr = EXEC
                    return
                endif
            enddo
        else
            read(luin,err=1000,end=2000,iostat=io) &
            ((SURF(i,n),i=1,NSD),n=1,NSURF)
        endif
        if (iecho.eq.2) write(stderr,170)
  170   format(1x,'surface uplifts have been read')
    endif

endif ! read mechanical results

return 

!.... ....................... R E A D   E R R O R S ............................

1000    write(stderr,*) 'femout: Error occurred while trying to read ', trim(datatype), &
                        ' of timestep ', IREQ, &
                        ' error: ', io
ierr = EXEC
return

2000    write(stderr,*) 'femout: Unexpected EOF while trying to read ', trim(datatype), &
                        ' of timestep ', IREQ, &
                        ' error: ', io

ierr = EXEC
return

end subroutine
!-------------------------------------------------------------------------------

#ifdef EXPERIMENTAL_NORATESTORAGE
subroutine fewrite (D,STN,STR,DXE,IDSLE,PRPMAT, &
    IEN,PRPLAS,SURF,T,FLUX,TX,IREQ,MECHRD,TEMPRD,LUFO,IERR)
#else
subroutine fewrite (D,STN,STR,DXE,IDSLE,PRPMAT, &
    IEN,PRPLAS,RATE,SURF,T,FLUX,TX,IREQ,MECHRD,TEMPRD,LUFO,IERR)
#endif
!
! routine to write FE output file in original TECTON format
!
USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif
USE MATERIALSMODULE, only: LMAT
use debugmodule,     only: iecho
use algebra,         only: idig
use iomodule,        only: ascout, stderr, luout
use errors,          only: exec, ok


implicit none
!-pass
logical MECHRD,TEMPRD
integer IREQ,LUFO

integer          :: idsle, ien


double precision :: d, stn, dxe, prpmat, t, flux, prplas, surf, str, tx

dimension D(NDOF,*),STN(NSTR,*),DXE(NDOF,NEN,*),PRPMAT(6,*), &
  T(*),FLUX(NSD,*),PRPLAS(9,*),IDSLE(*),SURF(NSD,*), &
    IEN(NEN,*),STR(NSTR,*),TX(*)

#ifndef EXPERIMENTAL_NORATESTORAGE
double precision :: rate
dimension RATE(NSTR,*)
#endif

integer          :: ierr


!-locl
character(len=40)frmt
integer :: i,j,m,n,io
!-init
ierr  = OK
!
if ((.not.temprd).and.(.not.mechrd)) return
!
!    write temperatures data block
if (TEMPRD) then
if (IREQ.gt.0 .and. IADVEC.ne.0) then
!              switch for advecting grid implemented, but not tested.
   if (ascout) then
       write(frmt,1) NDOF
1              format('(I5,1X,',I1,'(1PG15.6))')
       do n=1,NUMNP
           write(lufo,fmt=frmt,err=1000,iostat=io) &
            n,(D(i,n),i=1,NDOF)
       enddo
   else
       write(lufo,err=1000,iostat=io) &
        ((D(i,n),i=1,NDOF),n=1,NUMNP)
   endif
   if (iecho.eq.2) write(stderr,2)
2          format(1x,'displacements have been written')
endif
 if (ascout) then
 write(lufo,3,err=1000,iostat=io) &
 (n,T(n),n=1,NUMNP)
3        format(I5,1X,1PG15.6)
 else
 write(lufo,err=1000,iostat=io) &
 (T(n),n=1,NUMNP)
 endif
 if (iecho.eq.2) write(stderr,4)
4        format(1x,'temperatures have been written')
 if (ascout) then
 write(frmt,5) NSD
5        format('(I5,1X,',I1,'(1PG15.6))')
 do n=1,NUMEL
write(lufo,fmt=frmt,err=1000,iostat=io) &
       n,(FLUX(i,n),i=1,NSD)
  enddo
 else
 write(lufo,err=1000,iostat=io) &
 ((FLUX(i,n),i=1,NSD),n=1,NUMEL)
 endif
 if (iecho.eq.2) write(stderr,6)
6        format(1x,'heat fluxes have been written')
if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) then
!                 write (*,*) 'Going to write heat fluxes'
    if (ascout) then
        write(luout(2),23) (n,TX(n),n=1,NUMNP)
23               format(I5,1X,F11.5)
    else
        write(luout(2)) (TX(n),n=1,NUMNP)
    endif
 endif
endif

if (mechrd) then
 if (IREQ.gt.0) then
if (LMAT) then
  m = NUMEL
else
  m = NUMAT
endif
    if (IDIG(ICVIS,4,1).ne.0.and.MODE.eq.6) then
   if (ascout) then
       write(lufo,fmt=7,err=1000,iostat=io) &
            (n,PRPMAT(3,n),PRPMAT(4,n),n=1,m)
7                format(I5,1X,1PG15.6,1X,0PF4.2)
   else
       write(lufo,err=1000,iostat=io) &
            (PRPMAT(3,n),PRPMAT(4,n),n=1,m)
   endif
   if (iecho.eq.2) write(stderr,8)
8               format(1x,'viscosities have been written')
    endif
    if (IDIG(ICVIS,4,2).ne.0) then
   if (ascout) then
   write(lufo,fmt=9,err=1000,iostat=io) &
          (n,PRPLAS(3,n),PRPLAS(4,n),n=1,m)
9               format(I5,1X,1PG15.6,1X,1PG15.6)
   else
   write(lufo,err=1000,iostat=io) &
          (PRPLAS(3,n),PRPLAS(4,n),n=1,m)
   endif
   if (iecho.eq.2) write(stderr,10)
10               format(1x,'plastic parameters have been written')
    endif
 endif
 if (ascout) then
 write(frmt,11) NDOF
11        format('(I5,1X,',I1,'(1PG15.6))')
 do n=1,NUMNP
write(lufo,fmt=frmt,err=1000,iostat=io) &
       n,(D(i,n),i=1,NDOF)
    enddo
 else
write(lufo,err=1000,iostat=io) &
     ((D(i,n),i=1,NDOF),n=1,NUMNP)
 endif
 if (iecho.eq.2) write(stderr,12)
12        format(1x,'displacements have been written')
 if (NUMSN.gt.0) then
if (ascout) then
write(frmt,13) NDOF
13            format('(I5,1X,I5,1X,',I1,'(1PG15.6))')
#ifdef EXPERIMENTAL_ReducedSlipAlloc
do n=1,nElemsWithSlipperyNodes 
    do j=1,NEN
        write(lufo,fmt=frmt,err=1000,iostat=io) &
           elementSlipperyID(n,130),IEN(j,n),(DXE(i,j,n),i=1,NDOF)
    enddo
enddo 
#else
do n=1,NUMEL
    do j=1,NEN
        write(lufo,fmt=frmt,err=1000,iostat=io) &
           n,IEN(j,n),(DXE(i,j,n),i=1,NDOF)
    enddo
enddo
#endif
else
#ifdef EXPERIMENTAL_ReducedSlipAlloc
     write(lufo,err=1000,iostat=io) (((DXE(i,j,n), &
       i=1,NDOF),j=1,NEN),n=1,nElemsWithSlipperyNodes)
#else
     write(lufo,err=1000,iostat=io) (((DXE(i,j,n), &
       i=1,NDOF),j=1,NEN),n=1,NUMEL)
#endif
endif

if (iecho.eq.2) write(stderr,14)

14            format(1x,'differential displacements written')
 endif

 if (ascout) then
 write(frmt,16) NSTR
16        format('(I5,1X,',I1,'(1PG15.6))')
 do n=1,NUMEL
write(lufo,fmt=frmt,err=1000,iostat=io) &
       n,(STN(i,n),i=1,NSTR)
    enddo
 else
write(lufo,err=1000,iostat=io) &
     ((STN(i,n),i=1,NSTR),n=1,NUMEL)
 endif
 if (iecho.eq.2) write(stderr,17)
17        format(1x,'stresses have been written')
 do n=1,NUMEL
 STR(3,n) = STR(3,n)*2d0
 enddo
 if (ascout) then
 do n=1,NUMEL
write(lufo,fmt=frmt,err=1000,iostat=io) &
       n,(STR(i,n),i=1,NSTR)
    enddo
 else
write(lufo,err=1000,iostat=io) &
     ((STR(i,n),i=1,NSTR),n=1,NUMEL)
 endif
 do n=1,NUMEL
STR(3,n) = STR(3,n)*5d-1
 enddo
 if (iecho.eq.2) write(stderr,18)
18        format(1x,'strains have been written')



 if (IREQ.gt.0) then

#ifndef EXPERIMENTAL_NORATESTORAGE
do n=1,NUMEL
RATE(3,n) = RATE(3,n)*2d0
 enddo
if (ascout) then
write(frmt,19) NSTR
19            format('(I5,1X,',I1,'(1PG15.6))')
do n=1,NUMEL
write(lufo,fmt=frmt,err=1000,iostat=io) &
      n,(RATE(i,n),i=1,NSTR)
        enddo
else
write(lufo,err=1000,iostat=io) &
       ((RATE(i,n),i=1,NSTR),n=1,NUMEL)
endif
do n=1,NUMEL
RATE(3,n) = RATE(3,n)*5d-1
 enddo
#endif

if (iecho.eq.2) write(stderr,20)
20            format(1x,'total strain rates have been written')
    if (NSURF.gt.0) then
   if (ascout) then
   write(frmt,21) NSTR
21               format('(I5,1X,',I1,'(1PG15.6))')
   do n=1,NSURF
      write(lufo,fmt=frmt,err=1000,iostat=io) &
         n,(SURF(i,n),i=1,NSD)
           enddo
   else
   write(lufo,err=1000,iostat=io) &
          ((SURF(i,n),i=1,NSD),n=1,NSURF)
   endif
   if (iecho.eq.2) write(stderr,22)
22               format(1x,'surface uplifts have been written')
    endif
 endif
endif

900    return 
!.... ................... W R I T E   E R R O R S .........................
1000    write(stderr,901) io
  901    format(1x,'fewrite: write error ',I3)
call prerr('fewrite')
ierr = EXEC
goto 900
end
!-------------------------------------------------------------------------------
subroutine fewritalt (X,D,STN,IEN,IREQ,MECHRD,TEMPRD, &
                      LUFO,IERR)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use debugmodule,     only: iecho
use iomodule,        only: ascout, stderr
use errors,          only: exec, ok


! routine to write FE output file in alternative format
!
implicit none
!-pass
logical MECHRD,TEMPRD
integer IREQ,LUFO
integer :: ien
double precision :: d, stn, x
dimension D(NDOF,*),STN(NSTR,*),IEN(NEN,*),X(NSD,*)



!-locl
character(len=40)frmt
integer :: i, n, k
integer :: ierr, io
!-init
ierr  = OK
!
if ((.not.temprd).and.(.not.mechrd)) return
!
if (mechrd) then
!            if (IREQ.eq.0) then
if (ascout) then
!                   Write header including number of elms and time.
        write(lufo,1,err=1000,iostat=io) NUMEL,IREQ
1            format(2I5)
write(frmt,2) (NEN*NSD)+(NEN*NDOF)+NSTR
2            format('(I5,',I2,'(1PG15.6))')
do n=1,NUMEL
        write(lufo,fmt=frmt,err=1000,iostat=io) &
        n,((X(i,IEN(k,n)),i=1,NSD),k=1,NEN), &
                ((D(i,IEN(k,n)),i=1,NDOF),k=1,NEN), &
                (STN(i,n),i=1,NSTR)
enddo
else
     write(lufo,err=1000,iostat=io) NUMEL,IREQ
do n=1,NUMEL
   write(lufo,err=1000,iostat=io) n, &
           ((X(i,IEN(k,n)),i=1,NSD),k=1,NEN), &
           ((D(i,IEN(k,n)),i=1,NDOF),k=1,NEN), &
           ((STN(i,n)),i=1,NSTR)
        enddo
endif
if (iecho.eq.2) write(stderr,3)
3            format(1x,'displacements written')
 endif
!    endif

900    return 
!.... ................... W R I T E   E R R O R S .........................
1000    write(stderr,901) io
  901    format(1x,'fewritalt: write error ',I3)
call prerr('fewritalt')
ierr = EXEC
goto 900
end
!-------------------------------------------------------------------------------
subroutine FOUTCP (IREQ,NREQ,IFLAG,IERR)

use iomodule,        only: lufo, stderr, luin
USE MODELDEFINITION, only: MODE
use errors,          only: exec, ok

!
! Routine to copy FE output at requested time IREQ.
! iflag = 0     called by xcp - write in original TECTON output format
! iflag = 1     called by xdu - write FE output in alternative format
! iflag not really used in foutcp, but in subroutine fcp
!
implicit none
!-pass
integer iflag
integer IREQ,NREQ,IERR



!
if (MODE.ne.2 .and. MODE.ne.5 .and. MODE.ne.6) then
 write(stderr,1) MODE
1        format(1x,'FoutCp: MODE = ',I1,' no data on output file')
 ierr = EXEC
 return
endif
!
if (IREQ.eq.0) then
 rewind(luin)
 rewind(lufo)
 if (ierr.ne.OK) return
 call OPNDSK(.false.,ierr)
 if (ierr.ne.0) ierr = EXEC
 if (ierr.ne.OK) return
endif
!
!    copy data block
call fcp (IREQ,IFLAG,IERR)
if (IERR.ne.OK) return
!
return
end
!-------------------------------------------------------------------------------
subroutine writfe (IERR,lu)

USE CONSTANTS, only: year
USE MESHDATAMODULE
USE MODELDATAMODULE
use materialsmodule, only: matdatactx, grav
use modeldefinition, only: NUMROT, numfn, icvis, incomp, &
                           iresdu, ivlim, nlink, nocompr, &
                           igrav, lgdef, npre, nsurf, &
                           numnp, numel, numwnk, mode, &
                           nsed, numslp, numat, numstr, &
                           numpr, nwink
use modeltopology,   only: nen, ndof, nstr
use timestepmodule,  only: nmprt, nprmat, timestepdatactx, &
                           outputcontroldatactx, maxit, nintg
use algebra,         only: idig
use iomodule,        only: stderr
use plotcontrol,     only: numz
use errors,          only: exec, ok, prgerr


!
! Routine to write an f3d (TECTON) input file
! (if that is so, why is it here, and not in plt3d?...)
!
implicit none
!-pass
integer IERR,lu



!-locl
integer, parameter :: MXNUMZ=100   ! hard coded maximum? This can only fail...
character(len=256) mal
double precision :: Z
dimension Z(MXNUMZ)
integer :: i, ios, k
!
write(lu,10,err=1000,iostat=ios)
10    format('automatically generated by OPNPLT')
write(lu,20,err=1000,iostat=ios) NUMZ*NUMNP
20    format(I5)
write(lu,30,err=1000,iostat=ios) MODE,NINTG,NMPRT,MAXIT,NLINK, &
    ICVIS,NPRMAT
30    format(8I5)
mal = '(51I5)'

call wrint(outputcontroldatactx%imprint,MAX(1,NMPRT),mal,lu,IERR)

if (IERR.ne.OK) return
if (NPRMAT.gt.0) then
! todo
! call wrint(A(IPMPRM),NPRMAT,mal,lu,IERR)
 if (IERR.ne.OK) return
endif
mal = '(2I5,1PG14.6,1PG14.6)'

call xwrite (meshdatactx%IEN, &
        meshdatactx%X,  &
       modeldatactx%D,  &
  modeldatactx%LMF,modeldatactx%tfault,modeldatactx%DXE, &
    modeldatactx%SKEW,Z,lu,IERR)


if (IERR.ne.0) return
write(lu,40,err=1000,iostat=ios) 'nodal coordinates'
40    format('end ',40a)
mal = '(5I5)'

do k=1,NUMZ
 call wrigen(modeldatactx%IBOND,NDOF,NUMNP,mal,(k-1)*NUMNP,lu,IERR)
 if (IERR.ne.0) return

enddo

write(lu,40,err=1000,iostat=ios) 'boundary condition codes'
mal = '(2I5,1PG14.6,1PG14.6,1PG14.6)'

do k=1,NUMZ
 call wrfgen(modeldatactx%BOND,NDOF,NUMNP,mal,(k-1)*NUMNP,lu,IERR)
 if (IERR.ne.0) return
enddo

if (IERR.ne.0) return
write(lu,40,err=1000,iostat=ios) 'boundary condition magnitudes'

if (NWINK.gt.0) then
 mal = '(5I5)'
 call wrigen(modeldatactx%IWINK,NDOF,NUMNP,mal,0,lu,IERR)
 if (IERR.ne.0) return
endif

write(lu,40,err=1000,iostat=ios) 'nodal Winkler codes'

if (NWINK.gt.0) then
 mal = '(2I5,1PG14.6,1PG14.6,1PG14.6)'
 call wrfgen(modeldatactx%WINK,NDOF,NUMNP,mal,0,lu,IERR)
 if (IERR.ne.0) return
endif

write(lu,40,err=1000,iostat=ios) 'nodal Winkler magnitudes'

if (NUMROT.gt.0) then
 mal = '(2I5,1PG14.6,1PG14.6)'
 call wrfgen(modeldatactx%SKEW,NDOF,NUMNP,mal,0,lu,IERR)
 if (IERR.ne.0) return
endif

write(lu,40,err=1000,iostat=ios) 'skew angles'
mal = '(50i5)'

if (IERR.ne.0) then
 return
endif

call wrflt(timestepdatactx%DELT,MAX(1,NINTG),1d0/YEAR,mal,lu,IERR,1)

if (IERR.ne.0) then 
 return
endif

write(lu,50,err=1000,iostat=ios) ('year',i=1,MAX(1,NINTG))
50    format(50a5)
mal = '(50f5.1)'
call wrflt(timestepdatactx%alfa,MAX(1,NINTG),1d0,mal,lu,IERR,0)
if (IERR.ne.0) return

write(lu,60,err=1000,iostat=ios) (NUMZ-1)*NUMEL,NUMAT,0,8,NUMPR, &
    NUMSTR,NUMWNK,NUMSLP*2*(NUMZ-1),NUMFN*NUMZ,NPRE,LGDEF,IRESDU, &
    IGRAV,IVLIM,NSURF,NSED,INCOMP,NOCOMPR
60   format(22I5)
! those are 18 numbers printed in a 22 number format...


mal = '(2I5,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6)'
call wrfgen(matdatactx%prpmat,6,NUMAT,mal,0,lu,IERR)

if (IERR.ne.0) return

write(lu,40,err=1000,iostat=ios) 'mechanical property data'

if (NINTG.gt.0.and.IDIG(ICVIS,4,1).ne.0) then
 mal = '(2I5,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6)'
 call wrfgen(matdatactx%PRPCRP,5,NUMAT,mal,0,lu,IERR)
 if (IERR.ne.0) return
 write(lu,40,err=1000,iostat=ios) 'power law creep data'
endif

if (NINTG.gt.0.and.IDIG(ICVIS,4,3).ne.0) then
 mal = '(2I5,9PG14.6)'
!            mal = '(2I5,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6)'                                 '
 call wrfgen(matdatactx%PRPDIF,9,NUMAT,mal,0,lu,IERR)
 if (IERR.ne.0) return
 write(lu,40,err=1000,iostat=ios) 'creep weakening data'
endif

if (NINTG.gt.0.and.IDIG(ICVIS,4,2).ne.0) then
 mal = '(2I5,9PG14.6)'
!            mal = '(2I5,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6)'                                 '
 call wrfgen(matdatactx%PRPLAS,9,NUMAT,mal,0,lu,IERR)
 if (IERR.ne.0) return
 write(lu,40,err=1000,iostat=ios) 'plasticity data'
endif

mal = '(1PG14.6,1PG14.6)'

call wrflt(grav,2,1d0,mal,lu,IERR,0)

if (IERR.ne.0) return

call elmwrt(meshdatactx%ien,meshdatactx%mat,NUMEL,NUMNP,lu,IERR)

if (IERR.ne.0) return

if (NPRE.gt.0) then
    call prewrt(modeldatactx%STN0, &
                modeldatactx%iselm, &
                modeldatactx%ISTIME,NPRE,lu,IERR)

 if (IERR.ne.0) return
endif

if (NUMPR.ne.0) then
 call preswt(modeldatactx%ielno,modeldatactx%iside,modeldatactx%pres,IABS(NUMPR),lu, &
  IERR)
 if (IERR.ne.0) return
endif

if (NUMSTR.gt.0) then
    call strswt(modeldatactx%IELSTR, &
                modeldatactx%iside, &
                modeldatactx%STRS, &
                modeldatactx%ISTR, &
                NUMSTR,lu,IERR)
    if (IERR.ne.0) return
endif

if (NUMWNK.ne.0) then
 call winkwt(modeldatactx%iwelm,modeldatactx%IWSIDE, &
             modeldatactx%IWTIME,modeldatactx%WPRES, &
   NUMWNK,lu,IERR)
 if (IERR.ne.0) return
endif

if (NUMFN.gt.0) then
write(stderr,90)
90       format(1x,'writfe: NUMFN<>0 not implemented')
ierr = PRGERR
return
endif

if (NUMSLP.gt.0) then
 call slipwt(modeldatactx%NSLIP,modeldatactx%DIFORC,NUMSLP,NUMNP,NUMEL, &
   NUMZ,lu,IERR)
 if (IERR.ne.0) return
endif

if (NSURF.gt.0 .and. NINTG.gt.0) then
write(stderr,100)
  100       format(1x,'writfe: NSURF<>0 not implemented')
ierr = PRGERR
return
endif
!
return
!
1000    write(stderr,1001) ios
 1001    format(1x,'writfe: write error ',i4)
IERR = EXEC
return
end
!-------------------------------------------------------------------------------
subroutine wrint (IARRAY,NA,FMT,lu,ierr)

use iomodule,        only: stderr
use errors,          only: exec

implicit none

integer          :: NA,i, ierr, ios,lu
integer          :: iarray(na)
character(len=*) FMT

write(lu,fmt=FMT,err=100,iostat=ios) (int(IARRAY(i)),i=1,NA)
return
100    write(stderr,101) ios
  101    format(1x,'wrint: write error ',I5)
ierr = EXEC
return
end
!-------------------------------------------------------------------------------
subroutine wrflt (ARRAY,NA,factor,FMT,lu,ierr,iflag)

use iomodule,        only: stderr
use errors,          only: exec

implicit none

double precision ARRAY,factor
integer NA,i,iflag, lu, ierr
character(len=*) FMT
dimension ARRAY(NA)


integer :: ios

if (iflag.eq.0) then
 write(lu,fmt=FMT,err=100,iostat=ios)(ARRAY(i)*factor,i=1,NA)
else
 write(lu,fmt=FMT,err=100,iostat=ios) &
   (NINT(ARRAY(i)*factor),i=1,NA)
endif

return
100    write(stderr,101) ios
  101    format(1x,'wrflt: write error ',I5)
ierr = EXEC
return
end
!-------------------------------------------------------------------------------
subroutine wrfgen (ARRAY,N,M,FMT,IOFF,lu,ierr)

use iomodule,        only: stderr
use errors,          only: exec

implicit none

double precision ARRAY
integer N,M,i,j,ios,ierr,ioff,lu
character(len=*) FMT
dimension ARRAY(N,M)

write(lu,fmt=FMT,err=100,iostat=ios) &
    (j+IOFF,0,(ARRAY(i,j),i=1,N),j=1,M)
return
100    write(stderr,101) ios
  101    format(1x,'wrfgen: write error ',I5)
ierr = EXEC
return
end
!-------------------------------------------------------------------------------
subroutine xwrite (IEN,X,D,LMF,TFAULT,DXE,SKEW,Z,lu,ierr)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NEN
use constants,       only: dfmax, useDXE
use debugmodule,     only: iecho
use iomodule,        only: stderr
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use plotControl,     only: numz, iopnflag, opnaxnu, opnscale
use errors,          only: exec

!
implicit none

!-pass
 integer lu,ierr
integer :: ien, lmf
double precision :: x, d, tfault, dxe, skew, z
 dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
    TFAULT(NDOF,*),DXE(NDOF,NEN,*),SKEW(2,*),Z(NUMZ)

!-more comm


!-locl
double precision :: DL(12)
integer :: i, k, m, n
double precision :: dz, zmin, zmax
integer :: ios

!
if (iopnflag.eq.0) then
zmin = dfmax
zmax = -dfmax
 do n=1,NUMEL
do i=1,NSD
        if (X(i,n).lt.zmin) zmin = X(i,n)
        if (X(i,n).gt.zmax) zmax = X(i,n)
    enddo
enddo
else if (iopnflag.eq.1) then
zmin = dfmax
zmax = -dfmax
 do n=1,NUMEL
call LDISP (DL,D,IEN(1,n),NDOF,NEN)
 call ADDSNE (DL,n,NDOF,NEN,useDXE)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN, 182)
do i=1,NEN
k = 3*i
if (DL(k).gt.zmax) zmax = DL(k)
if (DL(k).lt.zmin) zmin = DL(k)
enddo
 enddo
endif
!
!    Expand boundaries
dz = (zmax-zmin)*opnaxnu
zmin = (zmin - dz)*opnscale
zmax = (zmax + dz)*opnscale
if (deformed) then
 zmin = zmin * deffac
 zmax = zmax * deffac
endif
if (iecho.eq.1) write(stderr,1) zmin,zmax
1    format(1x,'out-of-plane dimensions [',1PG12.3,',',1PG12.3,']')
!
!    Generate
dz = (zmax-zmin)/DBLE(NUMZ-1)
do n=1,NUMZ
 z(n) = zmin + DBLE(n-1)*dz
enddo
!
m = 0
do k=1,NUMZ
 do n=1,NUMNP
 m = m + 1
 write(lu,3,err=100,iostat=ios) m,0,(X(i,n),i=1,NSD), z(k)
3        format(2I5,1PG14.6,1PG14.6,1PG14.6)
 enddo
enddo
!
return
!
100    write(stderr,101) ios
  101    format(1x,'xwrite: write error ',I5)
ierr = EXEC
return
end
!-------------------------------------------------------------------------------
subroutine wrigen (IARRAY,N,M,FMT,IOFF,lu,ierr)

use iomodule,        only: stderr
use errors,          only: exec

implicit none

integer          :: N,M,i,j, ioff, ierr, ios,lu
integer          :: iarray(n,m)
character(len=*) FMT

write(lu,fmt=FMT,err=100,iostat=ios) &
    (j+IOFF,0,(int(IARRAY(i,j)),i=1,N),j=1,M)
return

100    write(stderr,101) ios
  101    format(1x,'wrigen: write error ',I5)
ierr = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine elmwrt(IEN,MAT,NUMEL,NUMNP,lu,IERR)

use modeltopology,   only: nen, nstr
use iomodule,        only: stderr
use plotControl,     only: numz
use errors,          only: exec

implicit none
!-pass
integer :: ien, mat
dimension IEN(NEN,*),MAT(*)
integer :: NUMEL,NUMNP
integer :: ierr, lu



!-local
integer :: k, m, n
integer :: node1, node2, node3, node4, node5, node6, node7, node8
integer :: ios

m = 0
do k=1,NUMZ-1
    do n=1,NUMEL
        m = m + 1
        node1 = IEN(1,n) + (k-1)*NUMNP
        node2 = IEN(2,n) + (k-1)*NUMNP
        node3 = IEN(3,n) + (k-1)*NUMNP
        node4 = IEN(4,n) + (k-1)*NUMNP
        node5 = IEN(1,n) + (k-0)*NUMNP
        node6 = IEN(2,n) + (k-0)*NUMNP
        node7 = IEN(3,n) + (k-0)*NUMNP
        node8 = IEN(4,n) + (k-0)*NUMNP

        write(lu,10,err=100,iostat=ios) m,MAT(n),node1,node2, &
            node3,node4,node5,node6,node7,node8
10        format(11I5)
    enddo
enddo

write(lu,20,err=100,iostat=ios)
20    format('end element definitions')
return

100    write(stderr,101) ios
  101    format(1x,'elmwrt: write error ',i4)

IERR = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine prewrt(STN0,ISELM,ISTIME,NPRE,lu,IERR)
!
use modeltopology,   only: nstr
use iomodule,        only: stderr
use errors,          only: exec

implicit none
!-pass
integer :: ISELM, ISTIME
double precision :: STN0
dimension STN0(NSTR,*),ISELM(*),ISTIME(*)
integer :: npre, lu, ierr


!-local
integer :: i, ios, j

write(lu,10,err=100,iostat=ios) &
    (ISELM(j),0,ISTIME(j),(STN0(i,j),i=1,NSTR),j=1,NPRE)
10    format(2I5,I6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6)
write(lu,20,err=100,iostat=ios)
20    format('end pre-stresses')
return
!
100    write(stderr,101) ios
  101    format(1x,'prewrt: write error ',i4)
IERR = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine preswt(IELNO,ISIDE,PRES,NUMPR,lu,IERR)

use iomodule,        only: stderr
use errors,          only: exec

implicit none
!-pass
integer          :: ielno, iside
double precision :: pres
dimension PRES(*),IELNO(*),ISIDE(*)
integer          :: lu, numpr, ierr


!-local
integer :: j, ios

write(lu,10,err=100,iostat=ios) &
    (IELNO(j),0,ISIDE(j),PRES(j),j=1,IABS(NUMPR))
10    format(3I5,1PG14.6)
write(lu,20,err=100,iostat=ios)
20    format('end pressure boundary conditions')
return
!
100    write(stderr,101) ios
  101    format(1x,'preswt: write error ',i4)
IERR = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine strswt(IELSTR,ISSIDE,STRS,ISTR,NUMSTR,lu,IERR)
!
use modeltopology,   only: nstr
use iomodule,        only: stderr
use errors,          only: exec

implicit none
!-pass
integer NUMSTR
integer :: ielstr, isside, istr
double precision :: strs
 dimension STRS(6,*),IELSTR(*),ISSIDE(*),ISTR(2,*)
integer :: ierr

integer :: lu

!-local
integer :: i, j, ios

! todo i=1,6   nstr hard coded?
write(lu,10,err=100,iostat=ios) &
    (IELSTR(j),0,ISSIDE(j),(ISTR(i,j),i=1,2),(STRS(i,j),i=1,6), &
    j=1,NUMSTR)
10    format(3I5,2I6,1PG14.6,1PG14.6,1PG14.6,1PG14.6,1PG14.6)
write(lu,20,err=100,iostat=ios)
20    format('end stress boundary conditions')
return
!
100    write(stderr,101) ios
  101    format(1x,'strswt: write error ',i4)
IERR = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine winkwt(IWELM,IWSIDE,IWTIME,WPRES,NUMWNK,lu,IERR)
!
use iomodule,        only: stderr
use errors,          only: exec

implicit none
!-pass
integer          :: iwelm, iwside, iwtime
double precision :: wpres
 dimension WPRES(*),IWELM(*),IWSIDE(*),IWTIME(*)
integer          :: lu, ierr, numwnk


!-local
integer :: j, ios

write(lu,10,err=100,iostat=ios) &
    (IWELM(j),0,IWSIDE(j),IWTIME(j),WPRES(j),j=1,ABS(NUMWNK))
10    format(3I5,I6,1PG14.6,1PG14.6)
write(lu,20,err=100,iostat=ios)
20    format('end Winkler pressures')
return
!
100    write(stderr,101) ios
  101    format(1x,'winkwt: write error ',i4)
IERR = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine slipwt(NSLIP,DIFORC,NUMSLP,NUMNP,NUMEL,NUMZ,lu, &
    ierr)

use modeltopology,   only: ndof
use iomodule,        only: stderr
use errors,          only: exec

implicit none
!-pass
integer :: nslip
double precision :: diforc
integer :: numslp, numnp, numel, numz, lu, ierr
 dimension NSLIP(5,*),DIFORC(NDOF,*)


!-local
integer :: ielm, node, i, ix, iy, iz, ios
integer :: n1, j, jelm, k

do j=1,NUMSLP
 ielm = NSLIP(1,j)
 node = NSLIP(2,j)
 ix   = NSLIP(3,j)
 iy   = NSLIP(4,j)
 iz   = NSLIP(5,j)
 do k=1,NUMZ-1
 jelm = ielm+(k-1)*NUMEL
 n1 = node + (k-1)*NUMNP

 write(lu,10,err=100,iostat=ios) jelm,n1,ix,iy,iz, &
   (DIFORC(i,j),i=1,NDOF)

 n1 = node + k*NUMNP

 write(lu,10,err=100,iostat=ios) jelm,n1,ix,iy,iz, &
   (DIFORC(i,j),i=1,NDOF)
10        format(5I5,1PG10.3,1PG10.3,1PG10.3)
 enddo
enddo
write(lu,20,err=100,iostat=ios)
20    format('end slippery node data'/'end differential Winkler data')
return
!
100    write(stderr,101) ios
  101    format(1x,'slipwt: write error ',i4)
IERR = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine FOUTRD (IREQ,IERR)

USE ALGEBRA
USE MODELDEFINITION
USE MODELDATAMODULE
USE MESHDATAMODULE
use TIMESTEPMODULE,  only: chktim,SetTimeStepSize,timdif,maxitime, &
                           timestepdatactx, nstep
use filemodule,      only: index_filename, &
                           BaseFedsk
use debugmodule,     only: iecho
use modeltopology,   only: nsd, ndof
use iomodule,        only: ascin, fedsk, stderr, fedsk, ifeout, luin, &
                           MOUT,TOUT
use postProcessData, only: setria,deformed
use errors,          only: exec, ok


! Routine to read FE output at requested time IREQ.
! ifeout = 0, (re-)opened output file - initialize pointers &
!     > 0, re-entering the routine.
! NSTEP = step that has been read in last reading operation

implicit none
!-pass
integer IREQ,IERR



!-locl
integer NSKIP
logical there

integer ncntr
integer, external :: lnblk, ifopen, ifbopen

if (ifeout.eq.0) then
    rewind(luin)
    NSTEP  = -1
endif
 
if (MODE.ne.2 .and. MODE.ne.5 .and. MODE.ne.6) then
    write(stderr,1) MODE
1   format(1x,'FoutRd: MODE = ',I1,' no data on output file')
    ierr = EXEC
    goto 100
endif

! determine the number of blocks that needs to be skipped
NSKIP = IREQ - NSTEP
if (NSKIP.eq.0 .and. ifeout.ne.0) return

call chktim (IREQ,MOUT,TOUT)

if ((.not.MOUT).and.(.not.TOUT)) then
    ierr = EXEC
    goto 100
endif

if (MAXITIME.ge.1) then
    ncntr = log10(dble(MAXITIME)) + 1
else
    ncntr = 1
endif
FEDSK = index_filename(BaseFedsk, maxitime, IREQ, .false.)
inquire(file=FEDSK, exist=there)
if (.not.there) then
    write(stderr, 3) '"', FEDSK(1:lnblk(FEDSK)), '"'
3   format(1x, 'FE-output file NOT found: ', 258a)
    goto 100
endif

if (ascin) then
    ierr = ifopen  (luin,FEDSK,'old')
    if (ierr.ne.0) then
        ierr = EXEC
        goto 100
    endif
else
    ierr = ifbopen (luin,FEDSK,'old')
    if (ierr.ne.0) then
        ierr = EXEC
        goto 100
    endif
endif

NSTEP = IREQ

call fread (ireq,ierr)

if (ierr.ne.OK) then
    ifeout = 0
    goto 100
endif

if (iecho.eq.1) then 
    write(stderr,2) NSTEP
2   format(1x,'FE-output file has been read successfully',1x, &
     '(itime=',i5,')')
endif

ifeout = 1

! set faulted node displacements and velocities
call setfaulteddispvelo (NSTEP,modeldatactx%NFAULT,modeldatactx%FAULT, &
  modeldatactx%TFAULT,modeldatactx%DFAULT, &
  timestepdatactx%MAXSTP,timestepdatactx%DELT, &
  modeldatactx%SKEW)

call extrme ()

if (deformed) then
    setria=.true.
endif

if (NSTEP.eq.0) then 
    call ADDSRF (modeldatactx%ISURF,modeldatactx%SURF, &
     meshdatactx%x, modeldatactx%d,NSD,NDOF,NSURF,0)
endif


100    return
end subroutine
!-------------------------------------------------------------------------------
subroutine fread (ireq, ierr)

USE MODELDEFINITION
USE ALGEBRA
USE MODELDATAMODULE
use timestepmodule,  only: chktim, &
                           SetTimeStepSize, &
                           nstep
use materialsmodule, only: matdatactx
use timestepmodule,  only: deltP
use plotcontrol,     only: irelt
use errors,          only: ok

!
! Routine to read a single fedsk file
!
implicit none



!-locl
logical            :: temprd, mechrd
integer            :: ireq, ierr
integer, parameter :: IPIDNT=1

ierr = OK
IREQ = NSTEP

call chktim (IREQ,mechrd,temprd)

if (temprd) then
    call EQUATE(modeldatactx%TLST,modeldatactx%T,NUMNP)
    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
        call EQUATE(modeldatactx%TLSTX,modeldatactx%TX,NUMNP)
    endif
endif

#ifdef EXPERIMENTAL_NORATESTORAGE
call femout (modeldatactx%d, &
             modeldatactx%STN, &
             modeldatactx%STR, &
             modeldatactx%DXE, &
               matdatactx%prpmat, &
               matdatactx%PRPLAS,&
            modeldatactx%T,modeldatactx%TX, modeldatactx%Thermalflux, &
    modeldatactx%SURF,modeldatactx%idsle,&
     modeldatactx%DELD, &
    modeldatactx%DX,IREQ,mechrd,temprd,ierr)


#else
!    read the output file
call femout (modeldatactx%d, &
     modeldatactx%STN, &
     modeldatactx%STR, &
    modeldatactx%DXE,matdatactx%prpmat, &
    matdatactx%PRPLAS,&
     modeldatactx%RATE,modeldatactx%T,modeldatactx%TX, modeldatactx%Thermalflux, &
    modeldatactx%SURF,modeldatactx%idsle,&
     modeldatactx%DELD, &
    modeldatactx%DX,IREQ,mechrd,temprd,ierr)
#endif

if (ierr.ne.OK) then
    write(*,*)  "Failed to read fedsk file"
    goto 100
endif

NSTEP = IREQ

call SetTimeStepSize ()

if (NSTEP.eq.irelt .and. temprd) then
    call EQUATE(modeldatactx%T0,modeldatactx%T, NUMNP)
endif

100    return
end subroutine
!-------------------------------------------------------------------------------
subroutine faultequalize (NFAULT,FAULT,NUMFN,NDOF)
!
! recomputes faulted displacement/velocity such that it is equally distributed
! on both sides of a fault. This is relevant for plotting displacements, which are
! relative to the normal/average node position
!
use iomodule,    only: stderr
use debugmodule, only: iecho

implicit none
!-pass
integer :: NFAULT(3,*),NUMFN,NDOF
double precision :: FAULT(NDOF,*)
!-locl
logical :: doecho
integer :: nf,node,idof,j
double precision :: minslp(3),maxslp(3),shift,slip,ffe

do nf=1,NUMFN

    ! already processed nodes are/will be labeled with a negative node number entry
    if (NFAULT(2,nf).gt.0) then

        ! find maximum and minimum faulted displacement for this node
        node = NFAULT(2,nf)
        do idof=1,NDOF
            minslp(idof) = FAULT(idof,nf)
            maxslp(idof) = FAULT(idof,nf)
        enddo
        NFAULT(2,nf) = -NFAULT(2,nf)
        do j=nf+1,NUMFN
            if (NFAULT(2,j).eq.node) then
                do idof=1,NDOF
                   if (FAULT(idof,j).lt.minslp(idof)) minslp(idof)=FAULT(idof,j)
                   if (FAULT(idof,j).gt.maxslp(idof)) maxslp(idof)=FAULT(idof,j)
                enddo
                NFAULT(2,j) = -NFAULT(2,j)
            endif
        enddo
        
!       test
        do idof=1,NDOF
            if (minslp(idof).ne.0d0 .and. maxslp(idof).ne.0d0) then
                if (minslp(idof).eq.maxslp(idof)) then
                    write(stderr,10) nf,node,idof
10                  format('Input error for NFAULT entry ',I6,':', &
                     'no differential slip for node ',I12,' dof=',I1,/ & 
                     'Probably, entries for one fault side are missing')
                    call exitp(1)
                endif
            endif
        enddo

        ! recompute all differential displacement entries for this node
        do j=nf,NUMFN
            if (IABS(NFAULT(2,j)).eq.node) then
                doecho = .false.
                do idof=1,NDOF
                    shift = ( maxslp(idof) + minslp(idof) ) * 0.5d0
                    slip  = ( maxslp(idof) - minslp(idof) ) * 0.5d0
                    doecho = (doecho .or. (shift.gt.0d0 .and. iecho.eq.2) )
                    ffe = ABS(FAULT(idof,j)-shift) ! = magnitude of expected new value
                    if (ABS(ffe-slip).gt.1d-25) then
                        write(stderr,20) j,node,idof
20                      format('Input error for NFAULT entry ',I6,':', &
                         'variable differential slip in FAULT for node ',I12,' dof=',I1)
                        call exitp(1)
                    endif
                    FAULT(idof,j) = SIGN(slip,FAULT(idof,j)-shift)
                enddo
                if (doecho) then
                    write(*,15) j,NFAULT(1,j),IABS(NFAULT(2,j)), &
                     (FAULT(idof,j),idof=1,NDOF)
15                  format('FAULT entry ',I6,' for element ',I12,' and node ',&
                     I12,': split redefined to',3(1x,1PE12.4))
                endif
            endif
        enddo

    endif

enddo

! make node numbers positive again
do nf=1,NUMFN
    NFAULT(2,nf) = IABS(NFAULT(2,nf))
enddo

return
end subroutine

!-------------------------------------------------------------------------------

subroutine setfaulteddispvelo (IREQ,NFAULT,FAULT,TFAULT,DFAULT,MAXSTP,DELT,SKEW)

! loads faulted node displacements and velocities at requested time step number IREQ
! TFAULT holds the total differential displacement of split nodes
! DFAULT hold the differential velocity between IREQ-1 and IREQ

use ALGEBRA,         only: EQUATE,VMPRD,FORMRT
use CONSTANTS,       only: ZERO,EPS
use MODELTOPOLOGY,   only: NDOF
use TIMESTEPMODULE,  only: NINTG,faultvelo
use MODELDEFINITION, only: NUMROT,NUMFN,MODE
implicit none
!-pass
integer :: IREQ,NFAULT,MAXSTP
double precision :: FAULT,TFAULT,DFAULT,DELT,SKEW
dimension NFAULT(3,*),FAULT(NDOF,*),TFAULT(NDOF,*),DFAULT(NDOF,*),MAXSTP(*),DELT(*),SKEW(2,*)
!-locl
logical :: VELBC,doit
integer :: iflt,idof,node,itime,igrp,grpindx
double precision :: rot(3,3),dl(3)

VELBC = (faultvelo.ne.0)

if ((MODE.ne.2 .and. MODE.ne.6) .or. NUMFN.eq.0) return

! initialize TFAULT with zero-time displacements, and  DFAULT=0 (large array,
! so don't do a wholesale zeroing)
do iflt=1,NUMFN
    do idof=1,NDOF
        TFAULT(idof,iflt) = ZERO
        DFAULT(idof,iflt) = ZERO
    enddo
    if (.not.VELBC) then
        if (NFAULT(3,iflt).le.0) then
            do idof=1,NDOF
                TFAULT(idof,iflt) = FAULT(idof,iflt)
            enddo
        endif
    endif
enddo

if (NINTG.eq.0 .or. IREQ.eq.0) then
    if (NUMROT.gt.0) then
        ! Euler rotations on TFAULT before we go
        do iflt=1,NUMFN
            node = NFAULT(2,iflt)
            if (ABS(SKEW(1,node)).gt.EPS .or. ABS(SKEW(2,node)).gt.EPS) then
                call FORMRT(SKEW(:,node),rot,NDOF)
                call EQUATE(dl,TFAULT(1:NDOF,iflt),NDOF)
                call VMPRD (rot,dl,TFAULT(1,iflt),NDOF,3)
            endif
        enddo
    endif
    return
endif

! step through timestep groups and increment TFAULT to the value at itime=IREQ
do iflt=1,NUMFN
    itime = 0
    igrp = 1
    loop_ireq: do while (itime.lt.IREQ)
        grpindx = 1 
        loop_grps: do while (grpindx.le.MAXSTP(igrp))
            itime = itime + 1
            if (VELBC) then
                do idof=1,NDOF
                    TFAULT(idof,iflt) = TFAULT(idof,iflt) + DELT(igrp)*FAULT(idof,iflt)
                enddo
            else
                if (NFAULT(3,iflt).ge.0) then
                    doit = (itime.eq.NFAULT(3,iflt))
                else
                    doit = (MOD(itime,-NFAULT(3,iflt)).eq.0)
                endif
                if (doit) then
                    do idof=1,NDOF
                        TFAULT(idof,iflt) = TFAULT(idof,iflt) + FAULT(idof,iflt)
                    enddo
                endif
            endif
            if (itime.eq.IREQ) then
                exit loop_ireq
            endif
            grpindx = grpindx + 1
        enddo loop_grps
        igrp = igrp + 1
    enddo loop_ireq
    if (VELBC) then
        do idof=1,NDOF
            DFAULT(idof,iflt) = FAULT(idof,iflt)
        enddo
    else
        ! there is no realistic velocity that can be attributed to incidental displacement
        ! therefore, only attribute a velocity for NFAULT(3,:) =  -1
        do idof=1,NDOF
            if (NFAULT(3,iflt).eq.-1 .and. DELT(igrp).gt.EPS) then
                DFAULT(idof,iflt) = FAULT(idof,iflt)/DELT(igrp)
            else
                DFAULT(idof,iflt) = ZERO
            endif
        enddo
    endif
    ! apply Euler rotations
    node = NFAULT(2,iflt)
    if (ABS(SKEW(1,node)).gt.EPS .or. ABS(SKEW(2,node)).gt.EPS) then
        call FORMRT(SKEW(:,node),rot,NDOF)
        call EQUATE(dl,TFAULT(1:NDOF,iflt),NDOF)
        call VMPRD (rot,dl,TFAULT(1,iflt),NDOF,3)
        call EQUATE(dl,DFAULT(1:NDOF,iflt),NDOF)
        call VMPRD (rot,dl,DFAULT(1,iflt),NDOF,3)
    endif
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine fcp (Ireq,iflag,ierr)

USE MODELDATAMODULE
USE MESHDATAMODULE
use materialsmodule, only: matdatactx
use modeldefinition, only: mode
use timestepmodule,  only: nmprt, ntprt, outputcontroldatactx
use iomodule,        only: lufo, stderr
use errors,          only: exec, ok

!
! Routine to copy a single data block (iflag=0) in original TECTON
! format or to write FE output in an alternative format (iflag=1).

implicit none
!-pass
integer :: ierr,Ireq
integer :: iflag



!-locl
logical temprd,mechrd,RDTOO
!-init
ierr = OK
!
!    determine if mechanical and thermal blocks need to be read
mechrd = RDTOO (outputcontroldatactx%imprint,IREQ,NMPRT, &
                MODE.eq.2.or.MODE.eq.6)
temprd = RDTOO (outputcontroldatactx%itprint,IREQ,NTPRT, &
                MODE.eq.5.or.MODE.eq.6)

!    read file
#ifdef EXPERIMENTAL_NORATESTORAGE
call femout (modeldatactx%d, &
             modeldatactx%stn,modeldatactx%str,modeldatactx%DXE,matdatactx%prpmat, &
             matdatactx%PRPLAS,modeldatactx%T,modeldatactx%TX,modeldatactx%Thermalflux, &
             modeldatactx%SURF,modeldatactx%idsle, &
             modeldatactx%DELD, & 
             modeldatactx%dx,Ireq,mechrd,temprd,ierr)
#else
call femout (modeldatactx%d, &
             modeldatactx%stn,modeldatactx%str,modeldatactx%DXE,matdatactx%prpmat, &
             matdatactx%PRPLAS,modeldatactx%rate,modeldatactx%T,modeldatactx%TX,modeldatactx%Thermalflux, &
             modeldatactx%SURF,modeldatactx%idsle, &
             modeldatactx%DELD, &
             modeldatactx%dx,Ireq,mechrd,temprd,ierr)
#endif

if (ierr.ne.OK) return
!
!    write file
 if (iflag.eq.0) then
!       Write in TECTON output format

#ifdef EXPERIMENTAL_NORATESTORAGE
call fewrite (modeldatactx%d, &
modeldatactx%stn,modeldatactx%str,modeldatactx%DXE,modeldatactx%idsle, &
  matdatactx%prpmat,meshdatactx%ien, &
  matdatactx%PRPLAS,modeldatactx%SURF,modeldatactx%T,modeldatactx%Thermalflux, &
  modeldatactx%TX, &
  Ireq,mechrd,temprd,lufo,ierr)
#else
call fewrite (modeldatactx%d, &
modeldatactx%stn,modeldatactx%str,modeldatactx%DXE,modeldatactx%idsle, &
  matdatactx%prpmat,meshdatactx%ien, &
  matdatactx%PRPLAS,modeldatactx%rate,modeldatactx%SURF,modeldatactx%T,modeldatactx%Thermalflux, &
  modeldatactx%TX, &
  Ireq,mechrd,temprd,lufo,ierr)
#endif

 elseif (iflag.eq.1) then
   call fewritalt(meshdatactx%X, &
               modeldatactx%d, &
               modeldatactx%stn,meshdatactx%ien,Ireq, &
    mechrd,temprd,lufo,ierr)
 else
    write(stderr,1) iflag
1      format(1x,'fcp: IFLAG= ',I1,' (must be 0 or 1)')
    ierr = EXEC
    return
 endif
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine fskip (ireq,NSKIP,ierr)

! Routine to skip NSKIP data blocks

use errors,          only: ok


implicit none

!-pass
integer :: ireq,NSKIP,ierr
!-locl
integer :: i

if (NSKIP.le.0) then
    return
endif

do i=1,NSKIP
     call fread (ireq, ierr)
    if (ierr.ne.OK) then
        return
    endif
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine extrme ()

use constants,       only: dfmax
USE MODELDEFINITION
USE MODELDATAMODULE
use modeltopology,   only: nstr, nsd, ndof
use iomodule,        only: MOUT,TOUT
use timestepmodule,  only: nstep
use plotcontrol,     only: strsmn, strsmx, &
                           dispmn, dispmx, &
                           ratemn, ratemx, &
                           fluxmn, fluxmx, &
                           tempmn, tempmx
!use modeltopoly,    only: nstr
!  determines extrema in arrays

implicit none

!-local
#ifdef EXPERIMENTAL_NORATESTORAGE
integer :: iElem
double precision :: strainRate(nstr)
#endif

if (Mout) then
    call vecxtr (modeldatactx%d, &
                 NDOF,NUMNP,dispmn,dispmx,0)
    call vecxtr (modeldatactx%stn,NSTR,NUMEL,strsmn,strsmx,1)
    if (NSTEP.gt.0) then

#ifdef EXPERIMENTAL_NORATESTORAGE
!            ratemn = dfmax
!            ratemx = -dfmax
!            do iElem=1,meshdatactx%neglobal
!                call StrainrateOfElem (meshdatactx%X,&
!                              modeldatactx%D,&
!                              modeldatactx%DELD, &
!                               strainRate,&
!                               meshdatactx%IEN,&
!                              modeldatactx%LMF, &
!                              modeldatactx%TFAULT,&
!                              modeldatactx%DFAULT,&
!                              modeldatactx%LMX, &
!                              modeldatactx%DX, &
!                              modeldatactx%SKEW,&
!                              modeldatactx%DXE, &
!                               iElem)

!                do iStrain =1,NSTR
!                    if (strainRate(iStrain).lt.ratemn) then
!                        ratemn = strainRate(iStrain)
!                    endif
!                    if (strainRate(iStrain).gt.ratemx) then
!                        ratemx = strainRate(iStrain)
!                    endif
!                enddo

!            enddo
#else
            call vecxtr (modeldatactx%rate,NSTR,NUMEL,ratemn,ratemx,0)
#endif

    endif
endif

if (Tout) then
    call sclxtr (modeldatactx%T,NUMNP,tempmn,tempmx)
    call vecxtr (modeldatactx%Thermalflux, &
         NSD,NUMEL,fluxmn,fluxmx,0)
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine vecxtr (A,N,M,AMIN,AMAX,IDIF)

use constants, only: zero, dfmax

! Routine for determining bounds on vector data in array A.
!         AMIN(1)      <= positive As   <=   AMAX(1)
!         AMIN(2)      <= negative As   <=   AMAX(2)

implicit none
!-pass
integer N,M,IDIF
double precision :: A(N,M),AMAX(2),AMIN(2)
!-locl
logical pos,neg
double precision :: B(6)
integer i, j
!-init

AMIN(1) =  dfmax
AMAX(1) =  ZERO

AMAX(2) = -dfmax
AMIN(2) =  ZERO

pos = .false.
neg = .false.

do j=1,M
 do i=1,N
B(i) = A(i,j)
 enddo

 do i=1,N
if (B(i).ge.ZERO) then
    pos = .true.
    if (B(i).gt.AMAX(1)) then
        AMAX(1) = B(i)
    endif
    if (B(i).lt.AMIN(1)) then
        AMIN(1) = B(i)
    endif
else
    neg = .true.
    if (B(i).gt.AMAX(2)) then
        AMAX(2) = B(i)
    endif
    if (B(i).lt.AMIN(2)) then
        AMIN(2) = B(i)
    endif
endif
 enddo
enddo

if (.not.pos) then
 AMIN(1) = ZERO
endif

if (.not.neg) then
 AMAX(2) = ZERO
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine sclxtr (A,N,AMIN,AMAX)

use constants, only: dfmax

! Routine for determining bounds on scalar data in array A.

implicit none

!-pass
integer          :: N
integer          :: i
double precision :: A(N)
double precision :: amin, amax

!-init
AMIN =  dfmax
AMAX = -dfmax

do i=1,N
 if (A(i).gt.AMAX) then
AMAX = A(i)
 endif
 if (A(i).lt.AMIN) then
AMIN = A(i)
 endif
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine RdCore (T,FLUX,PRPMAT,PRPLAS,D,DX,STN, &
                   lu,lue,ierr)

! routine to read FE core file

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact
use errors,          only: exec


implicit none
!-pass
integer lu,lue,ierr
double precision :: t, flux, prpmat, prplas, d, dx, stn
dimension T(*),FLUX(NSD,*),PRPMAT(6,*),PRPLAS(9,*),D(NDOF,*), &
    DX(NDOF,*),STN(NSTR,*)
!-locl
integer, parameter :: blockSize = 51
character(len=1)c
integer n,nn,m,i

read (lu,2,err=1000,end=2000) NSTEP
2    format(///13X,I6)

read (lu,4,err=1000,end=2000) c
4    format(/A)

if (MODE.eq.5 .or. MODE.eq.6) then
!        Temperatures block
nn = blockSize

    do n=1,NUMNP
        if (nn.gt.50) then
            read(lu,6,err=1000,end=2000) c
6            format(//////A)
            nn = 1
        endif
        read(lu,8,err=1000,end=2000) m,T(n)
8        format(5X,I6,17X,E15.0)
        if (n.ne.m) then
            write(lue,10) m,n
10               format(1x,'RdCore: synchronization error'/ &
             9X,'temperatures block, read node=',I5, &
             3X,'expected node=',I5)
            ierr = EXEC
            return
        endif
        nn = nn + 1
    enddo
endif

if (MODE.eq.5 .or. MODE.eq.6) then
!        Heat fluxes block
    nn = blockSize
    do n=1,NUMEL
        if (nn.gt.50) then
            read(lu,12,err=1000,end=2000) c
12            format(//////A)
            nn = 1
        endif
        read(lu,14,err=1000,end=2000) m,(FLUX(i,n),i=1,NSD)
14        format(1X,I6,8X,3(E15.0,4X))
        if (n.ne.m) then
            write(lue,16) m,n
16               format(1x,'RdCore: synchronization error'/ &
             9X,'heat fluxes block, read element=',I5, &
             3X,'expected element=',I5)
            ierr = EXEC
            return
         endif
         nn = nn + 1
    enddo
endif

if (MODE.eq.5 .or. MODE.eq.6) then
!        Creep properties block
    nn = blockSize
    do n=1,NUMAT
        if (nn.gt.50) then
            read(lu,18,err=1000,end=2000) c
18            format(///////A)
            nn = 0
        endif
        read(lu,20,err=1000,end=2000) m,(PRPMAT(i,n),i=3,4)
20        format(I5,5X,E10.0,5X,F4.0)
        if (n.ne.m) then
            write(lue,22) m,n
22               format(1x,'RdCore: synchronization error'/ &
             9X,'creep properties block, read material=',I5, &
             3X,'expected material=',I5)
            ierr = EXEC
            return
        endif
        nn = nn + 1
    enddo
endif

!    Plastic properties block
nn = blockSize
do n=1,NUMAT
    if (nn.gt.50) then
        read(lu,24,err=1000,end=2000) c
24        format(///////A)
        nn = 0
    endif
    read(lu,26,err=1000,end=2000) m,(PRPLAS(i,n),i=3,4)
26        format(1X,I5,5X,E10.0,4X,E10.0)
    if (n.ne.m) then
        write(lue,28) m,n
28           format(1x,'RdCore: synchronization error'/ &
         9X,'plastic properties block, read material=',I5, &
         3X,'expected material=',I5)
         ierr = EXEC
         return
     endif
    nn = nn + 1
enddo

if (MODE.eq.2 .or. MODE.eq.6) then
!        Displacements block
    nn = blockSize
    do n=1,NUMNP
        if (nn.gt.50) then
            read(lu,30,err=1000,end=2000) c
30            format(////////A)
            nn = 1
        endif
        read(lu,32,err=1000,end=2000) m,(D(i,n),i=1,NDOF)
32        format(6X,I5,10X,6E18.0)
        if (n.ne.m) then
            write(lue,34) m,n
34               format(1x,'RdCore: synchronization error'/ &
             9X,'displacements block, read node=',I5, &
             3X,'expected node=',I5)
            ierr = EXEC
            return
        endif
        nn = nn + 1
    enddo
endif

if ((MODE.eq.2 .or. MODE.eq.6).and.NUMSLP.gt.0) then
!        Differential displacements block
    read(lu,36,err=1000,end=2000) c
36        format(///////A)
    do n=1,NUMSN
        read(lu,38,err=1000,end=2000) m,(DX(i,n),i=1,NDOF)
38        format(I5,1X,3E15.0)
    enddo
endif

if (MODE.eq.2 .or. MODE.eq.6) then
!        Stresses block
    nn = blockSize
    do n=1,NUMEL
        if (nn.gt.50) then
            read(lu,40,err=1000,end=2000) c
40            format(////////A)
            nn = 1
        endif
        read(lu,42,err=1000,end=2000) m,(STN(i,n),i=1,NSTR)
42        format(1X,I5,1X,6E14.0)
        if (n.ne.m) then
            write(lue,44) m,n
44               format(1x,'RdCore: synchronization error'/ &
            9X,'stresses block, read element=',I5, &
            3x,'expected element=',I5)
            ierr = EXEC
            return
        endif
        nn = nn + 1
    enddo
endif

return

1000    call prerr('RdCore')
ierr = EXEC
return
2000    write(lue,2001)
2001    format(1x,'RdCore: EOF error')
ierr = EXEC
return
end subroutine
!-----------------------------------------------------------------------
subroutine tecout (T,lufo,luMon)

USE MODELDEFINITION, only: numnp
use monitorVars,     only: npmon, nnp, monnp

! writes "tecton" output file (for initial temperatures)

implicit none
!-pass
integer lufo,luMon
double precision :: T
dimension T(*)
!-local
integer :: n, nn

1    format(2I5,1PG14.6)
2    format('end initial temperatures')

if (npmon) then
    do nn=1,nnp
        n = monnp(nn)
        write(luMon,1) n,0,T(n)
    enddo
    write(luMon,2)
else
    do n=1,NUMNP
        write(lufo,1) n,0,T(n)
    enddo
    write(lufo,2)
endif

return
end subroutine
!-----------------------------------------------------------------------
subroutine SVSRF (ISURF,SURF,SURF1,X,D,NSURF)

use modeltopology, only: ndof, nsd

! Saves SURF1 = SURF - X - D

implicit none
!-pass
integer NSURF, isurf
double precision :: x, d, surf, surf1
dimension ISURF(*),SURF(NSD,*),SURF1(NSD,*),X(NSD,*),D(NDOF,*)
!-local
integer :: n, i, j

if (NSURF.le.0) then
    return
endif

do j=1,NSURF

    n = ISURF(j)

    do i=1,NSD
        SURF1(i,j) = SURF(i,j) - X(i,n)
    enddo

    do i=1,NSD
        SURF1(i,j) = SURF1(i,j) - D(i,n)
    enddo

enddo

return
end subroutine
!------------------------------------------------------------------------------
subroutine RdXinp (FILE,XINP,NSD,NUMNP,ierr)

use debugmodule,     only: iecho
use iomodule,        only: stderr, luin
use errors,          only: ok, syntax


implicit none
!-pass
integer NSD,NUMNP
double precision :: xinp
dimension XINP(NSD,*)
character(len=*) FILE
integer :: ierr


integer, external :: lnblk, ifopen, nextlu
!-locl
character(len=12)fmt
integer :: lusave, nread

!
ierr = OK

if (NUMNP.le.0) then
    return
endif

lusave = luin
luin = nextlu(0)
ierr = ifopen (luin,FILE,'old')

if (ierr.ne.0) then
 ierr = SYNTAX
 goto 100
endif

fmt = '(2I12,3G14.0)'

call FGEN (XINP,NSD,NUMNP,fmt,nread)
call closef (luin)
if (nread.ne.NUMNP) then
 write(stderr,1) nread,NUMNP
1       format(1x,'"read": read ',I5,' nodal point coordinates,',1x, &
 'expected ',I5)
ierr = SYNTAX
 else
 if (iecho.eq.1) write(stderr,2) 'EOF "',FILE(1:lnblk(FILE)), &
   '"'
2        format(1x,80A)
endif
100    luin = lusave
return
end subroutine
!-----------------------------------------------------------------------
subroutine CURFOR (MAT,IEN,X,LMF,TFAULT, &
                   PRPMAT,STN,D,SKEW,DXE,monside,pcur)
USE ALGEBRA
use materialsmodule, only: MaterialIDOfElement
use modeltopology,   only: nsd, ndof, nstr, nen
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel
use errors,          only: exec


! Routine for calculating current force level PCUR from stress deviator

implicit none
!-pass

integer :: mat, ien, lmf
double precision :: x, prpmat, stn, pcur, d, tfault, skew, dxe
dimension MAT(*),IEN(NEN,*),X(NSD,*),PRPMAT(6,*), &
    STN(NSTR,*),PCUR(NDOF),D(NDOF,*),LMF(NDOF,NEN,*), &
    TFAULT(NDOF,*),SKEW(2,*), DXE(NDOF,NEN,*)
!-locl
double precision :: PL(NDOF)
integer          :: i,n,nn,m,iside, monside
double precision :: thick
integer          :: ierr

if (.not.elmon) then
    write(stderr,10)
10        format(1x,'force: no monitor elements')
    ierr = EXEC
    return
endif

call CLEAR (PCUR,NDOF,"PCUR")

do nn=1,nnel
    n     = monel(nn)
    iside = monside
    m = MaterialIDOfElement(n)
    THICK = PRPMAT(6,m)

!   Calculate force working on this side
    call ELFORC (IEN(1,n),X,D,LMF(1,1,n),TFAULT, &
                 STN(1,n),SKEW, n,THICK,ISIDE,PL,ierr)

    if (ierr.ne.0) then
        write(stderr,1) n
1        format(1x,'curfor: error in element ',I5)
        call stoper()
    endif

    do i=1,NDOF
        PCUR(i) = PCUR(i) + PL(i)
    enddo

enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine ELFORC (IEN,X,D,LMF,TFAULT, &
                   STN,SKEW,elemID,THICK,NSIDE,PCUR,IERR)

USE ALGEBRA
use constants,       only: useDXE
use modeltopology,   only: nsd, ndof, nstr, nen
USE MODELDEFINITION, only: lgdef
!
!    Calculates force PCUR working on element side NSIDE
!
implicit none

integer, parameter :: NEE=8 ! nen * ndof
!-pass
integer IERR

integer :: ien,lmf
double precision :: x,d,stn,pcur,tfault, skew
dimension IEN(NEN),X(NSD,*),D(NDOF,*),STN(NSTR),PCUR(NDOF), &
    LMF(NDOF,NEN),TFAULT(NDOF,*), SKEW(2,*)
integer :: elemID, nside
double precision :: thick
!-locl
logical QUAD
double precision :: XL(NSD,NEN),PL(NDOF,NEN),DL(NDOF,NEN),SDEV(NSTR)
integer          :: i,j,k
!
call CLEAR (PCUR,NDOF,"PCUR")
!
QUAD = (IEN(3).ne.IEN(4))

i = NSIDE
j = i + 1

if (QUAD) then
    if (j.gt.NEN) then
        j = 1
    endif
else
    if (j.ge.NEN) then
        j = 1
    endif
endif

call LCOORD (X,XL,IEN)

if (LGDEF.ne.0) then
    call LDISP (DL,D,IEN,NDOF,NEN)
    ! TODO, add proper element ID. Should probably be 1
    call AddFaultDisplacement (DL,1,TFAULT,NDOF,NEN,183)
    call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
    call REZONE (XL,DL,1d0)
endif

call CLEAR (PL,NEE,"PL")
! use SDEV or STN ?
!    call DEVSTN (STN,SDEV)
 call EFORCE (STN,PL,XL,THICK,QUAD,ierr)
if (ierr.ne.0) then
    return
endif

do k=1,NDOF
    PCUR(k) = PCUR(k) + PL(k,i) + PL(k,j)
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine DETVELO (D,DELD,DXE,DELDXE,IMPRINT,ITPRINT,MAXSTP, &
                    DELT)

!TODO, remove imprint etc from call arguments

! Routine for determining velocities from output (IVELOUT=0)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDEFINITION, only: numel, numnp, numsn, nElemsWithSlipperyNodes
#else
USE MODELDEFINITION, only: numel, numnp, numsn
#endif
use timestepmodule,  only: timdif
use modeltopology,   only: NDOF,NEN
use constants,       only: zero

implicit none

!-pass

integer          :: imprint, itprint, maxstp
double precision :: d,deld,dxe,deldxe,delt

dimension D(NDOF,*),DELD(NDOF,*), &
    DXE(NDOF,NEN,*),DELDXE(NDOF,NEN,*), &
    IMPRINT(*),ITPRINT(*),MAXSTP(*),DELT(*)

!-locl
double precision :: dtime, tdummy
integer          :: i,j,n
!

!       get time since previous output (in DELD and DELDXE)
dtime=1d0
call timdif (dtime,tdummy)

if (dtime.le.ZERO) then
    return
endif

do n=1,NUMNP
    do i=1,NDOF
           DELD(i,n) = (D(i,n)-DELD(i,n))/dtime
    enddo
enddo

if (NUMSN.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
    do n=1,nElemsWithSlipperyNodes
        do j=1,NEN
            do i=1,NDOF
                DELDXE(i,j,n) = (DXE(i,j,n)-DELDXE(i,j,n))/dtime
            enddo
        enddo
    enddo
#else
    do n=1,NUMEL
        do j=1,NEN
            do i=1,NDOF
                DELDXE(i,j,n) = (DXE(i,j,n)-DELDXE(i,j,n))/dtime
            enddo
        enddo
    enddo
#endif
endif


return
end subroutine
