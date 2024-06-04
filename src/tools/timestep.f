	program main
*
	implicit none
*
	integer MAXGRP,MAXPRT
	parameter (MAXGRP=500, MAXPRT=500)
	double precision YEAR
	parameter (YEAR=3.1556736D+07)
c-locl
	integer nout,i,igrp,nstep,NINT,istep,nprt,IMPRINT,MAXSTP,
     >	 lu,NEXTLU,nmax,INMAX,MAXINT
	character UNIT*5
	double precision E,v,eta,mhu,maxwell,dt0,time,stepsize,tmax,out,
     >	 DELT,effe
	dimension IMPRINT(MAXPRT),MAXSTP(MAXGRP),DELT(MAXGRP),
     >	 UNIT(MAXGRP)
	external NINT,NEXTLU,INMAX
*
	MAXINT = INMAX()
*
	write(*,1)
    1	format(
     >   '     ************************************************'/,
     >   '     * Generate time step scheme for stress flow,   *'/,
     >   '     * based on (approximately) uniform stress drop *'/,
     >   '     * during each time step.                       *'/,
     >   '     ************************************************'//
     >   ' Young''s modulus (Pa) > ',$)
	read(*,*,end=100,err=100) E
	write(*,2)
    2	format(1x,'Poisson''s ratio > ',$)
	read(*,*,end=100,err=100) v
	write(*,3)
    3	format(1x,'Viscosity (Pa s) > ',$)
	read(*,*,end=100,err=100) eta
*
	mhu = E/(2d0*(1d0+v))
	maxwell = eta/mhu
	write(*,4) maxwell,maxwell/year
    4	format(1x,'Maxwell time = ',1PG14.6,' sec (',
     >   1PG14.6,' year)')
*
c	default time step size at oscillation limit
	dt0 = eta/E*2d0*(1-v)*(1-2d0*v)
	write(*,5) dt0
    5	format(1x,'Initial time step size (',1PG8.2,' sec) > ',$)
	read(*,'(E20.0)',end=100,err=100) effe
	if (ABS(effe).gt.1e-6) dt0 = effe
	write(*,6)
    6	format(1x,'Total integration time (number of Maxwell ',
     >   'times) > ',$)
	read(*,*,end=100,err=100) nmax
	write(*,7)
    7	format(1x,'Output frequency (fraction of Maxwell time) > ',$)
	read(*,*,err=100,end=100) out
	write(*,'(1X)')
*
	time = 0d0				! total integration time
	nout = 0				! total number of time steps
	igrp = 1				! time step group number
	nprt = 1				! number of outputs
50	tmax = maxwell * LOG(DBLE(igrp+1))	! integration time contribution by current time step group
	stepsize = DBLE(igrp)*dt0		! time step size in current time step group
	nstep = MAX(1,NINT(tmax/stepsize))	! number of time steps in current time step group
	if (time+nstep*stepsize.gt.nmax*maxwell)
     >   nstep = MAX(1,NINT((nmax*maxwell-time)/stepsize))
	if (time+nstep*stepsize.gt.nmax*maxwell)
     >   stepsize = (nmax*maxwell-time)/DBLE(nstep)
	if (igrp.gt.MAXGRP) stop 'MAXGRP dimensioned too small'
	MAXSTP(igrp) = nstep
	DELT(igrp)   = stepsize
	write(*,8) time/maxwell,(time+nstep*stepsize)/maxwell,nstep,
     >   stepsize/dt0
    8	format(1x,'From time = ',F7.3,' to ',F7.3,
     >   ' Maxwell times: ',I5,' steps of ',F7.3,
     >   ' initial stepsizes')
	istep = NINT((nprt*out*maxwell-time)/stepsize)
	do while (istep.ge.1 .and. istep.le.nstep)
	    write(*,9) nout+istep,(time+istep*stepsize)/maxwell
    9	    format(1x,'Output at istep = ',I5,' @ ',F7.3,
     >       ' Maxwell times')
	    if (nprt.gt.MAXPRT) stop 'MAXPRT dimensioned too small'
	    IMPRINT(nprt) = nout+istep
	    nprt = nprt + 1
	    istep = NINT((nprt*out*maxwell-time)/stepsize)
	enddo
	nout = nout + nstep
	time = time + DBLE(nstep)*stepsize
	igrp = igrp + 1
	if (time.lt.nmax*maxwell) goto 50
*
	nprt = nprt - 1
	igrp = igrp - 1
	write(*,10) nout,igrp
   10	format(1x,I5,' time steps, ',I5,' time step groups')
	write(*,11) nprt
   11	format(1x,I5,' outputs')
*
	lu = nextlu(0)
	call openf(lu,'tecin.dat.tstp','unknown')
	write(lu,'(50i5)') (MAXSTP(i),i=1,igrp)
	do i=1,igrp
	    stepsize = DELT(i)/YEAR
C	    write(*,*) 'ORG: ',stepsize,' year, ',DELT(i),' s'
c	    try to fit into integer format (works up to inmax)
	    if (DELT(i).le.MAXINT .and. NINT(DELT(i)).le.99999) then
		write(lu,'(I5,$)') NINT(DELT(i))
		UNIT(i) = '  sec'
C	write(*,*) NINT(DELT(i)),UNIT(i) 
	    elseif (stepsize.lt.1.0) then
		write(lu,'(F5.4,$)') stepsize
		UNIT(i) = ' year'
C	write(*,*) stepsize,UNIT(i)
	    else if (stepsize.lt.10.0) then
		write(lu,'(F5.3,$)') stepsize
		UNIT(i) = ' year'
C	write(*,*) stepsize,UNIT(i)
	    else if (stepsize.lt.100.0) then
		write(lu,'(F5.2,$)') stepsize
		UNIT(i) = ' year'
C	write(*,*) stepsize,UNIT(i)
	    else if (stepsize.lt.1000.0) then
		write(lu,'(F5.1,$)') stepsize
		UNIT(i) = ' year'
C	write(*,*) stepsize,UNIT(i)
	    else if (stepsize.lt.1.0D5) then
		write(lu,'(I5,$)') NINT(stepsize)
		UNIT(i) = ' year'
C	write(*,*) NINT(stepsize),UNIT(i)
	    else if (stepsize.lt.1.0D6) then
		write(lu,'(F5.4,$)') stepsize*1D-6
		UNIT(i) = '   Ma'
C	write(*,*) stepsize*1D-6,UNIT(i)
	    else if (stepsize.ge.1.0D6) then
	        stop 'error'
	    endif
	    if (MOD(i,50).eq.0) write(lu,'(1X)')
	enddo
	if (MOD(igrp,50).ne.0) write(lu,'(1X)')
	write(lu,'(50A5)') (UNIT(i),i=1,igrp)
	write(lu,'(50F5.1)') (0.5,i=1,igrp)
	call closef(lu)
	write(*,12)
   12	format(1x,'Written time step groups to ',
     >   '"tecin.dat.tstp"')
	call openf(lu,'tecin.dat.mprt','unknown')
	write(lu,'(51i5)') (IMPRINT(i),i=1,nprt)
	call closef(lu)
	write(*,13)
   13	format(1x,'Written output times to ',
     >   '"tecin.dat.mprt"')
*
100	end
