!c The program computes the analytical solution of a bending clamped plate
!c The program requires a file argument (f.e. input.dat) which look like
!c
!c    # used option for flexure
!c    2
!c    # inputdata for flexure, name of the output file
!c    output.dat
!c    # variables for the model
!c    # Young Modulus (in Pa)
!c    7.5E10
!c    # Poisson ratio
!c    0.30
!c    # model thickness (in km)
!c    1
!c    # length (in km)
!c    500
!c    # Applied force (in N/m)
!c    1E5
!c    # xmin (in km)
!c    0
!c    # xmax (in km)
!c    500
!c    # dx (in km)
!c    10
 


        program flexure
        implicit none
        character :: ch
        character (len=15)::filename,filename1

        call getarg(1,filename1)

   	    open(unit = 11,file=filename1)
        read (11,*);read(11,*) ch
        read(11,*);read(11,'(a15)') filename

        if (ch=='h') then
                call help
        elseif (ch=='1') then
                call uno(filename)
        elseif (ch=='2') then
                call due(filename)
        elseif (ch=='3') then
                call tre(filename) 
        elseif (ch=='4') then
                call qua(filename)
        elseif (ch=='5') then
                call cin(filename)
        else
                print *, 'OPTION NOT RECOGNIZED!!!!'
                call help
        endif

        end program        

!----------------------------------------------------------------
        subroutine help
        print *,'Use: flexure option filename'
        print *,'filename: output filename'
        print *,'h this file'
        print *,'1 Embedded plate at one end subject to an apply torque at the other'
        print *,'2 Embedded plate at one end and subject to concentrated linear force at the other end'
        print *,'3 Infinite slab with line load @ 0 and density contrast'
        print *,'4 As 3 but for broken slab'
        print *,'5 Half infinite slab with Moment and linear force applied at 0'
        end subroutine help
        
!----------------------------------------------------------------
        subroutine slabinput
        implicit none
        common /slab/ E,nu,h,L,D
        real::E,nu,h,L,D
!c        print *,'Input the values for the slab:'
!c        print *,'Input Young modulus (Pa)  :';read *,E
!c        print *,'Input Poisson ratio       :';read *,nu
!c        print *,'Input Plate thickness (km):';read *,h
!c        print *,'Input Slab lenght (km)    :';read *,L
	read(11,*);read(11,*);read(11,*) E
	read(11,*);read(11,*) nu
	read(11,*);read(11,*) h
	read(11,*);read(11,*) L

        h=h*1e3;L=L*1e3
        D=E*h*h*h/(12.*(1.-nu*nu))

        end subroutine slabinput

!----------------------------------------------------------------
        subroutine uno (filename)

        implicit none

        common /slab/ E,nu,h,L,D
        real::E,nu,h,L,D
        real::M,xmin,xmax,x,w,dx,D1
        integer :: xi
        character (len=15)::filename

        call slabinput

        D1=1./D
!c        print *;print *,'Input the applied torque (Nm)';read *,M
!c        print *;print *,'Input xmin (km)';read *,xmin;xmin=xmin*1e3
!c        print *,'Input xmax (km)';read *,xmax;xmax=xmax*1e3
!c        print *,'Input dx (km)  ';read *,dx;dx=dx*1e3
        read(11,*);read(11,*) M
        read(11,*);read(11,*) xmin;xmin=xmin*1e3
        read(11,*);read(11,*) xmax;xmax=xmax*1e3
        read(11,*);read(11,*) dx;dx=dx*1e3
        open (21,file=filename)
        write(21,'("Slab embedded at one side and subject to a torque at the other")')

1       format ("E=",e10.3,"Pa nu=",e10.3," h=",e10.3,"km L=",e10.3,"km")
        write(21,1) E,nu,h*1e-3,L*1e-3

2       format ("M= ",e10.3,"Nm")
        write(21,2) M


!        do x = xmin,xmax,dx
        do xi = 0,floor((xmax-xmin)/dx)
            x = xmin + dble(xi) * dx
						
            w=-M*x*x*0.5*D1
            write(21,'(2e18.8)') x,w
        enddo

        close(21)
        end subroutine uno

!----------------------------------------------------------------
        subroutine due (filename)
        implicit none
        common /slab/ E,nu,h,L,D
        real::E,nu,h,L,D
        real::V,xmin,xmax,x,w,dx,D1
        integer :: xi
        character (len=15)::filename
        call slabinput
        D1=1./D
!c        print *;print *,'Input the applied force (N/m)';read *,V
!c        print *;print *,'Input xmin (km)';read *,xmin;xmin=xmin*1e3
!c        print *,'Input xmax (km)';read *,xmax;xmax=xmax*1e3
!c        print *,'Input dx (km)  ';read *,dx;dx=dx*1e3
	read(11,*);read(11,*) V
	read(11,*);read(11,*) xmin;xmin=xmin*1e3
	read(11,*);read(11,*) xmax;xmax=xmax*1e3
	read(11,*);read(11,*) dx;dx=dx*1e3
        open (21,file=filename)
        write(21,'("Slab embedded at one side and subject to a linear force at the other")')
1       format ("E=",e10.3,"Pa nu=",e10.3," h=",e10.3,"km L=",e10.3,"km")
        write(21,1) E,nu,h*1e-3,L*1e-3
2       format ("V= ",e10.3,"N/m")
        write (21,2) V
!        do x = xmin,xmax,dx
        do xi = 0,floor((xmax-xmin)/dx)           
            x = xmin + dble(xi) * dx 
         w=-V*x*x*0.5*D1*(L-x/3.)
         write(21,'(2e18.8)') x,w
        enddo
        close(21)
        end subroutine due 

!----------------------------------------------------------------
        subroutine tre (filename)
        implicit none
        common /slab/ E,nu,h,L,D
        integer :: xi
        real::E,nu,h,L,D
        real::V,xmin,xmax,x,w,dx,rw,rm,g,alpha,w0,tmp
        character (len=15)::filename
        call slabinput
!c        print *;print *,'Input the applied force (N/m)';read *,V
!c        print *,'Input density material above (kg/m3)';read *,rw
!c        print *,'Input density material below (kg/m3)';read *,rm
!c        print *;print *,'Input xmin (km)';read *,xmin;xmin=xmin*1e3
!c        print *;print *,'Input xmax (km)';read *,xmax;xmax=xmax*1e3
!c        print *,'Input dx (km)  ';read *,dx;dx=dx*1e3
        read(11,*);read(11,*) V
        read(11,*);read(11,*) rw
        read(11,*);read(11,*) rm
        read(11,*);read(11,*) xmin;xmin=xmin*1e3
        read(11,*);read(11,*) xmax;xmax=xmax*1e3
        read(11,*);read(11,*) dx;dx=dx*1e3
        open (21,file=filename)
        write(21,'("Infinite slab with line load @ 0 and restoring forces")')
1       format ("E=",e10.3,"Pa nu=",e10.3," h=",e10.3,"km L=",e10.3,"km")
        write(21,1) E,nu,h*1e-3,L*1e-3
2       format("rhom=",f6.1," kg/m3 rhow=",f6.1," kg/m3 V= ",e10.3,"N/m")
        write(21,2) rm,rw,V
        g=9.81
        alpha=(4.*D/((rm-rw)*g))**0.25
        w0=-0.125*V*alpha**3./D
!       do x=xmin,xmax,dx
        do xi = 0,floor((xmax-xmin)/dx)           
            x = xmin + dble(xi) * dx 
         tmp=x/alpha
         w=w0*exp(-tmp)*(cos(tmp)+sin(tmp))
         write(21,'(2e18.8)') x,w
        enddo
        close (21)
        end subroutine tre

!----------------------------------------------------------------
        subroutine qua(filename)
        implicit none
        common /slab/ E,nu,h,L,D
        integer :: xi
        real::E,nu,h,L,D
        real::V,xmin,xmax,x,w,dx,rw,rm,g,alpha,w0,tmp
        character (len=15)::filename
        call slabinput
!c        print *;print *,'Input the applied force (N/m)';read *,V
!c        print *,'Input density material above (kg/m3)';read *,rw
!c        print *,'Input density material below (kg/m3)';read *,rm
!c        print *;print *,'Input xmin (km)';read *,xmin;xmin=xmin*1e3
!c        print *;print *,'Input xmax (km)';read *,xmax;xmax=xmax*1e3
!c        print *,'Input dx (km)  ';read *,dx;dx=dx*1e3
	read(11,*);read(11,*) V
	read(11,*);read(11,*) rw
	read(11,*);read(11,*) rm
	read(11,*);read(11,*) xmin;xmin=xmin*1e3
	read(11,*);read(11,*) xmax;xmax=xmax*1e3
	read(11,*);read(11,*) dx;dx=dx*1e3
        open (21,file=filename)
        write(21,'("Infinite broken slab with line load @ 0 and restoring forces")')
1       format ("E=",e10.3,"Pa nu=",e10.3," h=",e10.3,"km L=",e10.3,"km")
        write(21,1) E,nu,h*1e-3,L*1e-3
2       format("rhom=",f6.1," kg/m3 rhow=",f6.1," kg/m3 V= ",e10.3,"N/m")
        write(21,2) rm,rw,V
        g=9.81
        alpha=(4.*D/((rm-rw)*g))**0.25
        w0=-0.25*V*alpha**3./D
!        do x=xmin,xmax,dx
        do xi = 0,floor((xmax-xmin)/dx)           
            x = xmin + dble(xi) * dx 
         tmp=x/alpha
         w=w0*exp(-tmp)*cos(tmp)
         write(21,'(2e18.8)') x,w
        enddo
        close(21)
        end subroutine qua

!----------------------------------------------------------------
        subroutine cin(filename)
        implicit none
        common /slab/ E,nu,h,L,D
        integer :: xi
        real::E,nu,h,L,D
        real::V,M,xmin,xmax,x,w,dx,rw,rm,g,alpha,w0,tmp
        character (len=15)::filename
        call slabinput
!c        print *;print *,'Input the applied force (N/m)';read *,V
!c        print *;print *,'Input the applied torque (Nm)';read *,M
!c        print *,'Input density material above (kg/m3)';read *,rw
!c        print *,'Input density material below (kg/m3)';read *,rm
!c        print *;print *,'Input xmin (km)';read *,xmin;xmin=xmin*1e3
!c        print *;print *,'Input xmax (km)';read *,xmax;xmax=xmax*1e3
!c        print *,'Input dx (km)  ';read *,dx;dx=dx*1e3
	read(11,*);read(11,*) V
	read(11,*);read(11,*) M
	read(11,*);read(11,*) rw
	read(11,*);read(11,*) rm
	read(11,*);read(11,*) xmin;xmin=xmin*1e3
	read(11,*);read(11,*) xmax;xmax=xmax*1e3
	read(11,*);read(11,*) dx;dx=dx*1e3
        open (21,file=filename)
        write(21,'("Semi-infinite slab with line load and torque @ 0 and restoring forces")')
1       format ("E=",e10.3,"Pa nu=",e10.3," h=",e10.3,"km L=",e10.3,"km")
        write(21,1) E,nu,h*1e-3,L*1e-3
2       format("rhom=",f6.1," kg/m3 rhow=",f6.1," kg/m3 V= ",e10.3,"N/m M=",e10.3,"Nm")
        write(21,2) rm,rw,V,M
        g=9.81
        alpha=(4.*D/((rm-rw)*g))**0.25
        w0=-alpha*alpha*0.5/D
!        do x=xmin,xmax,dx
        do xi = 0,floor((xmax-xmin)/dx)           
            x = xmin + dble(xi) * dx 
         tmp=x/alpha
         w=w0*exp(-tmp)*(-M*sin(tmp)+(V*alpha+M)*cos(tmp))
         write(21,'(2e18.8)') x,w
        enddo
        close(21)
        end subroutine cin
        
