    program main
*
c reads coordinates and data from separate files
c resamples data onto a regular grid
*
    lue = iflu('stderr')
    luc = nextlu(0)
    call openf(luc,'coordinates.dat','old')
    lud = nextlu(lud,'data.dat','old')

