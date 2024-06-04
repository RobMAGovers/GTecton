
FC = gfortran
FFLAGS = -O -fPIC -DLINUX -fbounds-check -fno-automatic
FSOURCE = Maxwell.f timestep.f

all: $X/Maxwell $X/timestep

$X/Maxwell: Maxwell.f
	$(FC) $(FFLAGS) -o $X/Maxwell Maxwell.f

$X/timestep: timestep.f
	$(FC) $(FFLAGS) -o $X/timestep timestep.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/dtools.a

clean:
	rm -f $X/Maxwell $X/timestep
