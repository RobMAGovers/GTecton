include $(GTECTONHOME)/tools/src/bmake/$(PETSC_ARCH)/system_variables
CC = cc
FC = mpif90
CFLAGS = -O -fPIC
# several fortran program have extra options. Because they do not apply to all
# fortran programs, they are not included in the de fflags.
FFLAGS = -std=f2008 -ffree-form

SRC = ./
BIN = $X/

# CSWITCHES = 

# RM should be set to the name of your favorite rm (file deletion program).

RM = /bin/rm -f

# The action starts here.

all: $(BIN)tri2fe $(BIN)tri2felong $(BIN)tri2feoldlong $(BIN)picknps $(BIN)picknpsold $(BIN)elmside $(BIN)setarea $(BIN)getarea $(BIN)elmsidespecial $(BIN)elmsideJanneke $(BIN)polysort $(BIN)reorder $(BIN)partition $(BIN)elmside2D $(BIN)helloworld

#all: $(BIN)tri2fe $(BIN)tri2feold $(BIN)picknps $(BIN)picknpsold                $(BIN)elmsideold $(BIN)setarea $(BIN)getarea $(BIN)elmsidespecial $(BIN)elmsideJanneke $(BIN)polysort $(BIN)reorder 
#$(BIN)partition


$(BIN)picknps: $(SRC)picknps.c
	$(CC) $(CFLAGS) -o $(BIN)picknps $(SRC)picknps.c -lm

$(BIN)picknpsold: $(SRC)picknpsold.c
	$(CC) $(CFLAGS) -o $(BIN)picknpsold $(SRC)picknpsold.c -lm

$(BIN)elmsideC: $(SRC)elmside.c
	$(CC) $(CFLAGS) -o $(BIN)elmsideC $(SRC)elmside.c -lm

$(BIN)elmside: $(SRC)elmside.f
	$(FC) $(FFLAGS) -ffree-form -o $(BIN)elmside $(SRC)elmside.f

$(BIN)elmsideold: $(SRC)elmsideold.c
	$(CC) $(CFLAGS) -o $(BIN)elmsideold $(SRC)elmsideold.c -lm

$(BIN)elmsidespecial: $(SRC)elmsidespecial.c
	$(CC) $(CFLAGS) -o $(BIN)elmsidespecial $(SRC)elmsidespecial.c -lm

$(BIN)elmsideJanneke: $(SRC)elmsideJanneke.c
	$(CC) $(CFLAGS) -o $(BIN)elmsideJanneke $(SRC)elmsideJanneke.c -lm

$(BIN)elmside2D: $(SRC)elmside2D.c
	$(CC) $(CFLAGS) -o $(BIN)elmside2D $(SRC)elmside2D.c -lm

$(BIN)elmsidenew: $(SRC)elmsidenew.c
	$(CC) $(CFLAGS) -o $(BIN)elmsidenew $(SRC)elmsidenew.c -lm

$(BIN)reorder: $(SRC)reorder.c
	$(CC) $(CFLAGS) -o $(BIN)reorder $(SRC)reorder.c -lm -I$(GTECTONHOME)/tools/metis/Lib $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/libmetis.a

$(BIN)partition: $(SRC)partition.c
	$(CC) $(CFLAGS) -o $(BIN)partition $(SRC)partition.c -lm -I$(GTECTONHOME)/tools/metis/Lib $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/libmetis.a

$(BIN)colorize: $(SRC)colorize.c
	$(CC) $(CFLAGS) -o $(BIN)colorize $(SRC)colorize.c -lm

$(BIN)plotsol: $(SRC)plotsol.c
	$(CC) $(CFLAGS) -o $(BIN)plotsol $(SRC)plotsol.c -lm

$(BIN)unitsquare: $(SRC)unitsquare.c
	$(CC) $(CFLAGS) -o $(BIN)unitsquare $(SRC)unitsquare.c -lm

$(BIN)opt2fe: $(SRC)opt2fe.c
	$(CC) $(CFLAGS) -o $(BIN)opt2fe $(SRC)opt2fe.c -lm

$(BIN)tri2opt: $(SRC)tri2opt.f
	$(FC) $(FFLAGS) -o $(BIN)tri2opt $(SRC)tri2opt.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/stools.a

$(BIN)tri2fe: $(SRC)tri2fe.f
	$(FC) $(FFLAGS) -ffree-form -o $(BIN)tri2fe $(SRC)tri2fe.f $(GTECTONHOME)/lib/tools.a

$(BIN)tri2feoldlong: $(SRC)tri2feoldlong.f
	$(FC) $(FFLAGS) -ffree-form -o $(BIN)tri2feoldlong $(SRC)tri2fe.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/dtools.a

$(BIN)tri2felong: $(SRC)tri2felong.f
	$(FC) $(FFLAGS) -ffree-form -o $(BIN)tri2felong $(SRC)tri2fe.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/dtools.a

#$(BIN)tri2feold: $(SRC)tri2feold.f
#	$(FC) $(FFLAGS) -o $(BIN)tri2feold $(SRC)tri2feold.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/stools.a

$(BIN)tri2feNEW: $(SRC)tri2feNEW.f
	$(FC) $(FFLAGS) -o $(BIN)tri2feNEW $(SRC)tri2feNEW.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/stools.a

$(BIN)setarea: $(SRC)setarea.f
	$(FC) $(FFLAGS) -ffree-form -o $(BIN)setarea $(SRC)setarea.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/stools.a

$(BIN)getarea: $(SRC)getarea.f
	$(FC) $(FFLAGS) -o $(BIN)getarea $(SRC)getarea.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/dtools.a

$(BIN)polysort: $(SRC)polysort.f
	$(FC) $(FFLAGS) -o $(BIN)polysort $(SRC)polysort.f $(GTECTONHOME)/tools/lib/$(PETSC_ARCH)/dtools.a

$(BIN)helloworld:
	gcc -o $(BIN)helloworld test.c

dist:
	tar -cvzf triangle-helpers.tgz makefile *.f *.c

clean: 
	$(RM) $(BIN)tri2fe $(BIN)tri2feold $(BIN)picknps $(BIN)picknpsold $(BIN)elmside $(BIN)elmsideold $(BIN)setarea $(BIN)getarea $(BIN)elmsidespecial $(BIN)elmsideJanneke $(BIN)polysort $(BIN)reorder $(BIN)partition $(BIN)elmside2D

veryclean:
	$(RM) $(BIN)triangle $(BIN)showme $(BIN)tri2opt $(BIN)opt2fe $(BIN)tri2fe $(BIN)picknps $(BIN)elmside #$(BIN)elmsidenew
