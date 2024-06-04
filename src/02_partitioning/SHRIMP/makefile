sources = enumerates.F files.F MPIhandling.F geometry.F meshData.F readInput.F shrimp.F splitPointCloud.F findNeighbors.F writeOutput.F matchElements.F

all: clean compile run combine


compile: ${sources}
	sed -i 's/\t/    /g' shrimp.F
	mpif90 -o shrimp -g -ffree-form -fallow-argument-mismatch \
           -I/trinity/opt/apps/software/openmpi/openmpi-4.0.4/build/include \
           ${sources}

run:
	mpirun -np 2 ./shrimp -p /scratch/lukas/PRIMP/testMeshes/3D/kubus/1e5/tecin.dat.nps -e /scratch/lukas/PRIMP/testMeshes/3D/kubus/1e5/tecin.dat.elm -i gtecton -o gtecton -n 4 -d 3 &> shrimp.log

vtu:
	mpirun -np 2 ./shrimp -p /scratch/lukas/PRIMP/testMeshes/3D/kubus/1e5/tecin.dat.nps -e /scratch/lukas/PRIMP/testMeshes/3D/kubus/1e5/tecin.dat.elm -i gtecton -o vtu -n 1 -d 3 &> shrimp.log

combine:
	cat tecin.dat.partf.nps1 >  tecin.dat.partf.nps
	cat tecin.dat.partf.nps2 >> tecin.dat.partf.nps 
	cat tecin.dat.partf.elm1 >  tecin.dat.partf.elm
	cat tecin.dat.partf.elm2 >> tecin.dat.partf.elm



clean:
	rm -f shrimp.log compile.log
	rm -f job*
	rm -f *.mod

plot2D:
	./partitioning.gmt
	gmt psconvert partitionedPoints.ps
	gmt psconvert partitionedElements.ps


count: countBorders.f
	sed -i 's/\t/    /g' countBorders.f
	gfortran -g -o countBorders -fcheck=bounds -std=f2008 -ffree-form countBorders.f
	chmod +x countBorders
	./countBorders
