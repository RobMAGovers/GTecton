all: clean generate plnplt plot show

generate:
	sed -i 's/\t/    /g' trialMesh2D.f
	gfortran -std=f2008 -ffree-form -o trialMesh2D trialMesh2D.f
	chmod +xr trialMesh2D
	./trialMesh2D

plnplt:
	plnplt mesh.plnplt

plot:
	./mesh.gmt

show:
#	gv partitions.ps
	cat partition.info

clean:
	rm -rf plnplt
	rm -f tecin.dat*
	rm -f dummyTECIN.dat
