all: generate convert

generate:
	sed -i 's/\t/    /g' trialMesh2D.f
	gfortran -std=f2008 -ffree-form -o trialMesh2D trialMesh2D.f
	chmod +xr trialMesh2D
	./trialMesh2D
	triangle -pqACa10 domain.poly
	tri2fe -n domain.1.node -e domain.1.ele

convert:
	partition -d 2 -p 16 -f 1 -n tecin.dat.nps -e tecin.dat.elm
	partition2paraview

show:
	showme domain.poly
