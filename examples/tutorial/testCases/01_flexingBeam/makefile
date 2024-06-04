nProcs=1
elemSize=50000

test: clean grid partition bcs tecin soelim parsim merge extract analytic plot

extract:
	plnplt midline.plnplt > midline.dat
	tail -n 51 midline.dat >| midline_data.dat
	awk '{print $$3*0.001, $$7}' < midline_data.dat >| numericalResultPlot.dat

grid:
	# initial triangulation. max area influences hard coded element IDs in domain.1.area, do not change!
	triangle -pqACa$(elemSize) domain.poly
	tri2fe -n domain.1.node -e domain.1.ele

bcs:
	./setbcs > tecin.dat.bcs

partition:
	partition -n tecin.dat.nps -e tecin.dat.elm -p $(nProcs) -d 2 -f

tecin:
	./tecin.sh

soelim:
	soelim -r TECIN.DAT >| TECIN.FULL

parsim:
	@rm -f pln.log
	mpirun -np $(nProcs) pln bi FEOUT \
    workpath=$(PWD) \
    partinfo=partition.info noso fein=TECIN.FULL fedsk=fedsk.par &> pln.log

merge:
	@rm -f merge.log
	mergefiles partinfo=partition.info bi fedsk=fedsk.par &> merge.log

analytic:
	# compute the analytical solution
	gfortran -fno-automatic -ffree-form -Walign-commons -ffree-form -o flexure flexure.f
	chmod +x flexure
	./flexure input.dat
	tail -n 51 output.dat >| analyticalResult.dat
	awk '{print $$1*0.001, $$2}' < analyticalResult.dat >| analyticalResultPlot.dat

plot:
	./compareDisp.gmt
	ps2pdf compareDisp.ps compareDisplacement.pdf

clean:
	@rm -f sim.log
	@rm -f domain.[12].poly
	@rm -f domain.[12].node
	@rm -f domain.[12].ele
	@rm -f tecin.dat.*
	@rm -f fedsk*
	@rm -rf elastc
	@rm -f partition.info
	@rm -f FEDSK.DAT FEOUT.DAT
	@rm -f STATUS
	@rm -f GTecton.rc
	@rm -f fecore
	@rm -f modeldata.dat*
	@rm -f outputf.*
	@rm -f residu*
	@rm -rf results
	@rm -f output*dat midline*dat
	@rm -f flexure
	@rm -f 2D_bending_bar_compare.ps
	@rm -f Buf*
	@rm -f job* FEDSK*
	@rm -f conn.dat
	@rm -f GTectonLogFile
	@rm -f BTOT.DAT
	@rm -f merge.log
	@rm -rf FEM*
	@rm -f pln.log
	@rm -f np.dat
	@rm -f petsc_interface.log
	@rm -f TECIN.FULL
