LTARGET = opnlib.o
LPRJS   = $(LTARGET:.o=.prj)

${LTARGET} : ../io.i
cln: FORCE
	$(RM) *.prj a.out core
verycln: FORCE
	$(RM) ${LTARGET}

chk: $(LPRJS)

include ../Makefile
