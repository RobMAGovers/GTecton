DIR=../cat1
GROFF=/usr/bin/groff -spet -S
CP=cp

all: elmside.1c f3d.1c getarea.1c gmsh2tecton.1c gtecton.1c mergefiles.1c opn.1c opnplt.1c partition.1c picknps.1c pln.1c plnplt.1c plt3d.1c polysort.1c setarea.1c showme.1c tecinbuilder.1c tri2fe.1c triangle.1c setprem.1c axi.1c axiplt.1c

clean:
	$(RM) elmside.1c f3d.1c getarea.1c gmsh2tecton.1c gtecton.1c mergefiles.1c opn.1c opnplt.1c partition.1c picknps.1c pln.1c plnplt.1c plt3d.1c polysort.1c setarea.1c showme.1c tecinbuilder.1c tri2fe.1c triangle.1c setperm.1c axi.1c axiplt.1c
	$(RM) elmside.html f3d.html getarea.html gmsh2tecton.html gtecton.html mergefiles.html opn.html opnplt.html partition.html picknps.html pln.html plnplt.html plt3d.html polysort.html setarea.html showme.html tecinbuilder.html tri2fe.html triangle.html

.SUFFIXES: .1 .1c .html
.1.1c: 
	$(GROFF) -Tascii -mandoc $< > $@
	$(CP) $@ ../cat1/$<
.1.html: 
	$(MANHTML) $< > $@
	$(CP) $@ ../html1/$<.html
