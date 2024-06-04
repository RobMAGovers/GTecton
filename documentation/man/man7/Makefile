DIR=../cat7
GROFF=/usr/bin/groff -spet -S
CP=cp

all: f3d.7c opn.7c pln.7c axi.7c

clean:
	$(RM) f3d.7c  opn.7c pln.7c

.SUFFIXES: .7 .7c .html
.7.7c: 
	$(GROFF) -Tascii -mandoc $< > $@
	$(CP) $@ ../cat7/$<
.7.html: 
	$(MANHTML) $< > $@
	$(CP) $@ ../html7/$<.html
