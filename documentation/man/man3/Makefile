GROFF=/usr/bin/groff -spet -S
CP=cp

all: cstat.3c csystem.3c dtime.3c exitp.3c fdate.3c flush.3c fnblnk.3c hostname.3c howtouse.3c lnblk.3c rm.3c signalc.3c stat.3c traper.3c triangulate.3c cstat.html csystem.html dtime.html exitp.html fdate.html flush.html fnblnk.html hostname.html howtouse.html lnblk.html rm.html signalc.html stat.html traper.html triangulate.html

clean:
	$(RM) cstat.3c csystem.3c dtime.3c exitp.3c fdate.3c flush.3c fnblnk.3c hostname.3c howtouse.3c lnblk.3c rm.3c signalc.3c stat.3c traper.3c triangulate.3c cstat.html csystem.html dtime.html exitp.html fdate.html flush.html fnblnk.html hostname.html howtouse.html lnblk.html rm.html signalc.html stat.html traper.html triangulate.html

.SUFFIXES: .3 .3c .html
.3.3c: 
	$(GROFF) -Tascii -mandoc $< > $@
	$(CP) $@ ../cat3/$<
.3.html: 
	$(GROFF) -Thtml $< > $@
	$(CP) $@ ../html3/$<
