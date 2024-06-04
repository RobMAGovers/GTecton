'\" t
.TH PLNPLT 1 "August 26, 2013" "Utrecht University"
.UC 4
.SH NAME
feplt, plnplt \- finite element post-processor.
.SH SYNOPSIS
\fBplnplt\fP [options]
.SH DESCRIPTION
\fIplnplt\fP is an interactive program for examining finite element (FE)
input and output data. \fIplnplt\fP can also be viewed upon as an interface between FE program
\fIpln\fP and plot programs \fIGMT\fP and \fIP\fP. 
.PP
When called without a file argument, \fIplnplt\fP reads commands from 
standard input (stdin). On-line help is provided in this interactive mode by typing "?", or for example "? femi".
When a file argument is given, commands are read from this file until,
after the end-of-file (EOF) has been reached succesfully, commands
are read from stdin.
.SH VERSION
2011
.SH OPTIONS
.TP
.BI "noecho"
Silent mode, i.e. set echoing off (this is the default).
.TP
.BI "echo"
Echoing mode. Major stages of the program are reported.
.TP
.BI "file"
Commands file name.
.SH FILES
Files which are special to the program are described more extensively
in the \fIpln\fP(L) manual pages. More specific, FE input files are
by default assumed to be "TECIN.DAT" and "TEMPIN.DAT". FE output is
by default assumed to reside in "FEDSK.DAT".
.SH ENVIRONMENT
Environmental variable PU should be set to the subdirectory "putil" of
the directory where \fIp\fP resides.
.PP
If set, variable MYTMP is the directory
where temporary files are kept. If not set, temporary storage occurs in
/tmp.
.SH SEE ALSO
GMT(1), pln(1), P(1), pln(7)
.SH "BUGS"
Report bugs to Rob Govers (email: r.govers@uu.nl).
.SH AUTHORS
Written by Rob Govers.
