'\" t
.TH PLN 1 "August 26, 2013" "Utrecht University"
.UC 4
.SH NAME
pln \- thermo-mechanical finite element program for 2-D tectonic model
calculations.
.SH SYNOPSIS
\fBpln\fP workpath=directory containing the modelfiles [options]
.SH DESCRIPTION
\fIpln\fP performs finite element (FE) calculations to solve the mechanical
equilibrium as well as the heat equation. Temperature-dependent
material properties and non-linear rheologies are implemented.
\fIpln\fP is the 2D version of the \fBGTECTON\fP FE package.
.PP
When called with only the \fBworkpath\fP as argument arguments, \fIpln\fP assumes that an input file, named
TECIN.DAT, exists in the current working directory. In this file, 
nodal point coordinates, mechanical boundary conditions, material properties
and elements are defined. This file also contains switching variables
which control, for instance, the solution mode. If thermal solutions are
requested, a second file (TEMPIN.DAT) is assumed, which defines the thermal
problem.
Refer to textbooks on the finite element method and the GTECTON manual
(see LITERATURE) for an elaborate discussion of the program and input files.
.SH VERSION
2014
.SH OPTIONS
.TP
.BI "as[cii]"
Writes  ASCII data to
.IR workpath/FEDSK.DAT.
Same as \fBnobi[nary]\fP.
The default is binary.
.TP
.BI "bi[nary]"
Write binary FE-output data to
.IR workpath/FEDSK.DAT.
Same as \fBnoas[cii]\fP.
This is the default
.TP
.BI "echo=0"
Silent mode, i.e. set echoing off (this is the default). See also 'noecho
.TP
.BI "echo, echo=1"
Normal echoing mode. Major stages of the program are reported.
.TP
.BI "echo=2"
Print info on mechanical input file structure. Typically used
when "debugging" the mechanical input file.
.TP
.BI "echo=22"
Print very detailed info on mechanical input file structure. 
Typically used when "debugging" the mechanical input file. 
It is well advised to direct output to file, when using this echo value.
.TP
.BI "echo=3"
Print what is written to "FEDSK.DAT" (only in "ASCII"-mode,
otherwise equivalent to "echo=1")
.TP
.BI "echo=4"
Print info on thermal input file structure. Typically used
when "debugging" the thermal input file.
.TP
.BI "echo=5"
Print info on the interaction between GTecton and PETSc.
Warning: This can yield a lot of output
.TP
.BI "echo=6"
Print info on the time spent in various part of the code. Useful for profiling.
.TP
.BI "echo=8"
Print info on the health of the various threads. Due to the need to synchronise the threads in order to produce readable output, 
this echo value will slow the code down considerably.
.TP
.BI "echo=9"
Print info for general debugging. Here be dragons.
.TP
.BI "fedsk=[dirname/][basename]"
Writes FE-output according to how
.IR dirname
and
.IR basename
are supplied.

When only \fIbasename\fP is given, the output will be written
into \fIworkpath/basename.timestep.processorID\fP.

When only \fIdirname/\fP is given -- the '/' is required in this case
-- the output will be written into \fIworkpath/dirname/FEDSK.DAT.timestep.processorID\fP. (FEDSK
being the default output name.)

When both are given, the output will be written into \fIworkpath/dirname/basename.timestep.processorID\fP.

When this option is not used (\fBfedsk\fP should not be used without an argument)
the output will be written into \fIworkpath/FEDSK.DAT.timestep.processorID\fP.

After a run, the output of the various processors can be combined with \fImergefiles\fP, to obtain a single file for each timestep. See `man mergefiles'.
.TP
.BI "fein=file"
Reads input from
.IR workpath/file
instead of
.IR workpath/TECIN.DAT.
.TP
.BI "feout[=file]"
Writes extensive information on input and output of the FE
calculations to 
.IR workpath/FEOUT.DAT.
Optional =
.IR file 
replaces output file name
.IR workpath/FEOUT.DAT
by
.IR workpath/file. 
Note that this option does not work for runs on multiple processors. For detailed information on the input on multiple-processors runs, consider using echo=22 instead.
.TP
.BI noas[cii]
Write binary FE-output data to
.IR workpath/FEDSK.DAT.
Same as \fBbi[nary]\fP.
This is the default
.TP
.BI "nobi[nary]"
Writes ASCII data to
.IR workpath/FEDSK.DAT.
Same \fBas[cii]\fP.
The default is to write data as binary
.TP
.BI "noecho"
Silent mode, i.e. set echoing off (this is the default). Same as \fBecho=0\fP.
.TP
.BI "nofeout"
Does not write 
.IR workpath/FEOUT.DAT.
This is the default.
.TP
.BI "noso"
During a \fIpln\fP run, all processors generate a temporary file in \fI/tmp\fP containing the TECIN.DAT file with all the included files explicitly included.
These file file can easily become larger than the available space on \fI/tmp\fP. To prevent this, include the files explicitly before running pln,
using \fBsoelim\fP (see: `man soelim'). During the run, provide the result of \fBsoelim\fP as argument to the \fBfein\fP option and add the \fBnoso\fP option.
This will cause all the processes to read the complete file directly.
.TP
.BI "out=file"
Merge standard output and standard error into 
.IR workpath/file.
.TP
.BI "tempin=file"
Reads thermal FE-input from
.IR workpath/file
instead of 
.IR workpath/TEMPIN.DAT
(the default file name).
.TP
.BI "test"
Perform initialization and one single time step (if required) before stopping.
.TP
.BI "workpath"
Directory that contains the model files. This is particularly useful when submitting a pln run to a cluster. 
The nodes of the cluster do not necessarily know where the model is located. This argument provides that reference.
.TP
.BI "-pc_type TYPE"
When using the
.IR PETSc
iterative solvers, this flag allows to set the preconditioner that is used. Note that not all KSP's give stable results! TYPE should be one of   
"none", "jacobi" (Jacobi), "sor" (successive overrelaxation), "lu" (direct solver, based on LU factorization, DEFAULT for mechanical problem), "bjacobi" (block Jacobi), "eisenstat" (symmetric successive overrelaxation, symmetric Gauss-Seidel with Eisenstat's trick), "ilu" (incomplete factorization, DEFAULT for thermal problem),
"icc" (Incomplete Cholesky factorization),
"asm" (additive Schwarz),
"spai" (Sparse Approximate Inverse method),
"nn" (Balancing Neumann-Neumann),
"cholesky" (Cholesky factorization based direct solver),
"pbjacobi" (Point block Jacobi),
"prometheus" (diagonal scaling). More information can be found on the
.IR PETSc
website.
.TP
.BI "-ksp_type TYPE"
When using the
.IR PETSc
iterative solvers, this flag can be used to set the Krylov subspace method (KSP). Not all types work for the GTECTON matrix problem! TYPE should be one of
"richardson" (preconditioned Richardson),
"chebychev" (preconditioned Chebychev),
"cg" (conjugate gradient, DEFAULT for mechanical and thermal problems),
"cgne" (cheeper preconditioned conjugate gradient),
"gmres" (Generalized Minimal Residual, DEFAULT for thermal advection problem),
"fgmres" (Flexible Generalized Minimal Residual),
"tcqmr" (variant of quasi minimal residual),
"bcgs" (Stabilized version of BiConjugate Gradient Squared),
"bcgsl" (Enhanced BiCGStab),
"cgs" (Conjugate Gradient Squared),
"cr" (preconditioned conjugate residuals),
"lsqr" (least squares),
"bicg" (Biconjugate gradient),
"minres" (Minimum Residual),
"symmlq" (SYMMLQ).
More information can be found on the
.IR PETSc
website.
.TP
.BI "-ksp_rtol NUMBER"
Sets relative tolerance for iterative solver in
.IR PETSc
(default = 1e-7).
.SH FILES
The names of all files can be changed (see OPTIONS). Names below refer to
defaults.
.SS Input files
Input files can be built from other include-files. Include-files are indicated
by .so
.IR file,
\fBsoelim\fP can be used to examine the resulting file.
Input files all are in fixed format.
.TP
.BI petsc_options.txt
If the workpath contains this file, the content is read as options for the initialisation of PETSc. 
For available options, see PETSC documentation, such as:
 
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscInitialize.html and
 
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatView.html etc.
.TP
.BI TECIN.DAT
FE input file. In this file the geometry of the finite element mesh and the mechanical properties and boundary conditions are defined.
.TP
.BI TEMPIN.DAT
FE input file. In this file the thermal properties and boundary conditions
are defined.
.SS Output files
.TP
.BI FEOUT.DAT
FE output file. This file gives a verbose report of in- and output
of the finite element program. The file is in ASCII format and tends to
get very big for many applications. It is not used for post-processing
and is intended to be used for elaborate checking of input data only.
.TP
.BI FEDSK.DAT
FE output file. Compact output used for post-processing (see \fIplnplt\fP).
If not explicitly set to ASCII (see OPTIONS), this file has binary format.
.SH ENVIRONMENT
Environment variable MYTMP defines the directory where temporary files are
kept. If not set, temporary storage occurs in /tmp.
.PP
LD_LIBRARY_PATH should point to the PETSC home directory /lib in case that sparse
solvers are being used.
.SH SEE ALSO
pln(7), plnplt(1), triangle(1)
.SH "BUGS OR POSSIBLE FLAWS"
Report bugs to Rob Govers (email: r.govers@uu.nl).
.SH LITERATURE
Govers, R., GTECTON technical reference manual.
