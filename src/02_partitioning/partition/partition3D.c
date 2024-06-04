/* routine to compute partitioning of a mesh created in triangle format
 * input is a pair of element and nodal point information, and a 
 * number indicating the number of partitions
 * output is an adjecency list and properly formatted elements and nodal point files.
 * In essence, this program has three stages:
 * 1: Reading data
 * 2: Calling METSI to do the partitioning
 * 3: Write the partitioning to file.
 */

#define TRUE       1
#define FALSE      0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <metis.h>

main (int argc, char **argv)
{
  extern char *optarg;
  extern int optind;

  char *elementfile = NULL, *nodalpointfile = NULL;
  char elementfileout[256], nodalpointfileout[256];
  char partitioninfo[100];
  char record[BUFSIZ],error[BUFSIZ];

  int  parts;
  FILE *fel, *fnp;
  FILE *felo, *fnpo;

  int i,j,flag,NUMNP, NDOF, NATTRIB, BM, ID;
  int NUMEL, AT1, AT2, n1, n2, n3;

  int *elmnts, *npart, *epart, *matel;
  int *nxadj, *nadjncy;
  double x, y;
  int fortranstyle;


  fortranstyle = FALSE;
  while ((flag = getopt(argc, argv, "e:n:p:f")) != EOF)
    switch (flag) {
	case 'e':
	  elementfile = optarg; break;
    case 'n':
      nodalpointfile =  optarg; break;
	case 'p':
	  sscanf(optarg, "%d", &parts); break;
	case 'f':
      fortranstyle = TRUE; break;
    default:
      fprintf(stderr,"Usage: %s -n tecin.dat.nps -e tecin.dat.elm -p #parts [-f]\n",argv[0]); 
      exit(1);
    }

  if ( elementfile == NULL || nodalpointfile == NULL || parts <= 0 ){
      fprintf(stderr,"Usage: %s -n tecin.dat.nps -e tecin.dat.elm -p #parts [-f]\n",argv[0]); 
      exit(1);
  }

  if ( (fel = fopen(elementfile,"r")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],elementfile);
    perror(error); exit(1);
  }
  if ( (fnp = fopen(nodalpointfile,"r")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],nodalpointfile);
    perror(error); exit(1);
  }


//  while (EOF != (scanf("%*[^\n]"), scanf("%*c"))) 
//    ++lines;






  if ( fscanf(fnp,"%d %d %d %d",&NUMNP,&NDOF,&NATTRIB,&BM) == EOF ) {
    sprintf(error,"%s: EOF \"%s\"",argv[0],nodalpointfile);
    perror(error); exit(1);
  }
  if ( fscanf(fel,"%d %d %d",&NUMEL,&AT1,&AT2) == EOF ) {
    sprintf(error,"%s: EOF \"%s\"",argv[0],elementfile);
    perror(error); exit(1);
  }
  fprintf(stdout, "Number of nodal points %d \n", NUMNP);
  fprintf(stdout, "Number of elements %d \n", NUMEL);

  elmnts 	= (int*) malloc(sizeof(int)*(NUMEL*3));
  epart 	= (int*) malloc(sizeof(int)*(NUMEL));
  matel 	= (int*) malloc(sizeof(int)*(NUMEL));
  npart 	= (int*) malloc(sizeof(int)*(NUMNP));
  nxadj 	= (int*) malloc(sizeof(int)*(NUMNP+1));
  nadjncy 	= (int*) malloc(sizeof(int)*(NUMNP*6));

  for (i=0;i<NUMEL;i++) {
	  if (AT2 == 1) {
		  fscanf(fel,"%d %d %d %d %d", &ID, &n1, &n2, &n3, &matel[i]);
	  } else if (AT2 == 0) {
		  fscanf(fel,"%d %d %d %d", &ID, &n1, &n2, &n3);
	  }
	  if (fortranstyle) {
		  elmnts[3*i] = n1+1;  elmnts[3*i+1] = n2+1;  elmnts[3*i+2] = n3+1;
	  } else {
		  elmnts[3*i] = n1;  elmnts[3*i+1] = n2;  elmnts[3*i+2] = n3;
	  }
  }
  fclose(fel);

  int etype = 1; // triangles
  int numflag = 0; //C or FORTRAN style numbering
  if (fortranstyle) {
	  numflag = 1;
  } 
  int edgecut;

//*************************** METIS calls *********************

  // get mesh partitioning
  if (parts > 1) {
	  METIS_PartMeshDual(&NUMEL, &NUMNP, elmnts, &etype, &numflag, &parts, &edgecut, epart, npart);

  // in case of fortran style, epart and npart start at 1; change this for Petsc compatibility
  // and subtract one

      if (fortranstyle) {
		  for (i=0;i<NUMEL;i++) {
			  epart[i] = epart[i]-1;
		  }
		  for (i=0;i<NUMNP;i++) {
			  npart[i] = npart[i] - 1;
		  }
	  }

  } else {
	  // Partition into one single partition. There is effectively no partition.

	  // always keep partition numbering started from 0!!!! 
	  // this has to do with compatibility with Petsc rank and size 

	  if (fortranstyle) {
		  for (i=0;i<NUMEL;i++) epart[i] = 0;
		  for (i=0;i<NUMNP;i++) npart[i] = 0;
	  } else {
		  for (i=0;i<NUMEL;i++) epart[i] = 0;
		  for (i=0;i<NUMNP;i++) npart[i] = 0;
	  }

  }

  // get graph information
  METIS_MeshToNodal(&NUMEL, &NUMNP, elmnts, &etype, &numflag, nxadj, nadjncy);

//********************* Done with METIS ********************

//  for (i=0;i<NUMNP+1;i++) {
//	  fprintf(stdout,"nxadj[%d] = %d\n",i,nxadj[i]);
// }
//for (i=0;i<NUMNP;i++) {
//	  fprintf(stdout, "%d %d %d %d %d %d\n",nadjncy[6*i + 0],nadjncy[6*i + 1],nadjncy[6*i + 2],nadjncy[6*i + 3],nadjncy[6*i + 4],nadjncy[6*i + 5]);
//  }

  // create element file with partitioning info
  if (fortranstyle) {
	  sprintf(elementfileout,"tecin.dat.partf.elm");
  } else {
	  sprintf(elementfileout,"tecin.dat.part.elm");
  }
  if ( (felo = fopen(elementfileout,"w")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],elementfileout);
    perror(error); exit(1);
  }
  for (i=0;i<NUMEL;i++) {
	  if (fortranstyle) {
		  fprintf(felo, "%5d %8d %5d %8d %8d %8d %8d\n", epart[i],i+1, matel[i], elmnts[3*i], elmnts[3*i+1],elmnts[3*i+2],elmnts[3*i+2]);
	  } else {
		  fprintf(felo, "%5d %8d %8d %8d\n", epart[i],elmnts[3*i], elmnts[3*i+1],elmnts[3*i+2]);
	  }
  }
  if (fortranstyle) fprintf(felo, "end partitioned element file\n");
  fclose(felo);

  // create nodal point / connectivity file 
  int nnbrs = 0;
  if (fortranstyle) {
	  sprintf(nodalpointfileout,"tecin.dat.partf.nps");
  } else {
	  sprintf(nodalpointfileout,"tecin.dat.part.nps");
  }
  if ( (fnpo = fopen(nodalpointfileout,"w")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],nodalpointfileout);
    perror(error); exit(1);
  }

  for (i=0;i<NUMNP;i++) {
	  fscanf(fnp,"%d %lf %lf %d",&ID,&x, &y, &AT1);
	  if (fortranstyle) {
//		  fprintf(fnpo, "%5d %10d %24.18g %24.18g", npart[i],i+1,x,y);
		  fprintf(fnpo, "%5d %10d %14.6g %14.6g", npart[i],i+1,x,y);
	  } else {
		  fprintf(fnpo, "%5d  %24.18g %24.18g", npart[i],x,y);
	  }
	  nnbrs = nxadj[i+1]-nxadj[i];
//	  fprintf(stdout,"node %d nbrs %d\n", i, nnbrs);			 
	  fprintf(fnpo, "%5d ", nnbrs);
	  for (j=0;j<nnbrs;j++) {
		  if (fortranstyle) {
			  fprintf(fnpo, "%10d ", nadjncy[nxadj[i]-1+j]);
		  } else {
			  fprintf(fnpo, "%10d ", nadjncy[nxadj[i]+j]);
		  }

	  }
	  fprintf(fnpo, "\n");
  }
  if (fortranstyle) fprintf(fnpo, "end partitioned nodal point file\n");
  fclose(fnp);
  fclose(fnpo);
  sprintf(partitioninfo,"partition.info");
  int partinfo[parts][2];
  for (i=0;i<parts;i++) {
	  for (j=0;j<2;j++) {
		  partinfo[i][j] = 0;
	  }
  }

  if ( (fnp = fopen(partitioninfo,"w")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],partitioninfo);
    perror(error); exit(1);
  }
  for (i=0;i<NUMEL;i++) {	  
	  partinfo[epart[i]][1]++;
  }
  for (i=0;i<NUMNP;i++) {
	  partinfo[npart[i]][0]++;
  }
  fprintf(fnp,"%10d\n",parts);
  for (i=0;i<parts;i++) {
	  fprintf(fnp,"%10d%10d%10d\n", i, partinfo[i][0],partinfo[i][1]);
  }
  fclose(fnp);
}
