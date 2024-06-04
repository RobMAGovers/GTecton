/* routine to compute a fill reducing ordering of sparse matrices
 * Metis routines are used, and the workflow is as follows:
 * create a nodal and element file with triangle 
 * run this program with the domain.1.ele and domain.1.node as inputs
 * the programs returns two new files new_*** with the new ordering
 * proceed as before (i.e., use tri2fe, picknps, etc etc)
 * this routine works directly on the triangle outputs because
 * that saves the boundary markers otherwise lost
 *
 * Version 0.0.1 21-03-2011
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

//  char elementfile[256], nodalpointfile[256];
  char *elementfile = NULL, *nodalpointfile = NULL;
  char *graphfile = NULL;
  char elementfileout[256], nodalpointfileout[256];
  char tmp_line[256];
  char *numberingstyle = NULL;
  char buff[256]; 
  char record[BUFSIZ],error[BUFSIZ];
  FILE *fnp, *fel, *fpe;
  FILE *fnpo, *felo;
  FILE *ftopo;

  int i,j,flag,NUMNP, NDOF, NATTRIB, BM, ID;
  int NUMEL, AT1, AT2, n1, n2, n3;

  double *X,*Y;
  int *M, *MAT;
  int table[10][2];

  int *elmnts, *nxadj, *nadjncy;
  int options[8];
  int *perm, *iperm;

  while ((flag = getopt(argc, argv, "e:n:o:s:")) != EOF)
    switch (flag) {
	case 'e':
	  elementfile = optarg; break;
    case 'n':
      nodalpointfile =  optarg; break;
	case 'o':
	  graphfile = optarg; break;
	case 's':
	  numberingstyle = optarg; break;
    default:
      fprintf(stderr,"Usage: %s -n triangle.node -e triangle.ele [-o graphfile_out -s numbering_style] \n",argv[0]); 
      exit(1);
    }

  if ( elementfile == NULL || nodalpointfile == NULL ){
      fprintf(stderr,"Usage: %s -n triangle.node -e triangle.ele [-o graphfile_out -s numbering_style] \n",argv[0]); 
      exit(1);
  }
  if (numberingstyle == NULL) {
	  numberingstyle = "fortran";
  }

  if ( (fel = fopen(elementfile,"r")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],elementfile);
    perror(error); exit(1);
  }
  if ( (fnp = fopen(nodalpointfile,"r")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],nodalpointfile);
    perror(error); exit(1);
  }

  if ( fscanf(fnp,"%d %d %d %d",&NUMNP,&NDOF,&NATTRIB,&BM) == EOF ) {
    sprintf(error,"%s: EOF \"%s\"",argv[0],nodalpointfile);
    perror(error); exit(1);
  }
  if ( fscanf(fel,"%d %d %d",&NUMEL,&AT1,&AT2) == EOF ) {
    sprintf(error,"%s: EOF \"%s\"",argv[0],elementfile);
    perror(error); exit(1);
  }

  strcpy(nodalpointfileout,"new_"); strcat(nodalpointfileout, nodalpointfile);
  strcpy(elementfileout,"new_");  strcat(elementfileout, elementfile);

  fprintf(stdout, "Number of nodal points %d \n", NUMNP);
  fprintf(stdout, "Number of elements %d \n", NUMEL);

  elmnts 	= (int*) malloc(sizeof(int)*(NUMEL*3));
  nxadj 	= (int*) malloc(sizeof(int)*(NUMNP+1));
  nadjncy 	= (int*) malloc(sizeof(int)*(NUMNP*6));
  
  perm = (int*) malloc(sizeof(int)*NUMNP);
  iperm = (int*) malloc(sizeof(int)*NUMNP);
  
  MAT = (int*) malloc(sizeof(int)*NUMEL);

  for (i=0;i<NUMEL;i++) {
	  fscanf(fel,"%d %d %d %d %d", &ID, &n1, &n2, &n3, &MAT[i]);
	  elmnts[3*i] = n1;  elmnts[3*i+1] = n2;  elmnts[3*i+2] = n3;
  }

  int numflag = 0;
  int etype = 1;
  options[0]=0;
  METIS_MeshToNodal(&NUMEL, &NUMNP, elmnts, &etype, &numflag, nxadj, nadjncy);
  METIS_NodeND(&NUMNP, nxadj, nadjncy, &numflag, options, perm, iperm);


/*
  for(i=0;i<NUMNP+1;i++) {
	  fprintf(stdout,"xadj[%d]:%d\n", i,nxadj[i]);
  }
  for(i=0;i<nxadj[NUMNP];i++) {
	  fprintf(stdout,"adjncy[%d]:%d\n", i,nadjncy[i]);
  }
*/			  

  X = (double*) malloc(sizeof(double)*NUMNP);
  Y = (double*) malloc(sizeof(double)*NUMNP);
  M = (int*) malloc(sizeof(int)*NUMNP);

  if ( (fnpo = fopen(nodalpointfileout,"w")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],nodalpointfileout);
    perror(error); exit(1);
  }
  fprintf(fnpo,"%d %d %d %d\n",NUMNP,NDOF,NATTRIB,BM);
  for (i=0;i<NUMNP;i++) {
	  int ipermute = iperm[i];
	  fscanf(fnp,"%d %lf %lf %d",&ID,&X[ipermute],&Y[ipermute], &M[ipermute]);
  }

  for (i=0;i<NUMNP;i++) { 
	  fprintf(fnpo,"%8d %20.14g %20.14g %3d\n", i,X[i],Y[i],M[i]);
  }


  while(fgets(tmp_line,sizeof(tmp_line),fnp) != NULL)
   {
      // strip trailing '\n' if it exists
      int len = strlen(tmp_line)-1;
      if(tmp_line[len] == '\n') 
         tmp_line[len] = 0;
   }
  fprintf(fnpo,"%s\n", tmp_line);

  fclose(fnp);
  if (fclose(fnpo) == 0) {
      fprintf(stdout, "Succesfully created new nodal point file \"%s\" \n", nodalpointfileout);
  }

  if ( (felo = fopen(elementfileout,"w")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],elementfileout);
    perror(error); exit(1);
  }
  fprintf(felo,"%d %d %d\n",NUMEL, AT1, AT2);
  for (i=0;i<NUMEL;i++) {
	  n1 = iperm[elmnts[3*i]]; n2 = iperm[elmnts[3*i+1]]; n3 = iperm[elmnts[3*i+2]];
	  fprintf(felo, "%10d %8d %8d %8d %3d\n", i, n1, n2, n3, MAT[i]);
  }
  
  while(fgets(tmp_line,sizeof(tmp_line),fel) != NULL)
   {
      // strip trailing '\n' if it exists
      int len = strlen(tmp_line)-1;
      if(tmp_line[len] == '\n') 
         tmp_line[len] = 0;
   }
  fprintf(felo,"%s\n", tmp_line);


  fclose(fel);
  if (fclose(felo) == 0) {
      fprintf(stdout, "Succesfully created new element definition file \"%s\" \n", elementfileout);
  }

  if ( (felo = fopen(elementfileout,"r")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],elementfileout);
    perror(error); exit(1);
  }
  if ( fscanf(felo,"%d %d %d",&NUMEL,&AT1,&AT2) == EOF ) {
    sprintf(error,"%s: EOF \"%s\"",argv[0],elementfile);
    perror(error); exit(1);
  }

  for (i=0;i<NUMEL;i++) {
	  fscanf(felo,"%d %d %d %d %d", &ID, &n1, &n2, &n3, &MAT[i]);
	  elmnts[3*i] = n1;  elmnts[3*i+1] = n2;  elmnts[3*i+2] = n3;
  }

  // finally, compute the connectivity graph for the new domain
//  METIS_MeshToNodal(&NUMEL, &NUMNP, elmnts, &etype, &numflag, nxadj, nadjncy);
  fclose(felo);

  if (graphfile == NULL ) {
	  graphfile = "sparsity.dat";
  }
  if ( (ftopo = fopen(graphfile,"w")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],graphfile);
    perror(error); exit(1);
  }

  int shift = 0;
  if (strcmp(numberingstyle,"fortran")==0 ) {
	  shift = 1;
  }
  int rows = (NUMNP+1) / 10;
  int remn = (NUMNP+1) - 10*rows;
  for (i=0;i<rows;i++) {
	  for (j=0;j<10;j++) {
		  fprintf(ftopo,"%10d", nxadj[i*10+j]+ shift);
	  }
	  fprintf(ftopo, "\n");
  }
  for (j=0;j<remn;j++) {
	  fprintf(ftopo,"%10d", nxadj[rows*10+j]+shift);
  }
  fprintf(ftopo, "\n");

  rows=nxadj[NUMNP]/10;
  remn=nxadj[NUMNP] - 10*rows;

  for (i=0;i<rows;i++) {
	  for (j=0;j<10;j++) {
		  fprintf(ftopo,"%10d", nadjncy[i*10+j]+shift);
	  }
	  fprintf(ftopo, "\n");
  }
  for (j=0;j<remn;j++) {
	  fprintf(ftopo,"%10d", nadjncy[rows*10+j]+shift);
  }
  fprintf(ftopo, "\n");
   if (fclose(ftopo) == 0) {
      fprintf(stdout, "Succesfully written graph file \"%s\" \n", graphfile);
   }
  exit(0);

}

