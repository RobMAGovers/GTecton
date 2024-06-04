/* reads TOAST (optimize) file and writes TECTON nps and elm files,
   making sure that orientation of elements is OK
*/
#define TRUE       1
#define FALSE      0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

main (int argc, char **argv)
{
  extern char *optarg;
  extern int optind;

  FILE *fp;
  char *toastfile = NULL, *npfile = NULL,  *elmfile = NULL;
  char record[BUFSIZ],error[BUFSIZ];
  int i, j, flag, NUMNP=0, NUMEL=0, NDOF, BOUNDARY, *EL1, *EL2, *EL3, *MAT;
  float *X, *Y;

  while ((flag = getopt(argc, argv, "t:n:e:")) != EOF)
    switch (flag) {
    case 't':
      toastfile = optarg; break;
    case 'n':
      npfile = optarg; break;
    case 'e':
      elmfile = optarg; break;
    default:
      fprintf(stderr,"Usage: %s -t TOAST [-n npfile] [-e elmfile]\n",argv[0]); 
      exit(1);
    }
   if ( toastfile == NULL ) {
     fprintf(stderr,"Usage: %s -t TOAST [-n npfile] [-e elmfile]\n",argv[0]); 
     exit(1);
   }

  if ( (fp = fopen(toastfile,"r")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],toastfile);
    perror(error); exit(1);
  }
  while ( NUMNP == 0) {
    fscanf(fp,"%s",record);
    if ( strncmp(record,"NodeList",8) == 0 ) {
	fscanf(fp,"%d %d",&NUMNP,&NDOF);
	if (NDOF != 2) {
	    fprintf(stderr,"\"%s\": %s specifies NDOF != 2\n",toastfile,record);
	    exit(1);
        }
        X = (float *) malloc(sizeof(float)*NUMNP);
        Y = (float *) malloc(sizeof(float)*NUMNP);
        for (i=0; i<NUMNP; i++) {
          if (fscanf(fp,"%s",record) == EOF) {
	    sprintf(error,"%s: EOF error on \"%s\"",argv[0],toastfile);
    	    perror(error); exit(1);
          }
          if ( strncmp(record,"2B",2) == 0 ) {
            if (fscanf(fp,"%f %f %d",&X[i],&Y[i],&BOUNDARY) == EOF) {
	      sprintf(error,"%s: EOF error on \"%s\"",argv[0],toastfile);
    	      perror(error); exit(1);
            }
          } else if ( strncmp(record,"2N",2) == 0 ) {
            if (fscanf(fp,"%f %f",&X[i],&Y[i]) == EOF) {
	      sprintf(error,"%s: EOF error on \"%s\"",argv[0],toastfile);
    	      perror(error); exit(1);
            }
 	  }
	}
    }
  }

  while ( NUMEL == 0 ) {
    fscanf(fp,"%s",record);
    if ( strncmp(record,"ElementList",11) == 0 ) {
      fscanf(fp,"%d",&NUMEL);
      EL1 = (int *) malloc(sizeof(int)*NUMEL);
      EL2 = (int *) malloc(sizeof(int)*NUMEL);
      EL3 = (int *) malloc(sizeof(int)*NUMEL);
      MAT = (int *) malloc(sizeof(int)*NUMEL);
      for (i=0; i<NUMEL; i++) {
        if (fscanf(fp,"%s",record) == EOF) {
	  sprintf(error,"%s: EOF error on \"%s\"",argv[0],toastfile);
    	  perror(error); exit(1);
        }
        if (fscanf(fp,"%d %d %d %d %d %d",&MAT[i],&j,&j,&EL1[i],&EL2[i],&EL3[i]) == EOF) {
	  sprintf(error,"%s: EOF error on \"%s\"",argv[0],toastfile);
    	  perror(error); exit(1);
 	}
        if ( MAT[i] <= 0 ) MAT[i] = 1;
      }
    }
  }
  fclose(fp);

  if ( npfile != NULL ) {
    if ( (fp = fopen(npfile,"w")) == NULL ) {
      sprintf(error,"%s: open error on \"%s\"",argv[0],npfile);
      perror(error); exit(1);
    }
    for (i=0; i<NUMNP; i++) fprintf(fp,"%5d%5d%14.6e%14.6e\n",i+1,0,X[i],Y[i]);
    fprintf(fp,"end nodal point coordinates\n");
    fclose(fp);
  } else {
    for (i=0; i<NUMNP; i++) fprintf(stdout,"%5d%5d%14.6e%14.6e\n",i+1,0,X[i],Y[i]);
    fprintf(stdout,"end nodal point coordinates\n");
  }

  if ( elmfile != NULL ) {
    if ( (fp = fopen(elmfile,"w")) == NULL ) {
      sprintf(error,"%s: open error on \"%s\"",argv[0],elmfile);
      perror(error); exit(1);
    }
  }
  for (i=0; i<NUMEL; i++) {
    int i1,i2,i3;
    double x1, y1, x2, y2, x3, y3, det;
    i1 = EL1[i]-1; i2 = EL2[i]-1; i3 = EL3[i]-1;
    x1 = X[i1]; y1 = Y[i1]; x2 = X[i2]; y2 = Y[i2]; x3 = X[i3]; y3 = Y[i3];
    det = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
    if (fabs(det) < 1e-6) {
      fprintf(stderr,"determinant = %lf in element %d\n",det,i+1);
      fprintf(stderr,"%5d%5d%5d%5d\n",i+1,EL1[i],EL2[i],EL3[i]);
      fprintf(stderr,"%5d%14.6e%14.6e\n",EL1[i],x1,y1);
      fprintf(stderr,"%5d%14.6e%14.6e\n",EL2[i],x2,y2);
      fprintf(stderr,"%5d%14.6e%14.6e\n",EL3[i],x3,y3);
      exit(1);
    }
    if (det < 0) {
      if ( elmfile != NULL )
        fprintf(fp,"%5d%5d%5d%5d%5d%5d\n",i+1,MAT[i],EL1[i],EL3[i],EL2[i],EL2[i]);
      else
        fprintf(stdout,"%5d%5d%5d%5d%5d%5d\n",i+1,MAT[i],EL1[i],EL3[i],EL2[i],EL2[i]);
    } else {
      if ( elmfile != NULL )
        fprintf(fp,"%5d%5d%5d%5d%5d%5d\n",i+1,MAT[i],EL1[i],EL2[i],EL3[i],EL3[i]);
      else
        fprintf(stdout,"%5d%5d%5d%5d%5d%5d\n",i+1,MAT[i],EL1[i],EL2[i],EL3[i],EL3[i]);
    }
  }
  if (elmfile != NULL)
    fprintf(fp,"end element definitions\n");
  else
    fprintf(stdout,"end element definitions\n");
  fclose(fp);
  exit(0);
}
