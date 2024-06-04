// This new version includes the possibily to select points in 3D

/* routine to select specificly marked nodes */
#define TRUE       1
#define FALSE      0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

void output2D(int printcoords, int ID, float x, float y);
void output3D(int printcoords, int ID, float x, float y, float z);
int coordmatch(float coord, float min, float max, float EPS, int range);

void main (int argc, char **argv)
{
  extern char *optarg;
  extern int optind;

  
  FILE *fp, *fs, *pipe;
  char *tectonfile = NULL;
  char *selectionfile = NULL;
  char record[BUFSIZ],error[BUFSIZ],wordcount[BUFSIZ];
  char dummy[200];
  int dimensions = 2;
  int flag, marker = -1, i, j, k,NUMNP, NDOF, NATTRIB, BM, ID, *M, INCR, MARKER, selected, selected_o; 
  float *X,*Y,*Z,*AT1,*AT2, x, y, z;
  float selx = -999.0, sely=-999.0, selz=-999.0;
  float selmaxx = -999.0, selmaxy=-999.0, selmaxz=-999.0;
  int *selm, nmarkers=0, nx=0, ny=0, nz=0, nxy=0, nxz=0, nyz=0, nxyz=0;
  float *selectx, *selecty, *selectz;
  int   xrange =FALSE, yrange=FALSE, zrange=FALSE;  // false when a single value is given, true when range is given
  float EPS=0.1;
  int VERBOSE = FALSE;
  int printcoords = FALSE;
  int selectByRadius=FALSE;


  while ((flag = getopt(argc, argv, "re:vn:m:cd:x:y:z:s:")) != EOF)
    switch (flag) {
    case 's':
      selectionfile = optarg; break;
    case 'n':
      tectonfile = optarg; break;
    case 'm':
      marker = atoi(optarg); break;
    case 'v':
      VERBOSE = TRUE; break;
    case 'd':
      dimensions = atoi(optarg); 
      if (dimensions != 2 && dimensions != 3) {printf("picknps: dimension should be 2 or 3\n"); exit(1); }
      break;
    case 'x':
      sscanf(optarg,"%f:%f",&selx,&selmaxx); break;
    case 'y':
      sscanf(optarg,"%f:%f",&sely,&selmaxy); break;
    case 'z':
      sscanf(optarg,"%f:%f",&selz,&selmaxz); break;
    case 'e':
      sscanf(optarg,"%f",&EPS); break;
    case 'c':
      printcoords = TRUE; break;
    case 'r':
      selectByRadius = TRUE; break;
    default:
      printf("default\n");
      fprintf(stderr,"Usage: %s -n tecton.np [-d 2|3] [-s selectionfile] [-m marker] [-x xvalue[:xmax]] [-y yvalue[:ymax]] [-z yvalue[:zmax]] [-e epsilon]\n",argv[0]); 
      exit(1);
    }

	if (selmaxx != -999)
	{
		xrange = TRUE;
	}
    if (selmaxy != -999)
    {
        yrange = TRUE;
    }
    if (selmaxz != -999)
    {
        zrange = TRUE;
    }



  	if ( tectonfile == NULL || (dimensions != 2 && dimensions != 3) )
	{
        printf("Tectonfile null or no dimensions\n");
      	fprintf(stderr,"Usage: %s -n tecton.np [-d 2|3] [-s selectionfile] [-m marker] [-x xvalue[:xmax]] [-y yvalue[:ymax]] [-z yvalue[:zmax]] [-e epsilon]\n",argv[0]); 
    	exit(1);
  	}

  	if ( selectionfile == NULL && marker == -1 && selx == -999.0 && sely == -999.0 && selz == -999.0) 
	{
    	fprintf(stderr,"%s: no selection criterion specified\n",argv[0]); exit(1);
  	}

    // Check whether the ranges or points are given.
	



	// either select points from the command line...

	if (selectionfile == NULL) 
	{
	    if (dimensions==3) 
    	{
			if (marker != -1) 
			{
		  		nmarkers = 1;
				selm = (int *) malloc(sizeof(int)*nmarkers);
		  		selm[0] = marker;
	  		}
	  		if (selx != -999.0 && sely != -999.0 && selz != -999.0) 
			{
		  		nxyz = 1;
		  		selectx = (float *) malloc(sizeof(float)*nxyz);
		  		selecty = (float *) malloc(sizeof(float)*nxyz);
          		selectz = (float *) malloc(sizeof(float)*nxyz);
		  		selectx[0] = selx;
		  		selecty[0] = sely;
          		selectz[0] = selz;
	  		}

      		if (selx == -999.0 && sely != -999.0 && selz != -999.0) 
			{
          		nyz = 1;
          		selecty = (float *) malloc(sizeof(float)*nyz);
          		selectz = (float *) malloc(sizeof(float)*nyz);
          		selecty[0] = sely;
          		selectz[0] = selz;
      		}
      		if (selx != -999.0 && sely == -999.0 && selz != -999.0) 
			{
          		nxz = 1;
          		selectx = (float *) malloc(sizeof(float)*nxz);
          		selectz = (float *) malloc(sizeof(float)*nxz);
          		selectx[0] = selx;
          		selectz[0] = selz;
      		}
      		if (selx != -999.0 && sely != -999.0 && selz == -999.0) 
			{
          		nxy = 1;
          		selectx = (float *) malloc(sizeof(float)*nxy);
          		selecty = (float *) malloc(sizeof(float)*nxy);
          		selectx[0] = selx;
          		selecty[0] = sely;
      		}
	  		if ( selx != -999.0 && sely == -999.0 && selz == -999.0) 
			{
		  		nx = 1;
		  		selectx = (float *) malloc(sizeof(float)*nx);
		  		selectx[0] = selx;
	  		}
	  		if ( selx == -999.0 && sely != -999.0 && selz == -999.0) 
			{
		  		ny = 1;
		  		selecty = (float *) malloc(sizeof(float)*ny);
		  		selecty[0] = sely;
	  		}
      		if ( selx == -999.0 && sely == -999.0 && selz != -999.0) {
          		nz = 1;
		        selectz = (float *) malloc(sizeof(float)*nz);
          		selectz[0] = selz;
      		}

		}
		if (dimensions==2)
        {
            if (marker != -1)
            {
                nmarkers = 1;
                selm = (int *) malloc(sizeof(int)*nmarkers);
                selm[0] = marker;
            }
            if (selx != -999.0 && sely != -999.0)
            {
                nxy = 1;
                selectx = (float *) malloc(sizeof(float)*nxy);
                selecty = (float *) malloc(sizeof(float)*nxy);
                selectx[0] = selx;
                selecty[0] = sely;
            }
            if ( selx != -999.0 && sely == -999.0)
            {
                nx = 1;
                selectx = (float *) malloc(sizeof(float)*nx);
                selectx[0] = selx;
            }
            if ( selx == -999.0 && sely != -999.0)
            {
                ny = 1;
                selecty = (float *) malloc(sizeof(float)*ny);
                selecty[0] = sely;
            }

		}
	}




  if ( (fp = fopen(tectonfile,"r")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],tectonfile);
    perror(error); exit(1);
  }

  if ( tectonfile != NULL ) {
    sprintf(wordcount,"wc -l %s",tectonfile); pipe = popen(wordcount,"r");
    fgets(record,BUFSIZ,pipe); (void) pclose(pipe);
    sscanf(record,"%d",&NUMNP); NUMNP--;
  }

  if ( selectionfile != NULL ) {
	 if ((fs = fopen(selectionfile,"r")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],selectionfile);
    perror(error); exit(1);
	 }
  }


	// read the input from a selection file

  	if ( selectionfile != NULL ) 
	{
		if (dimensions==2)
        {
			fscanf(fs, "%d %d %d %d", &nmarkers, &nx, &ny,&nxy);
	  		if (nmarkers > 0) 
			{
		  		selm = (int *) malloc(sizeof(int)*nmarkers);
		  		fscanf(fs, "\n%[^\n]", dummy);
		  		for (k=0; k< nmarkers; k++){
		      		fscanf(fs, "\n%[^\n]", dummy);
			  		sscanf(dummy, "%d", &selm[k]);
		  		}
	  		}
	  		if (nx > 0|| nxy > 0) 
			{
		  		selectx = (float *) malloc(sizeof(float)*(nx+nxy));
		  		if (nx > 0) {
		  			fscanf(fs, "\n%[^\n]", dummy);
		  			for (k=0; k< nx; k++){
		      			fscanf(fs, "\n%[^\n]", dummy);
			  			sscanf(dummy, "%f", &selectx[nxy+k]);
		  			}
		  		}
	  		}
	  		if (ny > 0|| nxy > 0) 
			{
		  		selecty = (float *) malloc(sizeof(float)*(ny+nxy));
		  		if (ny > 0) {
		  			fscanf(fs, "\n%[^\n]", dummy);
 		  			for (k=0; k< ny; k++){
		      			fscanf(fs, "\n%[^\n]", dummy);
			  			sscanf(dummy, "%f", &selecty[nxy+k]);
		  			}
		  		}
	  		}
	  		if (nxy > 0) 
			{
		  		fscanf(fs, "\n%[^\n]", dummy);
		  		for (k=0; k< nxy; k++){
		      		fscanf(fs, "\n%[^\n]", dummy);
			  		sscanf(dummy,"%f %f",&selectx[k], &selecty[k]);
		  		}
	  		}
	  		fclose(fs);
	
		}
		if (dimensions==3)
        {
			fscanf(fs, "%d %d %d %d", &nmarkers, &nx, &ny, &nz, &nxy, &nxz, &nyz, &nxyz);
            if (nxyz > 0)
            {
             	fscanf(fs, "\n%[^\n]", dummy);
                for (k=0; k< nxyz; k++){
                    fscanf(fs, "\n%[^\n]", dummy);
                    sscanf(dummy,"%f %f %f",&selectx[k], &selecty[k], &selectz[k]);
                }
            }
            else if (nmarkers > 0 || nx > 0 || nxy > 0 || nxz > 0 || ny > 0 || nz > 0 || nyz > 0)
			{
				printf("Sampling via selection file in 3D is only implemented via nxyz\n"); 
                printf("and triplets of coordinates");

			}
            fclose(fs);

		}
  	}



  for (i=0;i<nxy;i++) {
	 if (VERBOSE) fprintf(stdout, "select xy %d: %f %f\n", i+1, selectx[i], selecty[i]);
  }

  /* actual selection starts here 
   * for each nodal point in the input file (tecin.dat.nps for example)
   * check if it satisfies one of the selection criteria
   * 1. marker ("coloring")
   * 2. (x,y) coordinate  (optional: (x,y,z))
   * 3. single x value
   * 4. single y value
   * (5. single z value)
   *
   * */

	if (dimensions==2)
    {

		for (i=0; i<NUMNP; i++) 
		{

    		fscanf(fp,"%d %d %f %f",&ID,&MARKER,&x,&y);
			if (VERBOSE) {
				fprintf(stdout,"read from tecin %d in xy %d %f %f\n", ID, MARKER,x, y );
			}

			selected = FALSE;

        // ********** testing markers ***************
			if (!selected) {
	    		if (VERBOSE) fprintf(stdout,"testing markers\n");
				for (k=0;k< nmarkers; k++) {
					selected = (MARKER == selm[k]);
					if (selected)
					{
                        output2D(printcoords, ID, x, y);
						break;
					}
				}
			}

        // ********** testing for x,y ***************
			if (!selected) {
	    		if (VERBOSE) fprintf(stdout,"testing xy combinations\n");
        		for (k=0; k < nxy; k++) {

            		selected = (sqrt((x-selectx[k])*(x-selectx[k]) + (y-selecty[k])*(y-selecty[k])) < EPS);

                    if (selected)
					{
                        output2D(printcoords, ID, x, y);
						break;
					}
    			} 
			}

        // ********** testing for x ***************
			if (!selected) {
	    		if (VERBOSE) fprintf(stdout,"testing x combination\n");
				for (k=0;k< nx; k++) {
                    selected = coordmatch(x, selx, selmaxx, EPS, xrange);
//					selected = (fabs(x-selectx[nxy + k]) < EPS);
                    if (selected)
					{ 
                        output2D(printcoords, ID, x, y);
						break;
					}
				}
			}

        // ********** testing for y ***************
			if (!selected) {
	    		if (VERBOSE) fprintf(stdout,"testing y combination\n");
				for (k=0;k< ny; k++) {
                    selected = coordmatch(y, sely, selmaxy, EPS, yrange);
//					selected = (fabs(y-selecty[nxy + k]) < EPS);
                    if (selected)
					{
                        output2D(printcoords, ID, x, y);
						break;
					}
				}
			}
  		}
	}

    if (dimensions==3)
    {
		
		for (i=0; i<NUMNP; i++)
        {
            fscanf(fp,"%d %d %f %f %f",&ID,&MARKER,&x,&y,&z);
            if (VERBOSE) 
			{
				fprintf(stdout,"read from tectin %d in xyz %d %f %f %f\n", ID, MARKER,x, y,z );
			}
//            printf("xyz %f %f %f\n",x,y,z);

            selected = FALSE;

        // ********** testing for markers ***************

            if (!selected) {
                if (VERBOSE) fprintf(stdout,"testing markers\n");
                for (k=0;k< nmarkers; k++) {
                    selected = (MARKER == selm[k]) ;
                    if (selected)
                    {
                        output3D(printcoords, ID, x, y, z);
                        break;
                    }
                }
            }

        // ********** testing for x,y,z ***************


            if (!selected) {
                if (VERBOSE) fprintf(stdout,"testing xyz combinations\n");
                for (k=0; k < nxyz; k++) {
                    if(selectByRadius) {
                       selected = ( (x-selectx[k])*(x-selectx[k]) + (y-selecty[k])*(y-selecty[k]) + (z-selectz[k])*(z-selectz[k]) < EPS*EPS);
                    }
                    else {
                        // select by Manhattan distance.
  	 					selected = coordmatch(x, selx, selmaxx, EPS, xrange) * coordmatch(y, sely, selmaxy, EPS, yrange) * coordmatch(z, selz, selmaxz, EPS, zrange); 
					}
                    if (selected)
                    {
                        output3D(printcoords, ID, x, y, z);
                        break;
                    }
                }
            }

        // ********** testing for x,y ***************

			if (!selected) {
                if (VERBOSE) fprintf(stdout,"testing xy combinations\n");
                for (k=0; k < nxy; k++) {
					selected = coordmatch(x, selx, selmaxx, EPS, xrange) * coordmatch(y, sely, selmaxy, EPS, yrange);
//                    selected = (sqrt((x-selectx[k])*(x-selectx[k]) + (y-selecty[k])*(y-selecty[k])) < EPS);
                    if (selected)
                    {
                        output3D(printcoords, ID, x, y, z);
                        break;
                    }
                }
            }

        // ********** testing for x,z ***************
            if (!selected) {
                if (VERBOSE) fprintf(stdout,"testing xz combinations\n");
                for (k=0; k < nxz; k++) {
                    selected = coordmatch(x, selx, selmaxx, EPS, xrange) * coordmatch(z, selz, selmaxz, EPS, zrange);
//                    selected = (sqrt((x-selectx[k])*(x-selectx[k]) + (z-selectz[k])*(z-selectz[k])) < EPS);
                    if (selected)
                    { 
                        output3D(printcoords, ID, x, y, z);
                        break;
                    }
                }
            }

        // ********** testing for y,z ***************
            if (!selected) {
                if (VERBOSE) fprintf(stdout,"testing yz combinations\n");
                for (k=0; k < nyz; k++) {
                    selected = coordmatch(y, sely, selmaxy, EPS, yrange) * coordmatch(z, selz, selmaxz, EPS, zrange);
//                    selected = (sqrt((y-selecty[k])*(y-selecty[k]) + (z-selectz[k])*(z-selectz[k])) < EPS);
                    if (selected)
                    { 
                        output3D(printcoords, ID, x, y, z);
                        break;
                    }
                }
            }

        // ********** testing for x ***************
            if (!selected) {
                if (VERBOSE) fprintf(stdout,"testing x combination\n");
                for (k=0;k< nx; k++) {
                    selected = coordmatch(x, selx, selmaxx, EPS, xrange);
//                    selected = (fabs(x-selectx[k]) < EPS);
                    if (selected)
					{
                        output3D(printcoords, ID, x, y, z);
						break;
                    }
                }
            }

        // ********** testing for y ***************
            if (!selected) {
                if (VERBOSE) fprintf(stdout,"testing y combination\n");
                for (k=0;k< ny; k++) {
                    selected = coordmatch(y, sely, selmaxy, EPS, yrange);
//                    selected = (fabs(y-selecty[k]) < EPS);
                    if (selected)
                    {
                        output3D(printcoords, ID, x, y, z);
						break;
                    }
                }
            }

        // ********** testing for z ***************
            if (!selected) {
                if (VERBOSE) fprintf(stdout,"testing z combination\n");
                for (k=0;k< nz; k++) {
                    selected = coordmatch(z, selz, selmaxz, EPS, zrange);
//                    selected = (fabs(z-selectz[k]) < EPS);
                    if (selected)
                    {
                        output3D(printcoords, ID, x, y, z);
                        break;
                    }
                }
            }
		}
	}

  fclose(fp); 
  exit(0);
}


int coordmatch(float coord, float min, float max, float EPS, int range)
{
	if (range==TRUE)
	{
		if (coord > min - EPS && coord < max + EPS) 
		{
			return TRUE;
		}
		else
		{
			return FALSE;
		}
	}
	else
	{
		if (coord > min - EPS && coord < min + EPS) 
        {
            return TRUE;
        }
        else
        {
            return FALSE;
        }
	}
}


void output2D(int printcoords, int ID, float x, float y)
{
    if (printcoords)
    {
        printf("%d %f %f\n",ID,x,y);
    }
    else
    {
        printf("%d\n",ID);
    }
}

void output3D(int printcoords, int ID, float x, float y, float z)
{
    if (printcoords)
    {
        printf("%d %f %f %f\n",ID,x,y,z);
    }
    else
    {
        printf("%d\n",ID);
    }
}

