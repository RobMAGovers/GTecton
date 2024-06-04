/* routine to compute partitioning of a mesh created in triangle format
 * input is a pair of element and nodal point information, and a 
 * number indicating the number of partitions
 * output is an adjecency list and properly formatted elements and nodal point files.
 * In essence, this program has three stages:
 * 1: Reading data
 * 2: Calling METIS to do the partitioning
 * 3: Write the partitioning to file.
 */

#define TRUE  1
#define FALSE 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <metis.h>
#include <getopt.h>


void main (int argc, char **argv)
{
	extern char *optarg;

	char *elementfile = NULL, *nodalpointfile = NULL;
	char elementfileout[256], nodalpointfileout[256];
	char partitioninfo[100];
	char record[BUFSIZ], error[BUFSIZ];
	int  verbose;

//	partition file
	FILE    *FP_partition;

// move to tecin:
	FILE    *FP_elems_input,  *FP_points_input;
	FILE    *FP_elems_output, *FP_points_output;

	int      nElems, nPoints;
    char     c;

  	int *connectivity, *nodePartitioning, *elemPartitioning, *material;
  	int *adjacencyOffset, *adjacency;
    int nCommon;
    int *connectivityOffset;

    int      dimensions;
	int      parts, nPartitions;

// *********

    int countFrom0or1;

	int i, j, flag, NUMNP, NDOF, NATTRIB, BM, ID;
	int NUMEL, AT1, AT2, n1, n2, n3, n4, ch;

	int *elmnts, *npart, *epart, *matel, *dimension;
	int *nxadj, *nadjncy;
	double x, y, z;

  	while ((flag = getopt (argc, argv, "e:n:p:d:fzv")) != EOF)
    	switch (flag)
    	{
     	case 'e':
			elementfile = optarg;
			break;
      	case 'n':
			nodalpointfile = optarg;
			break;
        case 'd':
            sscanf (optarg, "%d", &dimensions);
            break;
      	case 'p':
	  		sscanf (optarg, "%d", &parts);
      		nPartitions = parts;
			break;
      	case 'z':
			break;
		case 'f':
            break;
        case 'v':
			verbose = TRUE;
            break;
      	default:
			fprintf (stderr,
		 	"Usage: %s -n tecin.dat.nps -e tecin.dat.elm -p #parts -d #dimensions (2 or 3) [-f]\n",
		 	argv[0]);
			exit (1);
      	}


	if (verbose == TRUE)
	{
		printf("Validating command line arguments\n");
	}


  	// validate command line arguments
	if (elementfile == NULL || nodalpointfile == NULL || parts <= 0 || (dimensions != 2 && dimensions != 3))
    {
      	fprintf (stderr,
	       "Usage: %s -n tecin.dat.nps -e tecin.dat.elm -p #parts -d #dimensions (2 or 3) [-f]\n",
	       argv[0]);
      	exit (1);
    }


    if (verbose	== TRUE) 
    {
    	printf("Command line arguments OK\n");
    }



// *************************************************
//  test; read tecin.dat.nps and ...elm in stead of domain files
// *************************************************

    if (verbose	== TRUE) 
    {
    	printf("Opening elements file\n");
    }

	if ((FP_elems_input = fopen (elementfile, "r")) == NULL)
	{
		sprintf (error, "%s: open error on \"%s\"", argv[0], elementfile);
		perror (error);
		exit (1);
	}

    if (verbose == TRUE)
    {
     	printf("Opening nodal points file\n");
    }

    if ((FP_points_input = fopen (nodalpointfile, "r")) == NULL)
    {
        sprintf (error, "%s: open error on \"%s\"", argv[0], nodalpointfile);
        perror (error);
        exit (1);
    }

	nElems = 0;
	nPoints = 0;

    if (verbose == TRUE)
    {
     	printf("Count elements\n");
    }


	while ( (c=fgetc(FP_elems_input)) != EOF ) {
        if ( c == '\n' )
            nElems++;
    }

    if (verbose == TRUE)
    {
     	printf("Count nodal points\n");
    }


    while ( (c=fgetc(FP_points_input)) != EOF ) {
        if ( c == '\n' )
            nPoints++;  
    }

	// remove 1 to compensate for 'end' statement
	nElems--;
	nPoints--;

    if (verbose == TRUE)
    {
     	printf("Counted %i nodal points\n", nPoints);
     	printf("Counted %i elements\n", nElems);
    }

    rewind(FP_points_input);
    rewind(FP_elems_input);


//////////////////////////////////////////////////////////
// allocate room for the elements 

	if (dimensions == 2)	
	{
		connectivity  = (int *) malloc (sizeof (int) * (nElems * 3));
	}
	if (dimensions == 3)
    {
        connectivity  = (int *) malloc (sizeof (int) * (nElems * 4));
    }

  	elemPartitioning   = (int *) malloc (sizeof (int) * (nElems + 1));
  	material   = (int *) malloc (sizeof (int) * (nElems));
  	nodePartitioning   = (int *) malloc (sizeof (int) * (nPoints + 1));

    connectivityOffset = (int *) malloc (sizeof (int) * (nElems + 1));

	// build index that is needed by PartMeshDual
    for (i=0;i<nElems+1;i++){
	    if (dimensions == 2){
			connectivityOffset[i] = 3 * i;
		}
        if (dimensions == 3){
    	    connectivityOffset[i] =	4 * i;
        }

	}



////////////////////////////////////////////////////////
// Read the elements input data.

    if (verbose == TRUE)
    {
     	printf("Read elements\n");
    }


	for (i=0;i<nElems;i++) 
	{
        fscanf(FP_elems_input,"%d %d %d %d %d %d", &ID, &material[i], &n1, &n2, &n3, &n4);

		// -1 to switch to GTECTON numbering to C numbering.
		if (dimensions==2)
		{
        	connectivity[3*i]   = n1 - 1;  
			connectivity[3*i+1] = n2 - 1;  
			connectivity[3*i+2] = n3 - 1;
		}
		if (dimensions==3) 
        {
            connectivity[4*i]   = n1 - 1;  
			connectivity[4*i+1] = n2 - 1;  
			connectivity[4*i+2] = n3 - 1; 
            connectivity[4*i+3] = n4 - 1;  
        }
  	}

  	// number of nodes an element must share to share an graph edge
    if (dimensions==2)
    {
        nCommon = 2; // triangles
    }
    if (dimensions==3)
    {
        nCommon = 3; // tetrahedrons
    }


    int numflag = 0; //C or FORTRAN style numbering


  	int edgecut;



// ********************************************************** 
// *************************** METIS calls ******************
// **********************************************************
// (note that only the element input file has been read. 
// The NPS file not yet, as it is not necessary to compute the
// partitioning.


    if (verbose == TRUE)
    {
     	printf("Calling METIS to create the partitioning\n");
    }

  	// get mesh partitioning; expressed in epart and npart
  	if (nPartitions > 1) 
	{

	    if (verbose == TRUE)
    	{
	        printf("needs multiple partitions, calling METIS PartMeshDual\n");
	    }

/*
		printf("nElems %i \n", nElems);
    	printf("nPoints %i \n", nPoints);

	    for(int i=0; i<nElems+1; i++)
	    {
	     	printf("connectivityOffset[%i] : %i\n", i, connectivityOffset[i]);
	    }
	    for(int i=0; i<nElems; i++)
	    {
	     	printf("connectivity[%i] : %i %i %i\n", i, connectivity[3*i], \
                                                       connectivity[3*i+1], \
                                                       connectivity[3*i+2]);
	    }
*/

      	METIS_PartMeshDual(&nElems, \
                           &nPoints, \
                           connectivityOffset, \
                           connectivity, \
                           NULL, \
       	       	       	   NULL, \
                           &nCommon, \
                           &nPartitions, \
       	       	       	   NULL, \
       	       	       	   NULL, \
                           &edgecut, \
                           elemPartitioning, \
                           nodePartitioning);


        if (verbose == TRUE)
        {
            printf("finished METIS PartMeshDual\n");
        }

  	} 
	else 
	{

        if (verbose == TRUE)
        {
            printf("needs only 1 partition, adding partition number\n");
        }

      // Only 1 partition, no need to call expensive METIS call.
      // always keep partition numbering started from 0!!!!
      // this has to do with compatibility with Petsc rank and size.
      for (i=0;i<nElems;i++) elemPartitioning[i] = 0;
      for (i=0;i<nPoints;i++) nodePartitioning[i] = 0;
    }

    if (verbose == TRUE)
    {
     	printf("Calling METIS MeshToNodal\n");
    }

	countFrom0or1 = 0;

    // get graph information
    METIS_MeshToNodal(&nElems, \
                      &nPoints, \
                      connectivityOffset, \
                      connectivity, \
                      &countFrom0or1, \
                      &adjacencyOffset, \
                      &adjacency);

    // adjacencyOffset[i] contains the offset in adjacency where the neighbors of vertex i begin
    // adjacency contains the IDs of neighboring vertices of each vertex.

    // example mesh from METIS manual:

    // 0 -- 1 -- 2 -- 3 -- 4
    // |    |    |    |    |
    // 5 -- 6 -- 7 -- 8 -- 9
    // |    |    |    |    |
    // 10 - 11 - 12 - 13 - 14

    // would give:

    // adjacencyOffset: 0 2 5 8 11 13 16 20 24 28 31 33 36 39 42 44

    // and adjacency: 1 5 0 2 6 1 3 7 2 4 8 3 9 0 6 10 1 5 7 11 2 6 8 12 3 
    //                  7 9 13 4 8 14 5 11 6 10 12 7 11 13 8 12 14 9 13

	// from offset 0 up to but not including 2 are the neighbours of point 0, which are 1 and 5
    // from offset 2 up to but not including 5 are the neighbours of point 1, which are 0, 2 and 6
    // from offset 5 up to but not including 8 are the neighbours of point 2, which are 1, 3 and 7
	// etc

    // number of neighbors of one point can be determined by substracting offsets

    if (verbose == TRUE)
    {
     	printf("Finished with METIS\n");
    }

// **********************************************************
// ********************* Done with METIS ********************
// ********************************************************** 

  	// create element file with partitioning info
  	if (countFrom0or1 == 1) 
	{
    	sprintf(elementfileout,"tecin.dat.partf.elm");
  	} 
	else 
	{
      	sprintf(elementfileout,"tecin.dat.partf.elm");
  	}

  	if ( (FP_elems_output = fopen(elementfileout,"w")) == NULL ) 
	{
    	sprintf(error,"%s: open error on \"%s\"",argv[0],elementfileout);
	    perror(error); exit(1);
  	}

//////////////////////////////////////////////////////////

    if (verbose == TRUE)
    {
     	printf("Write partitioned elements\n");
    }


    if (dimensions==2)
    {
		for (i=0;i<3*nElems;i++)
		{
			connectivity[i]++;
		}
	}
    if (dimensions==3)
    {
    	for (i=0;i<4*nElems;i++)
    	{
    	    connectivity[i]++;
    	}
    }




	if (dimensions==2)
	{
  		for (i=0;i<nElems;i++) 
		{
       		fprintf(FP_elems_output, "%5d %12d %12d %12d %12d %12d %12d\n", \
                    elemPartitioning[i], \
                    i+1, \
                    material[i], \
                    connectivity[3*i], \
                    connectivity[3*i+1],\
                    connectivity[3*i+2],\
                    connectivity[3*i+2]);
  		}
	}
    if (dimensions==3)
    {
        for (i=0;i<nElems;i++)  
        {
            fprintf(FP_elems_output, "%5d %12d %12d %12d %12d %12d %12d\n", \
                    elemPartitioning[i], \
                    i+1, \
                    material[i], \
                    connectivity[4*i], \
                    connectivity[4*i+1],\
                    connectivity[4*i+2],\
                    connectivity[4*i+3]);
        }
    }

	fprintf(FP_elems_output, "end partitioned element file\n");


	fclose(FP_elems_output);





	// create nodal point / connectivity file
	int nnbrs = 0;

	if (countFrom0or1 == 1) 
	{
    	sprintf(nodalpointfileout,"tecin.dat.partf.nps");
	} 
	else 
	{
    	sprintf(nodalpointfileout,"tecin.dat.partf.nps");
	}

	if ( (FP_points_output = fopen(nodalpointfileout,"w")) == NULL ) 
	{
    	sprintf(error,"%s: open error on \"%s\"",argv[0],nodalpointfileout);
    	perror(error); exit(1);
	}

    if (verbose == TRUE)
    {
     	printf("Write partitioned nodes\n");
    }




	if (dimensions==2)
	{
		for (i=0;i<nPoints;i++) 
		{
      		fscanf(FP_points_input,"%d %d %lf %lf ",&ID, &AT1, &x, &y);

            fprintf(FP_points_output, "%5d %12d %12d %25.17E %25.17E", \
                 nodePartitioning[i], i+1, AT1, x,y);

            // compute numbers of neighbours

			nnbrs = adjacencyOffset[i+1]-adjacencyOffset[i];
//    		fprintf(stdout,"node %d nbrs %d\n", i, nnbrs);
      		fprintf(FP_points_output, "%5d ", nnbrs);

      		for (j=0;j<nnbrs;j++) 
			{
          		fprintf(FP_points_output, "%12d ", adjacency[adjacencyOffset[i]+j]+1);
      		}
      		fprintf(FP_points_output, "\n");
  		}
	}
	if (dimensions==3)
	{
        for (i=0;i<nPoints;i++)
        {
            fscanf(FP_points_input,"%d %d %lf %lf %lf",&ID, &AT1, &x, &y, &z);

            fprintf(FP_points_output, "%5d %12d %12d %25.17E %25.17E %25.17E", nodePartitioning[i],i+1,AT1,x,y,z);

            nnbrs = adjacencyOffset[i+1]-adjacencyOffset[i];
            fprintf(FP_points_output, "%5d ", nnbrs);

            for (j=0;j<nnbrs;j++)  
            {
                fprintf(FP_points_output, "%12d ", adjacency[adjacencyOffset[i]+j]+1);
            }
            fprintf(FP_points_output, "\n");
        }
    }

	fprintf(FP_points_output, "end partitioned nodal point file\n");

  fclose(FP_points_input);
  fclose(FP_elems_input);
  fclose(FP_points_output);




  // finally write the partition info

  sprintf(partitioninfo,"partition.info");
  int partitionInfo[parts][2];

  for (i=0;i<nPartitions;i++) {
      for (j=0;j<2;j++) {
          partitionInfo[i][j] = 0;
      }
  }



  if ( (FP_partition = fopen(partitioninfo,"w")) == NULL ) {
    sprintf(error,"%s: open error on \"%s\"",argv[0],partitioninfo);
    perror(error); exit(1);
  }
  for (i=0;i<nElems;i++) {
      partitionInfo[elemPartitioning[i]][1]++;
  }
  for (i=0;i<nPoints;i++) {
      partitionInfo[nodePartitioning[i]][0]++;
  }


  fprintf(FP_partition,"%5d\n",nPartitions);

  for (i=0;i<nPartitions;i++) {
      fprintf(FP_partition,"%5d%12d%12d\n", i, partitionInfo[i][0],partitionInfo[i][1]);
  }

	fclose(FP_partition);

}


