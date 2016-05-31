/* -------------------------------------------
 * PARAVT - PARALLEL VORONOI TESSELLATION
 *
 * COMPUTE VT ON LARGE PARTICLE DISTRIBUTIOINS 
 * USING QHULL LIBRARY AND MPI.
 * ALSO COMPUTE SEVERAL QUANTITIES AND OUTPUTS
 * FOR ASTROPHYSICAL PURPOSES.
 * Roberto Gonzalez -- regonzar@astro.puc.cl
 * https://github.com/regonzar/paravt
--------------------------------------------*/
#include "qhullsrc_r/qhull_ra.h"
#include "config.h"		/*Configuration file */
#include <mpi.h>
#include <stdio.h>
#include "vars.h"
#include "func.h"
#include "vt.h"

/*Main Program*/

int
main (int argc, char *argv[])
{

  double t1, t2, t3, t4, t5, t6, t7;	/*timers */

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size (MPI_COMM_WORLD, &NTask);

  t1 = MPI_Wtime ();

  initvars ();

  check_params (argc, argv);

  splitdomain ();		/*select mesh size for domain decomposition based on NTask, produce spx,spy,spz */

  read_header ();		/*read file header and set file format */

  MPI_Barrier (MPI_COMM_WORLD);

  read_part ();			/*read particles, allocate each task >> locind, locbuffer, NumThis */

  selectposition ();		/*set cell position this task > ithis,jthis,kthis */

  setbordersize ();		/*set border lenght bsize */

  MPI_Barrier (MPI_COMM_WORLD);
  t2 = MPI_Wtime ();

  /*Transfer temporary boundary particles in the 6 directions */
  selectborder (0);/*-x*/
  selectborder (1);		/*+x */
  selectborder (2);/*-y*/
  selectborder (3);		/*+y */
  selectborder (4);/*-z*/
  selectborder (5);		/*+z */

  MPI_Barrier (MPI_COMM_WORLD);

  t3 = MPI_Wtime ();
  printf ("Task=%d NumThis=%d NumThis_boundary=%d\n", ThisTask, NumThis,
	  NumThisb);

  if (NumThisb > NumThis) printf ("%s",warning2);

  /*run vt */
#ifndef LOWMEMORY
  runvt ();
#else
  if ((ThisTask % 2) == 0)
    runvt ();			/*run half jobs first to save some memory, slower */
  MPI_Barrier (MPI_COMM_WORLD);
  if ((ThisTask % 2) == 1)
    runvt ();
#endif
  t4 = MPI_Wtime ();

  computeneig ();		/* Neighbors computation */
#ifdef VERBOSE
  if (ThisTask == 0)
    printf ("Neighbors complete.\n");
#endif
  computedens ();		/* Density computation */
#ifdef VERBOSE
  if (ThisTask == 0)
    printf ("Density and Volume complete.\n");
#endif

#ifdef COMPUTEGRAD
 computegrad ();   	        /*VT gradient computation*/
  if (ThisTask == 0)
    printf ("Gradient complete.\n");
#endif
  t5 = MPI_Wtime ();
  MPI_Barrier (MPI_COMM_WORLD);
  t6 = MPI_Wtime ();


#ifdef WRITEDENSITY
  writedensopen ();
  if (ThisTask == 0)
    printf ("write density...");
  savedata (0);			/*gather and save density */
  if (ThisTask == 0)
    printf ("volumes...");
  savedata (1);			/*gather and save volumes */
#endif

#ifdef COMPUTEGRAD
  if (ThisTask == 0)
    printf ("gradient...");
  savedata (3);                 /*gather and save gradient x */
  savedata (4);
  savedata (5);
#endif

#ifdef WRITEDENSITY
  writedensclose ();
#endif

#ifdef WRITENEIGHBORS
  if (ThisTask == 0)
    printf ("neighbors...\n");
  writeneigopen ();
  savedata (2);			/*gather and save neighbors */
  writeneigclose ();
#else
  printf ("\n");
#endif

#ifdef COMPUTEGRID
  if (checkgrid ())
    {
      writegridopen ();
      computegrid ();
      writegridclose ();
    }
#endif
  t7 = MPI_Wtime ();
  MPI_Barrier (MPI_COMM_WORLD);
  if (ThisTask == 0)
    printf
      ("Time Stats: Read=%.2f Bdry=%.2f Comp=%.2f Write=%.2f Tot-I/O=%.2F Tot=%.2f\n",
       t2 - t1, t3 - t2, t6 - t3, t7 - t6, t6 - t2, t7 - t1);
#ifdef VERBOSE
  MPI_Barrier (MPI_COMM_WORLD);
  printf ("Task=%d Time VT_only=%.2f\n", ThisTask, t4 - t3);
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  finalize ();			/*free memory and finalize code */
  return 0;
}				/* End Main */
