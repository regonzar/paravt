#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "config.h"
#include "vars.h"
#include "vt.h"

//inline static float sqr(float x) {
//        return x*x;
//    }

void
stopcode (void)
{
  stopit = 1;
  if (ThisTask == 0)
    printf ("Stop code...\n");
  MPI_Abort (MPI_COMM_WORLD, 0);

}

void
check_stop (void)
{
  MPI_Barrier (MPI_COMM_WORLD);	/*wait root task completion before check */
  MPI_Bcast (&stopit, 1, MPI_INT, root, MPI_COMM_WORLD);
  if (stopit)
    MPI_Abort (MPI_COMM_WORLD, 0);
}


void
check_params (int argc, char *argv[])
{
  int condition;

  if (NTask <= 1)
    {
      if (ThisTask == 0)
	printf
	  ("This code is intended for parallel run, using 1 processor may run slower than any other voronoi tessellation code due MPI overhead.\n");
    }

  if (argc < 3)
    {
      if (ThisTask == 0)
	printf
	  ("Input filename and format required: mpirun -np <Nproc> paravt <InputFilename> <fileformat>\n");
      stopcode ();
    }

  if (argc >= 3)
    {
      strcpy (InputFile, argv[1]);
      fileformat = 0;
      fileformat = atoi (argv[2]);
      if (ThisTask == 0)
	printf ("Filename=%s fileformat=%d\n", InputFile, fileformat);
    }


#ifdef BOX
  if (ThisTask == 0)
    printf ("BOX enabled.\n");
#else
  if (ThisTask == 0)
    printf ("BOX disabled.\n");
#endif

#ifdef PERIODIC
  if (ThisTask == 0)
    printf ("PERIODIC enabled.\n");
  if (condition == 1)
    condition = 2;
#else
  if (ThisTask == 0)
    printf ("PERIODIC disabled.\n");
#endif

#ifdef MULTIMASS
  if (ThisTask == 0)
    printf ("MULTIMASS enabled.\n");
#else
  if (ThisTask == 0)
    printf ("MULTIMASS disabled.\n");
#endif

#ifdef WRITEASCII
  if (ThisTask == 0)
    printf ("WRITEASCII enabled.\n");
#else
  if (ThisTask == 0)
    printf ("WRITEASCII disabled.\n");
#endif

#if !defined(BOX) && defined(PERIODIC)
  if (ThisTask == 0)
    printf ("Disabling PERIODIC because BOX is disabled.\n");
#undef PERIODIC
#endif

#if defined(COMPUTEGRAD) && !defined(WRITEDENSITY)
  if (ThisTask == 0)
    printf ("WRITEDENSITY MUST BE ENABLED IN COMPUTEGRAD OPTION, PLEASE CHANGE IT.\n");
#endif

}

void
splitdomain (void)
{
  int num, res, in, i;
  int split[3];
  div_t divres;

  if ((ThisTask == root) && (DOMAINTYPE == 0))
    {
      num = NTask;
      res = 0;
      in = 0;
      for (i = 0; i < 3; i++)
	split[i] = 1;
      while ((res == 0) && (num > 1))
	{
	  divres = div (num, 2);
	  num = divres.quot;
	  res = divres.rem;
	  if ((res == 0) && (num >= 1))
	    {
	      split[in] = split[in] * 2;
	      in++;
	      if (in > 2)
		in = 0;
	    }
	}
      if (split[0] * split[1] * split[2] != NTask)
	{
	  printf ("Not allowed number of tasks (allowed Ntask=2^N)\n");
	  stopcode ();
	}
      spx = split[0];
      spy = split[1];
      spz = split[2];
      printf ("indices = %d %d %d\n", spx, spy, spz);
    }


  if ((ThisTask == root) && (DOMAINTYPE == 1))	/*linear split along x (for irregular volumes) */
    {
      spx = NTask;
      spy = 1;
      spz = 1;
    }

  MPI_Bcast (&spx, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast (&spy, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast (&spz, 1, MPI_INT, root, MPI_COMM_WORLD);

}

void
initvars (void)
{
  locbuffer = NULL;
  locind = NULL;
  root = 0;			/*root process */
  NumThisb = 0;			/*boundary parts */
  stopit = 0;
  dx = 0;
  dy = 0;
  dz = 0;
}

void
finalize (void)
{

#ifdef WRITEVTSTRUCT
  if (ThisTask == 0)
    printf ("Finish.\n");
  MPI_Finalize ();
  exit (0);
#endif

  free (locind);
  free (pos);
#ifdef MULTIMASS
  free (mass);
#endif
  free (vol);
  free (rho);
#ifdef COMPUTEGRAD
  free (grad);
  free (neigborsg);
#endif
  free (neignum);
  free (neigbors);
  free (neigind);

  MPI_Finalize ();


}

void
selectposition (void)
{
  int ii, jj, kk, mm;

  for (ii = 0; ii < spx; ii++)
    {
      for (jj = 0; jj < spy; jj++)
	{
	  for (kk = 0; kk < spz; kk++)
	    {
	      mm = ii * spy * spz + jj * spz + kk;
	      if (mm == ThisTask)
		{
		  ithis = ii;
		  jthis = jj;
		  kthis = kk;
		}
	    }
	}
    }
}

void
setbordersize (void)
{
  float partdist;
#ifdef BOX
  partdist = pow (lbox * lbox * lbox / NumPart, 0.333);
  pdensity = NumPart/(lbox * lbox * lbox);
#else
  partdist = pow (lboxx * lboxy * lboxz / NumPart, 0.333);
  pdensity = NumPart/(lboxx * lboxy * lboxz);
#endif
  bsize = partdist * BORDERFACTOR;
#ifdef VERBOSE
  if (ThisTask == root)
    printf ("Border size=%f , density=%f\n", bsize,pdensity);
#endif
}

void
selectborder (int direction)
{
  float xmin1, ymin1, zmin1, xmax1, ymax1, zmax1;
  float thr,bvolume,bbsize;
  int i, count, TaskSend, Nsend, Nrec, tag, TaskRec;
  int ii, jj, kk, ii2, jj2, kk2;
  float offx, offy, offz;
  float *sendbuf, *recbuf;
  int *sendind, *recind;
  int bufrank;
  MPI_Request request, request2;
  MPI_Status status;

  xmin1 = ithis * dx;
  ymin1 = jthis * dy;
  zmin1 = kthis * dz;
  xmax1 = xmin1 + dx;
  ymax1 = ymin1 + dy;
  zmax1 = zmin1 + dz;

  if (direction == 0)
    thr = xmin1 + bsize;
  if (direction == 1)
    thr = xmax1 - bsize;
  if (direction == 2)
    thr = ymin1 + bsize;
  if (direction == 3)
    thr = ymax1 - bsize;
  if (direction == 4)
    thr = zmin1 + bsize;
  if (direction == 5)
    thr = zmax1 - bsize;
  count = 0;

  if ((direction == 0)||(direction == 1))
    bvolume=dy*dz*bsize;
  if ((direction == 2)||(direction == 3))
    bvolume=dx*dz*bsize;
  if ((direction == 4)||(direction == 5))
    bvolume=dx*dy*bsize;

#ifdef MULTIMASS
  bufrank = 4;
#else
  bufrank = 3;
#endif

  member = malloc (NumThis * sizeof (int));
  for (i = 0; i < NumThis; i++)
    {
      member[i] = 0;
      if (direction == 0)
	{
	  if (locbuffer[i * bufrank] < thr)
	    {
	      member[i] = 1;
	      count++;
	    }
	}
      if (direction == 1)
	{
	  if (locbuffer[i * bufrank] > thr)
	    {
	      member[i] = 1;
	      count++;
	    }
	}
      if (direction == 2)
	{
	  if (locbuffer[i * bufrank + 1] < thr)
	    {
	      member[i] = 1;
	      count++;
	    }
	}
      if (direction == 3)
	{
	  if (locbuffer[i * bufrank + 1] > thr)
	    {
	      member[i] = 1;
	      count++;
	    }
	}
      if (direction == 4)
	{
	  if (locbuffer[i * bufrank + 2] < thr)
	    {
	      member[i] = 1;
	      count++;
	    }
	}
      if (direction == 5)
	{
	  if (locbuffer[i * bufrank + 2] > thr)
	    {
	      member[i] = 1;
	      count++;
	    }
	}
    }

/*check if boundary region has enough particles*/
#ifndef AUTOBORDER
if (count/bvolume < MIN_BDENSITY*pdensity) 
  printf("WARNING! TOO FEW BORDER PARTICLES, increase BORDERFACTOR or change domain decomposition scheme.\n Task=%3d Direction=%1d Npthis=%d Nborder=%d density=%f pdensity=%f\n",ThisTask,direction,NumThis,count,count/bvolume,pdensity);
#endif

#ifdef AUTOBORDER
if (count/bvolume < MIN_BDENSITY*pdensity ) {
printf("AUTOBODER: increasing boundary for task %d, direction=%d\n boundary particles found=%d --- minimum expected=%.0f\n",ThisTask,direction,count,MIN_BDENSITY*pdensity*bvolume);
if (ThisTask==0){
  if (BORDERFACTOR < 1.0) printf("Warning: BORDERFACTOR in AUTOBORDER mode should be larger than 1.0!!!\n");
  if (MIN_BDENSITY > 1.0) printf("Warning: MIN_BDENSITY must be lower than 1.0!!!\n");
  }

bbsize=bsize;
while ((count < MIN_BDENSITY*pdensity*bvolume)&&(bbsize < 0.5*dx)) {
  bbsize=bbsize*1.2;
  if (direction == 0)
    thr = xmin1 + bbsize;
  if (direction == 1)
    thr = xmax1 - bbsize;
  if (direction == 2)
    thr = ymin1 + bbsize;
  if (direction == 3)
    thr = ymax1 - bbsize;
  if (direction == 4)
    thr = zmin1 + bbsize;
  if (direction == 5)
    thr = zmax1 - bbsize;
  count = 0;
  for (i = 0; i < NumThis; i++)
    {
      member[i] = 0;
      if (direction == 0)
        {
          if (locbuffer[i * bufrank] < thr)
            {
              member[i] = 1;
              count++;
            }
        }
      if (direction == 1)
        {
          if (locbuffer[i * bufrank] > thr)
            {
              member[i] = 1;
              count++;
            }
        }
      if (direction == 2)
        {
          if (locbuffer[i * bufrank + 1] < thr)
            {
              member[i] = 1;
              count++;
            }
        }
      if (direction == 3)
        {
          if (locbuffer[i * bufrank + 1] > thr)
            {
              member[i] = 1;
              count++;
            }
        }
      if (direction == 4)
        {
          if (locbuffer[i * bufrank + 2] < thr)
            {
              member[i] = 1;
              count++;
            }
        }
      if (direction == 5)
        {
          if (locbuffer[i * bufrank + 2] > thr)
            {
              member[i] = 1;
              count++;
            }
        }
    }

}/* end while*/
printf("New number of particles task %d: %d -- boundary expanded a factor of %.1f \n",ThisTask,count,bbsize/bsize);

} /*autoborder*/
#endif

/*Define who send/receive and Nparticles send/receive*/
  Nsend = count;		/*particles to send this task */
  ii = ithis;
  jj = jthis;
  kk = kthis;
  if (direction == 0)
    ii = ii - 1;
  if (direction == 1)
    ii = ii + 1;
  if (direction == 2)
    jj = jj - 1;
  if (direction == 3)
    jj = jj + 1;
  if (direction == 4)
    kk = kk - 1;
  if (direction == 5)
    kk = kk + 1;
  ii2 = ithis;
  jj2 = jthis;
  kk2 = kthis;
  if (direction == 0)
    ii2 = ii2 + 1;
  if (direction == 1)
    ii2 = ii2 - 1;
  if (direction == 2)
    jj2 = jj2 + 1;
  if (direction == 3)
    jj2 = jj2 - 1;
  if (direction == 4)
    kk2 = kk2 + 1;
  if (direction == 5)
    kk2 = kk2 - 1;
#if !defined(PERIODIC)
  if (ii < 0)
    Nsend = 0;
  if (jj < 0)
    Nsend = 0;
  if (kk < 0)
    Nsend = 0;
  if (ii == spx)
    Nsend = 0;
  if (jj == spy)
    Nsend = 0;
  if (kk == spz)
    Nsend = 0;
#endif
  offx = 0;
  offy = 0;
  offz = 0;

#ifdef BOX
  if (ii < 0)
    {
      ii = spx - 1;
      offx = lbox;
    }
  if (jj < 0)
    {
      jj = spy - 1;
      offy = lbox;
    }
  if (kk < 0)
    {
      kk = spz - 1;
      offz = lbox;
    }
  if (ii == spx)
    {
      ii = 0;
      offx = -lbox;
    }
  if (jj == spy)
    {
      jj = 0;
      offy = -lbox;
    }
  if (kk == spz)
    {
      kk = 0;
      offz = -lbox;
    }
#else
  if (ii < 0)
    {
      ii = spx - 1;
      offx = lboxx;
    }
  if (jj < 0)
    {
      jj = spy - 1;
      offy = lboxy;
    }
  if (kk < 0)
    {
      kk = spz - 1;
      offz = lboxz;
    }
  if (ii == spx)
    {
      ii = 0;
      offx = -lboxx;
    }
  if (jj == spy)
    {
      jj = 0;
      offy = -lboxy;
    }
  if (kk == spz)
    {
      kk = 0;
      offz = -lboxz;
    }
#endif

  if (ii2 < 0)
    ii2 = spx - 1;
  if (jj2 < 0)
    jj2 = spy - 1;
  if (kk2 < 0)
    kk2 = spz - 1;
  if (ii2 == spx)
    ii2 = 0;
  if (jj2 == spy)
    jj2 = 0;
  if (kk2 == spz)
    kk2 = 0;

  TaskSend = ii * spy * spz + jj * spz + kk;	/*task to send particles */
  TaskRec = ii2 * spy * spz + jj2 * spz + kk2;	/* task from which this task will receive parts */
  tag = 0;			/*packet identifier */
  Nrec = 0;			/*reset Nrec=0 */
  MPI_Send (&Nsend, 1, MPI_INT, TaskSend, tag, MPI_COMM_WORLD);
  MPI_Recv (&Nrec, 1, MPI_INT, TaskRec, tag, MPI_COMM_WORLD,
	    MPI_STATUS_IGNORE);

/*printf("Task=%3d Npthis=%d Nsend=%d Nrec=%d Send=%d Rec=%d\n",ThisTask,NumThis,Nsend,Nrec,TaskSend,TaskRec);*/

/*create send buffer, correct offset*/

  if (Nsend > 0)
    {
      sendbuf = malloc (Nsend * sizeof (float) * bufrank);
      sendind = malloc (Nsend * sizeof (int));
      ii = 0;
      for (i = 0; i < NumThis; i++)
	{
	  if (member[i] == 1)
	    {
	      sendbuf[ii * bufrank] = locbuffer[i * bufrank] + offx;
	      sendbuf[ii * bufrank + 1] = locbuffer[i * bufrank + 1] + offy;
	      sendbuf[ii * bufrank + 2] = locbuffer[i * bufrank + 2] + offz;
#ifdef MULTIMASS
	      sendbuf[ii * bufrank + 3] = locbuffer[i * bufrank + 3];
#endif
	      sendind[ii] = locind[i];
	      ii++;
	    }
	}
    }

/*send particles asynchronous. using mpi_send block the buffer if large or many send requests*/

  if (Nrec > 0)
    recbuf = malloc (Nrec * sizeof (float) * bufrank);
  if (Nrec > 0)
    MPI_Irecv (recbuf, Nrec * bufrank, MPI_FLOAT, TaskRec, 0, MPI_COMM_WORLD,
	       &request);
  if (Nsend > 0)
    MPI_Isend (sendbuf, Nsend * bufrank, MPI_FLOAT, TaskSend, 0,
	       MPI_COMM_WORLD, &request2);
  if (Nrec > 0)
    MPI_Wait (&request, &status);
  if (Nsend > 0)
    MPI_Wait (&request2, &status);
  if (Nrec > 0)
    recind = malloc (Nrec * sizeof (int));
  if (Nrec > 0)
    MPI_Irecv (recind, Nrec, MPI_INT, TaskRec, 0, MPI_COMM_WORLD, &request);
  if (Nsend > 0)
    MPI_Isend (sendind, Nsend, MPI_INT, TaskSend, 0, MPI_COMM_WORLD,
	       &request2);
  if (Nrec > 0)
    MPI_Wait (&request, &status);
  if (Nsend > 0)
    MPI_Wait (&request2, &status);

/*append received particles*/
  if (Nrec > 0)
    {
      locbuffer =
	realloc (locbuffer,
		 (NumThis + NumThisb + Nrec) * sizeof (float) * bufrank);
      locind = realloc (locind, (NumThis + NumThisb + Nrec) * sizeof (int));
      ii = NumThis + NumThisb;
      for (i = 0; i < Nrec; i++)
	{
	  locbuffer[(ii + i) * bufrank] = recbuf[i * bufrank];
	  locbuffer[(ii + i) * bufrank + 1] = recbuf[i * bufrank + 1];
	  locbuffer[(ii + i) * bufrank + 2] = recbuf[i * bufrank + 2];
#ifdef MULTIMASS
	  locbuffer[(ii + i) * bufrank + 3] = recbuf[i * bufrank + 3];
#endif
	  locind[(ii + i)] = recind[i];
	}
    }

  NumThisb = NumThisb + Nrec;
  MPI_Reduce (&Nrec, &i, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/*only to be sure all complete before another send */

  if (Nrec > 0)
    {
      free (recbuf);
      free (recind);
    }
  if (Nsend > 0)
    {
      free (sendbuf);
      free (sendind);
    }
  free (member);
}

/*GRID FUNCTIONS*/
int
checkgrid (void)		/*check conditions to compute grid */
{
  int rr;
  rr = 1;
  if ((GRIDSIZE % spx) != 0)
    rr = 0;
  if ((GRIDSIZE % spy) != 0)
    rr = 0;
  if ((GRIDSIZE % spz) != 0)
    rr = 0;
  if (rr == 0)
    printf
      ("Density Grid computation skipped. It requires: \n i)a rectangular volume.\n ii)The number of tasks can be decomposed in NTask=spx*spy*spz where spx,spy,spz are integers.\n iii) (GRIDSIZE mod spx=0) (GRIDSIZE mod spy=0) (GRIDSIZE mod spz=0)\n");
  return rr;
}


void
computegrid (void)
{
  int nx, ny, nz, nn, i, nzero;
  int i1, j1, k1, ind;
  float offx, offy, offz;
  float xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
  float *masst, *volt;
  float dxx, dyy, dzz, ftemp;
  float *page, *gpage, *rpage;
  float *page2, *gpage2, *rpage2;
  int k, kz, npages, j, tag, rtask, in2, in1, ii, jj;
  MPI_Request request, request2;
  MPI_Status status;


  nx = (int) (GRIDSIZE / spx);
  ny = (int) (GRIDSIZE / spy);
  nz = (int) (GRIDSIZE / spz);
  if (ThisTask == root)
    printf ("Grid computation GRIDSIZE=%d nx=%d ny=%d nz=%d\n", GRIDSIZE, nx,
	    ny, nz);

  offx = ithis * dx;
  offy = jthis * dy;
  offz = kthis * dz;
  /*printf ("Task=%d ithis=%d jthis=%d kthis=%d offx=%f offy=%f offz=%f dx=%f\n",
     ThisTask, ithis, jthis, kthis, offx, offy, offz,dx); */
  dxx = dx / nx;
  dyy = dy / ny;
  dzz = dz / nz;

  xmin2 = ithis * dx;
  ymin2 = jthis * dy;
  zmin2 = kthis * dz;
  xmax2 = xmin2 + dx;
  ymax2 = ymin2 + dy;
  zmax2 = zmin2 + dz;
#ifdef VERBOSE
  printf ("Grid. Task=%d [%f %f] [%f %f] [%f %f]\n", ThisTask, xmin2, xmax2,
	  ymin2, ymax2, zmin2, zmax2);
#endif

  nn = nx * ny * nz;
  masst = malloc (nn * sizeof (float));
  for (i = 0; i < nn; i++)
    masst[i] = 0;
  volt = malloc (nn * sizeof (float));
  for (i = 0; i < nn; i++)
    volt[i] = 0;

  for (i = 0; i < NumThis; i++)
    {
      i1 = (int) floor ((pos[3 * i] - offx) / dxx);
      j1 = (int) floor ((pos[3 * i + 1] - offy) / dyy);
      k1 = (int) floor ((pos[3 * i + 2] - offz) / dzz);
      if (i1 == nx)
	i1 = nx - 1;
      if (j1 == ny)
	j1 = ny - 1;
      if (k1 == nz)
	k1 = nz - 1;
      ind = i1 * ny * nz + j1 * nz + k1;
      volt[ind] += vol[i];
#ifdef MULTIMASS
      masst[ind] += mass[i];
#else
      masst[ind] += 1.0;
#endif
    }

/*gather pages in z direction ordered and write*/
  rpage = malloc (nx * ny * sizeof (float));
  page = malloc (nx * ny * sizeof (float));
  if (ThisTask == root)
    gpage = malloc (GRIDSIZE * GRIDSIZE * sizeof (float));
  rpage2 = malloc (nx * ny * sizeof (float));
  page2 = malloc (nx * ny * sizeof (float));
  if (ThisTask == root)
    gpage2 = malloc (GRIDSIZE * GRIDSIZE * sizeof (float));

  for (k = 0; k < GRIDSIZE; k++)
    {
      kz = (int) floor (k / nz);	/*if kthis==kz this task send */
      npages = spx * spy;	/*number pages should be send */
      if (kthis == kz)
	{
	  i1 = 0;
	  for (i = 0; i < nx; i++)
	    {
	      for (j = 0; j < ny; j++)
		{
		  ind = i * ny * nz + j * nz + k;
		  page[i1] = volt[ind];
		  page2[i1] = masst[ind];
		  i1++;
		}		/*j */
	    }			/*i */
	  tag = ThisTask;
	  MPI_Isend (page, nx * ny, MPI_FLOAT, root, tag, MPI_COMM_WORLD,
		     &request);
	  MPI_Isend (page2, nx * ny, MPI_FLOAT, root, tag, MPI_COMM_WORLD,
		     &request2);
	}			/*kthis */

      if (ThisTask == root)
	{
	  for (i = 0; i < spx; i++)
	    {
	      for (j = 0; j < spy; j++)
		{
		  rtask = i * spy * spz + j * spz + kz;
		  MPI_Irecv (rpage, nx * ny, MPI_FLOAT, rtask, rtask,
			     MPI_COMM_WORLD, &request);
		  MPI_Wait (&request, &status);
		  MPI_Irecv (rpage2, nx * ny, MPI_FLOAT, rtask, rtask,
			     MPI_COMM_WORLD, &request2);
		  MPI_Wait (&request2, &status);
		  in2 = 0;
		  for (ii = 0; ii < nx; ii++)
		    {
		      for (jj = 0; jj < ny; jj++)
			{
			  in1 = (ii + nx * i) * GRIDSIZE + (jj + ny * j);
			  if (in1 > GRIDSIZE * GRIDSIZE - 1)
			    printf
			      ("FATAL GRID Error...(disable grid output)\n");
			  gpage[in1] = rpage[in2];
			  gpage2[in1] = rpage2[in2];
			  in2++;
			}
		    }
		}		/*j */
	    }			/*i */

	  /*write gpage - x swapped z */

	  for (i = 0; i < GRIDSIZE * GRIDSIZE; i++)
	    {
	      ftemp = log10 (minrho);	/*Empty cells are assigned with the lowest density */
	      if (gpage[i] > 0)
		ftemp = log10 (gpage2[i] / gpage[i]);
	      fprintf (filegrid2, "%g\n", ftemp);
	      fwrite (&ftemp, sizeof (float), 1, filegrid);
	    }
	  /*end write gpage=vol gpage2=mass */
	}			/*root */

      MPI_Barrier (MPI_COMM_WORLD);	/*sync at each page to avoid errors */

    }				/*k */

  free (page);
  free (rpage);
  if (ThisTask == root)
    free (gpage);
  free (page2);
  free (rpage2);
  if (ThisTask == root)
    free (gpage2);
  free (masst);
  free (volt);
#ifdef VERBOSE
  if (ThisTask == root)
    printf ("Grid complete.\n");
#endif
}

int linreg(int n, const float x[], const float y[], double* m, double* b, double* r)
{
 int i;
 double denom;

 double   sumx = 0.0;                        /* sum of x                      */
 double   sumx2 = 0.0;                       /* sum of x**2                   */
 double   sumxy = 0.0;                       /* sum of x * y                  */
 double   sumy = 0.0;                        /* sum of y                      */
 double   sumy2 = 0.0;                       /* sum of y**2                   */

       for (i=0;i<n;i++)   
          { 
          sumx  += x[i];       
          sumx2 += sqr(x[i]);  
          sumxy += x[i] * y[i];
          sumy  += y[i];      
          sumy2 += sqr(y[i]); 
          } 

       denom = (n * sumx2 - sqr(sumx));
       if (denom == 0) {
            *m = 0;
            *b = 0;
            if (r) *r = 0;
            return 0;
          }
          *m = (n * sumxy  -  sumx * sumy) / denom;
          *b = (sumy * sumx2  -  sumx * sumxy) / denom;
          if (r!=NULL) {
             *r = (sumxy - sumx * sumy / n) / sqrt((sumx2 - sqr(sumx)/n) * (sumy2 - sqr(sumy)/n)); 
             }
          return 1;         
             
}                                                
