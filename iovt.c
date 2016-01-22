#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "config.h"
#include "vars.h"
#include "func.h"
#include "iovt.h"
#include "vt.h"
#include <float.h>


/*Generic read functions*/
void
read_header (void)
{
  if (fileformat == 0)
    read_header_boxascii ();

  if (fileformat == 1)
    read_header_bin ();

  if (fileformat == 2)
    read_header_gadget ();

  if (fileformat == 3)
    read_header_rockstar ();
}

void
read_part (void)
{
  if (fileformat == 0)
    read_part_boxascii ();	/*ascii file */

  if (fileformat == 1)
    read_part_bin ();		/*binary format */

  if (fileformat == 2)
    read_part_gadget ();

  if (fileformat == 3)
    read_part_rockstar ();

}

void
read_header_boxascii (void)
{
  if (ThisTask == root)
    {
      if ((filein = fopen (InputFile, "r")) == NULL)
	{
	  printf ("Fail open file...\n");
	  stopcode ();
	}
#ifdef BOX
      fscanf (filein, "%d %f\n", &NumPart, &lbox);
      printf ("NumPart=%d  Lbox=%8.3f\n", NumPart, lbox);
      lboxx = 0;
      lboxy = 0;
      lboxz = 0;
#else
      fscanf (filein, "%d %f %f %f\n", &NumPart, &lboxx, &lboxy, &lboxz);
      printf ("NumPart=%d  lbox_x=%8.3f lbox_y=%8.3f lboxz=%8.3f\n", NumPart,
	      lboxx, lboxy, lboxz);
      lbox = 0;
#endif
    }
  check_stop ();
  MPI_Bcast (&NumPart, 1, MPI_INT, root, MPI_COMM_WORLD);
#ifdef BOX
  MPI_Bcast (&lbox, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
#else
  MPI_Bcast (&lboxx, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
  MPI_Bcast (&lboxy, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
  MPI_Bcast (&lboxz, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
#endif
}

void
read_header_bin (void)
{
  if (ThisTask == root)
    {
      if ((filein = fopen (InputFile, "rb")) == NULL)
	{
	  printf ("Fail open file...\n");
	  stopcode ();
	}
#ifdef BOX
      fread (&NumPart, sizeof (int), 1, filein);
      fread (&lbox, sizeof (float), 1, filein);
      printf ("NumPart=%d  Lbox=%8.3f\n", NumPart, lbox);
      lboxx = 0;
      lboxy = 0;
      lboxz = 0;
#else
      printf ("only BOX shaped binary input file supported.\n");
      stopcode ();
#endif
    }
  check_stop ();
  MPI_Bcast (&NumPart, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast (&lbox, 1, MPI_FLOAT, root, MPI_COMM_WORLD);

}


void
read_header_gadget (void)
{
  char filenew[300];
  int blk, blk2, i;

#define SKIP {fread(&blk,sizeof(int),1,filein);}
#define SKIP2 {fread(&blk2,sizeof(int),1,filein);}

#ifdef MULTIMASS
  if (ThisTask == root)
    printf ("MULTIMASS NOT AVAILABLE FOR GADGET FILES.\n");
  stopcode ();
#endif

  if (ThisTask == root)
    {
#ifdef GADGET_PTYPE
      printf ("Reading Gadget file with particle type %d ...\n",
	      GADGET_PTYPE);
#else
      printf ("Gadget format required define GADGET_PTYPE in config.h\n");
      stopcode ();
#endif

#ifndef BOX
      printf ("Gadget format only compatible with BOX flag...\n");
      stopcode ();
#endif

      if ((filein = fopen (InputFile, "r")) == NULL)
	{
	  printf ("Fail open file: %s\n", InputFile);
	  sprintf (filenew, "%s.0", InputFile);
	  if ((filein = fopen (filenew, "r")) == NULL)
	    {
	      printf ("Try Multiple file.. Fail open file: %s\n", filenew);
	      stopcode ();
	    }
	}

      SKIP;
      fread (&header, sizeof (header), 1, filein);
      SKIP2;

      if (blk != 256 || blk2 != 256)
	{
	  printf ("incorrect header format\n");
	  stopcode ();
	}

      NumPart = header.npartTotal[GADGET_PTYPE];
      lbox = header.BoxSize;

      printf ("NpartTotal=");
      for (i = 0; i < 6; i++)
	printf ("%d ", header.npartTotal[i]);
      printf ("\n");

      printf ("Numpart=%d  lbox=%8.3f num_files=%d\n", NumPart, lbox,
	      header.num_files);
      lboxx = 0;
      lboxy = 0;
      lboxz = 0;
      fclose (filein);
    }
  check_stop ();
  MPI_Bcast (&NumPart, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast (&lbox, 1, MPI_FLOAT, root, MPI_COMM_WORLD);

}

void
read_header_rockstar (void)
{
#define SKIPLINE { ch=' '; while (ch !='\n') ch=getc(filein);}

  int ngatos, nlines, i;
  float ssc;
  char ch, chlin[300];
#ifndef BOX
  printf ("Rockstar format only available for BOX condition.\n");
  stopcode ();
#endif
#ifndef MULTIMASS
  printf ("Rockstar format requires MULTIMASS.\n");
  stopcode ();
#endif

  if (ThisTask == root)
    {
      if ((filein = fopen (InputFile, "r")) == NULL)
	{
	  printf ("Fail open file...\n");
	  stopcode ();
	}
      nlines = 0;
      ngatos = 0;
      while ((ch = getc (filein)) != EOF)
	{
	  if (ch == '\n')
	    nlines++;
	  if (ch == '#')
	    ngatos++;
	}
      NumPart = nlines - ngatos;
      rewind (filein);
      SKIPLINE;
      SKIPLINE;
      i = 0;
      ch = ' ';
      while (ch != '\n')
	{
	  ch = getc (filein);
	  if (i >= 16)
	    chlin[i - 16] = ch;
	  i++;
	}
      for (i = 0; i < ngatos - 1 - 2; i++)
	SKIPLINE;
      sscanf (chlin, "%f", &lbox);
      printf
	("Rockstar File, nlines=%d headerlines=%d NumHalos=%d Lbox=%f\n",
	 nlines, ngatos, NumPart, lbox);

      printf ("Numpart=%d  lbox=%8.3f\n", NumPart, lbox);

    }
  lboxx = 0;
  lboxy = 0;
  lboxz = 0;

  check_stop ();
  MPI_Bcast (&NumPart, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast (&lbox, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
}



/*read particles and assign to each task  */
void
read_part_boxascii (void)
{
  int i, count, ii, jj, kk;
  float t0, t1, t2, t3;
  int bufrank;
#ifdef MULTIMASS
  bufrank = 4;
#else
  bufrank = 3;
#endif
  globalcount = 0;		/*read particles */
  globalindex = 0;		/*current assigned particles this */
  loopcount = 0;		/*number blocks reads */
  while (globalcount < NumPart)
    {
      localindex = 0;		/*part in this task in this buffer loop */
      if ((NumPart - globalcount) > BUFFERSIZE)
	currsize = BUFFERSIZE;
      if ((NumPart - globalcount) <= BUFFERSIZE)
	currsize = NumPart - globalcount;
      buffer = malloc (currsize * sizeof (float) * bufrank);
      member = malloc (currsize * sizeof (int));
      if (ThisTask == root)
	{
	  i = 0;
	  while (i < currsize)
	    {
#ifdef MULTIMASS
	      fscanf (filein, "%f %f %f %g", &t0, &t1, &t2, &t3);	/*HERE you can change ascii read format */
	      buffer[i * 4 + 3] = t3;
	      buffer[i * 4 + 0] = t0;
	      buffer[i * 4 + 1] = t1;
	      buffer[i * 4 + 2] = t2;
#else
	      fscanf (filein, "%f %f %f%*[^\n]", &t0, &t1, &t2);	/*skip remainder of lines.  (3 cols) */
	      buffer[i * 3 + 0] = t0;
	      buffer[i * 3 + 1] = t1;
	      buffer[i * 3 + 2] = t2;
#endif
	      i++;
	    }
	}
      count = currsize * bufrank;
      MPI_Bcast (buffer, count, MPI_FLOAT, root, MPI_COMM_WORLD);
#ifdef BOX
      dx = lbox / spx;
      dy = lbox / spy;
      dz = lbox / spz;
#else
      dx = lboxx / spx;
      dy = lboxy / spy;
      dz = lboxz / spz;
#endif
      for (i = 0; i < currsize; i++)
	{
	  if (DOMAINTYPE == 0)
	    {
	      ii = (int) floor (buffer[i * bufrank] / dx);
	      jj = (int) floor (buffer[i * bufrank + 1] / dy);
	      kk = (int) floor (buffer[i * bufrank + 2] / dz);
	      if (ii == spx)
		ii = spx - 1;
	      if (jj == spy)
		jj = spy - 1;
	      if (kk == spz)
		kk = spz - 1;
	      member[i] = ii * spy * spz + jj * spz + kk;	/*grid assign */
	      if (member[i] == ThisTask)
		localindex++;
	    }
	  if (DOMAINTYPE == 1)
	    {
	      ii = (int) floor (buffer[i * bufrank] / dx);
	      if (ii == spx)
		ii = spx - 1;
	      member[i] = ii;
	      if (member[i] == ThisTask)
		localindex++;
	    }
	}
      if (localindex > 0)
	{
	  if (locbuffer == NULL)
	    {
	      locbuffer = malloc (localindex * sizeof (float) * bufrank);
	    }
	  else
	    {
	      locbuffer =
		realloc (locbuffer,
			 (globalindex +
			  localindex) * sizeof (float) * bufrank);
	    }
	  if (locind == NULL)
	    {
	      locind = malloc (localindex * sizeof (int));
	    }
	  else
	    {
	      locind =
		realloc (locind, (globalindex + localindex) * sizeof (int));
	    }
	  count = 0;
	  for (i = 0; i < currsize; i++)
	    {
	      if (member[i] == ThisTask)
		{
		  locind[globalindex + count] = globalcount + i;

		  locbuffer[(globalindex + count) * bufrank] =
		    buffer[i * bufrank];
		  locbuffer[(globalindex + count) * bufrank + 1] =
		    buffer[i * bufrank + 1];
		  locbuffer[(globalindex + count) * bufrank + 2] =
		    buffer[i * bufrank + 2];
#ifdef MULTIMASS
		  locbuffer[(globalindex + count) * bufrank + 3] =
		    buffer[i * bufrank + 3];
#endif
		  count++;
		}
	    }
	}
      free (buffer);
      free (member);
      globalindex = globalindex + localindex;
      globalcount = globalcount + currsize;
      NumThis = globalindex;
      loopcount++;
    }				/*end loop */
  MPI_Reduce (&NumThis, &i, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
#ifdef VERBOSE
  printf ("Task=%3d Npthistot=%9d\n", ThisTask, NumThis);
#endif
  MPI_Barrier (MPI_COMM_WORLD);	/*ensure all particles are allocated */
  if (ThisTask == root)
    {
      printf ("Read Complete: Nread=%d in %d blocks of %d bytes.\n", i,
	      loopcount, BUFFERSIZE);
      fclose (filein);
    }
}



/*read particles */
void
read_part_rockstar (void)
{
  int i, j, count, ii, jj, kk;
  float t0, t1, t2, t3, ssc;
  int bufrank;
  char ch;
  bufrank = 4;			/*multimass for rockstar */
  globalcount = 0;		/*read particles */
  globalindex = 0;		/*current assigned particles this */
  loopcount = 0;		/*number blocks reads */
  while (globalcount < NumPart)
    {
      localindex = 0;		/*part in this task in this buffer loop */
      if ((NumPart - globalcount) > BUFFERSIZE)
	currsize = BUFFERSIZE;
      if ((NumPart - globalcount) <= BUFFERSIZE)
	currsize = NumPart - globalcount;
      buffer = malloc (currsize * sizeof (float) * bufrank);
      member = malloc (currsize * sizeof (int));
      if (ThisTask == root)
	{
	  i = 0;
	  while (i < currsize)
	    {
	      for (j = 0; j < 20; j++)
		{
		  fscanf (filein, "%f ", &ssc);
		  if (j == 10)
		    t3 = ssc;
		  if (j == 17)
		    t0 = ssc;
		  if (j == 18)
		    t1 = ssc;
		  if (j == 19)
		    t2 = ssc;
		}
	      SKIPLINE;
	      buffer[i * 4 + 3] = t3;
	      buffer[i * 4 + 0] = t0;
	      buffer[i * 4 + 1] = t1;
	      buffer[i * 4 + 2] = t2;
	      i++;
	    }
	}
      count = currsize * bufrank;
      MPI_Bcast (buffer, count, MPI_FLOAT, root, MPI_COMM_WORLD);
      dx = lbox / spx;
      dy = lbox / spy;
      dz = lbox / spz;
      for (i = 0; i < currsize; i++)
	{
	  if (DOMAINTYPE == 0)
	    {
	      ii = (int) floor (buffer[i * bufrank] / dx);
	      jj = (int) floor (buffer[i * bufrank + 1] / dy);
	      kk = (int) floor (buffer[i * bufrank + 2] / dz);
	      if (ii == spx)
		ii = spx - 1;
	      if (jj == spy)
		jj = spy - 1;
	      if (kk == spz)
		kk = spz - 1;
	      member[i] = ii * spy * spz + jj * spz + kk;	/*grid assign */
	      if (member[i] == ThisTask)
		localindex++;
	    }
	  if (DOMAINTYPE == 1)
	    {
	      ii = (int) floor (buffer[i * bufrank] / dx);
	      if (ii == spx)
		ii = spx - 1;
	      member[i] = ii;
	      if (member[i] == ThisTask)
		localindex++;
	    }
	}
      if (localindex > 0)
	{
	  if (locbuffer == NULL)
	    {
	      locbuffer = malloc (localindex * sizeof (float) * bufrank);
	    }
	  else
	    {
	      locbuffer =
		realloc (locbuffer,
			 (globalindex +
			  localindex) * sizeof (float) * bufrank);
	    }
	  if (locind == NULL)
	    {
	      locind = malloc (localindex * sizeof (int));
	    }
	  else
	    {
	      locind =
		realloc (locind, (globalindex + localindex) * sizeof (int));
	    }
	  count = 0;
	  for (i = 0; i < currsize; i++)
	    {
	      if (member[i] == ThisTask)
		{
		  locind[globalindex + count] = globalcount + i;

		  locbuffer[(globalindex + count) * bufrank] =
		    buffer[i * bufrank];
		  locbuffer[(globalindex + count) * bufrank + 1] =
		    buffer[i * bufrank + 1];
		  locbuffer[(globalindex + count) * bufrank + 2] =
		    buffer[i * bufrank + 2];
		  locbuffer[(globalindex + count) * bufrank + 3] =
		    buffer[i * bufrank + 3];
		  count++;
		}
	    }
	}
      free (buffer);
      free (member);
      globalindex = globalindex + localindex;
      globalcount = globalcount + currsize;
      NumThis = globalindex;
      loopcount++;
      /*if (ThisTask==root)printf("gi=%d gc=%d loop=%d currsize=%d dx=%f dy=%f dz=%f\n",globalindex,globalcount,loopcount,currsize,dx,dy,dz); */
    }				/*end loop */
  MPI_Reduce (&NumThis, &i, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
#ifdef VERBOSE
  printf ("Task=%d Npthistot=%d\n", ThisTask, NumThis);
#endif
  MPI_Barrier (MPI_COMM_WORLD);	/*ensure all particles are allocated */
  if (ThisTask == root)
    {
      printf ("Read Complete. Nread=%d in %d blocks of %d bytes.\n", i,
	      loopcount, BUFFERSIZE);
      fclose (filein);
    }
}





/*read binary*/
void
read_part_bin (void)
{
  int i, count, ii, jj, kk;
  float t0, t1, t2, t3;
  int bufrank;
#ifdef MULTIMASS
  bufrank = 4;
#else
  bufrank = 3;
#endif
  globalcount = 0;		/*read particles */
  globalindex = 0;		/*current assigned particles this */
  loopcount = 0;		/*number blocks reads */
  while (globalcount < NumPart)
    {
      localindex = 0;		/*part in this task in this buffer loop */
      if ((NumPart - globalcount) > BUFFERSIZE)
	currsize = BUFFERSIZE;
      if ((NumPart - globalcount) <= BUFFERSIZE)
	currsize = NumPart - globalcount;
      buffer = malloc (currsize * sizeof (float) * bufrank);
      member = malloc (currsize * sizeof (int));
      if (ThisTask == root)
	{
	  i = 0;
	  while (i < currsize)
	    {
	      fread (&t0, sizeof (float), 1, filein);
	      fread (&t1, sizeof (float), 1, filein);
	      fread (&t2, sizeof (float), 1, filein);
	      buffer[i * bufrank + 0] = t0;
	      buffer[i * bufrank + 1] = t1;
	      buffer[i * bufrank + 2] = t2;
#ifdef MULTIMASS
	      fread (&t3, sizeof (float), 1, filein);
	      buffer[i * bufrank + 3] = t3;
#endif
	      i++;
	    }
	}
      count = currsize * bufrank;
      MPI_Bcast (buffer, count, MPI_FLOAT, root, MPI_COMM_WORLD);
      dx = lbox / spx;
      dy = lbox / spy;
      dz = lbox / spz;
      for (i = 0; i < currsize; i++)
	{
	  if (DOMAINTYPE == 0)
	    {
	      ii = (int) floor (buffer[i * bufrank] / dx);
	      jj = (int) floor (buffer[i * bufrank + 1] / dy);
	      kk = (int) floor (buffer[i * bufrank + 2] / dz);
	      if (ii == spx)
		ii = spx - 1;
	      if (jj == spy)
		jj = spy - 1;
	      if (kk == spz)
		kk = spz - 1;
	      member[i] = ii * spy * spz + jj * spz + kk;	/*grid assign */
	      if (member[i] == ThisTask)
		localindex++;
	    }
	  if (DOMAINTYPE == 1)
	    {
	      ii = (int) floor (buffer[i * bufrank] / dx);
	      if (ii == spx)
		ii = spx - 1;
	      member[i] = ii;
	      if (member[i] == ThisTask)
		localindex++;
	    }
	}
      if (localindex > 0)
	{
	  if (locbuffer == NULL)
	    {
	      locbuffer = malloc (localindex * sizeof (float) * bufrank);
	    }
	  else
	    {
	      locbuffer =
		realloc (locbuffer,
			 (globalindex +
			  localindex) * sizeof (float) * bufrank);
	    }
	  if (locind == NULL)
	    {
	      locind = malloc (localindex * sizeof (int));
	    }
	  else
	    {
	      locind =
		realloc (locind, (globalindex + localindex) * sizeof (int));
	    }
	  count = 0;
	  for (i = 0; i < currsize; i++)
	    {
	      if (member[i] == ThisTask)
		{
		  locind[globalindex + count] = globalcount + i;

		  locbuffer[(globalindex + count) * bufrank] =
		    buffer[i * bufrank];
		  locbuffer[(globalindex + count) * bufrank + 1] =
		    buffer[i * bufrank + 1];
		  locbuffer[(globalindex + count) * bufrank + 2] =
		    buffer[i * bufrank + 2];
#ifdef MULTIMASS
		  locbuffer[(globalindex + count) * bufrank + 3] =
		    buffer[i * bufrank + 3];
#endif
		  count++;
		}
	    }
	}
      free (buffer);
      free (member);
      globalindex = globalindex + localindex;
      globalcount = globalcount + currsize;
      NumThis = globalindex;
      loopcount++;
    }				/*end loop */
  MPI_Reduce (&NumThis, &i, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
#ifdef VERBOSE
  printf ("Task=%3d Npthistot=%9d\n", ThisTask, NumThis);
#endif
  MPI_Barrier (MPI_COMM_WORLD);	/*ensure all particles are allocated */
  if (ThisTask == root)
    {
      printf ("Read Complete. Nread=%d in %d blocks of %d bytes.\n", i,
	      loopcount, BUFFERSIZE);
      fclose (filein);
    }
}

void
read_part_gadget (void)
{
  int i, count, ii, jj, kk, j;
  float t0, t1, t2, t3;
  int blk, blk2, hsize, loopsize;
  float *rbuf;
  int nfiles, numthisfile, typeskip, globalskip;
  unsigned int numtotal;
  char filenew[300];


  if (ThisTask == root)
    {
      nfiles = header.num_files;
      numtotal = header.npartTotal[GADGET_PTYPE];
      printf ("Number files=%d numtotal=%d\n", nfiles, numtotal);
    }
  MPI_Bcast (&nfiles, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast (&numtotal, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);


  globalcount = 0;		/*read particles */
  globalindex = 0;		/*current assigned particles this */


  for (j = 0; j < nfiles; j++)
    {
      globalskip = 0;
      sprintf (filenew, "%s", InputFile);
      if (nfiles > 1)
	sprintf (filenew, "%s.%d", InputFile, j);

      if (ThisTask == root)
	{
	  if ((filein = fopen (filenew, "r")) == NULL)
	    {
	      printf ("Fail open file: %s\n", filenew);
	      stopcode ();
	    }

	  SKIP;
	  fread (&header, sizeof (header), 1, filein);
	  SKIP2;

	  if (blk != 256 || blk2 != 256)
	    {
	      printf ("incorrect header format\n");
	      stopcode ();
	    }
	  numthisfile = header.npart[GADGET_PTYPE];
	  typeskip = 0;		/*seek due particle types */
	  for (i = 0; i < GADGET_PTYPE; i++)
	    typeskip += header.npart[i];	/*particles to skip */
	  printf ("Reading file=%d numthisfile=%d skiptype=%d\n", j,
		  numthisfile, typeskip);
	}			/*root */

      MPI_Bcast (&numthisfile, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast (&typeskip, 1, MPI_INT, root, MPI_COMM_WORLD);

      loopcount = 0;		/*number blocks reads */
      while (globalskip < numthisfile)
	{
	  localindex = 0;	/*part in this task in this buffer loop */
	  if ((numthisfile - globalskip) > BUFFERSIZE)
	    currsize = BUFFERSIZE;
	  if ((numthisfile - globalskip) <= BUFFERSIZE)
	    currsize = numthisfile - globalskip;
	  buffer = malloc (currsize * sizeof (float) * 3);
	  member = malloc (currsize * sizeof (int));
	  if (ThisTask == root)
	    {
	      hsize = 4 + 256 + 4;
	      loopsize =
		(typeskip + globalskip) * sizeof (float) * 3 + sizeof (int);
	      if (globalskip == 0)
		{
		  fseek (filein, hsize, SEEK_SET);
		  SKIP;
		  /*printf ("blk count=%d bytes_seek=%d\n",
		     blk / sizeof (float) / 3, hsize + loopsize); */
		}
	      fseek (filein, hsize + loopsize, SEEK_SET);
	      fread (buffer, sizeof (float), currsize * 3, filein);
	    }
	  count = currsize * 3;
	  MPI_Bcast (buffer, count, MPI_FLOAT, root, MPI_COMM_WORLD);
	  dx = lbox / spx;
	  dy = lbox / spy;
	  dz = lbox / spz;
	  for (i = 0; i < currsize; i++)
	    {
	      if (DOMAINTYPE == 0)
		{
		  ii = (int) floor (buffer[i * 3] / dx);
		  jj = (int) floor (buffer[i * 3 + 1] / dy);
		  kk = (int) floor (buffer[i * 3 + 2] / dz);
		  if (ii == spx)
		    ii = spx - 1;
		  if (jj == spy)
		    jj = spy - 1;
		  if (kk == spz)
		    kk = spz - 1;
		  member[i] = ii * spy * spz + jj * spz + kk;	/*grid assign */
		  if (member[i] == ThisTask)
		    localindex++;
		}
	      if (DOMAINTYPE == 1)
		{
		  ii = (int) floor (buffer[i * 3] / dx);
		  if (ii == spx)
		    ii = spx - 1;
		  member[i] = ii;
		  if (member[i] == ThisTask)
		    localindex++;
		}
	    }
	  if (localindex > 0)
	    {
	      if (locbuffer == NULL)
		{
		  locbuffer = malloc (localindex * sizeof (float) * 3);
		}
	      else
		{
		  locbuffer =
		    realloc (locbuffer,
			     (globalindex + localindex) * sizeof (float) * 3);
		}
	      if (locind == NULL)
		{
		  locind = malloc (localindex * sizeof (int));
		}
	      else
		{
		  locind =
		    realloc (locind,
			     (globalindex + localindex) * sizeof (int));
		}
	      count = 0;
	      for (i = 0; i < currsize; i++)
		{
		  if (member[i] == ThisTask)
		    {
		      locind[globalindex + count] = globalcount + i;
		      locbuffer[(globalindex + count) * 3] = buffer[i * 3];
		      locbuffer[(globalindex + count) * 3 + 1] =
			buffer[i * 3 + 1];
		      locbuffer[(globalindex + count) * 3 + 2] =
			buffer[i * 3 + 2];
		      count++;
		    }
		}
	    }
	  free (buffer);
	  free (member);
	  globalindex = globalindex + localindex;
	  globalcount = globalcount + currsize;
	  globalskip = globalskip + currsize;
	  NumThis = globalindex;
	  loopcount++;
	}			/*end loop */
      MPI_Reduce (&NumThis, &i, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
#ifdef VERBOSE
      printf ("Task=%3d Npthistot=%9d\n", ThisTask, NumThis);
#endif
      MPI_Barrier (MPI_COMM_WORLD);	/*ensure all particles are allocated */
      if (ThisTask == root)
	{
	  printf ("Read Complete. Nread=%d in %d blocks of %d bytes.\n", i,
		  loopcount, BUFFERSIZE);
	  fclose (filein);
	}

    }				/*end for j */

}


void
savedata (int dtype)
{
  int i, j, k, istart, iend, ii, irem, num, iloop, isize, nthis, nro, nro2;
  float *sdata, *srec, *stemp;
  int *stemp2, *stemp3, *stemp4, *srec2, *srec3, *srec4;
  int *sidata, *sirec, *sitemp, *sdata2, *sdata3, *sdata4;
  int in1, in2, tag, nrosum, id, indt, nrosum2, ind;
  int nroot[NTask], nroot2[NTask];
  MPI_Request request, request2, request3, request4;
  MPI_Status status;
  int nn, indtotal;
  int lastrec;
  int nneigtotal;

  MPI_Reduce (&nneig, &nneigtotal, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

  indtotal = 0;			/*index neigbors */
  lastrec = 0;
  istart = 0;
  iend = BUFFERSIZE;
  iloop = 0;			/*curr loop */
  ii = 0;			/*current index */
  if (NumPart < BUFFERSIZE)
    iend = NumPart;
  isize = iend - istart;
  while (ii < NumPart)		/*main loop */
    {
      nn = 0;			/*num neighbors this */
      nthis = 0;		/*number of particles to send this task in this loop */
      ind = 0;
      for (i = 0l; i < NumThis; i++)
	{
	  if ((locind[i] >= istart) && (locind[i] < iend))
	    nthis++;
	}
      if (nthis > 0)
	{
	  if (dtype == 0)
	    sdata = malloc (nthis * sizeof (float));
	  if (dtype == 1)
	    sdata = malloc (nthis * sizeof (float));
	  if (dtype == 2)
	    sdata2 = malloc (nthis * sizeof (int));
	  sidata = malloc (nthis * sizeof (int));
	  in1 = 0;
	  for (i = 0; i < NumThis; i++)
	    {
	      if ((locind[i] >= istart) && (locind[i] < iend))
		{
		  if (dtype == 0)
		    sdata[in1] = rho[i];
		  if (dtype == 1)
		    sdata[in1] = vol[i];
		  if (dtype == 2)
		    sdata2[in1] = neignum[i];
		  if (dtype == 2)
		    nn += neignum[i];
		  sidata[in1] = locind[i];
		  in1++;
		}
	    }

	  if (dtype == 2)
	    {
	      in1 = 0;
	      sdata3 = malloc (nn * sizeof (int));	/*neigbors */
	      sdata4 = malloc (nthis * sizeof (int));	/*indices */
	      for (i = 0; i < NumThis; i++)
		{
		  if ((locind[i] >= istart) && (locind[i] < iend))
		    {
		      id = neigind[i];
		      sdata4[in1] = ind;
		      for (j = 0; j < neignum[i]; j++)
			{
			  sdata3[ind] = neigbors[id + j];
			  ind++;
			}
		      in1++;
		    }
		}		/*for i */

	    }			/*endif dtype 2 */

	}			/*nthis>0 */
      /*send number particles thistask to root */
      tag = ThisTask;

      MPI_Isend (&nthis, 1, MPI_INT, root, tag, MPI_COMM_WORLD, &request);
      if (dtype == 2)
	MPI_Isend (&nn, 1, MPI_INT, root, tag, MPI_COMM_WORLD, &request2);

      if (ThisTask == root)
	{
	  nrosum = 0;
	  nrosum2 = 0;
	  for (i = 0; i < NTask; i++)
	    {
	      tag = i;
	      MPI_Irecv (&nro, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &request);
	      MPI_Wait (&request, &status);
	      nroot[i] = nro;
	      nrosum += nro;
	      if (dtype == 2)
		{
		  MPI_Irecv (&nro2, 1, MPI_INT, i, tag, MPI_COMM_WORLD,
			     &request2);
		  MPI_Wait (&request2, &status);
		  nroot2[i] = nro2;
		  nrosum2 += nro2;
		}

	    }
	  if (dtype == 0)
	    srec = malloc (isize * sizeof (float));
	  if (dtype == 1)
	    srec = malloc (isize * sizeof (float));
	  if (dtype == 2)
	    srec2 = malloc (isize * sizeof (int));
	  sirec = malloc (isize * sizeof (int));
	  if (dtype == 2)
	    srec3 = malloc (nrosum2 * sizeof (int));	/*neigbors */
	  if (dtype == 2)
	    srec4 = malloc (isize * sizeof (int));	/*neigind */
	  if (dtype == 2)
	    for (i = 0; i < isize; i++)
	      srec4[i] = 0;
	  if (dtype == 2)
	    for (i = 0; i < nrosum2; i++)
	      srec3[i] = -1;
	}


      MPI_Barrier (MPI_COMM_WORLD);

      /*send data */
      if (nthis > 0)
	{
	  tag = ThisTask;
	  MPI_Isend (sidata, nthis, MPI_INT, root, tag, MPI_COMM_WORLD,
		     &request);
	  if (dtype == 0)
	    MPI_Isend (sdata, nthis, MPI_FLOAT, root, tag, MPI_COMM_WORLD,
		       &request2);
	  if (dtype == 1)
	    MPI_Isend (sdata, nthis, MPI_FLOAT, root, tag, MPI_COMM_WORLD,
		       &request2);
	  if (dtype == 2)
	    MPI_Isend (sdata2, nthis, MPI_INT, root, tag, MPI_COMM_WORLD,
		       &request2);
	  if (dtype == 2)
	    MPI_Isend (sdata3, nn, MPI_INT, root, tag, MPI_COMM_WORLD,
		       &request3);
	  if (dtype == 2)
	    MPI_Isend (sdata4, nthis, MPI_INT, root, tag, MPI_COMM_WORLD,
		       &request4);
	}

      if (ThisTask == root)
	{
	  indt = 0;
	  for (i = 0; i < NTask; i++)
	    {

	      nro = nroot[i];
	      nro2 = nroot2[i];
	      if (nro > 0)
		{
		  if (dtype == 0)
		    stemp = malloc (nro * sizeof (float));
		  if (dtype == 1)
		    stemp = malloc (nro * sizeof (float));
		  if (dtype == 2)
		    stemp2 = malloc (nro * sizeof (int));
		  if (dtype == 2)
		    stemp3 = malloc (nro2 * sizeof (int));
		  if (dtype == 2)
		    stemp4 = malloc (nro * sizeof (int));
		  sitemp = malloc (nro * sizeof (int));
		  tag = i;
		  MPI_Irecv (sitemp, nro, MPI_INT, i, tag, MPI_COMM_WORLD,
			     &request);
		  if (dtype == 0)
		    MPI_Irecv (stemp, nro, MPI_FLOAT, i, tag, MPI_COMM_WORLD,
			       &request2);
		  if (dtype == 1)
		    MPI_Irecv (stemp, nro, MPI_FLOAT, i, tag, MPI_COMM_WORLD,
			       &request2);
		  if (dtype == 2)
		    MPI_Irecv (stemp2, nro, MPI_INT, i, tag, MPI_COMM_WORLD,
			       &request2);
		  if (dtype == 2)
		    MPI_Irecv (stemp3, nro2, MPI_INT, i, tag, MPI_COMM_WORLD,
			       &request3);
		  if (dtype == 2)
		    MPI_Irecv (stemp4, nro, MPI_INT, i, tag, MPI_COMM_WORLD,
			       &request4);

		  MPI_Wait (&request, &status);
		  MPI_Wait (&request2, &status);
		  if (dtype == 2)
		    {
		      MPI_Wait (&request3, &status);
		      MPI_Wait (&request4, &status);
		    }

		  /*gather stemp data into srec */

		  for (j = 0; j < nro; j++)
		    {
		      id = sitemp[j];	/*local index current part range: istart to iend */

		      sirec[id - istart] = id;	/*redundant - keep for debug */
		      if (dtype < 2)
			srec[id - istart] = stemp[j];
		      if (dtype == 2)
			{
			  srec2[id - istart] = stemp2[j];	/*nneig */
			  in1 = stemp4[j];	/*index in rec buffer(stemp3) */
			  in2 = stemp2[j];	/*number neigbors */
			  srec4[id - istart] = indtotal + indt;	/*new ind */
			  for (k = 0; k < in2; k++)
			    {
			      srec3[indt] = stemp3[in1 + k];
			      indt++;

			    }

			}	/*dtype2 */
		    }
		  if (dtype == 2)
		    {
		      free (stemp2);
		      free (stemp3);
		      free (stemp4);
		    }
		  else
		    {
		      free (stemp);
		    }
		  free (sitemp);
		}		/*nro */
	    }			/*for i */
	}			/* if root */

      MPI_Barrier (MPI_COMM_WORLD);

      /*DEBUG CODE -- REMOVE */
      if (ThisTask == root)
	{
	  for (j = 0; j < 3; j++)
	    {			/*isize */
	      if (sirec[j] != (j + istart))
		printf ("fail %d %d %d\n", j, sirec[j], j + istart);
	    }
	}
      /*END DEBUG */

      if (ThisTask == root)
	{
/*=======write srec */
	  if (dtype == 0)
	    {
#ifdef WRITEASCII
	      if (iloop == 0)
		fprintf (fileden, "%d\n", NumPart);
	      for (i = 0; i < isize; i++)
		fprintf (fileden, "%10.4g\n", srec[i]);
#else
	      if (iloop == 0)
		fwrite (&NumPart, sizeof (int), 1, fileden);
	      if (fwrite (srec, sizeof (float), isize, fileden) != isize)
		{
		  printf ("Error writing density.\n");
		}
#endif
	    }

	  if (dtype == 1)
	    {
#ifdef WRITEASCII
	      if (iloop == 0)
		fprintf (filevol, "%d\n", NumPart);
	      for (i = 0; i < isize; i++)
		fprintf (filevol, "%10.4g\n", srec[i]);
#else
	      if (iloop == 0)
		fwrite (&NumPart, sizeof (int), 1, filevol);
	      if (fwrite (srec, sizeof (float), isize, filevol) != isize)
		{
		  printf ("Error writing volume.\n");
		}
#endif
	    }

	  if (dtype == 2)
	    {
#ifdef WRITEASCII
	      if (iloop == 0)
		{
		  /*printf ("nb2 array lenght=%d\n", nneigtotal); */
		  fprintf (fileneig, "%d\n", NumPart);
		  fprintf (fileneig2, "%d\n", nneigtotal);
		}
	      for (i = 0; i < isize; i++)
		fprintf (fileneig, "%d %d\n", srec2[i], srec4[i]);
	      for (i = 0; i < nrosum2; i++)
		fprintf (fileneig2, "%d\n", srec3[i]);
#else
	      if (iloop == 0)
		{
		  /*printf ("nb2 array lenght=%d\n", nneigtotal); */
		  fwrite (&NumPart, sizeof (int), 1, fileneig);
		  fwrite (&nneigtotal, sizeof (int), 1, fileneig2);
		}
	      if (fwrite (srec2, sizeof (int), isize, fileneig) != isize)
		printf ("Error writing nb.\n");
	      if (fwrite (srec4, sizeof (int), isize, fileneig) != isize)
		printf ("Error writing nb.\n");
	      if (fwrite (srec3, sizeof (int), nrosum2, fileneig2) != nrosum2)
		printf ("Error writing nb2.\n");
#endif
	    }

/*===========end write section*/

	  if (dtype == 2)
	    {
	      free (srec2);
	      free (srec3);
	      free (srec4);
	    }
	  else
	    {
	      free (srec);
	    }

	  free (sirec);
	}			/*root */

      if (nthis > 0)
	{
	  if (dtype == 2)
	    {
	      free (sdata2);
	      free (sdata3);
	      free (sdata4);
	    }
	  else
	    {
	      free (sdata);
	    }
	  free (sidata);
	}

      ii = ii + isize;
      istart = iend;
      irem = NumPart - ii;	/*remaining particles */
      isize = irem;
      if (irem > BUFFERSIZE)
	isize = BUFFERSIZE;
      iend = istart + isize;
      iloop++;
      if (ThisTask == root)
	indtotal += indt;	/*continue neighbor indices */
    }				/*end main loop */
}

/*fileden filevol*/
void
writedensopen (void)
{
  char filename[250];
  char wmode[4];

#ifdef WRITEASCII
  sprintf (wmode, "w");
#else
  sprintf (wmode, "wb");
#endif

  if (ThisTask == root)
    {
      sprintf (filename, "%s.den", InputFile);
      if ((fileden = fopen (filename, wmode)) == NULL)
	{
	  printf ("Fail open density file...\n");
	  stopcode ();
	}
#ifdef VERBOSE
      printf ("Density write in file:%s\n", filename);
#endif
      sprintf (filename, "%s.vol", InputFile);
      if ((filevol = fopen (filename, wmode)) == NULL)
	{
	  printf ("Fail open volume file...\n");
	  stopcode ();
	}
#ifdef VERBOSE
      printf ("Volume write in file:%s\n", filename);
#endif
    }
  check_stop ();

}

void
writedensclose (void)
{

  if (ThisTask == root)
    fclose (fileden);
  if (ThisTask == root)
    fclose (filevol);
}

/*neigbors files*/
void
writeneigopen (void)
{
  char filename[250];
  char wmode[4];

#ifdef WRITEASCII
  sprintf (wmode, "w");
#else
  sprintf (wmode, "wb");
#endif

  if (ThisTask == root)
    {
      sprintf (filename, "%s.nb", InputFile);
      if ((fileneig = fopen (filename, wmode)) == NULL)
	{
	  printf ("Fail open nb file...\n");
	  stopcode ();
	}
#ifdef VERBOSE
      printf ("Neigbors write in file:%s\n", filename);
#endif
      sprintf (filename, "%s.nb2", InputFile);
      if ((fileneig2 = fopen (filename, wmode)) == NULL)
	{
	  printf ("Fail open nb2 file...\n");
	  stopcode ();
	}
#ifdef VERBOSE
      printf ("Neigbors indices write in file:%s\n", filename);
#endif
    }
  check_stop ();
}

void
writeneigclose (void)
{

  if (ThisTask == root)
    fclose (fileneig);
  if (ThisTask == root)
    fclose (fileneig2);
}


/*grid files*/
void
writegridopen (void)
{
  char filename[250];
  char wmode[4];
  int itemp;

  if (ThisTask == root)
    {
      sprintf (filename, "%s.grid", InputFile);
      if ((filegrid = fopen (filename, "wb")) == NULL)
	{
	  printf ("Fail open grid file...\n");
	  stopcode ();
	}
#ifdef VERBOSE
      printf ("Grid write in file:%s\n", filename);
#endif
      sprintf (filename, "%s.vtk", InputFile);
      if ((filegrid2 = fopen (filename, "w")) == NULL)
	{
	  printf ("Fail open grid vtk file...\n");
	  stopcode ();
	}
#ifdef VERBOSE
      printf ("Grid VTK write in file:%s\n", filename);
#endif
      itemp = GRIDSIZE;
      fwrite (&itemp, sizeof (int), 1, filegrid);

      fprintf (filegrid2, "# vtk DataFile Version 2.0\n");
      fprintf (filegrid2, "PARAVT density grid\n");
      fprintf (filegrid2, "ASCII\n");
      fprintf (filegrid2, "DATASET STRUCTURED_POINTS\n");
      fprintf (filegrid2, "DIMENSIONS %d %d %d\n", GRIDSIZE, GRIDSIZE,
	       GRIDSIZE);
      fprintf (filegrid2, "ORIGIN 0 0 0\n");
      fprintf (filegrid2, "SPACING 1 1 1\n");
      fprintf (filegrid2, "POINT_DATA %d\n", GRIDSIZE * GRIDSIZE * GRIDSIZE);
      fprintf (filegrid2, "SCALARS density float 1\n");
      fprintf (filegrid2, "LOOKUP_TABLE default\n");
    }
  check_stop ();
}


void
writegridclose (void)
{

  if (ThisTask == root)
    fclose (filegrid);
  if (ThisTask == root)
    fclose (filegrid2);
}
