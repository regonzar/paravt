#include "qhullsrc_r/qhull_ra.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "func.h"
#include "vars.h"
#include "config.h"
#include <mpi.h>
#include <float.h>

#define MAXGNB 300 /*used in gradient function*/


int nadj, nadjpos;		/*number of facets */
int *adj, *adji;		/*facets and facet index */
coordT *pos;
float *mass;
int nvertex, nfacets;
float *xv, *yv, *zv;		/*vertices */
int *neigbors;			/*array ids of neighbors */
int *neignum;			/*number of neighbors each particles */
int *neigind;			/*position in neighbors array where particle's neighbors begin */
float *vol;
float *rho;
float *grad;
int nbrank;
int *neigborsg;               /*array ids of neighbors for gradient computation*/


/*Uses Qhull to compute Voronoi Tessellation*/
void
runvt (void)
{

  vertexT *vertex;
  setT *vertices;
  facetT *facet;

  char flags[MAX_FILELEN];
  int exitcode;
  FILE *errfile = stdout;
  FILE *outfil;
  boolT ismalloc = False;
  int i, ntot;
  boolT islower;
  int numfacets, numsimplicial, numridges, totneighbors, numcoplanars,
    numtricoplanars, k, num, numcenters, totcount, vertex_i, vertex_n;
  qh_RIDGE innerouter;
  int curlong, totlong;
  int bufrank;
  float x1, x2, y1, y2, z1, z2;
  char fname[MAX_FILELEN];
  printvridgeT printvridge = qh_printvridge;
  qhT qh_qh;			/*qhull data structure */
  qhT *qh = &qh_qh;

  QHULL_LIB_CHECK qh_zero (qh, errfile);

#ifdef MULTIMASS
  bufrank = 4;
#else
  bufrank = 3;
#endif

  ntot = NumThis + NumThisb;
  pos = malloc (3 * ntot * sizeof (coordT));

#ifdef MULTIMASS
  mass = malloc (ntot * sizeof (float));
#endif

  x1 = FLT_MAX;
  y1 = FLT_MAX;
  z1 = FLT_MAX;
  x2 = FLT_MIN;
  y2 = FLT_MIN;
  z2 = FLT_MIN;
  for (i = 0; i < ntot; i++)
    {
      pos[i * 3] = locbuffer[i * bufrank];
      pos[i * 3 + 1] = locbuffer[i * bufrank + 1];
      pos[i * 3 + 2] = locbuffer[i * bufrank + 2];
#ifdef MULTIMASS
      mass[i] = locbuffer[i * bufrank + 3];
#endif

/*Debug line to check limits for each task*/
/*      if (pos[i * 3] < x1)
	x1 = pos[i * 3];
      if (pos[i * 3 + 1] < y1)
	y1 = pos[i * 3 + 1];
      if (pos[i * 3 + 2] < z1)
	z1 = pos[i * 3 + 2];
      if (pos[i * 3] > x2)
	x2 = pos[i * 3];
      if (pos[i * 3 + 1] > y2)
	y2 = pos[i * 3 + 1];
      if (pos[i * 3 + 2] > z2)
	z2 = pos[i * 3 + 2];*/
    }
/*printf("vt Task=%d x1=%f x2=%f y1=%f y2=%f z1=%f z2=%f\n",ThisTask,x1,x2,y1,y2,z1,z2);*/
  free (locbuffer);

  sprintf (flags, "qhull v %s p Fv", QHULLOPTIONS);
  if (ThisTask == 0)
    printf ("Running VT with options: %s\n", flags);

#ifdef WRITEVTSTRUCT
  sprintf (fname, "%s.VT.%d", InputFile, ThisTask);
  if ((outfil = fopen (fname, "w")) == NULL)
    {
      printf ("Fail open file: %s\n", fname);
      stopcode ();
    }
  /*We write local indices in VTSTRUCT file */
  fprintf (outfil, "%d %d\n", NumThis, NumThisb);
  for (i = 0; i < ntot; i++)
    fprintf (outfil, "%d\n", locind[i]);

  exitcode =
    qh_new_qhull (qh, DIM, ntot, pos, ismalloc, flags, outfil, errfile);
  printf ("Written VTSTRUCTURE: %s\n", fname);
  fflush (outfil);
  fflush (errfile);
  fclose (outfil);
#else
  exitcode =
    qh_new_qhull (qh, DIM, ntot, pos, ismalloc, flags, NULL, errfile);
#endif

  if (exitcode)
    {
      printf ("ERROR RUNNING QHULL!!!!!\n\n");
      free (locind);
      free (pos);
#ifdef MULTIMASS
      free (mass);
#endif
      MPI_Finalize ();
      exit (0);
    }
  if (!exitcode)
    {
      /*extract num of vertices */
      qh_setvoronoi_all (qh);
      qh_countfacets (qh, qh->facet_list, NULL, !qh_ALL, &numfacets,
		      &numsimplicial, &totneighbors, &numridges,
		      &numcoplanars, &numtricoplanars);
      nvertex = numfacets;

      /*Allocate mem for vertices */
      if ((xv = malloc (nvertex * sizeof (float))) == NULL)
	printf ("xv not allocated..\n");
      if ((yv = malloc (nvertex * sizeof (float))) == NULL)
	printf ("yv not allocated..\n");
      if ((zv = malloc (nvertex * sizeof (float))) == NULL)
	printf ("zv not allocated..\n");

      /*Loop over voronoi vertices */
      k = 0;
      i = 0;
      FORALLfacet_ (qh->facet_list)
      {
	num = qh->hull_dim - 1;	/*dims, not used yet, only 3D data */
	if (!facet->normal || !facet->upperdelaunay || !qh->ATinfinity)
	  {
	    //if (!facet->center)
	    //  {
	    //  facet->center = qh_facetcenter (qh, facet->vertices);
	    //  i++;
	    //  }
	    if (facet->center)
	      {
		xv[k] = facet->center[0];
		yv[k] = facet->center[1];
		zv[k] = facet->center[2];
		k++;
	      }
	  }
	else
	  {
	    printf ("qhINFINITE CASE \n");
	  }
      }
      if (k > nvertex)
	printf
	  ("Task=%d WARNING! check bounds vertices memory. k=%d nvertex=%d\n",
	   ThisTask, k, nvertex);


      innerouter = qh_RIDGEall;
      vertices = qh_markvoronoi (qh, qh->facet_list, NULL, 0, &islower, &numcenters);	/*arg 3,0->null */
      totcount =
	qh_printvdiagram2 (qh, NULL, NULL, vertices, innerouter, False);
      nfacets = totcount;


      /*here we use NFAC which is a guess of the mean number of vertices in a facet */
      if ((adj = malloc (nfacets * NFAC * sizeof (int))) == NULL)
	{			/*space for facets */
	  printf ("Not enought memory for facets (%d)\n",
		  nfacets * NFAC * sizeof (int));
	  finalize ();
	}
      if ((adji = malloc (nfacets * sizeof (int))) == NULL)
	printf ("adji not allocated..\n");
      for (i = 0; i < nfacets * NFAC; i++)
	adj[i] = 0;
      for (i = 0; i < nfacets; i++)
	adji[i] = 0;

      nadj = 0;
      nadjpos = 0;

      /*Make a loop over all facets & look for id1, id2 and vertices ids */
      totcount = 0;
      vertex_i = 0;
      vertex_n = 0;
      FORALLvertices vertex->seen = False;

      FOREACHvertex_i_ (qh, vertices)
      {
	if (vertex)
	  {
	    if (qh->GOODvertex > 0
		&& qh_pointid (qh, vertex->point) + 1 != qh->GOODvertex)
	      continue;
	    totcount +=
	      qh_eachvoronoi_mod (qh, errfile, printvridge, vertex, !qh_ALL,
				  innerouter, True);
	  }
      }


      if (nadjpos > nfacets*NFAC) printf("WARNING, you must increase NFAC value\n");

      /*printf ("Task=%d size adj=%d efective size adj=%d\n", ThisTask, nfacets*NFAC, nadjpos);*/

      qh_settempfree (qh, &vertices);

      qh->NOerrexit = True;
      qh_freeqhull (qh, !qh_ALL);
      qh_memfreeshort (qh, &curlong, &totlong);
      if (curlong || totlong)
	fprintf (errfile,
		 "warning: did not free %d bytes of long memory (%d pieces)\n",
		 totlong, curlong);
    }
}

/*Visit all facets and vertices & fills adj & adji arrays*/
int
qh_eachvoronoi_mod (qhT * qh, FILE * fp, printvridgeT printvridge,
		    vertexT * atvertex, boolT visitall, qh_RIDGE innerouter,
		    boolT inorder)
{
  boolT unbounded;
  int count, ind, k;
  facetT *neighbor, **neighborp, *neighborA, **neighborAp, *facet, **facetp;
  setT *centers;
  setT *tricenters = qh_settemp (qh, qh->TEMPsize);

  vertexT *vertex, **vertexp;
  boolT firstinf;
  unsigned int numfacets = (unsigned int) qh->num_facets;
  int totridges = 0;

  qh->vertex_visit++;
  atvertex->seen = True;
  if (visitall)
    {
      FORALLvertices vertex->seen = False;
    }
  FOREACHneighbor_ (atvertex)
  {
    if (neighbor->visitid < numfacets)
      neighbor->seen = True;
  }
  FOREACHneighbor_ (atvertex)
  {
    if (neighbor->seen)
      {
	FOREACHvertex_ (neighbor->vertices)
	{
	  if (vertex->visitid != qh->vertex_visit && !vertex->seen)
	    {
	      vertex->visitid = qh->vertex_visit;
	      count = 0;
	      firstinf = True;
	      qh_settruncate (qh, tricenters, 0);
	      FOREACHneighborA_ (vertex)
	      {
		if (neighborA->seen)
		  {
		    if (neighborA->visitid)
		      {
			if (!neighborA->tricoplanar
			    || qh_setunique (qh, &tricenters,
					     neighborA->center))
			  count++;
		      }
		    else if (firstinf)
		      {
			count++;
			firstinf = False;
		      }
		  }
	      }
	      if (count >= qh->hull_dim - 1)
		{		/* e.g., 3 for 3-d Voronoi */
		  if (firstinf)
		    {
		      if (innerouter == qh_RIDGEouter)
			continue;
		      unbounded = False;
		    }
		  else
		    {
		      if (innerouter == qh_RIDGEinner)
			continue;
		      unbounded = True;
		    }
		  totridges++;

		  trace4 ((qh, qh->ferr, 4017,
			   "qh_eachvoronoi: Voronoi ridge of %d vertices between sites %d and %d\n",
			   count, qh_pointid (qh, atvertex->point),
			   qh_pointid (qh, vertex->point)));

		  if (printvridge)
		    {
		      if (inorder && qh->hull_dim == 3 + 1)	/* 3-d Voronoi diagram */
			centers = qh_detvridge3 (qh, atvertex, vertex);
		      else
			centers = qh_detvridge (qh, vertex);
		      /*Loads Adj array */
		      adj[nadjpos] = qh_setsize (qh, centers) + 2;
		      adj[nadjpos + 1] = qh_pointid (qh, atvertex->point);
		      adj[nadjpos + 2] = qh_pointid (qh, vertex->point);
		      adji[nadj] = nadjpos;	/*starting index of current facet */
		      k = 0;
		      FOREACHfacet_ (centers)
		      {
			adj[nadjpos + 3 + k] = facet->visitid;
			k++;
		      }
		      nadj++;	/*next facet */
		      nadjpos += 3 + k;
		      qh_settempfree (qh, &centers);
		    }
		}
	    }
	}
      }
  }
  FOREACHneighbor_ (atvertex) neighbor->seen = False;
  qh_settempfree (qh, &tricenters);
  return totridges;
}				/* eachvoronoi end */

void
computeneig (void)
{
  int i, num, nmax, id1, id2, nsum, in;

  num = NumThis + NumThisb;
  neignum = malloc (num * sizeof (int));
  neigind = malloc (num * sizeof (int));
  for (i = 0; i < num; i++)
    neignum[i] = 0;
  for (i = 0; i < nfacets; i++)
    {
      id1 = adj[adji[i] + 1];
      id2 = adj[adji[i] + 2];
      /*if ((pos[3*id1]-pos[3*id2])>50) printf("VT I: %d %d %f %f %f %f %f %f\n",ThisTask,i,pos[3*id1],pos[3*id2],pos[3*id1+1],pos[3*id2+1],pos[3*id1+2],pos[3*id2+2]); *//*esta bien, solo hay particulas asi en bordes */
      if ((id1 >= 0) && (id2 >= 0))
	{
	  neignum[id1] = neignum[id1] + 1;
	  neignum[id2] = neignum[id2] + 1;
	}
    }
  nmax = 0;
  nsum = 0;
  neigind[0] = 0;
  for (i = 0; i < num; i++)
    {
      if (neignum[i] > nmax)
	nmax = neignum[i];
      nsum = nsum + neignum[i];
      if (i > 0)
	neigind[i] = neigind[i - 1] + neignum[i - 1];
    }
  /*printf ("Task=%d N.facets=%d N.vertices=%d N.pointfacets=%d\n", ThisTask, nfacets, nvertex, nadjpos); */
  /*printf ("Task=%d Max_neigbors=%d neigbors_lenght=%d lastne=%d lastnn=%d\n",
     ThisTask, nmax, nsum, neigind[num - 1], neignum[num - 1]); */

  printf
    ("Task=%d N.facets=%d N.vertices=%d N.pointfacets=%d Max_neighbors=%d\n",
     ThisTask, nfacets, nvertex, nadjpos, nmax);

  neigbors = malloc (nsum * sizeof (int));

#ifdef COMPUTEGRAD
  neigborsg = malloc (nsum * sizeof (int));
#endif
  for (i = 0; i < num; i++)
    neignum[i] = 0;
  for (i = 0; i < nfacets; i++)
    {
      id1 = adj[adji[i] + 1];
      id2 = adj[adji[i] + 2];
      if ((id1 >= 0) && (id2 >= 0))
	{
	  in = neigind[id1];
	  neigbors[in + neignum[id1]] = id2;
	  in = neigind[id2];
	  neigbors[in + neignum[id2]] = id1;
	  neignum[id1] = neignum[id1] + 1;
	  neignum[id2] = neignum[id2] + 1;
	}
    }

  for (i = 0; i < nsum; i++)	/*reordenamos indices vecinos a valores absolutos */
    {
#ifdef COMPUTEGRAD
      neigborsg[i]= neigbors[i];
#endif    
      in = neigbors[i];
      neigbors[i] = locind[in];
    }

  /*calculate nneig each task */
  nneig = 0;
  for (i = 0; i < NumThis; i++)
    {
      nneig += neignum[i];
    }
}

void
computedens (void)
{
  int i, j, num, nv, id1, id2, i1, i2, i3, vocount;
  float h, area, d1, d2, d3, s, vo, minvol, maxvol, maxrho, minmass,
    maxmass, vomax;
  num = NumThis + NumThisb;
  vol = malloc (num * sizeof (float));
  rho = malloc (num * sizeof (float));
  for (i = 0; i < num; i++)
    {
      vol[i] = 0.0;
      rho[i] = 0.0;
    }
  vocount = 0;
  vomax = bsize * bsize * bsize * 8.0;	/*Some vertices are located at larger distances in boundaries */
  /*then we truncate the maximum volume to a cube with side of interparticle distance */

  /*Calculate volume of voronoi cells
   * -loop over all facet
   * -calculate area of each triangle formed by 3 cotiguous vertices
   *  in a facet then computes the volume of each simplex formed by 
   *  this triangle and the particle position  
   * -sum all simplex in common to a particle to obtain the cell volume
   *                         */
  for (i = 0; i < nfacets; i++)
    {
      nv = adj[adji[i]] - 2;	/*number of vertices in facet */
      id1 = adj[adji[i] + 1];	/*IDs of the 2 particles in common with the facet */
      id2 = adj[adji[i] + 2];	/*(+1 in fortran) */
      if ((id1 > num - 1) || (id2 > num - 1))
	printf ("error id1 id2");
      if ((id1 > 0) && (id2 > 0))
	{			/*if bounded then continue */
	  h =
	    0.5 * sqrtf (sqr (pos[3 * id1] - pos[3 * id2]) +
			 sqr (pos[3 * id1 + 1] - pos[3 * id2 + 1]) +
			 sqr (pos[3 * id1 + 2] - pos[3 * id2 + 2]));
	  area = 0.0;
	  if (nv > 2)
	    {			/*loop over all vertices, must have more than 2 points to form a surface */
	      for (j = 1; j < nv - 1; j++)
		{
		  i1 = adj[adji[i] + 3] - 1;	/*1st vertex of facet, fixed common to all simplex */
		  i2 = adj[adji[i] + 3 + j] - 1;	/*other 2 contoguous vertices */
		  i3 = adj[adji[i] + 4 + j] - 1;
		  /*distances between vertices */
		  d1 =
		    sqrtf (sqr (xv[i1] - xv[i2]) + sqr (yv[i1] - yv[i2]) +
			   sqr (zv[i1] - zv[i2]));
		  d2 =
		    sqrtf (sqr (xv[i2] - xv[i3]) + sqr (yv[i2] - yv[i3]) +
			   sqr (zv[i2] - zv[i3]));
		  d3 =
		    sqrtf (sqr (xv[i1] - xv[i3]) + sqr (yv[i1] - yv[i3]) +
			   sqr (zv[i1] - zv[i3]));
		  s = 0.5 * (d1 + d2 + d3);	/*semiperimeter */
		  area += sqrtf (fabsf (s * (s - d1) * (s - d2) * (s - d3)));	/*heron's formula */
		}		/*endfor */
	    }
	  else
	    printf ("FAIL nv=2\n");
	  vo = 0.3333333333 * area * h;

	  if (vo > vomax)
	    {
	      vo = vomax;	/*very unfrequent(5 in 3million for only 
				   boundary cells where some coplanar configuration returns infinite volume) */
	      vocount++;
	    }
	  vol[id1] = vol[id1] + vo;	/*add polygon volume to cell 1 */
	  vol[id2] = vol[id2] + vo;	/*add polygon volume to cell 2, it is symetric w/r facet */
	}
    }				/*endfor nfacets */
  minvol = FLT_MAX;
  maxvol = FLT_MIN;
  for (i = 0; i < num; i++)
    {
      if ((vol[i] > 1.0e-20) && (vol[i] < minvol))
	minvol = vol[i];
      if ((vol[i] > maxvol))
	maxvol = vol[i];
    }
  for (i = 0; i < num; i++)
    {
      if ((vol[i] < minvol))
	vol[i] = minvol;
    }


  minrho = FLT_MAX;
  maxrho = FLT_MIN;
  minmass = FLT_MAX;
  maxmass = FLT_MIN;
  for (i = 0; i < num; i++)
    {
#ifdef MULTIMASS
      rho[i] = mass[i] / vol[i];
      if (mass[i] > maxmass)
	maxmass = mass[i];
      if (mass[i] < minmass)
	minmass = mass[i];
#else
      rho[i] = 1.0 / vol[i];
#endif
      if (rho[i] > maxrho)
	maxrho = rho[i];
      if (rho[i] < minrho)
	minrho = rho[i];

    }
#ifdef VERBOSE
  /*vocount is the number of incomplete particles with large volumes located in boundaries */
  printf
    ("Dens Task=%d mass[%8.2g %8.2g] vol=[%8.2g %8.2g] den=[%8.2g %8.2g]\n",
     ThisTask, minmass, maxmass, minvol, maxvol, minrho, maxrho);
#endif
  free (xv);			/* vertices no longer used */
  free (yv);
  free (zv);

}


void 
computegrad (void)        /* Gradient computation*/
{
int i,j,num,nne,ik,ij,dimin,dimax;
float xk[MAXGNB],yk[MAXGNB],zk[MAXGNB],denv[MAXGNB];
float divmin,divmax;
double mx,my,mz,b,r;

num = NumThis + NumThisb;
grad = malloc (3 * NumThis * sizeof(float));
//nbrank = malloc (num * sizeof(int));

xk[0]=0.;
yk[0]=0.;
zk[0]=0.;
for (i=0;i< NumThis;i++)
  {
    nne = neignum[i];
    if (nne < 1) printf("Particle with no neighbors!!!\n",i);
    denv[0]=rho[i];
    dimin=0; 
    dimax=0;
    divmin=rho[i];
    divmax=rho[i];
    for (j=0;j< nne;j++)
    {
     ik = neigind[i]+j;
     ij = neigborsg[ik]; 
     xk[j+1]=pos[3*ij]-pos[3*i];
     yk[j+1]=pos[3*ij+1]-pos[3*i +1];
     zk[j+1]=pos[3*ij+2]-pos[3*i +2];
     denv[j+1]=rho[ij];
     if (rho[ij] > divmax) {
        dimax=j+1;
        divmax=rho[ij];
        }
     if (rho[ij] < divmin) {
        dimin=j+1;
        divmin=rho[ij];
        }
    /*printf("-- %d %d %f %f %f %g\n",j,nne,xk[j+1],yk[j+1],zk[j+1],denv[j+1]);*/
    }

  //nbrank[i]=dimax+1; /*largest density neighbor index nbrank[i]-2 */
  //if (dimax == 0) nbrank[i]=1; /* maxima */
  //if (dimin == 0) nbrank[i]=0; /* minimum */
  /*tenemos pos y dens de vecinos, ahora fit lineal en cada direccion */

  if (linreg(nne+1,xk,denv,&mx,&b,&r)==0) mx=0;
  grad[3*i]=(float) mx;
  if (linreg(nne+1,yk,denv,&my,&b,&r)==0) my=0;
  grad[3*i+1]=(float) my;
  if (linreg(nne+1,zk,denv,&mz,&b,&r)==0) mz=0;
  grad[3*i+2]=(float) mz;
  /*printf(">> %d %d %g %g %g %g\n",i,NumThis,mx,my,mz,b);*/

  }/*i*/

}










