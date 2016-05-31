/*Configuration file*/
/*Check userguide.pdf for additional information*/

#define DOMAINTYPE 0
#define BOX
#define PERIODIC
#define MULTIMASS
#define BUFFERSIZE 100000
#define BORDERFACTOR 1.5
#define DIM 3
#define NFAC 15
#define WRITEASCII
#define WRITEDENSITY
#define WRITENEIGHBORS
/*#define WRITEVTSTRUCT*/
/*#define COMPUTEGRID*/
#define GRIDSIZE 64
#define GADGET_PTYPE 1
#define QHULLOPTIONS "QJ"
/*#define LOWMEMORY*/
#define VERBOSE
/*#define COMPUTEGRAD*/
/*#define AUTOBORDER*/


/* DOMAINTYPE: Domain decomposition strategy
 * For box volumes and periodic 0 is recommended,
 *   0 == rectangular blocks where Ntask must be a power of 2.
 * For irregular volumes 1 is recommended,
 *   1 == split volume along the larger axis in Ntask blocks, 
 *   Ntask can be any number. This disable PERIODIC.
 *
 * BOX: volume is a box
 * Set if volume is a box, it is automatically set for some
 * formats (gadget, rockstar). 
 * For boxes, DOMAINTYPE 0 is recommended.
 *
 * PERIODIC: Periodic conditions
 * Only used if BOX volume is defined
 *
 * MULTIMASS: if set, read an additional column of data with 
 * particle masses. (for GADGET it is automatically set, and
 * for Rockstar it is enabled)
 *
 * BUFFERSIZE: number of particles to read per slab
 * Notices a larger value requires more temporary memory to 
 * preallocate particles. Default(10000)
 *
 * BORDERFACTOR: Thickness boundary particles around each cell
 * in terms of the inter-particle distance. Default (3.0)
 *
 * DIM: NUmber of dimensions (only works for 3 dims for the moment)
 *
 * NFAC: Guess of the mean number of vertices of a facet (default 20)
 *
 * WRITE...: If defined write additional output
 *
 * COMPUTEGRID: Compute density in a grid using two methods, count
 * particles in each cell, and using average voronoi density. Grid
 * size is defined by GRIDSIZE parameter and must be divisible by 
 * domain decomposition.
 *  
 * GADGET_PTYPE: If read gadget file, it set the particle type to use. 
 *
 * QHULLOPTIONS: Qhull additional parameters, default "QJ", the option
 * "Qbb" is more precise but fail if repeated coordinates found.
 *
 * LOWMEMORY: Reduces total memory consumption to ~ 2/3 factor but 
 * doubles execution time.
 *
 * VERBOSE: write additional info on screen
 * 
 * COMPUTEGRAD: Compute density gradient vectors on VT, write 
 * if WRITEDENSITY defined.
 *
 * AUTOBORDER: Automatically increase boundary particles thickness (given 
 * initially by BORDERFACTOR) between two tasks, in the case there are too 
 * few boundary particles for correct VT computation. This paremeter may 
 * be helpful for very anisotropic particle distribution, but it takes more
 * computation time. If your data is well behaved with a given BORDERFACTOR, 
 * set this parameter off. It requires BORDERFACTOR > 1.
 *  
 * */
