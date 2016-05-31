/* Global variables*/

#include "vars.h"
#include <stdio.h>

int ThisTask;			/*local task */
int NTask;			/*total number of tasks */
int NumPart;			/*number of particles */
int NumThis;			/*number particles this */
int NumThisb;			/*number boundary particles this */
int fileformat;			/*input file format */

char InputFile[MAX_FILELEN];

float lbox;			/*box size */
float lboxx, lboxy, lboxz;	/*for rectangular shapes */
float dx, dy, dz;		/*cell spacing */
int ithis, jthis, kthis;	/*cell postion this task */
float bsize;			/*border size */

int localindex;			/*part in this task in this buffer loop */
int globalindex;		/*total particles allocated in this task */
int globalcount;		/*read particles all tasks */
int currsize;			/*read particles in current loop */
int loopcount;			/*current loop number */
int root;			/* root process (default 0) */
int spx, spy, spz;		/*domain divisions */
int nneig;			/*number neigbors of non boundary particles this task */

FILE *filein;			/*Input File */
FILE *fileden;			/*output density File */
FILE *filevol;			/*output volume File */
FILE *fileneig;			/*output neigbors File */
FILE *fileneig2;		/*output neigbors indices File */
FILE *filegrid;			/*output grid File */
FILE *filegrid2;		/*output grid file 2 */
FILE *filegrad;                 /*output gradient file*/

float *buffer;			/*Read buffer */
int *member;			/*save cell location for each particle */
float *locbuffer;		/*buffer allocated each task */
int *locind;			/*index part. for local buffer each task */

int stopit;

float minrho;

float pdensity;          /*particle density*/

/*! Header for the standard file format.
 *  */
struct io_header header;	/*!< holds header for snapshot files */
