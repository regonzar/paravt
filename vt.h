#include "qhullsrc_r/qhull_ra.h"

extern int nadj, nadjpos;	/*number of facets */
extern int *adj, *adji;		/*facets and facet index */
extern coordT *pos;
extern float *mass;
extern int nvertex, nfacets;
extern float *xv, *yv, *zv;	/*vertices */
extern vertexT *vertex;
extern setT *vertices;
extern facetT *facet;		/* set by FORALLfacets */
extern int *neigbors;		/*array ids of neighbors */
extern int *neignum;		/*number of neighbors each particles */
extern int *neigind;		/*position in neighbors array where particle's neighbors begin */
extern float *vol;
extern float *rho;
extern float *grad;
extern int nbrank;

extern int *neigborsg;           /*array ids of neighbors for gradient computation only */

void runvt (void);

int qh_eachvoronoi_mod (qhT * qh, FILE * fp, printvridgeT printvridge,
			vertexT * atvertex, boolT visitall,
			qh_RIDGE innerouter, boolT inorder);

void computeneig (void);

void computedens (void);

void computegrad (void);
