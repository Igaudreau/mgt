/* model.h:  model definitions */
#include "const.h"
#define MAXD 3
#define MAXLEV 31
#ifdef OLD
#define MAXT 20000 /* was 4096 */
#else
#define MAXT 2000
#endif
#define MAXSIDE 2
#define LEFT 0
#define RIGHT 1
#define MAXQ 128
#define SZCO MAXD*SZNAM
#define SZCQ MAXQ*SZNAM
#define MAXHIST 64

typedef struct grid_frame_type {
   float t;
   int nlev;
   int topn;
   int *sa;
   int *subn;
   float *qa;
   int topj;
   int *sb;
   int *subj;
   float *qb;
} grid_frame_struct;

typedef struct grid_model_type {
   float version;
   char date[SZDATE];
   char title[SZLIN];
   int nhist;
   char hist[MAXHIST][SZLIN];
   char geom[SZNAM];
   int nd;
   char coords[MAXD][SZNAM];
   char clabel[MAXD][SZNAM];
   char units[MAXD][SZNAM];
   int ndo;
   float xo[MAXD];
   char physics[SZNAM];
   int nq;
   char vlabels[MAXQ][SZNAM];
   char vunits[MAXQ][SZNAM];
   int vmlt[MAXQ][3];
   int limlev;
   int nc[MAXLEV][MAXD];
   float *xe[MAXLEV][MAXD];
   float *xa[MAXLEV][MAXD];
   int nt;
   grid_frame_struct frame[MAXT];
#ifdef CURV
   bool curv;
   int kad[MAXD];
#endif
} grid_model_struct;

typedef struct scatter_frame_type {
   float t;
   int *spt;
   float *qpt;
} scatter_frame_struct;

typedef struct scatter_model_type {
   bool empty;
   float version;
   char date[SZDATE];
   char title[SZLIN];
   int nhist;
   char hist[MAXHIST][SZLIN];
   char geom[SZNAM];
   int nd;
   char coords[MAXD][SZNAM];
   char clabel[MAXD][SZNAM];
   char units[MAXD][SZNAM];
   int ndo;
   float xo[MAXD];
   char physics[SZNAM];
   int nq;
   char vlabels[MAXQ][SZNAM];
   char vunits[MAXQ][SZNAM];
   int vmlt[MAXQ][3];
   int npt;
   float *xpt;
   int nt;
   scatter_frame_struct frame[MAXT];
#ifdef CURV
   bool curv;
   int kad[MAXD];
#endif
} scatter_model_struct;
