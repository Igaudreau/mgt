/* vis.h:  vis module definitions (assume const.h,model.h) */
#define MAXDT MAXD+1
#define MAXF 1024

typedef struct grid_vis_type {
   bool crunched;
   char title[SZLIN];
   char geom[SZNAM];
   int nd;
   int ndc;
   char coords[MAXDT][SZNAM];
   char clabel[MAXDT][SZNAM];
   char units[MAXDT][SZNAM];
   int ndo;
   float xo[MAXDT];
   int nc[MAXDT];
   bool nodal[MAXDT];
   float *xe[MAXDT];
   float *xa[MAXDT];
   bool nogrid;
   bool noscalar;
   char qlabel[SZNAM];
   char qunits[SZNAM];
   bool novector;
   int ndv;
   char vlabels[MAXD][SZLIN];
   int vector_index[MAXD];
   char vunits[SZNAM];
   int *sa;
   float *qa;
   float *va;
   int *sb;
   float *qb;
   float *vb;
   float xmin[MAXDT];
   float xmax[MAXDT];
   bool noxr[MAXDT];
   float qmin;
   float qmax;
   bool noqr;
   float vmax[MAXD];
   bool novmax[MAXD];
#ifdef CURV
   bool curv;
   int kad[MAXDT];
#endif
} grid_vis_struct;

typedef struct scatter_vis_type {
   char title[SZLIN];
   char geom[SZNAM];
   int nd;
   int ndc;
   char coords[MAXD][SZNAM];
   char clabel[MAXD][SZNAM];
   char units[MAXD][SZNAM];
   int ndo;
   float xo[MAXD];
   bool noscalar;
   char qlabel[SZNAM];
   char qunits[SZNAM];
   bool novector;
   int ndv;
   char vlabels[MAXD][SZLIN];
   int vector_index[MAXD];
   char vunits[SZNAM];
   int npt;
   float *xpt;
   int *spt;
   float *qpt;
   float *vpt;
   float xmin[MAXD];
   float xmax[MAXD];
   bool noxr[MAXD];
   float qmin;
   float qmax;
   bool noqr;
   float vmax[MAXD];
   bool novmax[MAXD];
#ifdef CURV
   bool curv;
   int kad[MAXDT];
#endif
} scatter_vis_struct;
