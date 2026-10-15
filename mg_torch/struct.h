#define MAXPE 256
#define MAXLEV 10
#define MAXQ 100
#define MAXPROP 10
#define MAXD 3
#define IMAX 211
#define JMAX 71
#define MAXSPECIES 100
#define MAXREACT 1000
#ifdef DOUBLE
#define float double
#endif
struct join {
  int ijoin, bound, done;
  bool empty;
  float *pa0, *dpdt, *dpdx, *cor, xj[MAXD], area;
  float dudx[MAXD][MAXD], divv, fvisc;	// Added by BM for ADI method.
  struct cell *left, *right;
  struct join  *prev, *nxt, *par;
};
struct cell {
  int icell, ichild, mat, ilev, nsub, solid;
  bool refined, empty, halo, gotdiff;
  float qsol, dqsoldt, qa[MAXQ], pa[MAXQ], dqdt[MAXQ], diff[MAXQ],
    xc[MAXD], vol, dudx[MAXD][MAXD], dtrat;
  float *x, *b, *Ac, *Al, *Ar, *dpdt, *tmp;			// Added by BM for ADI method. 
  struct cell *child[1<<MAXD], *par, *nxt, *prev;
  struct join *ljoin[MAXD], *rjoin[MAXD];
  struct mcell *mcptr;
  struct part *pptr;
};
struct part {
  int ilev;
  float t, mass, en, xp[MAXD], vp[MAXD], props[MAXPROP]; 
  struct part *nxt, *prev;
  struct cell *cptr;
};
struct level {
  int nc[MAXD];
  long long int nt, ndt;
  float dxbar[MAXD], xl[MAXD], xr[MAXD];
  struct cell *fcell, *lcell;
  struct join **fjoin, **ljoin;
};
struct cell_list {
  struct cell *cptr;
  struct cell_list *prev, *nxt;
};
struct join_list {
  struct join *jptr;
  struct join_list *nxt;
};
struct part_list {
  int pproc;
  struct part *pptr;
  struct part_list *nxt;
};
struct mcell {
  bool gotdiff, doit, ghost;
  float qa[MAXQ], qt[MAXQ], mso[MAXQ], pso[MAXQ], diff[MAXQ],
    mcofp[MAXD][MAXQ], mcofm[MAXD][MAXQ], vol, xc[MAXD];
  struct mcell *child[1<<MAXD], *lcell[MAXD], *rcell[MAXD], *par, *nxt;
  struct cell *cptr;
};

