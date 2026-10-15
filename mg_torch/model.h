using namespace std;
#undef MUG
#ifdef MAIN
string bcond[MAXD];
char  geom[SZNAM], title[SZNAM], coord[MAXD][SZITEM], item[SZNAM],
 kind[SZNAM], cas[SZNAM];
int nd, nq,  nal, scad, nlevs, alev, levc, levb, order, plmin, plmax, 
  nc[MAXD], bcl[MAXD], bcr[MAXD], ndiffstep, carrtop, celltop, jarrtop,
  jointop, nrcell, nrel, mcarrtop;
long long int nupdat, tnupdat;
bool admin, debug, dump, nom, diffbad, remhole, correct, defrag, 
  regrid, repart, keps, viscous, fast, grav, nomug, mugprt, boundpot,
  prof, timpart, particle;
float t, dt0, target, xl[MAXD], xr[MAXD], cn, cnf, cnv, cng, rtol, 
  refcrit,derefcrit, normrat, tcpu, cputime, mdcpu, initdcpu;
float wall_time, ostep_time;
struct level **lev;
struct cell **carr;
struct join **jarr;
struct cell_list *freecell;
struct join_list *freejoin;
struct mcell **mcarr;
#else
extern string bcond[MAXD];
extern char geom[SZNAM], title[SZNAM], coord[MAXD][SZITEM],
  item[SZNAM], kind[SZNAM], cas[SZNAM];
extern int nd, nq, nal, scad,  nlevs,alev,  levc, levb, order, plmin,
  plmax, nc[MAXD], bcl[MAXD], bcr[MAXD],  ndiffstep, carrtop, celltop, 
  jarrtop, jointop, nrcell, nrel, mcarrtop;
extern long long int nupdat, tnupdat;
extern bool admin, debug, dump, nom, diffbad, remhole, correct, 
  defrag, regrid, repart, keps, viscous, fast, grav, nomug, mugprt,
  boundpot, prof, timpart, particle;
extern float t, dt0, target, xl[MAXD], xr[MAXD], dxbar[MAXD], cn, cnf, 
  cnv, cng, rtol, refcrit, derefcrit, normrat, tcpu, cputime, mdcpu,
  initdcpu;
extern float wall_time, ostep_time;
extern struct level **lev;
extern struct cell **carr;
extern struct join **jarr;
extern struct cell_list *freecell;
extern struct join_list *freejoin;
extern struct mcell **mcarr;
#endif
