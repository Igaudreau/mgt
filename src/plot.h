#ifdef PLOT
int nnq, nsa, nqa, nxed, nxad, nl[MAXD], ncd[1][MAXD], *ssa, mergnsa,
  mergnqa, mid, ndm;
float dx[MAXD], dxt[MAXD], xlg[MAXD], xrg[MAXD], lxl[MAXD], *xed, *xad,
  *qqa, mval[MAXD], mvalt[MAXD],xlgt[MAXD], xrgt[MAXD];
char expname[SZNAM], cname[SZNAM];
FILE *expfile;
bool wrgraph, mcoord[MAXD];
#else
extern int nnq, nsa, nqa, nxed, nxad, nl[MAXD], ncd[1][MAXD], *ssa, 
  mergnsa, mergnqa, mid, ndm;
extern float dx[MAXD], dxt[MAXD], xlg[MAXD], xrg[MAXD], lxl[MAXD], *xed,
  *xad, *qqa, mval[MAXD], mvalt[MAXD],xlgt[MAXD], xrgt[MAXD];
extern char expname[SZNAM], cname[SZNAM];
extern FILE *expfile;
extern bool wrgraph, mcoord[MAXD];
#endif

extern "C"
{
  extern void model_flush();
  extern void vis_begin();
  extern void vis_set0r(char name[], float value, int *ios);
  extern void model_set0s(char name[], char value[], int *ios);
  extern void model_seta(char name[], int *ios);
  extern void model_prepare();
  extern void vis_end();
  extern void model_vis(int *);
  extern void model_put_ugrid(char name[], char title[], char geom[],
    float t, int  nd, char coords[], char units[], int nc[][MAXD],
    float xe[], float xa[], char physics[], int nq, char vlabels[],
    char vunits[],  int vmlt[][3], int topn, int sa[], float qa[],
    int topj, int sb[], float qb[], int *ios);
  extern void model_put_scatter(char name[], char title[], float t, 
    char geom[], int nd, char coords[], char units[], char physics[],
    int nq, char vlabel[], char vunits[], int vmlt[][3], int npt, 
    float xpt[], int spt[], float qpt[], int *ios);
  extern void model_remove(char name[], int *ios);
}
