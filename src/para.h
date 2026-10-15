#define SAFETY 10
#define com1 1
#define com2 2
#define com3 3
#define com4 4
#define com5 5
#define com6 6
#define com7 7
#define com8 8
#define com9 9
#define com10 10
#ifdef MAIN
int mype, limpe, *outs, *ins, nla[MAXD];
float pbd[MAXPE+1], levpbd[MAXLEV][MAXPE], pcpu[MAXPE], procwork[MAXPE],
  workwt, *wkda, *wka;
struct cell **locarr;
struct plevel **plev;
struct part_list *pstart;
#else
extern int mype, limpe, *outs, *ins, nla[MAXD];
extern float pbd[MAXPE+1], levpbd[MAXLEV][MAXPE], pcpu[MAXPE], 
  procwork[MAXPE], workwt, *wkda, *wka;
extern struct cell **locarr;
extern struct plevel **plev;
extern struct part_list *pstart;
#endif
struct proc_list{
  int proc;
  struct proc_list *nxt;
};
struct plevel{
  int maxcom, nlproc, lsproc, nrproc, rsproc;
  struct cell_list *lhalo, *lbound, *rbound, *rhalo, *lref, *rref;
  struct join **lfjoin, **rfjoin;
  struct proc_list *lproc, *rproc;
};
/* Parc templates */
extern "C"
{
  void parc_reset();
  void parc_size();
  void parc_args(int argc, char *argv[]);
  void parc_clock(float *t);
  void parc_rank(int *mype_r, int *limpe_r, int *ios);
  void parc_register(int *mype_r, int *limpe_r, int *ios);
  void parc_exit();
  void parc_broadcastr(float buf[], int *ndim, int dims[], int *ios);
  void parc_broadcast0c(char buf[], int *ios);
  void parc_xdatar(int *npe, int outs[], float outbuf[], int ins[],
    float inbuf[], int *ios);
}
