#undef MAIN
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <new>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "funcs.h"
extern bool usrqis(float xc[], float qa[]);
extern void usrqp(float pa[], float qa[]);
extern void usrpq(float qa[], float pa[]);
extern float usrqsol(int ilev, cell *cptr);
extern void usrbd(int id, int bound, join *jptr, float pl[], 
  float pr[]);
extern void usrfxq(struct cell *cptr);
extern void i_to_it(int ndim, int nt[], int icell, int it[]);
extern void i_to_it2(int i, int it[]);
extern int it_to_i(int ndim, int nt[], int it[]);
extern int it2_to_i(int it[]);
extern void clean_cell(struct cell *cptr);
extern void clean_join(struct join *jptr);
extern struct cell *find_cell(int ilev, float xc[]);
extern struct cell *find_fcell(int &ilev, float xc[]);
extern void mk_locarr();
extern void create_halo(int ilev);
extern void connect_halo(int ilev, struct cell_list *shalo);
extern void clear_halo(int ilev);
extern void com_halo(int com, int ilev, float tfrac);
extern void ref_proc_bound(int ilev);
extern void deref_proc_bound(int ilev, struct part_list *&start);
extern void com_cor(int ilev);
extern void set_proc_bound();
extern int re_part();
extern struct part_list *mov_part(bool all, struct part_list *start);
extern void get_part_proc(struct part_list *start);
extern float global_max(float value);
extern void clear_mug();
extern void clean_mug();
void refine(int ilev);
void set_alev();
void proj_init(int ilev, int flev, struct cell *cptr);
void inc_lev(bool doref);
void dec_lev();
bool newdom(float nxl[], float nxr[]);
void ext_bd(int side, int id, int dcell, float dx);
bool extend_bound(float nxl[], float nxr[]); 
void rem_g0();
void merge(int ilev);
void correct_flux(int ilev);
void re_grid(int ilev);
void test(int ilev);
void ref_cell(int ilev, struct cell *pcptr);
void connect(int ilev);
void connect_cell(int ilev, int id, struct join *&start,
  struct join *&end, struct cell *cptr);
void connect_par(int ilev, struct cell *cptr);
void rem_cell(int ilev, struct cell *cptr);
struct cell *new_cell();
struct join *new_join();
void del_cell(struct cell *cptr);
void del_join(struct join *jptr);
void map_up_cell(struct cell *pcptr);
void map_down_cell(bool test, int ilev, struct cell *pcptr);
bool cal_diffql(int id, struct cell *cptr, float diff[]);
bool cal_diffqr(int id, struct cell *cptr, float diff[]);
void mk_bound_join(struct join *jptr);
void umk_bound_join(struct join *jptr);
void set_bound_flag(int ilev);
void set_fine_bound(int ilev, int id, struct join *jptr);
void set_coarse_bound(int ilev, int id, struct join *jptr);
bool check_grid();
void diff_bad(int ilev, int nstep);
void rem_hole(int ilev);
void de_frag();
void num_cell(int ilev, int loff[], struct cell *pcptr);
void q_sort(cell *a[], int first, int last);
void move_cell(int oi, int ni);
void exchange_cell(int oi, int ni);
void del_cell_freelist();
void del_join_freelist();
float cell_vol(int ilev, float xc[]);
float join_area(int ilev, int id, float xj[]);
bool find_q(float xp[], float qa[]);
void mk_partlist(struct part_list *&start, struct part_list *&end);
struct part_list *assign_part(struct part_list *start);
void rm_part(struct part *pptr);
void del_allpart(struct cell *cptr);
//file refine.c
/*
+ refine - refine level ilev to ilev+1
     Description :
     Refines from  level ilev to  ilev+1 using ref_cell and  test.  It
     also calls  set_bound_flag to set the boundary  flags.  The cells
     have  to  be filled  with  the  initial  solution afterwards.  If
     particle=true, then it makes a list of particles in refined cells
     and reassigns them at the end.
*/   
void refine(int ilev)
{
  struct cell *cptr;
  struct part *pptr, *npptr;
  struct part_list *start, *plist, *pmov;
  /* Set processor boundaries */
  set_proc_bound();
  /* Ensure ref_proc_bound lists are NULL */
  plev[ilev+1]->lref = NULL;
  plev[ilev+1]->rref = NULL;
  /* Refine cells */
  cptr = lev[ilev]->fcell;
  /* Initialise particle list (in case) */
  start = NULL;
  if (ilev == 0){
    /* Base level -- refine all cells */
    while (cptr != NULL){
      ref_cell(ilev, cptr);
      if ((particle) && (cptr->pptr != NULL)){
	/* Remove particles and add to list */
	pptr = cptr->pptr;
	while (pptr != NULL){
	  npptr = pptr->nxt;
	  plist = new part_list;
	  plist->pptr = pptr;
	  plist->nxt = start;
	  start = plist;
	  rm_part(pptr);
	  pptr->cptr = NULL;
	  plist->pproc = -1;
	  pptr = npptr;
	}
	/* Set cell particle pointer to NULL */
	cptr->pptr = NULL;
      }
      cptr = cptr->nxt;
    }
  }
  else{
    /* Not base level -- check for refinement */
    test(ilev);
    while (cptr != NULL){
      if (cptr->qsol > refcrit){
	ref_cell(ilev, cptr);
	if ((particle) && (cptr->pptr != NULL)){
	  /* Remove particles and add to list */
	  pptr = cptr->pptr;
	  while (pptr != NULL){
	    npptr = pptr->nxt;
	    plist = new part_list;
	    plist->pptr = pptr;
	    plist->nxt = start;
	    start = plist;
	    rm_part(pptr);
	    pptr->cptr = NULL;
	    plist->pproc = -1;
	    pptr = npptr;
	  }
	  /* Set cell particle pointer to NULL */
	  cptr->pptr = NULL;
	}
      }
      cptr = cptr->nxt;
    }
  }
  if (limpe > 1){
    /* Refine cells on left/right processors */
    ref_proc_bound(ilev+1);
    /* Create ilev+1 halo */
    create_halo(ilev+1);
  }
  /* Set ilev+1 and ilev boundary flags */
  set_bound_flag(ilev+1);
  if (particle){
    /* Assign particles to cells if necessary */
    pmov = assign_part(start);
    if (limpe > 1){
      /* Move cells to other processors if necessary */
      get_part_proc(pmov);
      start = mov_part(false, pmov);
      /* Assign moved particles */
      if (assign_part(start) != NULL)
	cout << "refine pmov not NULL on proc " << mype << '\n';
    }
  }
  return;
}
//end refine.c

//file set_alev.c
/*
+ set_alev - Sets finest active level.
     Description :
        Finds the finest level that has cells. Used to ensure that
        gravity is calculated when nlevs - 1 has no cells.
*/
void set_alev()
{
  int ilev;
  alev = nlevs -1;
  for (ilev = 1; ilev < nlevs; ilev++){
    if (lev[ilev]->fcell != NULL)
      alev = ilev;      
  }
  if (limpe > 1)
    alev = (int)global_max((float)alev);
}
//end set_alev.c

//file proj_init.c
/*
+ proj_init - Projects initial solution onto ilev.
     Description :
        This uses usrqis to project the solution onto the cell cptr on
        grid ilev.  This is done by  sampling the cells  on the finest
        level (nlevs  - 1) and then  averaging over the  cell on level
        ilev.
*/
void proj_init(int ilev, int flev, struct cell *cptr)
{
  int ncell, id, i, k, icell, nt[MAXD], it[MAXD];
  float divncell, dxc[MAXD], xlb[MAXD], xc[MAXD], qt[MAXQ];
  /* No of child cells, cell position and mesh spacings */
  ncell = ipow(nd*(flev - ilev));
  k = ipow(flev - ilev);
  divncell = 1.0/(float)(ncell);
  for (id = 1; id <= nd; id++){
    xlb[id-1] = cptr->xc[id-1] - 0.5*lev[ilev]->dxbar[id-1];
    dxc[id-1] = lev[ilev]->dxbar[id-1]/(float)(k);
    nt[id-1] = k;
  }
  /* Average solution */
  for (i = 0; i < nq; i++)
    cptr->qa[i] = 0.0;
  for (icell = 1; icell <= ncell; icell++){
    i_to_it(nd, nt, icell, it);
    for (id = 1; id <= nd; id++)
      xc[id-1] = xlb[id-1] + (it[id-1] - 0.5)*dxc[id-1];
    usrqis(xc, qt);
    for (i = 0; i < nq; i++){
      cptr->qa[i] += qt[i];
    }
  }
  for (i = 0; i < nq; i++)
    cptr->qa[i] *= divncell;
  return;
}
//end proj_init.c 

//file inc_lev.c
/*
+ inc_lev - increase max level
     Description :
     Creates a new  finest level which is one  finer than the existing
     finest  level.  It uses  the  current  finest  level solution  to
     initialise this.
*/   
void inc_lev(bool doref)
{
  int imult, id;
  struct cell *cptr;
  /* Clear multigrid */
  if (grav || mugi) 
    clear_mug();
  /* Initialise new level */
  imult = ipow(nlevs);
  lev[nlevs] = new level;
  lev[nlevs]->fjoin = new join *[nd];
  lev[nlevs]->ljoin = new join *[nd];
  lev[nlevs]->fcell = NULL;
  lev[nlevs]->lcell = NULL;
  for (id = 1; id <= nd; id++){
    lev[nlevs]->nc[id-1] = imult*lev[0]->nc[id-1];
    lev[nlevs]->dxbar[id-1] = lev[0]->dxbar[id-1]/imult;
    lev[nlevs]->fjoin[id-1] = NULL;
    lev[nlevs]->ljoin[id-1] = NULL;
    lev[nlevs]->xl[id-1] = lev[nlevs-1]->xl[id-1];
    lev[nlevs]->xr[id-1] = lev[nlevs-1]->xr[id-1];
  }
  /* Initialise new plevel */
  plev[nlevs] = new plevel;
  plev[nlevs]->maxcom = 0;
  plev[nlevs]->lbound = NULL;
  plev[nlevs]->rbound = NULL;
  plev[nlevs]->lhalo = NULL;
  plev[nlevs]->rhalo = NULL;
  plev[nlevs]->lproc = NULL;
  plev[nlevs]->rproc = NULL;
  plev[nlevs]->lfjoin = new join *[nd-1];
  plev[nlevs]->rfjoin = new join *[nd-1];
  for (id = 1; id < nd; id++){
    plev[nlevs]->lfjoin[id-1] = NULL;
    plev[nlevs]->rfjoin[id-1] = NULL;
  }
  /* Increment nlevs */
  nlevs++;
  if (limpe > 1){
    /* Expand wka and wkda */
    delete [] wka;
    delete [] wkda;
    wka = new float[nc[nd-1]*ipow(nlevs-1)+SAFETY];
    wkda = new float[nc[nd-1]*ipow(nlevs-1)+SAFETY];
  }
  if (doref){
    /* Refine */
    refine(nlevs - 2);
    /* Check grid */
    if ((bool)global_max((float)!check_grid())){
      if (mype==0){
	cout << "Double refinements - failed to populate new level" 
          << '\n';
	cout << "Cells on new level will be generated by marching"
          << '\n';
      }
      dec_lev();
      inc_lev(false);
    }
    else{
      if (limpe > 1){
	/* Communicate nlevs - 2 halo qa, dqdx (com5, com6)*/
	com_halo(com5, nlevs - 2, 0.0);
	com_halo(com6, nlevs - 2, 0.0);
      }
      /* Initialise solution on finest level */
      cptr = lev[nlevs - 2]->fcell;
      while (cptr != NULL){
	if (cptr->refined)
	  map_down_cell(false, nlevs - 2, cptr);
	cptr = cptr->nxt;
      }
      if (limpe > 1)
	/* Repartition */
	re_part();
      else
	de_frag();
    }
  }
  /* Set plotting levels */
  plmin = 1;
  plmax = nlevs -1;
  return;
}
//end inc_lev.c

//file dec_lev.c
/*
+ dec_lev - decrease max level
     Description :
     Deletes the finest level.
*/   
void dec_lev()
{
  int icell;
  struct cell *cptr, *child;
  struct part *pptr, *npptr;
  struct part_list *start, *plist, *pmov;
  if (nlevs == 1){
    /* Only one level: return */
    if (mype == 0)
      cout << "Error: only one level" << '\n';
  }
  if (grav || mugi) 
    /* Clear multigrid */
    clear_mug();
  /* Set qsol=deref on level nlevs-2 */
  cptr = lev[nlevs-2]->fcell;
  while (cptr != NULL){
    if (cptr->refined)
      cptr->qsol = deref;
    else
      cptr->qsol = good;
    cptr = cptr->nxt;
  }
  start = NULL;
  if (limpe > 1){
    /* Communicate nlevs-2  halo qsol (com1) */
    com_halo(com1, nlevs-2, 0.0);
    /* Remove cells on left/right processors */
    deref_proc_bound(nlevs-1, start);
    clear_halo(nlevs-1);
    /* Ensure ref_proc_bound lists are NULL */
    plev[nlevs-1]->lref = NULL;
    plev[nlevs-1]->rref = NULL;
  }
  /* Derefine cells on level nlevs-1 */
  cptr = lev[nlevs - 2]->fcell;
  while (cptr != NULL){
    if (cptr->qsol == deref){
      for (icell = 0; icell < nrcell; icell++){
	child = cptr->child[icell];
	if ((particle) && (child != NULL) && (child->pptr != NULL)){
	  /* Remove child particles and add to list */
	  pptr = child->pptr;
	  while (pptr != NULL){
	    npptr = pptr->nxt;
	    plist = new part_list;
	    plist->pptr = pptr;
	    plist->nxt = start;
	    start = plist;
	    rm_part(pptr);
	    pptr->cptr = NULL;
	    plist->pproc = -1;
	    pptr = npptr;
	  }
	  /* Set child cell particle pointer to NULL */
	  child->pptr = NULL;
	}
	rem_cell(nlevs - 1, cptr->child[icell]);
      }
    }
    cptr = cptr->nxt;
  }
  nlevs--;
  /* Set pointers to NULL in level array */
  delete [] lev[nlevs]->fjoin;
  delete [] lev[nlevs]->ljoin;
  lev[nlevs]->fjoin = NULL;
  lev[nlevs]->ljoin = NULL;
  lev[nlevs]->fcell = NULL;
  lev[nlevs]->lcell = NULL;
  /* Set pointers to NULL in plevel array */
  plev[nlevs-1]->lbound = NULL;
  plev[nlevs]->rbound = NULL;
  plev[nlevs]->lhalo = NULL;
  plev[nlevs]->rhalo = NULL;
  delete [] plev[nlevs]->lfjoin;
  delete [] plev[nlevs]->rfjoin;
  plev[nlevs]->lfjoin = NULL;
  plev[nlevs]->rfjoin = NULL;
  if (particle){
    /* Assign particles to cells if necessary */
    pmov = assign_part(start);
    if (limpe > 1){
      /* Move particles to other processors if necessary */
      get_part_proc(pmov);
      start = mov_part(false, pmov);
      /* Assign moved particles */
      if (assign_part(start) != NULL)
	cout << "dec_lev pmov not NULL on proc " << mype << '\n';
    }
  }
  if (limpe > 1)
    /* Repartition */
    re_part();
  else
    de_frag();
  /* Set plotting levels */
  plmin = 1;
  plmax = nlevs - 1;
  return;
}
//end dec_lev.c

//file newdom.c
/*
+ newdom - change domain.
     Description :
        Changes the domain to that given  by nxl, nxr. The solution in
        the new domain is simply a  copy of the old solution where the
        new domain is a subset of the  old, but when it is larger, the
        solution is from usrqis.
*/   
bool newdom(float nxl[], float nxr[])
{
  int id, ilev, tnd, dcell_l[MAXD], dcell_r[MAXD];
  bool doit, fail;
  float dx;
  struct cell *cptr, *ncptr;
  /* Adjust new domain to integral number of level 0 cells */
  for (id = 1; id <= nd; id++){
    dx = lev[0]->dxbar[id-1];
    dcell_l[id-1] = (int)((nxl[id-1] - xl[id-1])/dx);
    dcell_r[id-1] = (int)((nxr[id-1] - xr[id-1])/dx);
    nxl[id-1] = xl[id-1] + dx*dcell_l[id-1];
    nxr[id-1] = xr[id-1] + dx*dcell_r[id-1];
    dcell_l[id-1] = abs(dcell_l[id-1]);
    dcell_r[id-1] = abs(dcell_r[id-1]);    
  }
  /* Check whether new domain has less than 2 level 0 cells */
  fail = false;
  for (id = 1; id <= nd; id++){
    dx = lev[0]->dxbar[id-1];
    if ((int)((nxr[id-1] - nxl[id-1] + 0.5*dx)/dx) < 2){
      fail = true;
      break;
    }
  }
  if (fail){
    if (mype ==0)
      cout << "Less than 2 level 0 cells\n";
    return(false);
  }
  /* Clear multigrid */
  if (grav || mugi) 
    clear_mug();
  /* Clear halos and remove cells outside domain */
  for (ilev = nlevs-1; ilev >=0 ; ilev--){
    clear_halo(ilev);
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      doit = false;
      for (id = 1; id <= nd; id++){
	if ((cptr->xc[id-1] < nxl[id-1])
	    || (cptr->xc[id-1] > nxr[id-1])){
	  doit = true;
	  break;
	}
      }
      if (doit){
	ncptr = cptr->nxt;
	rem_cell(ilev, cptr);
	cptr = ncptr;
      }
      else
	cptr = cptr->nxt;
    }
  }
  /* Create new cells where new domain boundary is outside old */
  for (id = 1; id <= nd; id++){
    dx = lev[0]->dxbar[id-1];			// BM
    if (nxl[id-1] < xl[id-1])
      ext_bd(-1, id, dcell_l[id-1], dx);
    else if (nxr[id-1] > xr[id-1])
      ext_bd(1, id, dcell_r[id-1], dx);
  }
  /* Set new boundaries and no of cells for level 0 in each direction  */
  for (id = 1; id <= nd; id++){
    xl[id-1] = nxl[id-1];
    xr[id-1] = nxr[id-1];
    dx = lev[0]->dxbar[id-1];
    nc[id-1] = (int)((xr[id-1] - xl[id-1] + 0.5*dx)/dx);
  }
  if (limpe > 1)
    tnd = nd - 1;
  else
    tnd = nd;  
  /* Set lev[ilev]->nc etc  */
  for (ilev = 0; ilev < nlevs; ilev++){
    for (id = 1; id <= tnd; id++){
      lev[ilev]->nc[id-1] = nc[id-1]*ipow(ilev);
      lev[ilev]->xl[id-1] = nxl[id-1];
      lev[ilev]->xr[id-1] = nxr[id-1];
    }
  }  
  /* Repartition and defrag */
  if (limpe > 1)
    re_part();
  else
    de_frag();    
  return(true);
}
//end newdom.c

//file ext_bd.c
/*
+ ext_bd - extend domain boundary.
     Description :
        Adds level 0 cells to extend the boundary given by id. side is
        -1 for  left boundary,  +1 for right  boundary. The  cells are
        filled with the solution from usrqis. and refined to level 1.
*/   
void ext_bd(int side, int id, int dcell, float dx)
{
  int tid, ncell, icell, nct[MAXD], it[MAXD];
  struct cell *cptr, *old_lcell;
  /* Check parallel for direction nd */
  if (limpe > 1){
    if (id == nd){
      if   (((side == -1) && (mype != 0)) ||
	    ((side ==  1) && (mype != limpe-1)))
	return;
    }
  }
  /* Add level 0 cells */
  cptr = lev[0]->lcell;
  if (cptr != NULL){
    ncell = 1;
    for (tid = 1; tid <= nd; tid++){
      if (tid == id)
	nct[tid-1] = dcell;
      else
	nct[tid-1] = lev[0]->nc[tid-1];
      ncell *= nct[tid-1];
    }
    for (icell = 1; icell <= ncell; icell++){
      cptr->nxt = new_cell();
      cptr = cptr->nxt;
      i_to_it(nd, nct, icell, it);
      for (tid = 1; tid <= nd; tid++){
	if (tid == id){
	  if (side == -1)
	    cptr->xc[tid-1] = lev[0]->xl[tid-1]
	      + (it[tid-1] - 0.5 - dcell)*dx;
	  else
	    cptr->xc[tid-1] = lev[0]->xr[tid-1]
	      + (it[tid-1] - 0.5)*dx;
	}
	else
	  cptr->xc[tid-1] = lev[0]->xl[tid-1]
	    + (it[tid-1] - 0.5)*dx;
      }
      cptr->vol = cell_vol(0, cptr->xc);
      /* Fill cells with solution */
      usrqis(cptr->xc, cptr->qa);
    }
    cptr->nxt = NULL;
    old_lcell = lev[0]->lcell;
    lev[0]->lcell = cptr;
    if (nlevs > 1){
      /* Ensure ref_proc_bound lists are NULL */
      plev[1]->lref = NULL;
      plev[1]->rref = NULL;
      /* Refine level new 0 cells */
      cptr = old_lcell->nxt;
      while (cptr != NULL){
	ref_cell(0, cptr);
	/* Fill child cells with solution */
	for (icell = 0; icell < nrcell; icell ++)
	  usrqis(cptr->child[icell]->xc, cptr->child[icell]->qa);
	cptr = cptr->nxt;
      }
    }
  }
}
//end ext_bd.c

//file rem_g0.c
/*
+ rem_g0 - remove level 0.
     Description :
        Removes level 0.
*/   
void rem_g0()
{
  int icell, ilev, id, olevc, lcell;
  struct cell *cptr;
  /* Ensure level 1 cells are refined */
  olevc = levc;
  levc = 2;
  re_grid(0);
  levc = olevc;
  /* Clear halos */
  for (ilev = 0; ilev < nlevs; ilev++)
    clear_halo(ilev);
  if (lev[0]->fcell != NULL){
    /* Remove level 0 cells, */
    lcell = lev[0]->lcell->icell;
    for (icell = 1; icell <= lcell; icell++)
      clean_cell(carr[icell]);
  }
  /* Reset lev structure */
  for (ilev = 1; ilev < nlevs; ilev ++){
    lev[ilev-1]->fcell = lev[ilev]->fcell;
    lev[ilev-1]->lcell = lev[ilev]->lcell;
    for (id = 1; id <= nd; id++){
      lev[ilev-1]->dxbar[id-1] = lev[ilev]->dxbar[id-1];
      lev[ilev-1]->xl[id-1] = lev[ilev]->xl[id-1];
      lev[ilev-1]->xr[id-1] = lev[ilev]->xr[id-1];
      lev[ilev-1]->nc[id-1] = lev[ilev]->nc[id-1];
    }
  }
  /* Ensure level 0 parent pointers are NULL */
  cptr = lev[0]->fcell;
  while (cptr != NULL){
    cptr->par = NULL;
    cptr = cptr->nxt;
  }
  /* Set new level 0 nc */
  for (id = 1; id <= nd; id++)
    nc[id-1] *= 2;
  /* New max level */
  nlevs--;
  /* Repartition and defrag */
  if (limpe > 1)
    re_part();
  else
    de_frag();
  /* Reset plmax */
  plmax = min(plmax, nlevs-1);
  return;
}
//end rem_g0.c

//file merge.c
/*
+ merge - merge cells from ilev+1 to ilev.
     Description :
     Merges solution from child cells  if these exist.
*/
void merge(int ilev)
{
  struct cell *cptr;
  if ((debug) && (mype == 0)){
    cout << "merging " << ilev << " with " << ilev + 1 << '\n';
  }
  if (limpe > 1)
    /* Communicate ilev+1 halo qa (com5)*/
    com_halo(com5, ilev+1, 0.0);
  /* Map up solution */
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    if (cptr->refined)
      /* Map up solution */
      map_up_cell(cptr);
    cptr = cptr->nxt;
  }
  return;
}
//end merge.c

//file correct_flux.c
/*
+ correct_flux - apply flux correction to  ilev.
     Description :
     Apply  the  flux correction  at  coarse/fine boundaries.
*/
void correct_flux(int ilev)
{ 
  int id, i;
  float dvol;
  struct cell *left, *right;
  struct join *jptr;
  if ((ilev < nlevs - 1) && (limpe > 1)) 
    /* Communicate halo coarse/fine correction */
    com_cor(ilev);
  /* Apply flux correction */
  for (id = 1; id <= nd; id++){
    jptr = lev[ilev]->fjoin[id-1];
    while (jptr != NULL){
      if (jptr->bound == cibound){
	/* Coarse/fine boundary */
	left = jptr->left;
	if (!left->refined){
	  /* Left cell unrefined -- apply correction */
	  dvol = 1.0/left->vol;
	  for (i = 0; i < nq; i++)
	    left->qa[i] -= dvol*jptr->cor[i];
	  usrfxq(left);
	}
	right = jptr->right;
	if (!right->refined){
	  /* Right cell unrefined -- apply correction */
	  dvol = 1.0/right->vol;
	  for (i = 0; i < nq; i++)
	    right->qa[i] += dvol*jptr->cor[i];
	  usrfxq(right);
	}
      }
      jptr = jptr->nxt;
    }
  }
  return;
}
//end correct_flux.c

//file re_grid.c
/*
+ re_grid - refine/derefine on level ilev+1
     Description :
     Checks  difference in  solution between  levels ilev  and ilev+1.
     The result is used to  refine/derefine cells on ilev+2. Cells are
     not destroyed  if they have children.  If  particle=true, then it
     makes a list  of particles in that need to  be moved to different
     cells and reassigns them at the end.
*/
void re_grid(int ilev)
{
  int id, icell, tilev;
  bool doit;
  //char temp[SZNAM];
  struct cell *cptr, *child;
  //struct cell_list *clist;
  struct part *pptr, *npptr;
  struct part_list *start, *plist, *pmov;
  /* Check level (safety) */
  if (ilev > nlevs - 3)
    return;
  /* Test ilev+1 */
  test(ilev+1);
  /* Mark ilev+1 cells for refinement/derefinement */
  cptr = lev[ilev+1]->fcell;
  while (cptr != NULL){
    if (cptr->qsol > refcrit){
      if (!cptr->refined){
	doit = true;
	/* Refinement must not create double refinements */
	for (id = 1; id <= nd; id++){
	  if ((cptr->ljoin[id-1]->left == NULL) && 
	    (cptr->ljoin[id-1]->bound != lebound))
	    doit = false;
	  if ((cptr->rjoin[id-1]->right == NULL) && 
	    (cptr->rjoin[id-1]->bound != rebound))
	    doit = false;
	}
	if (doit)
	  cptr->qsol = ref;
	else
	  cptr->qsol = good;
      }
    }
    else if (cptr->qsol < derefcrit)
      if (cptr->refined){
	doit = true;
	/* Derefinement must not create double refinements */
	for (icell = 0; icell < nrcell; icell++){
	  child = cptr->child[icell];
	  /* Check for grandchildren */
	  if (child->refined){
	    doit = false;
	    break;
	  }
	  else{
	    for (id = 1; id <= nd; id++){
	      if (child->ljoin[id-1]->bound == cibound)
		doit = false;
	      if (child->rjoin[id-1]->bound == cibound)
		doit = false;
	    }
	  }
	}
	if (doit)
	  cptr->qsol = deref;
	else
	  cptr->qsol = good;
      }
    cptr = cptr->nxt;
  }
  /* Initialise particle list (in case) */
  start = NULL;
  if (limpe > 1){
    /* Communicate ilev+1 halo qsol (com1)*/
    com_halo(com1, ilev+1, 0.0);
    /* Remove cells on left/right processors */
    deref_proc_bound(ilev+2, start);
    clear_halo(ilev+2);
  }
  /* Ensure ref_proc_bound lists are NULL */
  plev[ilev+2]->lref = NULL;
  plev[ilev+2]->rref = NULL;
  /* Refine/ derefine ilev+1 cells */
  cptr = lev[ilev+1]->fcell;
  while (cptr != NULL){
    if (cptr->qsol == ref){
      ref_cell(ilev+1, cptr);
      map_down_cell(false, ilev+1, cptr);
      if ((particle) && (cptr->pptr != NULL)){
	/* Remove particles and add to list */
	pptr = cptr->pptr;
	while (pptr != NULL){
	  npptr = pptr->nxt;
	  plist = new part_list;
	  plist->pptr = pptr;
	  plist->nxt = start;
	  start = plist;
	  rm_part(pptr);
	  pptr->cptr = NULL;
	  plist->pproc = -1;
	  pptr = npptr;
	}
	/* Set cell particle pointer to NULL */
	cptr->pptr = NULL;
      }
    }
    else if (cptr->qsol == deref){
      for (icell = 0; icell < nrcell; icell++){
	child = cptr->child[icell];
	if ((particle) && (child != NULL) && (child->pptr != NULL)){
	  /* Remove child particles and add to list */
	  pptr = child->pptr;
	  while (pptr != NULL){
	    npptr = pptr->nxt;
	    plist = new part_list;
	    plist->pptr = pptr;
	    plist->nxt = start;
	    start = plist;
	    rm_part(pptr);
	    pptr->cptr = NULL;
	    plist->pproc = -1;
	    pptr = npptr;
	  }
	  /* Set child cell particle pointer to NULL */
	  child->pptr = NULL;
	}
	rem_cell(ilev+2, child);
      }
    }
    cptr = cptr->nxt;
  }
  if (limpe > 1){
    /* Refine cells on left/right processors */
    ref_proc_bound(ilev+2);
    /* Recreate ilev+1 to nlevs-1 halos */
    for (tilev = ilev + 1; tilev < nlevs; tilev++)
      create_halo(tilev);
    /* Communicate halo refinement (com8) */
    com_halo(com8, ilev + 1, 0.0);
    /* Set boundary flags */
    for (tilev = ilev + 1; tilev < nlevs; tilev++)
      set_bound_flag(tilev);
  }
  else{
    /* Just set ilev+1 and ilev+2 boundary flags */
    set_bound_flag(ilev+2);
  }
  if (particle){
    /* Assign particles to cells if necessary */
    pmov = assign_part(start);
    if (limpe > 1){
      /* Move particles to other processors if necessary */
      get_part_proc(pmov);
      start = mov_part(false, pmov);
      /* Assign moved particles */
      if (assign_part(start) != NULL)
	cout << "re_grid pmov not NULL on proc " << mype << '\n';
    }
  }
  /* Reset alev */
  set_alev();
  return;
}
//end re_grid.c

//file test.c
/*
+ test - test cells for refinement.
     Description :
        This compares the solutions on ilev and ilev-1 and marks cells
        on ilev as not good  if the difference exceeds a tolerance. It
        then  optionally diffuses qsol (diffbad = true).
*/
void test(int ilev)
{
  int icell;
  struct cell *pcptr, *cptr;
  struct cell_list *clist;
  if (limpe > 1){
    /* Communicate ilev - 1 and ilev halo qa, dqdx (com5, com6) */
    com_halo(com5, ilev-1, 0.0);
    com_halo(com6, ilev-1, 0.0);
    com_halo(com5, ilev, 0.0);
    com_halo(com6, ilev, 0.0);
  }
  /* First pass -- compare solutions and mark bad cells */
  pcptr = lev[ilev-1]->fcell;
  while (pcptr != NULL){
    if (pcptr->refined){
      /* ilev cell exists - compare solutions */
      map_down_cell(true, ilev-1, pcptr);
      for (icell = 0; icell < nrcell; icell++){
	cptr = pcptr->child[icell];
	if (ilev < levc)
	  cptr->qsol = badcell;
	else
	  cptr->qsol = usrqsol(ilev, cptr);
      }
    }
    pcptr = pcptr->nxt;
  }
  if (limpe > 1){
    /* Do left halo cells */
    clist = plev[ilev-1]->lhalo;
    while (clist != NULL){
      pcptr = clist->cptr;
      if (pcptr->refined){
	/* ilev cell exists - compare solutions */
	map_down_cell(true, ilev-1, pcptr);
	for (icell = 0; icell < nrcell; icell++){
	  cptr = pcptr->child[icell];
	  if (cptr != NULL) {
	    if (ilev < levc)
	      cptr->qsol = badcell;
	    else
	      cptr->qsol = usrqsol(ilev, cptr);
	  }
	}
      }
      clist = clist->nxt;
    }
    /* Do right halo cells */
    clist = plev[ilev-1]->rhalo;
    while (clist != NULL){
      pcptr = clist->cptr;
      if (pcptr->refined){
	/* ilev cell exists - compare solutions */
	map_down_cell(true, ilev-1, pcptr);
	for (icell = 0; icell < nrcell; icell++){
	  cptr = pcptr->child[icell];
	  if (cptr != NULL)
	    cptr->qsol = usrqsol(ilev, cptr);
	}
      }
      clist = clist->nxt;
    }
  }
  /* Additional passes */
  if (remhole){
    if (limpe > 1)
      /* Communicate ilev halo qsol */
      com_halo(com1, ilev, 0.0);
    /* Remove holes */
    rem_hole(ilev);
    if (limpe > 1)
      /* Communicate ilev halo qsol */
      com_halo(com1, ilev, 0.0);
  }
  if (diffbad)
    /* Diffuse badness */
    diff_bad(ilev, ndiffstep);
  return;
}
//end test.c

//file ref_cell.c
/*
+ ref_cell - refine a cell on level ilev to ilev+1
     Description :
     Creates ilev+1 child cells for a cell on level ilev.  It connects
     these  to their  neighbours  and creates  the appropriate  joins.
     Cells   on   the   left   processor   are  put   in   list   with
     plev[ilev+1]->lref  pointing to  start of  list.  Cells  on right
     processor are pointed to  by plev[ilev+1]->rref.  Note that these
     must be set to NULL before  the first call to ref_cell. The cells
     then  have to  be  filled by  map_down_cell  and the  information
     communicated to the  neighbouring processors by ref_proc_bound.c,
     which then  creates the cells  on the neighbouring  processor and
     removes  them  from  this  processor.   It  does  not  deal  with
     particles.
*/
void ref_cell(int ilev, cell *pcptr)
{
  int id, icell, it[MAXD];
  float dx[MAXD], xlnd, xrnd;
  struct cell *cptr;
  struct cell_list *clist;
  /* Check whether cell is already refined */
  if (pcptr->refined)
    return;
  /* Cell size */
  for (id = 1; id <= nd; id++)
    dx[id-1] = lev[ilev+1]->dxbar[id-1];
  /* Check start of list */
  if (lev[ilev+1]->fcell == NULL){ 
    /* First cells on this level */
    cptr = new_cell();
    cptr->prev = NULL;
    lev[ilev+1]->fcell = cptr;
  }
  else{
    /* Level already has cells */
    cptr = lev[ilev+1]->lcell;
    cptr->nxt = new_cell();
    cptr->nxt->prev = cptr;
    cptr = cptr->nxt;
  }
  cptr->par = pcptr;
  pcptr->child[0] = cptr;
  /* Create and link cells to cell list */
  for (icell = 1; icell < nrcell; icell ++){
    cptr->nxt = new_cell();
    cptr->nxt->prev = cptr;
    cptr = cptr->nxt;
    cptr->par = pcptr;
    pcptr->child[icell] = cptr;
  }
  lev[ilev+1]->lcell = cptr;
  cptr->nxt = NULL;
  /* Parent is now refined */
  pcptr->refined = true;
  /* Cell coordinates and volume */
  for (icell = 0; icell < nrcell; icell++){
    cptr = pcptr->child[icell];
    cptr->ichild = icell;
    i_to_it2(icell+1, it);
    /* Cell coordinates */
    for (id = 1; id <= nd; id++)
      cptr->xc[id-1] = pcptr->xc[id-1] + (it[id-1] - 1.5)*dx[id-1];
    /* Set cell volume */
    cptr->vol = cell_vol(ilev+1, cptr->xc);
  }
  /* Connect cells */
  for (id = 1; id <= nd; id++)
    for (icell = 0; icell < nrcell; icell++)
      connect_cell(ilev+1, id, lev[ilev+1]->fjoin[id-1],
        lev[ilev+1]->ljoin[id-1], pcptr->child[icell]);
  if (limpe > 1){
    /* Check for cells on neighbouring processors */
    xlnd = lev[ilev+1]->xl[nd-1];
    xrnd = lev[ilev+1]->xr[nd-1];
    for (icell = 0; icell < nrcell; icell++){
      cptr = pcptr->child[icell];
      if (cptr->xc[nd-1] < xlnd){
	/* Left processor cells */
	clist = new cell_list;
	clist->cptr = cptr;
	clist->nxt = plev[ilev+1]->lref;
	plev[ilev+1]->lref = clist;
      }
      if (cptr->xc[nd-1] > xrnd){
	/* Right processor cells */
	clist = new cell_list;
	clist->cptr = cptr;
	clist->nxt = plev[ilev+1]->rref;
	plev[ilev+1]->rref = clist;
      }
    }
  }
  return;
}
//end ref_cell.c

//file connect.c
/*
+ connect - Connect level ilev.
     Description :
     Connects all cells on level  ilev to their neighbours on the same
     level and  to their  parents if ilev  > 0. Assumes that the joins
     on this level have been cleared and that the cell list exists.
*/
void connect(int ilev)
{
  int id;
  struct cell *cptr;
  /* Set join list pointers to NULL */
  for (id = 1; id <= nd; id++){
    lev[ilev]->fjoin[id-1] = NULL;
    lev[ilev]->ljoin[id-1] = NULL;
  }
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    if (ilev > 0)
      /* Connect to parent */
      connect_par(ilev, cptr);
    /* Connect cell */
    for (id = 1; id <= nd; id++)
      connect_cell(ilev, id, lev[ilev]->fjoin[id-1], 
        lev[ilev]->ljoin[id-1], cptr);
    cptr = cptr->nxt;
  }
  return;
}
//end connect.c

//file connect_cell.c
/*
+ connect_cell - Connect a cell.
     Description :
     Connects  the  cell on  level  ilev pointed  to  by  cptr to  its
     neighbours in  the id direction, creating joins  as necessary. It
     uses find_cell,  which assumes that ilev-1 exists.  The joins are
     added to the end  of the list pointed to by end.  start is set to
     the start of the list.
*/
void connect_cell(int ilev, int id, struct join *&start, 
  struct join *&end, struct cell *cptr)
{
  int tid;
  bool add;
  float dx, xc[MAXD];
  struct cell *left, *right;
  struct join *jptr;
  /* Mesh spacing */
  dx = lev[ilev]->dxbar[id-1];
  /* Check for left cell */
  for (tid = 1; tid <= nd; tid++)
    xc[tid-1] = cptr->xc[tid-1];
  xc[id-1] -= dx;
  /* Check whether join should be in list */
  add = !((id == nd) && (cptr->halo));
  left = find_cell(ilev, xc);
  jptr = NULL;
  if (left != NULL)
    /* Left cell exists */
    jptr = left->rjoin[id-1];
  if (jptr == NULL){
    /* No left join - create it */
    jptr = new_join();
    if (add){
      /* Not halo nd join - add to list */
      if (start == NULL){
	/* First join */
	start = jptr;
	end = jptr;
	jptr->prev = NULL;
      }
      else{
	/* Joins already exist */
	end->nxt = jptr;
	jptr->prev = end;
	end = jptr;
	jptr->nxt = NULL;
      }
    }
    /* Join coords and area */
    xc[id-1] += 0.5*dx;
    for (tid = 1; tid <= nd; tid++)
      jptr->xj[tid-1] = xc[tid-1];
    jptr->area = join_area(ilev, id, jptr->xj);
  }
  /* Connect left join */
  jptr->right = cptr;
  cptr->ljoin[id-1] = jptr; 
  jptr->left = left;

  /* Check for right cell */
  xc[id-1] = cptr->xc[id-1] + dx;
  right = find_cell(ilev, xc);
  jptr = NULL;
  if (right != NULL)
    /* Right cell exists */
    jptr = right->ljoin[id-1];
  if (jptr == NULL){
    /* No right join - create it */
    jptr = new_join();
    if (add){
      /* Not halo nd join - add to list */
      if (start == NULL){
	/* First join */
	start = jptr;
	end = jptr;
	jptr->prev = NULL;
      }
      else{
	/* Joins already exist */
	end->nxt = jptr;
	jptr->prev = end;
	end = jptr;
	jptr->nxt = NULL;
      }
    }
    /* Join coords and area */
    xc[id-1] -= 0.5*dx;
    for (tid = 1; tid <= nd; tid++)
      jptr->xj[tid-1] = xc[tid-1];
    jptr->area = join_area(ilev, id, jptr->xj);
  }
  /* Connect right join */
  jptr->left = cptr;
  cptr->rjoin[id-1] = jptr; 
  jptr->right = right;
  /* Fix halo join cell pointers at domain boundaries */
  if (bcl[nd-1] != bcper){
    if ((cptr->halo) && (id == nd)){
      if (cptr->xc[nd-1] < xl[nd-1])
	cptr->rjoin[nd-1]->left = NULL;
      else if (cptr->xc[nd-1] > xr[nd-1])
	cptr->ljoin[nd-1]->right = NULL;
    }
  }
  return;
}
//end connect_cell.c

//file connect_par.c
/*
+ connect_par - Connect a cell to its parent.
     Description :
     Connects  the  cell on  level  ilev pointed  to  by  cptr to  its
     parent.
*/
void connect_par(int ilev, struct cell *cptr)
{
  int id, it[MAXD];
  struct cell *pcptr;
  /* Check level (safety) */
  if (ilev == 0){
    cptr->par = NULL;
    return;
  }
  /* Find parent */
  pcptr = find_cell(ilev-1, cptr->xc);
  cptr->par = pcptr;
  /* Connect if parent exists */
  if (pcptr != NULL){
    for (id = 1; id <= nd; id++){
      if (cptr->xc[id-1] < pcptr->xc[id-1])
	it[id-1] = 1;
      else
	it[id-1] = 2;
    }
    cptr->ichild = it2_to_i(it);
    pcptr->child[cptr->ichild] = cptr;
    pcptr->refined = true;
  }
  return;
}
//end connect_par.c

//file rem_cell.c
/*
+ rem_cell - Remove a cell.
     Description :
     Removes  the  cell  on  level   ilev  pointed  to  by  cptr.   If
     neighbouring cells exist  on this level, then it  leaves the edge
     joins  intact and  sets the  appropriate pointers  to  NULL.  The
     parent pointers  are set appropriately. It updates  the cell list
     if the cell  is in the list i.e.  if this is not a  halo cell. It
     can be used when the cell is not connected.
*/
void rem_cell(int ilev, struct cell *cptr)
{
  int id, icell;
  bool refined;
  struct cell *lcptr, *ncptr, *pcptr;
  struct join *jptr, *njptr;
  /* Return if NULL or empty */
  if (cptr == NULL)
    return;
  else if (cptr->empty)
    return;
  /* Delete joins if necessary */
  for (id = 1; id <= nd; id++){
    /* Left join */
    jptr = cptr->ljoin[id-1];
    if (jptr != NULL){
      if ((jptr->left == NULL) || (jptr->left->halo)){
	/* No left cell or left is halo cell */
	if (!cptr->halo){
	  /* Not halo cell - update join list */
	  njptr = jptr->nxt;
	  if (jptr->prev != NULL){
	    if (njptr != NULL){
	      jptr->prev->nxt = njptr;
	      njptr->prev = jptr->prev;
	    }
	    else{
	      /* Last join */
	      jptr->prev->nxt = NULL;
	      lev[ilev]->ljoin[id-1] = jptr->prev;
	    }
	  }
	  else{
	    /* First join */
	    lev[ilev]->fjoin[id-1] = njptr;
	    if (njptr != NULL)
	      njptr->prev = NULL;
	  }
	}
	del_join(jptr);
      }
      else
	/* Left cell exists */
	jptr->right = NULL;
    }
    /* Right join */
    jptr = cptr->rjoin[id-1];
    if (jptr != NULL){
      if ((jptr->right == NULL) || (jptr->right->halo)){
	/* No right cell  or right is halo cell */
	if (!cptr->halo){
	  /* Not halo cell - update join list */
	  njptr = jptr->nxt;
	  if (jptr->prev != NULL){
	    if (njptr != NULL){
	      jptr->prev->nxt = njptr;
	      njptr->prev = jptr->prev;
	    }
	    else{
	      /* Last join */
	      jptr->prev->nxt = NULL;
	      lev[ilev]->ljoin[id-1] = jptr->prev;
	    }
	  }
	  else{
	    /* First join */
	    lev[ilev]->fjoin[id-1] = njptr;
	    if (njptr != NULL)
	      njptr->prev = NULL;
	  }
	}
	del_join(jptr);
      }
      else
	/* Right parent cell is refined */
	jptr->left = NULL;
    }
  }
  /* Set parent cell pointer */
  pcptr = cptr->par;
  if (pcptr != NULL){
    pcptr->child[cptr->ichild] = NULL;
    refined = false;
    for (icell = 0; icell < nrcell; icell++)
      if (pcptr->child[icell] != NULL)
	refined = true;
    pcptr->refined = refined;
  }
  /* Set child parent pointer */
  if (cptr->refined){
    for (icell = 0; icell < nrcell; icell++){
      if (cptr ->child[icell] != NULL)
	cptr->child[icell]->par = NULL;
    }
  }
  /* Prune cell list */
  lcptr = cptr->prev;
  ncptr = cptr->nxt;
  if (lcptr != NULL){
    if (ncptr != NULL){
      lcptr->nxt = ncptr;
      ncptr->prev = lcptr;
    }
    else if (cptr == lev[ilev]->lcell){
      /* Last cell */
      lcptr->nxt = NULL;
      lev[ilev]->lcell = lcptr;
    }
  }
  else if (cptr == lev[ilev]->fcell){
    /* First cell */
    lev[ilev]->fcell = ncptr;
    if (ncptr != NULL)
      ncptr->prev = NULL;
    else
      lev[ilev]->lcell = NULL;
  }
  /* Set particle pointer to NULL */
  cptr->pptr = NULL;
  /* Delete cell */
  del_cell(cptr);
  return;
}
//end rem_cell.c

//file new_cell.c
/*
+ new_cell - Create a new cell.
     Description :
     Grabs memory for  a new cell from the free list  if it exists and
     from  carr if  it does  not. Returns a pointer to the new cell.
*/
struct cell *new_cell()
{
  int icell;
  struct cell *cptr;
  struct cell_list *clist;
  /* Check free list */
  if (freecell == NULL){
    /* Grab memory from carr */
    if (carrtop < celltop + 100){
      /* Resize carr array */
      carr = (struct cell **)realloc((void *)carr, 
        (size_t)((carrtop + CINCR + SAFETY)*sizeof(cell *)));
      for (icell = carrtop + 1; icell <= carrtop + CINCR; icell++){
	carr[icell] = new cell;
	clean_cell(carr[icell]);
	carr[icell]->icell = icell;
      }
      carrtop += CINCR;
      if ((mype == 0) && (debug))
        cout << "increasing carrtop, old " <<  carrtop - CINCR <<
	  " new " << carrtop << '\n';
    }
    celltop++;
    cptr = carr[celltop];
    cptr->icell = celltop;
  }
  else{
    /* Grab memory from free list and update */
    cptr = freecell->cptr;
    clist = freecell;
    freecell = freecell->nxt;
    delete clist;
  }
  cptr->empty = false;
  return(cptr);
}
//end new_cell.c

//file new_join.c
/*
+ new_join - Create a new join.
     Description :
     Grabs memory for  a new join from the free list  if it exists and
     from  jarr if  it does  not. Returns a pointer to the new join.
*/
struct join *new_join()
{
  int ijoin;
  struct join *jptr;
  struct join_list *jlist;
  /* Check free list */
  if (freejoin == NULL){
    /* Grab memory from jarr */
    if (jarrtop < jointop + 100){
      /* Resize jarr array */
      jarr = (struct join **)realloc((void *)jarr, 
        (size_t)((jarrtop + JINCR + SAFETY)*sizeof(join *)));
      for (ijoin = jarrtop + 1; ijoin <= jarrtop + JINCR; ijoin++){
	jptr = new join;
	/* Ensure that joins are clean */
	jptr->pa0 = NULL;
	jptr->dpdt = NULL;
	jptr->dpdx = NULL;
	jptr->cor = NULL;
	jptr->bound = internal;
	jptr->prev = NULL;
	jptr->nxt = NULL;
	jptr->par = NULL;
	jptr->left = NULL;
	jptr->right = NULL;
	jptr->empty = true;
	jptr->ijoin = ijoin;
	jptr->empty = true;
	jarr[ijoin]= jptr;
      }
      jarrtop += JINCR;
      if ((mype == 0) && (debug))
        cout << "increasing jarrtop, old " <<  jarrtop - JINCR <<
	  " new " << jarrtop << '\n';
    }
    jointop++;
    jptr = jarr[jointop];
  }
  else{
    /* Grab memory from free list and update */
    jptr = freejoin->jptr;
    jlist = freejoin;
    freejoin = freejoin->nxt;
    delete jlist;
  }
  jptr->empty = false;
  return(jptr);
}
//end new_join.c

//file del_cell.c
/*
+ del_cell - Delete a cell.
     Description :
     Frees memory for cell cptr and puts it in the free list.
*/
void del_cell(struct cell *cptr)
{
  struct cell_list *tptr;
  if (cptr !=NULL){
    if (!cptr->empty){
      /* Add cell to free list */
      tptr = new cell_list;
      tptr->nxt = freecell;
      tptr->cptr = cptr;
      freecell = tptr;
    }
    /* Clean cell */
    clean_cell(cptr);
  }
  return;
}
//end del_cell.c

//file del_join.c
/*
+ del_join - Delete a join.
     Description :
     Frees memory for join jptr and puts it in the free list.
*/
void del_join(struct join *jptr)
{
  struct join_list *tptr;
  if (jptr != NULL){
    if (!jptr->empty){
      /* Add join to free list */
      tptr = new join_list;
      tptr->nxt = freejoin;
      tptr->jptr = jptr;
      freejoin = tptr;
    }
    /* Clean join */
    clean_join(jptr);
  }
  return;
}
//end del_join.c

//file map_up_cell.c
/*
+ map_up_cell - maps up solution from child cells.
     Description :
        Maps up  the solution from  the child cells given  the pointer
        pcptr to the parent cell.
*/
void map_up_cell(cell *pcptr)
{
  int icell, i;
  float divvolp, volr, qa[MAXQ];
  struct cell *cptr;
  /* Map up solution */
  divvolp = 1.0/pcptr->vol;
  for (i = 0; i < nq; i++)
    qa[i] = 0.0;
  for (icell = 0; icell < nrcell; icell++){
    cptr = pcptr->child[icell];
    /* Volume ratio */
    volr = cptr->vol*divvolp;
    /* Solution */
    for (i = 0; i < nq; i++)
      qa[i] += cptr->qa[i]*volr;
  }
  for (i = 0; i < nq; i++)
    pcptr->qa[i] = qa[i];
  return;
}
//end map_up_cell.c

//file map_down_cell.c
/*
+ map_down_cell - maps down solution from parent cell
     Description :
        Maps  down the  solution from  the parent  cell pointed  to by
        pcptr to its child cells. If test = true, then it does not use
        the averaging  function and puts  the mapped down  solution in
        pa.  Otherwise it  uses the  averaging function  and  puts the
        solution in qa.
*/
void map_down_cell(bool test, int ilev, cell *pcptr)
{
  int id, i, rind, icell;
  float pdivdx, pxc[MAXD], cxc[MAXD], dqdx[MAXD][MAXQ], diffl[MAXQ],
    diffr[MAXQ], qc[MAXQ], vc, vi, vo, tv;
  bool gotl, gotr;
  struct cell *cptr;
  /* Get gradients */
  for (id = 1; id <= nd; id++){
    /* Coarse cell coordinate */
    pxc[id-1] = pcptr->xc[id-1];
    /* Coarse level inverse mesh spacing */
    pdivdx = 1.0/lev[ilev]->dxbar[id-1];
    /* Left difference */
    gotl = cal_diffql(id, pcptr, diffl);
    /* Right difference */
    gotr = cal_diffqr(id, pcptr, diffr);
    if (gotl && gotr){
      if (test){
	for (i = 0; i < nq; i++)
	  dqdx[id-1][i] = 0.5*(diffl[i] + diffr[i])*pdivdx;
      }
      else{
	for (i = 0; i < nq; i++)
	  dqdx[id-1][i] = av(diffl[i], diffr[i])*pdivdx;
      }
    }
    else if (gotl){
      for (i = 0; i < nq; i++)
	dqdx[id-1][i] = diffl[i]*pdivdx;
    }
    else if (gotr){
      for (i = 0; i < nq; i++)
	dqdx[id-1][i] = diffr[i]*pdivdx;
    }
    else{
      for (i = 0; i < nq; i++)
	dqdx[id-1][i] = 0.0;
    }
  }
  /* Set central qa */
  for (i = 0; i < nq; i++)
    qc[i] = pcptr->qa[i];
  if (strcmp(geom, "RZP") == 0){
    /* Cylindrical polars - correct for cell volume */
    vc = pcptr->vol;
    /* Inner child cell volume */
    if (pcptr->child[0] != NULL)
      vi = pcptr->child[0]->vol;
    else{
      for (id = 1; id <= nd; id++)
	cxc[id-1] = pcptr->xc[id-1] - 0.5*lev[ilev+1]->dxbar[id-1];
      vi = cell_vol(ilev+1, cxc);
    }
    /* Outer child cell volume */
    if (strcmp(coord[0], "R") == 0){
      rind = 0;
      if (pcptr->child[1] != NULL)
	vo = pcptr->child[1]->vol;
      else{
	for (id = 1; id <= nd; id++)
	  cxc[id-1] = pcptr->xc[id-1] + 0.5*lev[ilev+1]->dxbar[id-1];
	vo = cell_vol(ilev+1, cxc);
      }
    }
    else{
      rind = 1;
      if (pcptr->child[2] != NULL)
	vo = pcptr->child[2]->vol;
      else{
	for (id = 1; id <= nd; id++)
	  cxc[id-1] = pcptr->xc[id-1] + 0.5*lev[ilev+1]->dxbar[id-1];
	vo = cell_vol(ilev+1, cxc);
      }
    }
    tv = lev[ilev+1]->dxbar[rind]*(vo - vi)/vc;
    for (i = 0; i < nq; i++)
      qc[i] -= tv*dqdx[rind][i];
  }
  /* Map down solution */
  for (icell = 0; icell < nrcell; icell++){
    cptr = pcptr->child[icell];
    if (cptr != NULL){
      if (test){
	/* For test put in pa */
	for (i = 0; i < nq; i++)
	  cptr->pa[i] = qc[i];
	for (id = 1; id <= nd; id++){
	  for (i = 0; i < nq; i++)
	    cptr->pa[i] += (cptr->xc[id-1] - pxc[id-1])*dqdx[id-1][i];
	}
      }
      else{
	/* For cell creation put in qa */
	for (i = 0; i < nq; i++)
	  cptr->qa[i] = qc[i];
	for (id = 1; id <= nd; id++){
	  for (i = 0; i < nq; i++)
	    cptr->qa[i] += (cptr->xc[id-1] - pxc[id-1])*dqdx[id-1][i];
	}
	/* Fix solution (qa) */
	usrfxq(cptr);
      }
    }
  }
  return;
}
//end map_down_cell.c

//file cal_diffql.c
/*
+ cal_diffql - calculates left conservative gradient.
     Description :
        Calculates  the left  gradient in  conservatives  in direction
        id,for the cell  pointed to by cptr. Returns  true if the left
        cell  exists on  this  level. The  difference  is returned  in
        diff).
*/
bool cal_diffql(int id, struct cell *cptr, float diff[])
{
  int i;
  float pl[MAXQ], pr[MAXQ], ql[MAXQ], qr[MAXQ];
  bool got;
  struct cell *left;
  struct join *ljoin;
  ljoin = cptr->ljoin[id-1];
  left = ljoin->left;
  if (left == NULL){
    if (ljoin->bound == lebound){
      /* Left external boundary */
      got = true;
      for (i = 0; i < nq; i++)
	qr[i] = cptr->qa[i];
      usrpq(qr, pr);
      usrbd(id, lebound, ljoin, pl, pr);
      usrqp(pl, ql);
      for (i = 0; i < nq; i++)
	diff[i] = qr[i] - ql[i];
    }
    else if ((id == nd) && (cptr->halo) && (cptr->gotdiff)){
      /* Left halo cell with difference - use cell diff */
      got = true;
      for (i = 0; i < nq; i++)
	diff[i] = cptr->diff[i];
    }
    else{
      got = false;
      for (i = 0; i < nq; i++)
	diff[i] = 0.0;
    }
  }
  else{
    got = true;
    for (i = 0; i < nq; i++)
      diff[i] = cptr->qa[i] - left->qa[i];
  }
  return(got);
}
//end cal_diffql.c

//file cal_diffqr.c
/*
+ cal_diffqr - calculates right conservative gradient.
     Description :
        Calculates  the right gradient  in conservatives  in direction
        id, for the cell pointed to by cptr. Returns true if the right
        cell  exists on  this  level. The  difference  is returned  in
        diff).
*/
bool cal_diffqr(int id, struct cell *cptr, float diff[])
{
  int i;
  float pl[MAXQ], pr[MAXQ], ql[MAXQ], qr[MAXQ];
  bool got;
  struct cell *right;
  struct join *rjoin;
  rjoin = cptr->rjoin[id-1];
  right = rjoin->right;
  if (right == NULL){
    if (rjoin->bound == rebound){
      /* Right external boundary */
      got = true;
      for (i = 0; i < nq; i++)
	ql[i] = cptr->qa[i];
      usrpq(ql, pl);
      usrbd(id, rebound, rjoin, pl, pr);
      usrqp(pr, qr);
      for (i = 0; i < nq; i++)
	diff[i] = qr[i] - ql[i];
    }
    else if ((id == nd) && (cptr->halo) && (cptr->gotdiff)){
      /* Right halo cell with difference - use cell diff */
      got = true;
      for (i = 0; i < nq; i++)
	diff[i] = cptr->diff[i];
    }
    else{
      got = false;
      for (i = 0; i < nq; i++)
	diff[i] = 0.0;
    }
  }
  else{
    got = true;
    for (i = 0; i < nq; i++)
      diff[i] = right->qa[i] - cptr->qa[i];
  }
  return(got);
}
//end cal_diffqr.c

//file mk_bound_join.c
/*
+ mk_bound_join - make a boundary join.
     Description :
        This  converts  an  ordinary   join  to  a  boundary  join  by
        allocating the child, pa0, dpdt, dpdx and cor arrays.
*/   
void mk_bound_join(join *jptr)
{
  if (jptr->pa0 != NULL)
      /* Already made -- return */
    return;
  else{
    /* Make interpolation arrays */
    jptr->pa0 = new float[MAXQ];
    jptr->dpdt = new float[MAXQ];
    jptr->dpdx = new float[MAXQ*(nd - 1)];
    jptr->cor = new float[MAXQ];
  }
  return;
}
//end mk_bound_join.c

//file umk_bound_join.c
/*
+ umk_bound_join - unmake a boundary join.
     Description :
        This  converts  a boundary  join  to  a  non-boundary join  by
        deallocating the pa0, dpdt, dpdx and cor arrays.
*/   
void umk_bound_join(join *jptr)
{
  /* Clear join arrays if any */
  if (jptr->pa0 != NULL){
    delete [] jptr->pa0;
    jptr->pa0 = NULL;
  }
  if (jptr->dpdt != NULL){
    delete [] jptr->dpdt;
    jptr->dpdt = NULL;
  }
  if (jptr->dpdx != NULL){
    delete [] jptr->dpdx;
    jptr->dpdx = NULL;
  }
  if (jptr->cor != NULL){
    delete [] jptr->cor;
    jptr->cor = NULL;
  }
  return;
}
//end umk_bound_join.c

//file set_bound_flag.c
/*
+ set_bound_flag - set boundary flags.
     Description :
        This  checks level  ilev  and ilev-1.   It  sets external  and
        internal  boundary  flags  and makes/unmakes  ilev-1  boundary
        joins.
*/
void set_bound_flag(int ilev)
{
  int id;
  float xt;
  struct join *jptr;
  struct cell_list *clist;
  for (id = 1; id <= nd; id++){
    /* Set level ilev flags for internal joins */
    jptr = lev[ilev]->fjoin[id-1];
    while (jptr != NULL){
      set_fine_bound(ilev, id, jptr);
      jptr = jptr->nxt;
    }
    if (limpe > 1)
      /* Set ilev halo outer join boundary flag */
      com_halo(com9, ilev, 0.0);
  }
  if (ilev > 0){
    for (id = 1; id <= nd; id++){
      /* Check for cibounds on level ilev-1 internal joins */
      jptr = lev[ilev-1]->fjoin[id-1];
      while (jptr != NULL){
	set_coarse_bound(ilev-1, id, jptr);
	jptr = jptr->nxt;
      }
    }
    if (limpe > 1)
      /* Set ilev-1 halo outer join boundary flag */
      com_halo(com9, ilev - 1, 0.0);
  }
  if (limpe > 1){
    for (id = 1; id < nd; id++){
      /* Set level ilev flags for left halo joins */
      jptr = plev[ilev]->lfjoin[id-1];
      while (jptr != NULL){
	set_fine_bound(ilev, id, jptr);
	jptr = jptr->nxt;
      }
      /* Set level ilev flags for right halo joins */
      jptr = plev[ilev]->rfjoin[id-1];
      while (jptr != NULL){
	set_fine_bound(ilev, id, jptr);
	jptr = jptr->nxt;
      }
    }
    if (ilev > 1){
      for (id = 1; id < nd; id++){
	/* Check for cibounds on level ilev-1 left halo joins */
	jptr = plev[ilev-1]->lfjoin[id-1];
	while (jptr != NULL){
	  set_coarse_bound(ilev-1, id, jptr);
	  jptr = jptr->nxt;
	}
	/* Check for cibounds on level ilev-1 right halo joins */
	jptr = plev[ilev-1]->rfjoin[id-1];
	while (jptr != NULL){
	  set_coarse_bound(ilev-1, id, jptr);
	  jptr = jptr->nxt;
	}
      }
      /* Check halo nd joins */
      if (lev[ilev-1]->fcell == NULL){
	/* No real cells -- check cibounds for halo nd joins */
	clist = plev[ilev-1]->lhalo;
	xt = xr[nd-1] - 0.5*lev[ilev-1]->dxbar[nd-1];
	while (clist != NULL){
	  jptr = clist->cptr->rjoin[nd-1];
	  /* Check for external boundary cell */
	  if ((jptr->xj[nd-1] > xt) && (bcl[nd-1] != bcper))
	    jptr->bound = rebound;
	  set_coarse_bound(ilev-1, nd, jptr);
	  clist = clist->nxt;
	}
      }
    }
  }
  return;
}
//end set_bound_flag.c

//file set_fine_bound.c
/*
+ set_fine-bound - set fine grid boundary flags.
     Description :
        This sets external and internal coarse/fine boundary flags for
        the join jptr on level ilev.
*/
void set_fine_bound(int ilev, int id, struct join *jptr)
{
  float dx;
  /* 1/2 Mesh spacing */
  dx = 0.5*lev[ilev]->dxbar[id-1];
  /* Default is internal unless cibound */
  if (jptr->bound != cibound)
    jptr->bound = internal;
  if (jptr->left == NULL){
    if (jptr->xj[id-1] - dx < xl[id-1]){
      /* Check for periodic boundary */
      if (bcl[id-1] != bcper)
	/* Non periodic left external boundary */
	jptr->bound = lebound;
      else{
	/* Periodic fine/coarse on left boundary */
	jptr->bound = lfibound;
	jptr->par = jptr->right->par->ljoin[id-1];
      }
    }
    else{
      /* Left fine/coarse internal boundary */
      jptr->bound = lfibound;
      jptr->par = jptr->right->par->ljoin[id-1];
    }
  }
  else if (jptr->right == NULL){
    if (jptr->xj[id-1] + dx > xr[id-1]){
      /* Check for periodic boundary */
      if (bcr[id-1] != bcper)
	/* Non periodic right external boundary */
	jptr->bound = rebound;
      else{
	/* Periodic fine/coarse on right boundary */
	jptr->bound = rfibound;
	jptr->par = jptr->left->par->rjoin[id-1];
      }
    }
    else{
      /* Right fine/coarse internal boundary */
      jptr->bound = rfibound;
      jptr->par = jptr->left->par->rjoin[id-1];
    }
  }
  return;
}
//end set_fine_bound.c

//file set_coarse_bound_flag.c
/*
+ set_coarse-bound_flag - set coarse grid boundary flags.
     Description :
        This  checks for  coarse/fine boundaries  on level  ilev+1 and
        sets cibound etc as required.
*/
void set_coarse_bound(int ilev, int id, struct join *jptr)
{
  bool lref, rref;
  struct cell *left, *right;
  if (jptr->bound < lebound){
    left = jptr->left;
    right = jptr->right;
    lref = false;
    rref = false;
    if (left != NULL){
      if (left->refined)
	lref = true;
    }
    if (right != NULL){
      if (right->refined)
	rref = true;
    }
    if (lref == rref){
      /* Remove internal boundary */
      umk_bound_join(jptr);
      jptr->bound = internal;
    }
    else{
      /* Create/confirm  internal boundary */
      mk_bound_join(jptr);
      jptr->bound = cibound;
    }
  }
  return;
}
//end set_coarse_bound.c

//file check_grid.c
/*
+ check_grid - check grid for double refinement.
     Description :
        This checks  all grid  levels for double  refinements. Returns
        false if it finds some.
*/
bool check_grid()
{

  int ilev, id, tid, icell;
  float dx;
  struct cell *cptr, *left, *right;
  struct join *jptr, *pjptr;
  bool fail;

  fail = false;
  for (ilev = 2; ilev < nlevs - 1; ilev++){
    for (id = 1; id <= nd; id++){
      pjptr = lev[ilev]->fjoin[id-1];
      dx = 0.5*lev[ilev]->dxbar[id-1];
      while (pjptr != NULL){
	left = pjptr->left;
	right = pjptr->right;
	if (left == NULL){
	  if (pjptr->xj[id-1] - dx > xl[id-1]){
	    /* Not left external boundary */
	    if (right->refined){
	      for (icell = 0; icell < nrcell; icell++){
		cptr = right->child[icell];
		jptr = cptr->ljoin[id-1];
		if (jptr->left == NULL){
		  fail = true;
		  if (debug){
		    cout << "double left refinement id " << id << '\n';
		    cout << " coarse level " << ilev << " join "
			 << pjptr->ijoin << " at ";
		    for (tid = 1; tid <=nd; tid++)
		      cout << " " << pjptr->xj[tid-1];
		    cout << '\n';
		    cout << " fine level " << ilev+1 << " join "
			 << jptr->ijoin << " at ";
		    for (tid = 1; tid <=nd; tid++)
		      cout << " " << jptr->xj[tid-1];
		    cout << '\n';
		  } 
		}
	      }
	    }
	  }
	}
	if (right == NULL){
	  if (pjptr->xj[id-1] + dx < xr[id-1]){
	    /* Not right external boundary */
	    if (left->refined){
	      for (icell = 0; icell < nrcell; icell++){
		cptr = left->child[icell];
		jptr = cptr->rjoin[id-1];
		if (jptr->right == NULL){
		  fail = true;
		  if(debug){
		    cout << "double right refinement id " << id <<'\n';
		    cout << " coarse level " << ilev << " join "
			 << pjptr->ijoin << " at ";
		    for (tid = 1; tid <=nd; tid++)
		      cout << " " << pjptr->xj[tid-1];
		    cout << '\n'; 
		    cout << " fine level " << ilev+1 << " join "
			 << jptr->ijoin << " at ";
		    for (tid = 1; tid <=nd; tid++)
		      cout << " " << jptr->xj[tid-1];
		    cout << '\n';
		  }  
		}
	      }
	    }
	  }
	}
	pjptr = pjptr->nxt;
      }
    }
  }

  return(!fail);
}
//end file check_grid.c

//file diff_bad.c
/*
+ diff_bad - diffuse cell badness
     Description :
        This  diffuses the qsol scalar on level ilev for nsteps.
*/
void diff_bad(int ilev, int nstep)
{
  int id, istep;
  float fqsol;
  struct join *jptr;
  struct cell *left, *right, *cptr;
  /* Set qsol = 1.0 for all bad cells, zero dqsoldt  */
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    cptr->qsol = min(cptr->qsol, 1.0);
    cptr->dqsoldt = 0.0;
    cptr = cptr->nxt;
  }
  if (limpe > 1)
      /* Communicate ilev halo qsol (com1) */
    com_halo(com1, ilev, 0.0);
  /* Calculate fluxes and increment dqsoldt */
  for (istep = 1; istep <= nstep; istep++){
    for (id = 1; id <= nd; id++){
      jptr = lev[ilev]->fjoin[id-1];
      while (jptr != NULL){
        left = jptr->left;
        right = jptr->right;
	fqsol = 0.0;
	if (left != NULL)
	  fqsol += 0.2*left->qsol;
	else if (jptr->bound == lebound)
	  fqsol += 0.2*right->qsol;
        if (right != NULL)
	  fqsol -= 0.2*right->qsol;
	else if (jptr->bound == rebound)
	  fqsol -= 0.2*left->qsol;
	if (left != NULL)
	  left->dqsoldt -= fqsol;
        if (right != NULL)
	  right->dqsoldt += fqsol;
        jptr = jptr->nxt;
      }
      /* Increment qsol */
      cptr = lev[ilev]->fcell;
      while (cptr != NULL){
        /* Always refine the initially bad-marked cells */
        if (cptr->qsol < 1.0)
          cptr->qsol += cptr->dqsoldt;
	cptr->dqsoldt = 0.0;
	cptr = cptr->nxt;
      }
      if (limpe > 1)
	/* Communicate ilev halo qsol (com1)*/
	com_halo(com1, ilev, 0.0);
    }
  }
  return;
}
//end diff_bad.c

//file rem_hole.c
/*
+ rem_hole - remove holes.propagate cell badness.
     Description :
        This sets qsol of a  hole to badneigh, which forces refinement
        of this cell.
*/
void rem_hole(int ilev)
{
  int id;
  bool doit;
  struct cell *cptr, *left, *right; 
  /* Run through cells */
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    if (cptr->qsol == good){
      /* Good cell -- check neighbours */
      doit = true;
      for (id = 1; id <= nd; id++){
	/* Check left cells */
	left = cptr->ljoin[id-1]->left;
	if (left == NULL)
	    doit = false;
	else{
	  if (left->qsol == good)
	    doit = false;
	}
	/* Check right cells */
	right = cptr->rjoin[id-1]->right;
	if (right == NULL)
	    doit = false;
	else{
	  if (right->qsol == good)
	    doit = false;
	}
      }
    }
    if (doit)
      cptr->qsol = badneigh;
    cptr = cptr->nxt;
  }
  return;
}
//end rem_hole.c

//file de_frag.c
/*
+ de_frag - Defragement cell and join memory.
     Description :
        This reorders the cells in carr on each level so that they are
        in the order  of their notional positions. The  joins are then
        reordered according to the cells to which the are connected.
*/
void de_frag()
{
  int i, id, ilev, icell, tncell, lncell[MAXLEV], loff[MAXLEV], ncell,
    il, ir, ijoin, pass, fail;
  bool done;
  float ocpu, ncpu;
  const float clocks = CLOCKS_PER_SEC;
  struct cell *cptr;
  struct part *pptr;
  struct part_list *start, *end, *plist, *pmov;
  ocpu = clock()/clocks;
  if (particle){
    /* Make particle list on this processor */
    mk_partlist(start, end);
    if ((limpe > 1) && (repart)){
      /* Add particles moved from other processors */
      if (end != NULL)
	end->nxt = pstart;
      else
	start = pstart;
    }
    /* Clear particle pointers */
    plist = start;
    while (plist != NULL){
      pptr = plist->pptr;
      pptr->prev = NULL;
      pptr->nxt = NULL;
      pptr->cptr = NULL;
      plist = plist->nxt;
    }
  }
  if (limpe > 1){
    /* Clear halos */
    for (ilev = 0; ilev < nlevs; ilev++)
      clear_halo(ilev);
  }
  /* Count cells on each level and total */
  tncell = 0;
  for (ilev = 0; ilev < nlevs; ilev++){
    lncell[ilev] = 0;
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      lncell[ilev]++;
      cptr = cptr->nxt;
    }
    /* Total number of cells */
    tncell += lncell[ilev];
  }
  /* Update locarr array (temporary) */
  mk_locarr();
  /* Connect to parents */
  for (ilev = 1; ilev < nlevs; ilev++){
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      connect_par(ilev, cptr);
      cptr = cptr->nxt;
    }
  }
  /* Size of locarr array */
  ncell = 1;
  for (id = 1; id <= nd; id++)
    ncell *= nla[id-1];
  /* Set level offsets */
  loff[0] = 0;
  for (ilev = 1; ilev < nlevs; ilev++)
    loff[ilev] = loff[ilev-1] + lncell[ilev-1];
  /* Put desired cell number in icell for each level */
  for (icell = 1; icell <= ncell; icell++){
    if (locarr[icell] != NULL){
      loff[0]++;
      locarr[icell]->icell = loff[0];
      /* Recursive cell numbering (sets halo to true if found) */
      num_cell(1, loff, locarr[icell]);
    }
  }
  if (limpe > 1){
    /* Orphan cells go at end of ilev section of carr */
    for (ilev = 1; ilev < nlevs; ilev++){
      cptr = lev[ilev]->fcell;
      while (cptr != NULL){
	if (!cptr->halo){
	  loff[ilev]++;
	  cptr->icell = loff[ilev];
	}
	/* Set halo to false (safety) */
	cptr->halo = false;
	cptr = cptr->nxt;
      }
    }
  }
  /* Exchange cells */
  done = false;
  pass = 0;
  while (!done){
    for (i = 1; i <= celltop; i++){
      if (!carr[i]->empty){
	icell = carr[i]->icell;
	if (icell != i){
	  if (carr[icell]->empty){
	    /* Target cell is empty -- put source in target */
	    move_cell(i, icell);
	    clean_cell(carr[i]);
	  }
	  else
	    /* Target cell is not empty -- exchange */
	    exchange_cell(icell, i);
	}
      }
    }
    pass++;
    done = true;
    fail = 0;
    for (i = 1; i <= tncell; i++){
      if ((carr[i]->icell != i) || (carr[i]->empty)){
	done = false;
	fail++;
      }
    }
  }
  /* Clean cells */
  for (icell = 1; icell <= celltop; icell++){
    clean_cell(carr[icell]);
    carr[icell]->pptr = NULL;
    carr[icell]->icell = icell;
    if (icell <= tncell)
      carr[icell]->empty = false;
  }
  /* Reset celltop */
  celltop = tncell;
  /* Delete cell freelist */
  del_cell_freelist();
  /* Recreate cell lists */
  il = 1;
  for (ilev = 0; ilev < nlevs; ilev++){
    if (lncell[ilev] != 0){
      ir = il + lncell[ilev] - 1;
      lev[ilev]->fcell = carr[il];
      lev[ilev]->lcell = carr[ir];
      for (icell = il; icell < ir; icell++)
	carr[icell]->nxt = carr[icell+1];
      carr[ir]->nxt = NULL;
      for (icell = ir; icell > il; icell--)
	carr[icell]->prev = carr[icell-1];
      carr[il]->prev = NULL;
      il = ir + 1;
    }
    else{
      lev[ilev]->fcell = NULL;
      lev[ilev]->lcell = NULL;
    }
  }
  /* Clean all joins */
  for (ijoin = 1; ijoin <= jointop; ijoin++)
    clean_join(jarr[ijoin]);
  jointop = 0;
  /* Delete join freelist */
  del_join_freelist();
  /* Create level 0 halo */
  create_halo(0);
  /* Update locarr array */
  mk_locarr();
  /* Reconnect level 0 */
  connect(0);
  /* Connect level 0 left halo cells */
  connect_halo(0, plev[0]->lhalo);
  /* Connect level 0 right halo cells */
  connect_halo(0, plev[0]->rhalo);
  /* Fix lev[0]->xl/xr */
  if (plev[0]->lhalo != NULL)
    lev[0]->xl[nd-1] = plev[0]->lhalo->cptr->rjoin[nd-1]->xj[nd-1];
  if (plev[0]->rhalo != NULL)
    lev[0]->xr[nd-1] = plev[0]->rhalo->cptr->ljoin[nd-1]->xj[nd-1];
  mk_locarr();
  /* Reconnect other levels */
  for (ilev = 1; ilev < nlevs; ilev++){
    connect(ilev);
    create_halo(ilev);
  }
  if (limpe > 1){
    for (ilev = 0; ilev < nlevs; ilev++)
      /* Communicate halo refinement (com8) */
      com_halo(com8, ilev, 0.0);
  }
  /* Set join boundary flags */
  for (ilev = 0; ilev < nlevs; ilev++)
    set_bound_flag(ilev);
  if (particle){
    /* Assign particles */
    pmov = assign_part(start);
    if (limpe > 1){
      /* Move particles if necessary */
      int temp;
      if (pmov != NULL){
	temp = 1;
	cout << "fixing de_frag pmov not NULL on proc " << mype << '\n';
      }
      else
	temp = -1;
      if ((int)global_max((float)temp) > 0){
	/* Move cells to other processors if necessary */
	get_part_proc(pmov);
	start = mov_part(false, pmov);
	/* Assign moved particles */
	pmov = assign_part(start);
	if (pmov != NULL){
	  cout << "failed to fix assign in de_frag on proc " << mype <<
	    '\n';
	  plist = pmov;
	  while (plist != NULL){
	    cout << " Bad particle at "<< '\n';
	    for (id = 1; id <= nd; id++)
	      cout << " " << plist->pptr->xp[id-1];
	    cout << '\n';
	    plist = plist->nxt;
	  }
	}
      }
    }
  }
  ncpu = clock()/clocks;
  if ((debug) && (mype == 0)){
    cout.setf(ios::scientific);
    cout << "defrag: that took " << setprecision(4)
	 << ncpu - ocpu << " secs\n";
  }
  if (grav || mugi)
    /* Clean multigrid */
    clean_mug();  
  /* Reset alev */
  set_alev();
  return;
}
//end de_frag.c

//file num_cell.c
/*
+ num_cell - Renumber cells.
     Description :
        This renumbers cells on level=ilev using the order of the cell
        list on  level=ilev-1. It is  used to determine  the defragged
        cell order. loff[ilev] holds the current cell number.
*/
void num_cell(int ilev, int loff[], struct cell *pcptr)
{
  int icell;
  struct cell *cptr;
  /* Return if no child cell */
  if (!pcptr->refined)
    return;
  /* Run through child cells putting desired order in icell */
  for (icell = 0; icell < nrcell; icell++){
    cptr = pcptr->child[icell];
    if (cptr != NULL){
      loff[ilev]++;
      cptr->icell = loff[ilev];
      cptr->halo = true;
      /* Recursive call to deal with sub-levels */
      if (cptr->refined)
	num_cell(ilev+1, loff, cptr);
    }
  }
  return;
}
//file num_cell.c

//file move_cell.c
/*
+ move_cell - Move a cell.
     Description :
        This moves  the data  in cell carr[oi]  to carr[ni].   It does
        nothing about the pointers, so these have to be reset.
*/
void move_cell(int oi, int ni)
{
  int i;
  struct cell *ocptr, *ncptr;
  /* Old and new cell pointers */
  ocptr =  carr[oi];
  ncptr = carr[ni];
  /* Set cell numbers etc */
  ocptr->icell = oi;
  ncptr->icell = ni;
  ncptr->ichild = ocptr->ichild;
  /* Move data */
  ncptr->refined = ocptr->refined;
  ncptr->empty = ocptr->empty;
  ncptr->halo = ocptr->halo;
  ncptr->vol = ocptr->vol;
  for (i = 0; i < nq; i++)
    ncptr->qa[i] = ocptr->qa[i];
  for (i = 0; i < nd; i++)
    ncptr->xc[i] = ocptr->xc[i];
  return;
}
//end move_cell.c

//file exchange_cell.c
/*
+ exchange_cell - Exchange cells.
     Description :
        This  exchanges the data  in cells  carr[oi] and  carr[ni]. It
        does  nothing about  the cell  pointers, so  these have  to be
        reset.
*/
void exchange_cell(int oi, int ni)
{
  int i, icell;
  bool state;
  float temp;
  struct cell *ocptr, *ncptr;
  /* Old and new cell pointers */
  ocptr =  carr[oi];
  ncptr = carr[ni];
  /* Exchange data */
  icell = ocptr->icell;
  ocptr->icell = ncptr->icell;
  ncptr->icell = icell;
  state = ocptr->refined;
  ocptr->refined = ncptr->refined;
  ncptr->refined = state;
  state = ocptr->empty;
  ocptr->empty = ncptr->empty;
  ncptr->empty = state;
  state = ocptr->halo;
  ocptr->halo = ncptr->halo;
  ncptr->halo = state;
  temp = ocptr->vol;
  ocptr->vol = ncptr->vol;
  ncptr->vol = temp;
  for (i = 0; i < nq; i++){
    temp = ocptr->qa[i];
    ocptr->qa[i] = ncptr->qa[i];
    ncptr->qa[i] = temp;
  }
  for (i = 0; i < nd; i++){
    temp = ocptr->xc[i];
    ocptr->xc[i] = ncptr->xc[i];
    ncptr->xc[i] = temp;
  }
  return;
}
//end exchange_cell.c

//file del_cell_freelist.c
/*
+ del_cell_freelist - Delete cell freelist.
     Description :
     Deletes the  cell freelist. Used after  defragmentation and shift
     of carr array.

*/
void del_cell_freelist()
{
  struct cell_list *clist, *tclist;
  /* Delete cell freelist */
  clist = freecell;
  while (clist != NULL){
    tclist = clist->nxt;
    delete clist;
    clist = tclist;
  }
  freecell = NULL;
  return;
}
//end del_cell_freelist.c

//file del_join_freelist.c
/*
+ del_join_freelist - Delete join freelist.
     Description :
     Deletes the  join freelist. Used after  defragmentation and shift
     of jarr array.

*/
void del_join_freelist()
{
  struct join_list *jlist, *tjlist;
  /* Delete join freelist */
  jlist = freejoin;
  while (jlist != NULL){
    tjlist = jlist->nxt;
    delete jlist;
    jlist = tjlist;
  }
  freejoin = NULL;
  return;
}
//end del_join_freelist.c

//file cell_vol.c
/*
+ cell_vol - set volume for cell.
     Description :
        This sets  the cell volume.   These include the factors  of pi
        etc for cylindrical and spherical polars.
*/
float cell_vol(int ilev, float xc[])
{
  int id;
  float vol, r1, r2;
  vol = 1.0;
  if (strcmp(geom, "XYZ") == 0){
    /* Cartesian */
    vol = 1.0;
    for (id = 1; id <= nd; id++)
      vol *= lev[ilev]->dxbar[id-1];
  }
  else if (strcmp(geom, "RZP") == 0){
    /* Cylindrical polars (nd <= 2) */
    if (strcmp(coord[0], "R") == 0){
      r2 = xc[0] + 0.5*lev[ilev]->dxbar[0];
      r1 = xc[0] - 0.5*lev[ilev]->dxbar[0];
      vol = gkpi*(r2 - r1)*(r2 + r1);
      if (nd == 2)
	vol *= lev[ilev]->dxbar[1];
    }
    else{
      r2 = xc[1] + 0.5*lev[ilev]->dxbar[1];
      r1 = xc[1] - 0.5*lev[ilev]->dxbar[1];
      vol = gkpi*(r2 - r1)*(r2 + r1);
      if (nd == 2)
	vol *= lev[ilev]->dxbar[0];
    }
  }
  else if (strcmp(geom, "RTP") == 0){
    /* Spherical polars (n = 1 only) */
    r2 = xc[0] + 0.5*lev[ilev]->dxbar[0];
    r1 = xc[0] - 0.5*lev[ilev]->dxbar[0];
    vol = 4.0*gkpi*(r2 - r1)*(r2*r2 + r1*r2 + r1*r1)/3.0;
  }
  return(vol);
}
//end cell_vol.c

//file join_area.c
/*
+ join_area - set area for join.
     Description :
        This sets the join area.   These include the factors of pi etc
        for cylindrical and spherical polars.
*/
float join_area(int ilev, int id, float xj[])
{
  int idd;
  float area, r1, r2;
  area = 1.0;
  if (strcmp(geom, "XYZ") == 0){
    /* Cartesian */
    area = 1.0;
    for (idd = 1; idd <= nd; idd++)
      if (idd != id)
	area *= lev[ilev]->dxbar[idd-1];
  }
  else if (strcmp(geom, "RZP") == 0){
    /* Cylindrical polars (nd <= 2) */
    if (id == 1){
      if (strcmp(coord[0], "R") == 0){
	r1 = xj[0];
	area = 2.0*gkpi*r1;
	if (nd == 2)
	  area *= lev[ilev]->dxbar[1];
      }
      else{
	r1 = xj[1] - 0.5*lev[ilev]->dxbar[1];
	r2 = xj[1] + 0.5*lev[ilev]->dxbar[1];
	area = gkpi*(r2 - r1)*(r2 + r1);
      }
    }
    else if (id == 2){
      if (strcmp(coord[0], "R") == 0){
	r1 = xj[0] - 0.5*lev[ilev]->dxbar[0];
	r2 = xj[0] + 0.5*lev[ilev]->dxbar[0];
	area = gkpi*(r2 - r1)*(r2 + r1);
      }
      else{
	r1 = xj[1];
	area = 2.0*gkpi*r1;
	if (nd == 2)
	  area *= lev[ilev]->dxbar[0];
      }
    }
  }
  else if (strcmp(geom, "RTP") == 0){
    /* Spherical polars (n = 1 only) */
    r1 = xj[0];
    area = 4.0*gkpi*r1*r1;
  }
  return(area);
}
//end join_area.c

//file find_q.c
/*
+ find_q - Get solution at given point.
     Description :
        This  finds  the conserved  solution  at  the point  given  by
        xp[]. It uses the finest  grid at that point. Implemented with
        find_fcell, cal_diffql, cal_diffrq.
*/
bool find_q(float xp[], float qa[])
{
  int i, id, ilev;
  float divdx, dxp, qc[MAXQ], diffl[MAXQ], diffr[MAXQ], dqdx[MAXQ];
  bool gotl, gotr;
  struct cell *cptr;
  /* Find finest cell */
  cptr = find_fcell(ilev, xp);
  if (cptr == NULL)
    return(false);
  else if (cptr->halo)
    return(false);
  /* Set qc to central solution */
  for (i = 0; i < nq; i++)
    qc[i] = cptr->qa[i];
  /* Get gradients */
  for (id = 1; id <= nd; id++){
    /* Inverse mesh spacing */
    divdx = 1.0/lev[ilev]->dxbar[id-1];
    /* Left difference */
    gotl = cal_diffql(id, cptr, diffl);
    /* Right difference */
    gotr = cal_diffqr(id, cptr, diffr);
    if (gotl && gotr){
      for (i = 0; i < nq; i++)
	dqdx[i] = av(diffl[i], diffr[i])*divdx;
    }
    else if (gotl){
      for (i = 0; i < nq; i++)
	dqdx[i] = diffl[i]*divdx;
    }
    else if (gotr){
      for (i = 0; i < nq; i++)
	dqdx[i] = diffr[i]*divdx;
    }
    else{
      for (i = 0; i < nq; i++)
	dqdx[i] = 0.0;
    }
    /* Calculate coordinate difference */
    dxp = xp[id-1] - cptr->xc[id-1];
    /* Increment solution at xp */
    for (i = 0; i < nq; i++)
      qc[i] += dqdx[i]*dxp; 
  }
  /* Put solution in qa */
  for (i = 0; i < nq; i++)
    qa[i] = qc[i];
  return(true);
}
//end find_q.c

//file mk_partlist.c
/*
+ mk_partlist - Make a list of particles on a processor.
     Description :
        This goes through the cells on a processor and makes a list of
        any particles. start and end are pointers to the start and end
        of list. It also sets  pptr->ilev. Any routine that calls this
        is responsible for deleting the list.
*/
void mk_partlist(struct part_list *&start, struct part_list *&end)
{
  int ilev;
  struct cell *cptr;
  struct part *pptr;
  struct part_list *plist;
  /* Initialise particle list */
  start = NULL;
  /* Find particles */
  for (ilev = 0; ilev < nlevs; ilev++){
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      if ((!cptr->refined) && (cptr->pptr != NULL)){
	/* Add particles */
	pptr = cptr->pptr;
	while (pptr != NULL){
	  plist = new part_list;
	  if (start == NULL)
	    end = plist;
	  plist->pptr = pptr;
	  plist->nxt = start;
	  start = plist;
	  find_fcell(pptr->ilev, pptr->xp);
	  pptr = pptr->nxt;
	}
      }
      cptr = cptr->nxt;
    }
  }
  if (start == NULL)
    end = NULL;
  return;
}
//end mk_partlist.c

//file assign_part.c
/*
+ assign_part - Assign list of particles to cells.
     Description :
        This assigns the particles in  the list pointed to by start to
        their appropiate cells. This list is deleted, but any particle
        not  on  this  processor  is  added to  the  move  list  whose
        beginning is returned by the function.
*/
struct part_list *assign_part(struct part_list *start)
{
  int ilev;
  bool move;
  struct cell *cptr;
  struct part *pptr;
  struct part_list *plist, *nplist, *pmov;
  /* Initialise pmov */
  pmov = NULL;
  /* Assign particles to cells and delete list */
  plist = start;
  while (plist != NULL){
    pptr = plist->pptr;
    /* Next element */
    nplist = plist->nxt;
    /* Initial NULL pointers (safety) */
    pptr->prev = NULL;
    pptr->nxt = NULL;
    cptr = pptr->cptr;
    /* Get cell */
    cptr = find_fcell(ilev, pptr->xp);
    move = ((cptr == NULL) || (cptr->halo));
    if (move){
      /* Add to move list */
      plist->pproc = -1;
      plist->nxt = pmov;
      pmov = plist;
    }
    else{
      /* Assign to cell */
      pptr->cptr = cptr;
      if (cptr->pptr == NULL)
	pptr->nxt = NULL;
      else{
	cptr->pptr->prev = pptr;
	pptr->nxt = cptr->pptr;
      }
      cptr->pptr = pptr;
      /* Delete list element */
      delete plist;
    }
    plist = nplist;
  }
  return(pmov);
}
//end assign_part.c

//file del_allpart.c
/*
+ del_allpart - Delete all particles in a cell.
     Description :
        This deletes all particles in  the cell pointed to by cptr and
        sets cptr->pptr=NULL;
*/
void del_allpart(struct cell *cptr)
{
  struct part *pptr, *tpptr;
  pptr = cptr->pptr;
  while (pptr != NULL){
    tpptr = pptr;
    pptr = pptr->nxt;
    delete tpptr;
  }
  cptr->pptr = NULL;
  return;
}
//end del_allpart.c

//file rm_part.c
/*
+ rm_part - Remove a particle from a cell.
     Description :
        This removes the particle pointed to by pptr from its cell. It
        sets the particle cell pointer to NULL for safety. It does not
        delete the particle.
*/
void rm_part(struct part *pptr)
{
  int id;
  struct cell *cptr;
  cptr = pptr->cptr;
  if (cptr == NULL){
    cout << "Error: no cell pointer, particle pos";
    for (id = 0; id < nd; id++)
      cout << " " << pptr->xp[id];
    cout << '\n';
    return;
  }
  if (pptr->prev == NULL){
    if (pptr->nxt == NULL)
      cptr->pptr = NULL;
    else{
      cptr->pptr = pptr->nxt;
      pptr->nxt->prev = NULL;
    }
  }
  else if (pptr->nxt == NULL)
    pptr->prev->nxt = NULL;
  else{
    pptr->prev->nxt = pptr->nxt;
    pptr->nxt->prev = pptr->prev;
  }
  pptr->cptr = NULL;
  pptr->nxt = NULL;
  pptr->prev = NULL;
  return;
}
//end rm_part.c
