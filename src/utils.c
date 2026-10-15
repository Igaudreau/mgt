#undef MAIN
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "funcs.h"
extern void clear_mug();
extern void usrquit();
extern void umk_bound_join(struct join *jptr);
extern void del_cell_freelist();
extern void del_join_freelist();
extern void del_allpart(struct cell *cptr);
extern void clear_halo(int ilev);
extern void merge_array(float buf[], int nvalues[]);
extern float global_sum(float value);
extern void usrpq(float qa[], float pa[]);	// BM
void mk_level();
void mk_plevel();
void mk_memory();
void check(int ilev);
void check_cell(int ilev, int icell);
void check_join(int ilev, int ijoin);
void peek_cell(int ilev, float xc[]);
struct cell *find_cell(int ilev, float xc[]);
struct cell *find_fcell(int &ilev, float xc[]);
void mk_locarr();
void clear();
void clean_cell(struct cell *cptr);
void clean_join(struct join *jptr);
void mem_calc(bool all);
void i_to_it(int ndim, int nc[], int i, int it[]);
void i_to_it2 (int i, int it[]);
int it_to_i(int ndim, int nt[], int it[]);
int it2_to_i(int it[]);
int dpdx_index(int jid, int id);

//file mk_level.c
/*
+ mk_level - make level structure.
     Description 
        This  creates the level structure for all levels.
*/
void mk_level()
{
  int imult, ilev, id;
  float dxbar[MAXD];
  /* No of child cells */
  nrcell = ipow(nd);
  /* Mesh spacing */
  for (id = 1; id <= nd; id++)
    dxbar[id-1] = (xr[id-1] - xl[id-1])/nc[id-1];
  /* Make level structure */
  lev = new level *[MAXLEV];
  imult = 1;
  for (ilev = 0; ilev < nlevs; ilev++){
    lev[ilev] = new level;
    lev[ilev]->fjoin = new join *[nd];
    lev[ilev]->ljoin = new join *[nd];
    lev[ilev]->fcell = NULL;
    lev[ilev]->lcell = NULL;
    for (id = 1; id <= nd; id++){
      lev[ilev]->nc[id-1] = imult*nc[id-1];
      lev[ilev]->dxbar[id-1] = dxbar[id-1]/imult;
      lev[ilev]->xl[id-1] = xl[id-1];
      lev[ilev]->xr[id-1] = xr[id-1];
      lev[ilev]->fjoin[id-1] = NULL;
      lev[ilev]->ljoin[id-1] = NULL;
    }
    imult *= 2;
  }
  return;
}
//end file mk_level.c

//file mk_plevel.c
/*
+ mk_plevel - make plevel structure.
     Description 
        This  creates the plev structure for all levels.
*/
void mk_plevel()
{
  int ilev, id;
  plev = new plevel *[MAXLEV];
  for (ilev = 0; ilev < nlevs; ilev++){
    plev[ilev] = new plevel;
    /* Set pointers to NULL */
    plev[ilev]->lbound = NULL;
    plev[ilev]->rbound = NULL;
    plev[ilev]->lhalo = NULL;
    plev[ilev]->rhalo = NULL;
    plev[ilev]->lproc = NULL;
    plev[ilev]->rproc = NULL;
    plev[ilev]->lfjoin = new join *[nd-1];
    plev[ilev]->rfjoin = new join *[nd-1];
    for (id = 1; id < nd; id++){
      plev[ilev]->lfjoin[id-1] = NULL;
      plev[ilev]->rfjoin[id-1] = NULL;
    }
  }
}
//end mk_plevel.c

//file mk_memory.c
/*
+ mk_memory - make memory.
     Description :
        This  uses new  to  allocate  cells and  joins  and put  their
        pointers in the carr and jarr arrays.
*/
void mk_memory()
{
  int icell, ijoin;
  struct join *jptr;
  /* Arrays of cell and join pointers (Offset 1) */
  carr = (cell **)malloc((size_t)((carrtop+SAFETY)*sizeof(cell *)));
  jarr = (join **)malloc((size_t)((jarrtop+SAFETY)*sizeof(join *)));
  /* Fill carr with clean cells */
  for (icell = 1; icell <= carrtop; icell++){
    carr[icell] = new cell;
    clean_cell(carr[icell]);
    carr[icell]->icell = icell;
  }
  /* Fill jarr with joins */
  for (ijoin = 1; ijoin <= jarrtop; ijoin++){
    jptr = new join;
    /* Ensure that join is clean */
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
    jarr[ijoin] = jptr;
  }
  /* Zero celltop, jointop (safety) */
  celltop = 0;
  jointop = 0;
  /* Set minicell array to NULL, zero mcarrtop (safety) */
  mcarr = NULL;
  mcarrtop = 0;
  return;
}
//end file mk_memory.c

//file peek_cell.c
/*
+ peek_cell - Look at contents of cell.
     Description :
        This  gives the  data  for cell  on  level ilev  at coords  in
        xc[MAXD].
*/
void peek_cell(int ilev, float xc[])
{
  int i, id, ipart, left, right;
  struct cell *cptr, *lptr, *rptr;
  struct part *pptr;
  /* Find cell */
  cptr = find_cell(ilev, xc);
  /* Print info */
  if (cptr == NULL){
    printf("No such cell \n");
    return;
  }
  printf("proc %d \n", mype);
  printf("cell %d, qsol %e, vol %e", cptr->icell, cptr->qsol,
	 cptr->vol);
  if (cptr->refined)
    printf(" refined ");
  else
    printf(" unrefined");
  if (cptr->par != NULL)
    printf(" parent %d \n", cptr->par->icell);
  else
    printf(" level 0 \n");
  for (id = 1; id <= nd; id++){
    lptr = cptr->ljoin[id-1]->left;
    rptr = cptr->rjoin[id-1]->right;
    if (lptr == NULL)
      left = -1;
    else
      left = lptr->icell;
    if (rptr == NULL)
      right = -1;
    else
      right = rptr->icell;
    printf("dir %d, x %e, l cell %d, r cell %d", id, 
	   cptr->xc[id-1], left, right);
    printf(" l join %d area %e, r join %d area %e\n",
	   cptr->ljoin[id-1]->ijoin, cptr->ljoin[id-1]->area,
	   cptr->rjoin[id-1]->ijoin, cptr->rjoin[id-1]->area) ;
  }
  if (cptr->halo)
    printf(" halo \n");
  printf("cell vol %e \n",cptr->vol);

  float pa[MAXQ];
  usrpq(cptr->qa, pa);
  printf("i %d, qa[i] %e, pa[i] %e \n", 0, cptr->qa[0], pa[0]);	// BM
  for (i = 1; i < nq; i++){
    printf("i %d, qa[i] %e, pa[i] %e \n", i, cptr->qa[i], pa[i]);  // BM
  }
  if (particle){
    /* Print particle info */
    pptr = cptr->pptr;
    ipart = 1;
    while (pptr != NULL){
      printf("particle %d, mass %e energy %e \n", ipart, pptr->mass,
        pptr->en);
      for (i = 0; i < nd; i++)
	printf("dir %d, xp %e vp %e \n", i, pptr->xp[i], pptr->vp[i]);
      pptr = pptr->nxt;
    }
  }
  return;
}
//end file peek_cell.c

//file find_cell.c
/*
+ find_cell - Find a cell.
     Description :
        Returns a pointer to the cell on level ilev which contains the
        point with coordinates xc. Done  by searching from the level 0
        cell, which  is found using  the locarr array. The  halos must
        have  been created  for levels  <= ilev.  Also works  for halo
        cells.  Returns NULL if the cell is not found.
*/
struct cell *find_cell(int ilev, float xc[])
{
  int id, tnd, tilev, it[MAXD];
  bool outside;
  float poff[MAXD], dx[MAXD], lxl[MAXD], lxr[MAXD];
  struct cell *cptr;
  /* Completely initialize arrays BM */
  for (id = 0; id < MAXD; id++) {
    it[id] = 0;
    poff[id] = 0.0;
    dx[id] = 0.0;
    lxl[id] = 0.0;
    lxr[id] = 0.0;
  }
  /* Check for legal level */
  if ((ilev < 0) || (ilev > nlevs - 1))
    return(NULL);
  /* Fix for periodic nd boundary */
  if (limpe > 1)
    tnd = nd;
  else
    tnd = nd + 1; 
  /* Mesh spacing, boundaries, periodic offset */
  for (id = 1; id <= nd; id++){
    dx[id-1] = lev[0]->dxbar[id-1];
    lxl[id-1] = lev[0]->xl[id-1];
    lxr[id-1] = lev[0]->xr[id-1];
    poff[id-1] = 0.0;
    /* Check for periodic */
    if (xc[id-1] < lxl[id-1]){
      if ((bcl[id-1] == bcper) && (id != tnd))
	/* Periodic - set poff */
	poff[id-1] = lxr[id-1] - lxl[id-1];
    }
    else if (xc[id-1] > lxr[id-1]){
      if ((bcl[id-1] == bcper) && (id != tnd))
	/* Periodic - set poff */
	poff[id-1] = lxl[id-1] - lxr[id-1];
    }
  }
  if (limpe > 1)
    /* Allow for halo cells */
    lxl[nd-1] -= dx[nd-1];
  /* Find coarse cell */
  outside = false;
  for (id = 1; id <= nd; id++){
    it[id-1] = (int)((xc[id-1] + poff[id-1] - lxl[id-1])/dx[id-1]
      + 1.0);
    if ((it[id-1] < 1) || (it[id-1] > nla[id-1])){
      outside = true;
      break;
    }
  }
  if (outside)
    return(NULL);
  cptr = locarr[it_to_i(nd, nla, it)];
  /* Search this cell down to level ilev */
  for (tilev = 1; tilev <= ilev; tilev++){
    if (cptr == NULL)
      break;
    if (!cptr->refined){
      cptr = NULL;
      break;
    }
    else{
      for (id = 1; id <= nd; id++){
	if (xc[id-1] + poff[id-1] < cptr->xc[id-1])
	  it[id-1] = 1;
	else
	  it[id-1] = 2;
      }
      cptr = cptr->child[it2_to_i(it)];
    }
  }
  return(cptr);
}
//end find_cell.c

//file find_fcell.c
/*
+ find_fcell - Find a fine grid cell.
     Description :
        Returns a  pointer to the  finest grid cell that  contains the
        point with coordinates xc. Done  by searching from the level 0
        cell, which  is found using  the locarr array. Also  works for
        halo cells. Returns NULL if the cell is not found. Returns the
        level  at which  the cell  is found  in ilev  (ilev=-1  if not
        found). Used for particles.
*/
struct cell *find_fcell(int &ilev, float xc[])
{
  int id, tnd, tilev, it[MAXD];
  bool outside;
  float poff[MAXD], dx[MAXD], lxl[MAXD], lxr[MAXD];
  struct cell *cptr;
  /* Completely initialize arrays BM */
  for (id = 0; id < MAXD; id++) {
    it[id] = 0;
    poff[id] = 0.0;
    dx[id] = 0.0;
    lxl[id] = 0.0;
    lxr[id] = 0.0;
  }
  /* Fix for periodic nd boundary */
  if (limpe > 1)
    tnd = nd;
  else
    tnd = nd + 1; 
  /* Mesh spacing, boundaries, periodic offset */
  for (id = 1; id <= nd; id++){
    dx[id-1] = lev[0]->dxbar[id-1];
    lxl[id-1] = lev[0]->xl[id-1];
    lxr[id-1] = lev[0]->xr[id-1];
    poff[id-1] = 0.0;
    /* Check for periodic */
    if (xc[id-1] < lxl[id-1]){
      if ((bcl[id-1] == bcper) && (id != tnd))
	/* Periodic - set poff */
	poff[id-1] = lxr[id-1] - lxl[id-1];
    }
    else if (xc[id-1] > lxr[id-1]){
      if ((bcl[id-1] == bcper) && (id != tnd))
	/* Periodic - set poff */
	poff[id-1] = lxl[id-1] - lxr[id-1];
    }
  }
  if (limpe > 1)
    /* Allow for halo cells */
    lxl[nd-1] -= dx[nd-1];
  /* Find coarse cell */
  outside = false;
  for (id = 1; id <= nd; id++){
    it[id-1] = (int)((xc[id-1] + poff[id-1] - lxl[id-1])/dx[id-1]
      + 1.0);
    if ((it[id-1] < 1) || (it[id-1] > nla[id-1])){
      outside = true;
      break;
    }
  }
  if (outside){
    ilev = -1;
    return(NULL);
  }
  cptr = locarr[it_to_i(nd, nla, it)];
  if (nlevs == 1){
    ilev = 0;
    return(cptr);
  }
  /* Search this cell down to level nlevs - 1 */
  for (tilev = 1; tilev < nlevs; tilev++){
    if (cptr == NULL)
      break;
    if (!cptr->refined)
      break;
    else{
      for (id = 1; id <= nd; id++){
	if (xc[id-1] + poff[id-1] < cptr->xc[id-1])
	  it[id-1] = 1;
	else
	  it[id-1] = 2;
      }
      cptr = cptr->child[it2_to_i(it)];
      ilev = tilev;
    }
  }
  if (cptr == NULL)
    ilev = -1;
  return(cptr);
}
//end find_fcell.c

//file mk_locarr.c
/*
+ mk_locarr - Make and fill location array.
     Description :
        This makes and fills the location array for level 0. The level
        0  halo  must  have  been  created. It  must  be  created  for
        find_cell and find_fcell to work.
*/
void mk_locarr()
{
  int icell, ncell, id, it[MAXD];
  float dx[MAXD], lxl[MAXD];
  struct cell *cptr; 
  struct cell_list *clist; 
  /* Completely initialize arrays BM */
  for (id = 0; id < MAXD; id++) {
    it[id] = 0;
    dx[id] = 0.0;
    lxl[id] = 0.0;
  }
  /* Boundaries of location array */
  for (id = 1; id <= nd; id++){
    dx[id-1] = lev[0]->dxbar[id-1];
    lxl[id-1] = lev[0]->xl[id-1];
    nla[id-1] = lev[0]->nc[id-1];
  }
  if (limpe > 1){
    /* Allow for halo cells */
    lxl[nd-1] -= dx[nd-1];
    nla[nd-1] += 2;
  }
  /* Size of array */
  ncell = 1;
  for (id = 1; id <= nd; id++)
    ncell *= nla[id-1];
  if (locarr != NULL){
    /* Delete old location array */
    // pragmas added by BM to instruct compiler to ignore warning
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wfree-nonheap-object"
    locarr++;
    delete [] locarr;
    locarr = NULL;		// BM
  }
  /* Create location array for coarsest level */
  locarr = new cell *[ncell+SAFETY];
  locarr--;
  /*Initialize location array */
  for (icell = 1; icell <= ncell; icell++)
    locarr[icell] = NULL;
  /* Put internal cells in location array */
  cptr = lev[0]->fcell;
  while (cptr != NULL){
    /* Position indices */
    for (id = 1; id <= nd; id++)
      it[id-1] = (int)((cptr->xc[id-1] - lxl[id-1])/dx[id-1]) + 1;
    locarr[it_to_i(nd, nla, it)] = cptr;
    cptr = cptr->nxt;
  }
  /* Put left halo cells in location array */
  clist = plev[0]->lhalo;
  while (clist != NULL){
    /* Position indices */
    cptr = clist->cptr;
    for (id = 1; id <= nd; id++)
      it[id-1] = (int)((cptr->xc[id-1] - lxl[id-1])/dx[id-1]) + 1;
    locarr[it_to_i(nd, nla, it)] = clist->cptr;
    clist = clist->nxt;
  }
  /* Put right halo cells in location array */
  clist = plev[0]->rhalo;
  while (clist != NULL){
    /* Position indices */
    cptr = clist->cptr;
    for (id = 1; id <= nd; id++)
      it[id-1] = (int)((cptr->xc[id-1] - lxl[id-1])/dx[id-1]) + 1;
    locarr[it_to_i(nd, nla, it)] = cptr;
    clist = clist->nxt;
  }
  return;
}
//end mk_locarr.c

//file clear.c
/*
+ clear - cleanup model memory
     Description :
        This deletes cells and joins etc from any previous model.
*/   
void clear()
{
  int icell, ijoin, ilev;
  /* Clear multigrid */
  clear_mug();  
  /* Clear parallel stuff */
  for (ilev = 0; ilev < nlevs; ilev++){
    /* Clear halos */
    clear_halo(ilev);
    /* Delete plev */
    delete [] plev[ilev]->lfjoin;
    delete [] plev[ilev]->rfjoin;
    delete plev[ilev];
  }
  delete [] plev;
  /* Delete ins/outs arrays */
  delete [] ins;
  delete [] outs;
  /* Delete work arrays */
  delete [] wka;
  delete [] wkda;
  if (carr != NULL){
    /* Delete cells */
    for (icell = 1; icell <= carrtop; icell++){
      del_allpart(carr[icell]);
      delete carr[icell];
    }
    free(carr);
    carr = NULL;
  }
  del_cell_freelist();
  if (jarr != NULL){
    /* Delete joins */
    for (ijoin = 1; ijoin <= jarrtop; ijoin++){
      clean_join(jarr[ijoin]);
      delete jarr[ijoin];
    }
    free(jarr);
    jarr = NULL;
  }
  del_join_freelist();
  /* Delete lev */
  for (ilev = 0; ilev < nlevs; ilev++){
    delete [] lev[ilev]->fjoin;
    delete [] lev[ilev]->ljoin;
    delete lev[ilev];
  }
  delete [] lev;
  if (locarr != NULL){
    /* Delete locarr array */
    locarr++;
    delete [] locarr;
    locarr = NULL;
  }
  /* Zero celltop, jointop, carrtop, jarrtop (safety) */
  celltop = 0;
  jointop = 0;
  carrtop = 0;
  jarrtop = 0;
  /* User clean up */
  usrquit();
  return;
}
//end clear.c

//file clean_cell.c
/*
+ clean_cell - clean cell
     Description :
        This cleans a cell so that it is in the same state as after it
        has just been  allocated. This is used before  a cell is added
        to the free list and after defragmentation.
*/
void clean_cell(struct cell *cptr)
{
  int id, icell;
  /* Clear pointers */
  cptr->prev = NULL;
  cptr->nxt = NULL;
  for (icell = 0; icell < nrcell; icell++)
    cptr->child[icell] = NULL;
  cptr->par = NULL;
  for (id = 1; id <= nd; id++){
    cptr->ljoin[id-1] = NULL;
    cptr->rjoin[id-1] = NULL;
  }
  cptr->refined = false;
  cptr->halo = false;
  cptr->empty = true;
  cptr->pptr = NULL;
  cptr->mcptr = NULL;
  return;
}
//end clean_cell.c

//file clean_join.c
/*
+ clean_join - clean join
     Description :
        This cleans a join so that it is in the same state as after it
        has just been  allocated. This is used before  a join is added
        to the free list and after defragmentation.
*/
void clean_join(struct join *jptr)
{
  /* Set bound to internal */
  jptr->bound = internal;
  /* Clear pointers */
  jptr->prev = NULL;
  jptr->nxt = NULL;
  jptr->par = NULL;
  jptr->left = NULL;
  jptr->right = NULL;
  /* Clear boundary arrays (if any) */
  umk_bound_join(jptr);
  jptr->empty = true;
  return;
}
//end clean_join.c

//file mem_calc.c
/*
+ mem_calc - Calculates memory use.
     Description :
        This calculates the  memory used by joins and  cells. It takes
        into  account  coarse/fine  boundary joins  and  particles  if
        any. If  all=true, then  it prints  stats for  all processors,
        otherwise it the grand total.
*/
void mem_calc(bool all)
{
  int sizej, sizebj, sizec, sizep, id, ilev, npe, nvalues[MAXPE];
  long long int njoin, ncell, pcell, ncellt, pcellt, nbjoin, njoint,
    nfcell, pnjoin[MAXLEV][MAXPE], pncell[MAXLEV][MAXPE], 
    pnfcell[MAXLEV][MAXPE], ppcell[MAXLEV][MAXPE];
  float memj, memc, memtj, memtc, memtp, memp, temp[8*MAXPE],
    pmemj[MAXLEV][MAXPE], pmemc[MAXLEV][MAXPE], pmemp[MAXLEV][MAXPE],
    carrmem[MAXPE], jarrmem[MAXPE], tcarrmem, tjarrmem, pworkwt[MAXPE];
  float twork;		// Added by BM
  struct cell *cptr;
  struct join *jptr;
  struct part *pptr;
  //struct cell_list *clist;
  //struct join_list *jlist;
   /* Join and cell sizes */
  sizej = sizeof(join);
  sizebj = sizej + nq*(nd + 2)*sizeof(float);
  sizec = sizeof(cell);
  if (particle)
    sizep = sizeof(part);
  if (mype == 0){
    cout << "cell size " << sizec << " join size " << sizej <<
      " boundary join size " << sizebj;
    if (particle)
        cout << " particle size " << sizep << '\n';
    else 
        cout << '\n';
  }
  /* Cell and join array sizes */
  for (npe = 0; npe < limpe; npe++)
    nvalues[npe] = 2;
  temp[0] = sizec*(float)carrtop;
  temp[1] = sizej*(float)jarrtop;
  merge_array(temp, nvalues);
  tcarrmem = 0;
  tjarrmem = 0;
  if (mype == 0){
    for (npe = 0; npe < limpe; npe++){
      carrmem[npe] = temp[2*npe];
      tcarrmem += carrmem[npe];
      jarrmem[npe] = temp[2*npe+1];
      tjarrmem += jarrmem[npe];
    }
  }
  /* Join Memory */
  for (ilev = 0; ilev < nlevs; ilev++){
    njoin = 0;
    nbjoin = 0;
    nfcell = 1;
    for (id = 1; id <= nd; id++){
      nfcell *= lev[ilev]->nc[id-1];
      jptr = lev[ilev]->fjoin[id-1];
      while(jptr != NULL){
	njoin++;
	if (jptr->bound == cibound)
	  nbjoin++;
	jptr = jptr->nxt;
      }
    }
    memj = sizej*(float)(njoin - nbjoin) + sizebj*(float)nbjoin;
    /* Count cells and particles */
    ncell = 0;
    pcell = 0;
    cptr = lev[ilev]->fcell;
    while(cptr != NULL){
      ncell++;
      if (particle){
	pptr = cptr->pptr;
	while (pptr != NULL){
	  pcell++; 
	  pptr = pptr->nxt;
	}
      }
      cptr = cptr->nxt;
    }
    /* Cell Memory */
    memc = sizec*(float)ncell;
    /* Particle memory */
    memp = sizep*(float)pcell;
    /* Merge join and cell info */
    for (npe = 0; npe < limpe; npe++)
      nvalues[npe] = 8;
    temp[0] = (float)njoin;
    temp[1] = memj;
    temp[2] = (float)nfcell;
    temp[3] = (float)ncell;
    temp[4] = memc;
    temp[5] = (float)pcell;
    temp[6] = memp;
    temp[7] = workwt;
    merge_array(temp, nvalues);
    if (mype == 0){
      for (npe = 0; npe < limpe; npe++){
	pnjoin[ilev][npe] = (long long int)temp[8*npe];
	pmemj[ilev][npe] = temp[8*npe+1];
	pnfcell[ilev][npe] = (long long int)temp[8*npe+2];
	pncell[ilev][npe] = (long long int)temp[8*npe+3];
	pmemc[ilev][npe] = temp[8*npe+4];
	ppcell[ilev][npe] = temp[8*npe+5];
	pmemp[ilev][npe] = temp[8*npe+6];
	pworkwt[npe] = temp[8*npe+7];
      }
    }
  }
  /* Master prints info for each processor */
  if (mype == 0){
    for (npe = 0; npe < limpe; npe++){
      if (all){
	cout << "Processor " << npe << '\n';
	if (debug){
	  cout << "Level" << "        No cells" << "   Memory (Mb)" << 
	    "      No Joins" << "   Memory (Mb)" << "      Frac";
	  if (particle)
	    cout << "    No part" << "      Memory(Mb)" << '\n';
	  else 
	    cout << '\n';
	}
      }
      njoint = 0;
      memtj = 0;
      ncellt = 0;
      memtc = 0;
      pcellt = 0;
      memtp = 0;
      twork = 0;
      for (ilev = 0; ilev < nlevs; ilev++){
	cout.setf(ios::scientific);
	cout << setprecision(2);
	if (all){
	  if (debug){
	    cout << setw(5) << ilev << "    " 
		 << setw(12) <<  pncell[ilev][npe] << "       " << setw(6) 
		 << pmemc[ilev][npe]*1.0e-6 << "  " << setw(12) 
		 << pnjoin[ilev][npe] << "     "
		 << pmemj[ilev][npe]*1.0e-6 << "  " 
		 << (float)pncell[ilev][npe]/
	      (max(1.0, (float)pnfcell[ilev][npe]));
	    if (particle)
	      cout << setw(12) << ppcell[ilev][npe] << "       " << 
		pmemp[ilev][npe]*1.0e-6 <<  '\n';
	    else 
	      cout << '\n';
	  }
	}
	njoint += pnjoin[ilev][npe];
	memtj += pmemj[ilev][npe];
	ncellt += pncell[ilev][npe];
	twork += pworkwt[npe]*(float)pncell[ilev][npe]*ipow(ilev);		// added by BM, total work on this proc
	memtc += pmemc[ilev][npe];
	pcellt += ppcell[ilev][npe];
	memtp += pmemp[ilev][npe];
	/* Accumulate grand totals in pnjoin[ilev][0] etc */
	if (npe != 0){
	  pnjoin[ilev][0] +=  pnjoin[ilev][npe];
	  pmemj[ilev][0] +=  pmemj[ilev][npe];
	  pncell[ilev][0] +=  pncell[ilev][npe];
	  pmemc[ilev][0] +=  pmemc[ilev][npe];
	  pnfcell[ilev][0] +=  pnfcell[ilev][npe];
	  ppcell[ilev][0] +=  ppcell[ilev][npe];
	  pmemp[ilev][0] +=  pmemp[ilev][npe];
	}
      }
      if (all){
	cout << "Total    " << setw(12) << ncellt << "       "
	     << setw(6) << (float)(memtc*1.0e-6) << "  " << setw(12) 
	     << njoint << "     " << memtj*1.0e-6 << "  " ;
	if (particle)
	  cout << "       " << setw(12) << pcellt << "      " << setw(6)
	       << (float)(memtp*1.0e-6) << "  Total memory " << 
	    (memtc + memtj + memtp)*1.0e-6 << " Mb" << '\n';
	else
	  cout << "  Total memory " << (memtc + memtj)*1.0e-6 << " Mb"
	    << '\n';
	cout << "Cell array " << (float)(carrmem[npe])*1.0e-6 << 
	  " Mb, Join array " << (float)(jarrmem[npe])*1.0e-6 << 
	  " Mb, Total " << (float)(carrmem[npe] + jarrmem[npe])*1.0e-6
	     << " Mb" << " Workweight " << pworkwt[npe] << '\n'; 
	cout << "Total work " << setw(12) << twork << '\n';		// added by BM
      }
    }
    cout << "Grand total " << '\n';
    cout << "Level" << "        No cells" << "    Memory (Mb)" << 
      "      No Joins" << "  Memory (Mb)" << "      Frac";
    if (particle)
       cout << "    No part" << "    Memory (Mb)"  <<'\n';
    else
       cout << '\n';
    njoint = 0;
    memtj = 0;
    ncellt = 0;
    memtc = 0;
    pcellt = 0;
    memtp = 0;
    for (ilev = 0; ilev < nlevs; ilev++){
      cout.setf(ios::scientific);
      cout << setprecision(2);
      cout << setw(5) << ilev << "    " 
	   << setw(12) <<  pncell[ilev][0] << "       " << setw(6) 
	   << pmemc[ilev][0]*1.0e-6 << "  " << setw(12) 
	   << pnjoin[ilev][0] << "     "
	   << pmemj[ilev][0]*1.0e-6 << "  " 
	   << (float)pncell[ilev][0]/(max(1.0,(float)pnfcell[ilev][0]));
      if (particle)
	cout << "   " << setw(12) << ppcell[ilev][0] << "   " << 
	  setw(6) << pmemp[ilev][0]*1.0e-6 << '\n';
      else
	cout << '\n';
      njoint += pnjoin[ilev][0];
      memtj += pmemj[ilev][0];
      ncellt += pncell[ilev][0];
      memtc += pmemc[ilev][0];
      pcellt += ppcell[ilev][0];
      memtp += pmemp[ilev][0];
    }
    cout << "Total    " << setw(12) << ncellt <<  "       "
	 << setw(6) << memtc*1.0e-6 << "  " << setw(12) 
	 << njoint << "     " << memtj*1.0e-6 << " " ; 
    if (particle )
       cout << setw(12) << pcellt << "      " << setw(6) << 
         (float)(memtp*1.0e-6) << "  Total memory " << 
         (memtc + memtj + memtp)*1.0e-6 << " Mb" << '\n';
     else
        cout <<	"  Total memory " << (memtc + memtj)*1.0e-6 << " Mb"
          << '\n';
    cout << "Cell array " << (float)(tcarrmem)*1.0e-6 <<
      " Mb, Join array " << (float)(tjarrmem)*1.0e-6 << 
      " Mb, Total " << (float)(tcarrmem + tjarrmem)*1.0e-6 << " Mb "
      <<'\n'; 
  }
  return;
}
//end file mem_calc.c

//file i_to_it.c
/*
  + i_to_it - convert one dimensional array index to multidimensions.
     Description :
        Converts  the one dimensional  index i  to the  nd dimensional
        index in  it[ndim].  The array  is nc[0]xnc[1]..nc[ndim-1].The
        array indices are assumed to start at 1;
*/
void i_to_it(int ndim, int nc[], int i, int it[])
{
  int id, j, k;
  k = 1;
  for (id = 1; id < ndim; id++)
    k = k*nc[id-1];
  for (id = ndim; id > 1; id--){
    j = (int)((i - 1)/k);
    it[id-1] = 1 + j;
    i -= j*k;
    k = k/nc[id-2];
  }
  it[0] = i;
  return;
}
//end i_to_it.c

//file i_to_it2.c
/*
+ i_to_it2 - nd array indices.
     Description :
        Gets nd array  indices corresponding to 1D array  index i. The
        nd array is 2x2x....
*/
void i_to_it2 (int i, int it[])
{
  int id, j, k;
  k = ipow(nd-1);
  for (id = nd; id >1; id--){
    j = (int)((i - 1)/k);
    it[id-1] = 1 + j;
    i -= j*k;
    k = k/2;
  }
  it[0] = i;
  return;
}
//end i_to_it2.c

//file it_to_i.c
/*
+ it_to_i - convert ndim  dimensional array index to one dimensions.
     Description :
        Returns the  one dimensional index  i corresponding to  the nd
        dimensional indices in it[ndim]. The array indices are assumed
        to start at 1. i starts at 1.
*/
int it_to_i(int ndim, int nt[], int it[])
{
  int i, j, id;
  if (ndim ==0)
    return(1);
  i = it[0];
  j = nt[0];
  for (id = 2; id <= ndim; id++){
    i += (it[id-1] - 1)*j;
    j *=nt[id-1];
  }
  return(i);
}
//end it_to_i.c

//file it2_to_i.c
/*
+ it2_to_i - convert ndim dimensional array index to one dimensions.
     Description :
        Returns the  one dimensional index  i corresponding to  the nd
        dimensional indices  in it[ndim]. The array is  2x2x.. and the
        array indices are assumed to start at 1. i starts at 0. 
*/
int it2_to_i(int it[])
{
  int i, k, id;
  i = 0;
  k = 1;
  for (id = 1; id <= nd; id++){
    i += (it[id-1] - 1)*k;
   k = k*2;
  }
  return(i);
}
//end it2_to_i.c

//file dpdx_index.c
/*
+ dpdx_index - index function for internal boundary join dpdx.
     Description :
        This returns  the variable indices in the  dpdx gradient array
        at  coarse/fine internal  boundaries. id  is the  direction of
        gradient, jid is the direction  of the join.  The gradients to
        be  put  in  dpdx[il]   ..  dpdx[iu].  This  just  does  nd=2,
        3. Generalise to nd later.
*/
int dpdx_index(int jid, int id)
{
  if (nd == 2){
    return(0);
  }
  else if (nd == 3){
    if (jid == 1){
      if (id == 2){
	return(0);
      }
      else if (id == 3){
	return(nq);
      }
    }
    else if (jid == 2){
      if (id == 1){
	return(0);
      }
      else if (id == 3){
	return(nq);
      }
    }
    else if (jid == 3){
      if (id == 1){
	return(0);
      }
      else if (id == 2){
	return(nq);
      }
    }
  }
  return(0);
}
//end dpdx_index.c
