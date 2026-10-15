#undef MAIN
#include <iostream>
#include <iomanip>
#include <new>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "funcs.h"

extern void clear_mug();
extern void cal_diffl(int ilev, int id, float tfrac, struct cell *cptr, 
  float diffl[]);
extern void cal_diffr(int ilev, int id, float tfrac, struct cell *cptr, 
  float diffr[]);
extern bool cal_diffql(int id, struct cell *cptr, float diff[]);
extern bool cal_diffqr(int id, struct cell *cptr, float diff[]);
extern struct cell *new_cell();
extern void del_cell(struct cell *cptr);
extern void del_join(struct join *cptr);
extern void clean_cell(struct cell *cptr);
extern void del_cell_freelist();
extern void del_join_freelist();
extern float cell_vol(int ilev, float xc[]);
extern void rem_cell(int ilev, struct cell *cptr);
extern void de_frag();
extern void connect_cell(int ilev, int id, struct join *&start, 
  struct join *&end, struct cell *cptr);
extern void connect_par(int ilev, struct cell *cptr);
extern void assign_part(struct part_list *start);
extern void mk_partlist(struct part_list *&start, 
  struct part_list *&end);
extern void rm_part(struct part *pptr);
extern void mk_locarr();
extern struct cell *find_fcell(int &ilev, float xc[]);
extern int it_to_i(int ndim, int nt[], int it[]);
void cal_part();
void work();
void create_halo(int ilev);
void put_halo(int &count, int boundn, struct cell_list *&sbound, 
  float outbufr[]);
void receive_halo(int ilev, struct cell_list *&shalo, float xnd, 
  float inbufr[]);
void connect_halo(int ilev, struct cell_list *shalo);
void clear_halo(int ilev);
void com_halo(int com, int ilev, float tfrac);
void put_com(int com, int ilev, bool left, float tfrac, 
  float outbufr[]);
void put_com1(int ilev, int &count, struct cell_list *clist, 
  float outbufr[]);
void put_com2(int ilev, int &count, struct cell_list *clist, 
  float outbufr[]);
void put_com3(int ilev, int &count, struct cell_list *clist, 
  float outbufr[]);
void put_com4(int ilev, bool left, int &count, struct cell_list *clist,
  float tfrac, float outbufr[]);
void put_com5(int ilev, int &count, struct cell_list *clist,
  float outbufr[]);
void put_com6(int ilev, bool left, int &count, struct cell_list *clist,
  float outbufr[]);
void put_com7(int ilev, int &count, struct cell_list *clist, 
  float outbufr[]);
void put_com8(int ilev, int &count, struct cell_list *clist,
  float outbufr[]);
void put_com9(int ilev, bool left, int &count, struct cell_list *clist,
  float outbufr[]);
void receive_com(int com, int ilev, bool left, float inbufr[]);
void receive_com1(int ilev, bool left, float inbufr[]);
void receive_com2(int ilev, bool left, float inbufr[]);
void receive_com3(int ilev, bool left, float inbufr[]);
void receive_com4(int ilev, bool left, float inbufr[]);
void receive_com5(int ilev, bool left, float inbufr[]);
void receive_com6(int ilev, bool left, float inbufr[]);
void receive_com7(int ilev, bool left, float inbufr[]);
void receive_com8(int ilev, bool left, float inbufr[]);
void receive_com9(int ilev, bool left, float inbufr[]);
void zero_halo_cor(int ilev);
void ref_proc_bound(int ilev);
void put_ref(int ilev, struct cell_list *clist, float outbufr[]);
void receive_ref(int ilev, float inbufr[]);
void deref_proc_bound(int ilev, struct part_list *&start);
void com_cor(int ilev);
void put_cor(int ilev, bool left, int &count, float outbufr[]);
void receive_cor(int ilev, bool left, int &count, float inbufr[]);
void re_part();
void set_proc_bound();
struct part_list *mov_part(bool all, struct part_list *start);
void get_part_proc(struct part_list *start);
void broadcast(int *value);
void broadcast(int buf[], int nvalues);
void broadcast(float *value);
void broadcast(float buf[], int nvalues);
void broadcast(bool *value);
void merge_array(int buf[], int nvalues[]);
void merge_array(float buf[], int nvalues[]);
void all2all(int *npe, int outs[], float outbuf[], int ins[],
  float inbuf[], int *ios); 
float global_max(float value);
float global_sum(float value);
//file cal_part.c
/*
+ cal_part - Calculate partition boundaries.
     Description :
        This uses  work() to get  the work array  and then uses  it to
        calculate  new partition  boundaries  between processors.
*/
void cal_part()
{
  int i, gncnd, npe;
  float wproc, wright, twork, dxf, dxc;
  (void)dxc;
  /* Get work array */
  work();
  /* Total no of cells in nd direction on finest grid */
  gncnd = nc[nd-1]*ipow(nlevs-1);
  /* Total work */
  twork = wka[gncnd-1];
  /* Work per processor */
  wproc = twork/limpe;
  /* Calculate new processor boundaries */
  dxf = lev[nlevs-1]->dxbar[nd-1];
  pbd[0] = xl[nd-1];
  for (npe = 0; npe < limpe - 1; npe++){
    /* Work at  right boundary */
    wright = (npe + 1)*wproc;
    for (i = 0; i < gncnd; i++){
      if (wka[i] > wright)
	break;
    }
    pbd[npe+1] = xl[nd-1] + dxf*((float)i + 0.01);
  }
  pbd[limpe] = xr[nd-1];
  return;
}
//end cal_part.c

//file work.c
/*
+ work - Calculate work array.
     Description :
        Calculates the  array wka[nc[nd-1]]  where nc[nd-1] is  the nc
        for the finest  level. wka[i] is the amount  of work for cells
        with centres  to the left of  (i+1)dx, where dx is  the z mesh
        spacing  on the  finest level.   The array  is then  merged on
        processor 0 and distributed to all processors
*/
void work()
{
  int ilev, iw, ncf[MAXPE], nvalues[MAXPE], gncnd, i, npe;
  float dxf, xc, xlf, cwork, tv;
  struct cell *cptr;
  /* Max no of cells in nd direction on finest grid */
  gncnd = nc[nd-1]*ipow(nlevs-1);
  /* Max no of cells in nd direction on finest grid on this processor */
  ncf[0]= (int)((lev[nlevs-1]->xr[nd-1] 
    - lev[nlevs-1]->xl[nd-1]+ 0.1*lev[nlevs-1]->dxbar[nd-1])/
    lev[nlevs-1]->dxbar[nd-1]);
  /* Zero work density array */
  for (i = 0; i < gncnd; i++)
    wkda[i] = 0.0;
  /* Calculate work density array */
  dxf = lev[nlevs-1]->dxbar[nd-1];
  xlf = lev[nlevs-1]->xl[nd-1];
  for (ilev = 0; ilev < nlevs; ilev++){
    cwork = workwt*ipow(ilev);
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      xc = cptr->xc[nd-1] - xlf;
      iw = (int)(xc/dxf);
      wkda[iw] += cwork;
      cptr = cptr->nxt;
    }
  }
  /* Get sizes from other processors */
  for (npe = 0; npe < limpe; npe++)
    nvalues[npe] = 1;
  merge_array(ncf, nvalues);
  /* Broadcast sizes */
  broadcast(ncf, limpe);
  /* Merge work density array */
  for (npe = 0; npe < limpe; npe++)
    nvalues[npe] = ncf[npe];
  merge_array(wkda, nvalues);
  /* Calculate cumulative work array */
  tv = 0.0;
  for (i = 0; i < gncnd; i++){
    tv += (float)wkda[i];
    wka[i] = tv;
  }
  /* Broadcast work array */
  broadcast(wka, gncnd);
  return;
}
//end work.c

//file create_halo.c
/*
+ create_halo - Creates ilev halo and processor boundary cell lists.
     Description :
        Creates  halo and  boundary  cell lists  for  level ilev.   It
        assumes  that  the plevel  structure  for  halo and  processor
        boundary cells has already  been created. It connects the halo
        to the corresponding boundary cells if ilev > 0. It also calls
        connect_par  for  boundary  cells  to  ensure  that  they  are
        connected to their parents.  This assumes that the ilev-1 halo
        has been created. The  processor boundaries are calculated for
        all levels and stored in levpbd[ilev][npe] (particles).
*/
void create_halo(int ilev)
{
  int maxcom, icell, count, lboundn, rboundn, ios, npe, nprocs, 
    maxprocs;
  (void)icell;
  float xlnd, xrnd, dx, *inbufr, *outbufr;
  bool info;
  (void)info;
  struct cell *cptr;
  struct cell_list *clist;
  struct proc_list *plist;
  /* Return if single processor (safety) */
  if (limpe == 1)
    return;
  /* Clear (safety) */
  clear_halo(ilev);
  /* Mesh spacing in direction nd */
  dx = lev[ilev]->dxbar[nd-1];
  /*Left/right processor boundaries */
  xlnd = lev[ilev]->xl[nd-1];
  xrnd = lev[ilev]->xr[nd-1];
  /* Calculate processor for this level (store in global levpbd) */
  levpbd[ilev][0] = xl[nd-1];
  for (npe = 1; npe <  limpe; npe++)
    levpbd[ilev][npe] = xl[nd-1] + 
      dx*(int)((pbd[npe] - xl[nd-1] + 0.5*dx)/dx);
  levpbd[ilev][limpe] = xr[nd-1];
  /* Assemble left boundary cell list */
  lboundn = 0;
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    /* Connect to parent for ilev > coarsest level */
    if (ilev > 0)
      connect_par(ilev, cptr);
    if (cptr->xc[nd-1] - xlnd < dx){
      /* Left boundary cell */
      lboundn++;
      clist = new cell_list;
      clist->cptr = cptr;
      clist->nxt = plev[ilev]->lbound;
      plev[ilev]->lbound = clist;
    }
    cptr = cptr->nxt;
  }
  /* Zero nlproc/nrproc for safety */
  plev[ilev]->nlproc = 0;
  plev[ilev]->nrproc = 0;
  /* Check if we have left boundary cells on this level */
  if (lboundn != 0){
    /* Create left recipient list (largest npe first) */
    nprocs = 0;
    /* Check for periodic */
    if ((bcl[nd-1] == bcper) && (xlnd - xl[nd-1] < 0.1*dx)){
      for (npe = limpe - 1; npe > mype; npe--){
	if (xr[nd-1] - levpbd[ilev][npe+1] < 0.5*dx){
	  /* Add to list */
	  nprocs++;
	  plist = new proc_list;
	  plist->proc = npe;
	  plist->nxt = plev[ilev]->lproc;
	  plev[ilev]->lproc = plist;
	}
      }
    }
    for (npe = mype - 1; npe >= 0; npe--){
      if ((xlnd + 0.5*dx > levpbd[ilev][npe+1])
	  && (xlnd - 0.5*dx < levpbd[ilev][npe+1])){
	/* Add to list */
	nprocs++;
	plist = new proc_list;
	plist->proc = npe;
	plist->nxt = plev[ilev]->lproc;
	plev[ilev]->lproc = plist;
      }
    }
    plev[ilev]->nlproc = nprocs;
  }
  /* Assemble right boundary cell list */
  rboundn = 0;
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    /* Connect to parent for ilev > coarsest level */
    if (ilev > 0)
      connect_par(ilev, cptr);
    if (xrnd - cptr->xc[nd-1] < dx){
      /* Right boundary cell */
      rboundn++;
      clist = new cell_list;
      clist->cptr = cptr;
      clist->nxt = plev[ilev]->rbound;
      plev[ilev]->rbound = clist;
    }
    cptr = cptr->nxt;
  }
  /* Check if we have right boundary cells on this level */
  if (rboundn != 0){
    /* Create right recipient list (largest npe first) */
    nprocs = 0;
    for (npe = limpe - 1; npe > mype; npe--){
      if ((xrnd + 0.5*dx > levpbd[ilev][npe])
	  && (xrnd - 0.5*dx < levpbd[ilev][npe])){
	/* Add to list */
	nprocs++;
	plist = new proc_list;
	plist->proc = npe;
	plist->nxt = plev[ilev]->rproc;
	plev[ilev]->rproc = plist;
      }
    }
    /* Check for periodic */
    if ((bcl[nd-1] == bcper) && (xr[nd-1] - xrnd < 0.1*dx)){
      for (npe = mype - 1; npe >= 0; npe--){
	if (levpbd[ilev][npe] - xl[nd-1] < 0.5*dx){
	  /* Add to to list */
	  nprocs++;
	  plist = new proc_list;
	  plist->proc = npe;
	  plist->nxt = plev[ilev]->rproc;
	  plev[ilev]->rproc = plist;
	}
      }
    }
    plev[ilev]->nrproc = nprocs;
  }
  maxprocs = max(plev[ilev]->nlproc, plev[ilev]->nrproc);
  maxprocs = (int)global_max((float)maxprocs);
  maxcom = max(lboundn*plev[ilev]->nlproc, rboundn*plev[ilev]->nrproc);
  maxcom = (int)global_max((float)maxcom);
  plev[ilev]->maxcom = maxcom;
  /* Create communications arrays */
  inbufr = new float[nd*maxcom+maxprocs+SAFETY];
  outbufr = new float[nd*maxcom+maxprocs+SAFETY];
  /* Zero ins/outs (safety) */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Check for left processors */
  if (plev[ilev]->lproc != NULL){
    /* Initialize count */
    count = -1;
    /* Fill outbufr for left processors */
    plist = plev[ilev]->lproc;
    while (plist != NULL){
      put_halo(count, lboundn, plev[ilev]->lbound, outbufr);
      /* Outs for this processor */
      outs[plist->proc] = nd*lboundn + 1;
      plist = plist->nxt;
    }
  }
  /* Send  cell data to left processors */
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  //  parc_xdatar(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Set right sending processor */
  plev[ilev]->rsproc = -1;
  for (npe = 0; npe < limpe; npe++){
    if (ins[npe] > 0){
      plev[ilev]->rsproc = npe;
      break;
    }
  }
  if (plev[ilev]->rsproc >= 0)
    /* Receive boundary cell information for right halo */
    receive_halo(ilev, plev[ilev]->rhalo, xrnd + 0.5*dx, inbufr);
  else
    plev[ilev]->rhalo = NULL;
  /* Zero ins/outs (safety) */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Check for valid right processor */
  if (plev[ilev]->rproc != NULL){
    /* Initialize count */
    count = -1;
    /* Fill outbufr for right processors */
    plist = plev[ilev]->rproc;
    while (plist != NULL){
      put_halo(count, rboundn, plev[ilev]->rbound, outbufr);
      /* Outs for this processor */
      outs[plist->proc] = nd*rboundn + 1;
      plist = plist->nxt;
    }
  }
  /* Send cell data to right processors */
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Set left sending processor */
  plev[ilev]->lsproc = -1;
  for (npe = 0; npe < limpe; npe++){
    if (ins[npe] > 0){
      plev[ilev]->lsproc = npe;
      break;
    }
  }
  if (plev[ilev]->lsproc >= 0)
     /* Receive boundary cell information for left halo */
     receive_halo(ilev, plev[ilev]->lhalo, xlnd - 0.5*dx, inbufr);
  else
    plev[ilev]->lhalo = NULL;
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  /* Connect halo for ilev > 0 */
  if (ilev > 0){
    /* Connect left halo cells */
    connect_halo(ilev, plev[ilev]->lhalo);
    /* Connect right halo cells */
    connect_halo(ilev, plev[ilev]->rhalo);
  }
  return;
}
//end create_halo.c

//file put_halo.c
/*
+ put_halo - Put halo information in outbufr.
     Description :
        Put  the  positions  in  nd  plane and  refinement  status  in
        outbufr. count  is incremented and must  have been initialized
        before  calling  this  function.   outs  is  set  outside  the
        function.
*/
void put_halo(int &count, int boundn, struct cell_list *&sbound, 
  float outbufr[])
{
  int id;
  struct cell * cptr;
  struct cell_list *clist;
  /* No of boundary cells */
  count++;
  outbufr[count] = (float)boundn;
  /* Cell boundary cell positions in nd plane */
  clist = sbound;
  while (clist != NULL){
    cptr = clist->cptr;
    count++;
    if (cptr->refined)
      outbufr[count] = 1.0;
    else
      outbufr[count] = 0.0;
    for (id = 1; id < nd; id++){
      count++;
      outbufr[count] = cptr->xc[id-1];
    }
    clist = clist->nxt;
  }
  return;
}
//file receive_halo.c
/*
+ receive_halo - Receive halo information.
     Description :
        Receives  the halo  information  and creates  halo cell  list.
        shalo is a pointer to the  start of the halo cell list, xnd is
        the nd  coordinate of the  halo cells and inbufr  contains the
        cell information.
*/
void receive_halo(int ilev, struct cell_list *&shalo, float xnd, 
  float inbufr[])
{
  int count, ncell, icell, id;
  struct cell *cptr;
  struct cell_list *clist, *tclist=0;
  count = 0;
  /* No of halo cells */
  ncell = (int)inbufr[count];
  if (ncell == 0){
    /* No halo list */
    shalo = NULL;
    return;
  }
  else{
    /* Create halo cell list */
    clist = new cell_list;
    shalo = clist;
    for (icell = 1; icell <= ncell; icell++){
      cptr = new_cell();
      cptr->halo = true;
      clist->cptr = cptr;
      /* Refined or not */
      count++;
      cptr->refined = ((int)inbufr[count] == 1);
      /* Left halo positions in nd plane */
      for (id = 1; id < nd; id++){
	count++;
	cptr->xc[id-1] = inbufr[count]; 
      }
      /* Cell nd position */
      cptr->xc[nd-1] = xnd;
      /* Cell volume */
      cptr->vol = cell_vol(ilev, cptr->xc);
      clist->nxt = new cell_list;
      tclist = clist;
      clist = clist->nxt;
    }
    /* Terminate halo cell list */
    delete clist;
    tclist->nxt= NULL;
  }
  return;
}
//end receive_halo.c 

//file connect_halo.c
/*
+ connect_halo - Connect halo cells.
     Description :
        Connects  the  halo cells  on  level  ilev  i.e.  creates  the
        appropriate joins and sets  their boundary flags. shalo is the
        start  of the list  (left or  right). Done  by using  the cell
        positions. Joins that  are not in the nd  direction are put in
        the halo join lists.
*/
void connect_halo(int ilev, struct cell_list *shalo)
{
  int id;
  struct cell_list *clist;
  struct cell *cptr;
  struct join *start, *end, *jptr;
  (void)jptr;
  /* Return if no halo cells */
  if (shalo == NULL)
    return;
  /* Connect cells to parents (assumes ilev - 1 halo exists) */
  clist = shalo;
  while (clist != NULL){
    cptr = clist->cptr;
    connect_par(ilev, cptr);
    clist = clist->nxt;
  }
  /* Connect cells */
  for (id = 1; id <= nd; id++){
    clist = shalo;
    start = NULL;
    end = NULL;
    while (clist != NULL){
      cptr = clist->cptr;
      if (id != nd)
	connect_cell(ilev, id, start, end, cptr);
      else
	connect_cell(ilev, id, lev[ilev]->fjoin[nd-1], 
	  lev[ilev]->ljoin[nd-1], cptr);
      clist = clist->nxt;
    }
    if (id < nd){
      if (shalo == plev[ilev]->lhalo)
	/* Set start of left halo join list */
	plev[ilev]->lfjoin[id-1] = start;
      else
	/* Set start of right halo join list */
	plev[ilev]->rfjoin[id-1] = start;
    }
  }
  return;
}
//end connect_halo.c

//file clear_halo.c
/*
+ clear_halo - Clear halo and boundary cell lists.
     Description :
        Removes the  halo and  boundary cell lists  for level  ilev so
        that they can be properly recreated.
*/
void clear_halo(int ilev)
{
  int id, icell;
  (void)icell;
  bool refined;
  (void)refined;
  struct cell_list *clist, *tclist;
  struct cell *cptr, *pcptr;
  (void)pcptr;
  struct proc_list *plist, *tplist;
  /* Clear Multigrid */
  if (grav || mugi)
    clear_mug();
  /* Remove left boundary list (if any) */
  clist = plev[ilev]->lbound;
  while (clist != NULL){
    tclist = clist;
    clist = clist->nxt;
    delete tclist;
  }
  plev[ilev]->lbound = NULL;
  /* Clear left halo cells (if any) */
  clist = plev[ilev]->lhalo;
  while (clist != NULL){
    cptr = clist->cptr;
    if (!cptr->empty)
      rem_cell(ilev, cptr);
    tclist = clist;
    clist = clist->nxt;
    delete tclist;
  }
  plev[ilev]->lhalo = NULL;
  /* Remove right boundary list (if any) */
  clist = plev[ilev]->rbound;
  while (clist != NULL){
    tclist = clist;
    clist = clist->nxt;
    delete tclist;
  }
  plev[ilev]->rbound = NULL;
  /* Clear right halo cells (if any) */
  clist = plev[ilev]->rhalo;
  while (clist != NULL){
    cptr = clist->cptr;
    if (!cptr->empty)
      rem_cell(ilev, cptr);
    tclist = clist;
    clist = clist->nxt;
    delete tclist;
  }
  plev[ilev]->rhalo = NULL;
  /* Set start of left/right halo join list to NULL */
  for (id = 1; id < nd; id++){
    plev[ilev]->lfjoin[id-1] = NULL;
    plev[ilev]->rfjoin[id-1] = NULL;
  }
  /* Remove left processor list */
  plist = plev[ilev]->lproc;
  while (plist != NULL){
    tplist = plist;
    plist = plist->nxt;
    delete tplist;
  }
  plev[ilev]->lproc = NULL;
  /* Remove right processor list */
  plist = plev[ilev]->rproc;
  while (plist != NULL){
    tplist = plist;
    plist = plist->nxt;
    delete tplist;
  }
  plev[ilev]->rproc = NULL;
  return;
}
//end clear_halo.c

//file com_halo.c
/*
+ com_halo - Communicate solution to halo cells.
     Description :
        This fills the halo cells  on level ilev with information from
        the boundary  cells on the neighbouring  processor. The nature
        of the information is determined by com as follows:
        com1  qsol
        com2  pa
        com3  qa, dqdt
        com4  dpdx
        com5  qa
        com6  dqdx
        com7  halo cell coarse/fine interpolant
        com8  tells  halo cell  whether corresponding boundary  cell is
              refined
        com9 outer halo join flag.
*/
void com_halo(int com, int ilev, float tfrac)
{
  int ios;
  float *inbufr, *outbufr;
  /* Return for limpe=1 (safety) */
  if (limpe == 1)
    return;
  /* Create communications arrays */
  inbufr = new float [(2*nq+nq*(nd-1)+2)*plev[ilev]->maxcom+SAFETY];
  outbufr = new float [(2*nq+nq*(nd-1)+2)*plev[ilev]->maxcom+SAFETY];
  /* Assemble information from left boundary cells */
  put_com(com, ilev, true, tfrac, outbufr); 
  /* Send information to left processor */
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Put received data in right halo cells */
  receive_com(com, ilev, false, inbufr);
  /* Assemble information from right boundary cells */
  put_com(com, ilev, false, tfrac, outbufr);
  /* Send information to right processor */
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Put received data in left halo cells */
  receive_com(com, ilev, true, inbufr);
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  return;
}
//end com_halo.c

//file put_com.c
/*
+ file put_com -  Put com values in outbufr.
     Description :
        Puts the values determined by com in outbufr.
*/
void put_com(int com, int ilev, bool left, float tfrac, float outbufr[])
{
  int npe, nprocs, count;
  struct cell_list *start;
  struct proc_list *plist, *splist;
  /* Zero ins/outs arrays for safety */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  if (left){
  /* Return if no data */
    if (plev[ilev]->lproc == NULL)
      return;
    /* Left boundary */
    start = plev[ilev]->lbound;
    /* Start of processor list */
    plist = plev[ilev]->lproc;
    splist = plist;
    /* No of left processors */
    nprocs = plev[ilev]->nlproc;
  }
  else{
  /* Return if no data */
    if (plev[ilev]->rproc == NULL)
      return;
    /* Right boundary */
    start = plev[ilev]->rbound;
    /* Start of processor list */
    plist = plev[ilev]->rproc;
    splist = plist;
    /* No of right processors */
    nprocs = plev[ilev]->nrproc;
  }
  /* Set count */
  count = -1;
  /* Put data in outbufr as appropriate */
  if (com == com1){
    while (plist != NULL){
      put_com1(ilev, count, start, outbufr);
      plist = plist->nxt;
    }
  }
  else if (com == com2){
    while (plist != NULL){
      put_com2(ilev, count, start, outbufr);
      plist = plist->nxt;
    }
  }
  else if (com == com3){
    while (plist != NULL){
      put_com3(ilev, count, start, outbufr);
      plist = plist->nxt;
    }
  }
  else if (com == com4){
    while (plist != NULL){
      put_com4(ilev, left, count, start, tfrac, outbufr);
      plist = plist->nxt;
    }
  }
  else if (com == com5){
    while (plist != NULL){
      put_com5(ilev, count, start, outbufr);
      plist = plist->nxt;
    }
  }
  else if (com == com6){
    while (plist != NULL){
      put_com6(ilev, left, count, start, outbufr);
      plist = plist->nxt;
    }
  }
  else if (com == com7){
    while (plist != NULL){
      put_com7(ilev, count, start, outbufr);
      plist = plist->nxt;
    }
  }
  else if (com == com8){
    while (plist != NULL){
      put_com8(ilev, count, start, outbufr);
      plist = plist->nxt;
    }
  }
  else if (com == com9){
    while (plist != NULL){
      put_com9(ilev, left, count, start, outbufr);
      plist = plist->nxt;
    }
  }
  /* Set outs */
  while (splist != NULL){
    outs[splist->proc] = (count + 1)/nprocs;
    splist = splist->nxt;
  }
  return;
}
//end put_com.c

//file put_com1.c
/*
+ file put_com1 -  Put qsol in outbufr.
     Description :
        Puts  qsol  from  the  boundary  cells  in  outbufr.  Side  is
        specified by left.
*/
void put_com1(int ilev, int &count, struct cell_list *clist, 
  float outbufr[])
{
  struct cell *cptr;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Solution quality */
    count++;
    outbufr[count] = cptr->qsol;
    clist = clist->nxt;
  }
  return;
}
//end put_com1.c

//file put_com2.c
/*
+ file put_com2 -  Put pa in outbufr.
     Description :
        Puts pa from the boundary  cells in outbufr. Side is specified
        by left.
*/
void put_com2(int ilev, int &count, struct cell_list *clist, 
  float outbufr[])
{
  int i;
  struct cell *cptr;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Primitives (pa) */
    for (i = 0; i < nq; i++){
      count++;
      outbufr[count] = cptr->pa[i];
    }
    clist = clist->nxt;
  }
  return;
}
//end put_com2.c

//file put_com3.c
/*
+ file put_com3 -  Put qa and dqdt  in outbufr.
     Description :
        Puts qa and  dqdt from the boundary cells  in outbufr. Side is
        specified by left.
*/
void put_com3(int ilev, int &count, struct cell_list *clist, 
  float outbufr[])
{
  int i;
  struct cell *cptr;
  while (clist != NULL){
    cptr = clist->cptr;
    for (i = 0; i < nq; i++){
      /* Conservatives (qa) */
      count++;
      outbufr[count] = cptr->qa[i];
      /* dqdt */
      count++;
      outbufr[count] = cptr->dqdt[i];
    }
    clist = clist->nxt;
  }
  return;
}
//end put_com3.c

//file put_com4.c
/*
+ file put_com4 -  Put dpdx  in outbufr.
     Description :
        Puts  dpdx  from  the  boundary  cells  in  outbufr.  Side  is
        specified by left.
*/
void put_com4(int ilev, bool left, int &count, struct cell_list *clist,
  float tfrac, float outbufr[])
{
  int i;
  float diff[MAXQ];
  struct cell *cptr;
  while (clist != NULL){
    cptr = clist->cptr;
    if (left)
      /* Right difference of primitives in direction nd */
      cal_diffr(ilev, nd, tfrac, cptr, diff);
    else
      /* Left difference of primitives in direction nd */
      cal_diffl(ilev, nd, tfrac, cptr, diff);
    for (i = 0; i < nq; i++){
      /* dpdx */
      count++;
      outbufr[count] = diff[i];
    }
    clist = clist->nxt;
  }
  return;
}
//end put_com4.c

//file put_com5.c
/*
+ file put_com5 -  Put qa in outbufr.
     Description :
        Puts qa from the boundary  cells in outbufr. Side is specified
        by left.
*/
void put_com5(int ilev, int &count, struct cell_list *clist,
  float outbufr[])
{
  int i;
  struct cell *cptr;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Conservatives (qa) */
    for (i = 0; i < nq; i++){
      count++;
      outbufr[count] = cptr->qa[i];
    }
    clist = clist->nxt;
  }
  return;
}
//end put_com5.c

//file put_com6.c
/*
+ file put_com6 -  Put qa and dqdx in outbufr.
     Description :
        Puts qa and  dqdx from the boundary cells  in outbufr. Side is
        specified by left.
*/
void put_com6(int ilev, bool left, int &count, struct cell_list *clist,
  float outbufr[])
{
  int i;
  float diff[MAXQ];
  bool got;
  struct cell *cptr;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Calculate differences */
    if (left)
      /* Left boundary */
      got = cal_diffqr(nd, cptr, diff);
    else
      /* right boundary */
      got = cal_diffql(nd, cptr, diff);
    /* Whether or not we have differences */
    count++;
    if (got)
      outbufr[count] = 1.0;
    else
      outbufr[count] = -1.0;
    /* Differences (dqdx) */
    for (i = 0; i < nq; i++){
      count++;
      outbufr[count] = diff[i];
    }
    clist = clist->nxt;
  }
  return;
}
//end put_com6.c

//file put_com7.c
/*
+ file put_com7 -  Put coarse/fine interpolant in outbufr.
     Description :
        Puts the  flux coarse/fine correction from  the boundary joins
        in outbufr. Side is specified by left.
*/
void put_com7(int ilev, int &count, struct cell_list *clist, 
  float outbufr[])
{
  int id, i;
  struct cell *cptr;
  struct join *jptr;
  while (clist != NULL){
    cptr = clist->cptr;
    for (id = 1; id < nd; id++){
      jptr = cptr->ljoin[id-1];
      if (jptr->bound == cibound){
	/* Coarse/fine boundary join */
	for (i = 0; i < nq; i++){
	  count++;
	  outbufr[count] = jptr->pa0[i];
	  count++;
	  outbufr[count] = jptr->dpdt[i];
	}
	for (i = 0; i < nq*(nd - 1); i++){
	  count++;
	  outbufr[count] = jptr->dpdx[i];
	}
      }
    }
    clist = clist->nxt;
  }
  return;
}
//end put_com7.c

//file put_com8.c
/*
+ file put_com8 -  Halo cell refinement flag.
     Description :
        Puts 1.0 in outbufr if boundary cell is refined, -1.0 if it is
        not.   This is to  ensure that  the halo  cells get  the right
        boundary flags. Side is specified by left.
*/
void put_com8(int ilev, int &count, struct cell_list *clist,
  float outbufr[])
{
  struct cell *cptr;
  while (clist != NULL){
    cptr = clist->cptr;
    count++;
    if (cptr->refined)
      outbufr[count] = 1.0;
    else
      outbufr[count] = -1.0;
    clist = clist->nxt;
  }
  return;
}
//end put_com8.c

//file put_com9.c
/*
+ file put_com9 -  Set halo outer join flag.
     Description :
        Puts a boundary  flag in outbufr.  This is  to ensure that the
        halo cells get the right  boundary flags. Side is specified by
        left.
*/
void put_com9(int ilev, bool left, int &count, struct cell_list *clist,
  float outbufr[])
{
  struct cell *cptr;
  while (clist != NULL){
    cptr = clist->cptr;
    count++;
    if (left)
      outbufr[count] = (float)cptr->rjoin[nd-1]->bound;
    else
      outbufr[count] = (float)cptr->ljoin[nd-1]->bound;
    clist = clist->nxt;
  }
  return;
}
//end put_com9.c

//file receive_com.c
/*
+ file receive_com -  Receive com values from inbufr.
     Description :
        Receives the values determined by  com in inbufr and puts them
        in the halo or boundary cells. Side is specified by left.
*/
void receive_com(int com, int ilev, bool left, float inbufr[])
{
  /* Return if no data */
  if (left && (plev[ilev]->lsproc < 0))
    return;
  if (!left && (plev[ilev]->rsproc < 0))
    return;
  if (com == com1)
    receive_com1(ilev, left, inbufr);
  else if (com == com2)
    receive_com2(ilev, left, inbufr);
  else if (com == com3)
    receive_com3(ilev, left, inbufr);
  else if (com == com4)
    receive_com4(ilev, left, inbufr);
  else if (com == com5)
    receive_com5(ilev, left, inbufr);
  else if (com == com6)
    receive_com6(ilev, left, inbufr);
  else if (com == com7)
    receive_com7(ilev, left, inbufr);
  else if (com == com8)
    receive_com8(ilev, left, inbufr);
  else if (com == com9)
    receive_com9(ilev, left, inbufr);
  return;
}
//end receive_com.c

//file receive_com1.c
/*
+ file receive_com1 -  Receive qsol.
     Description :
        Receives qsol from inbufr and  puts it in the halo cells. Side
        is specified by left.
*/
void receive_com1(int ilev, bool left, float inbufr[])
{
  int count;
  struct cell *cptr;
  struct cell_list *clist;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Solution quality */
    count++;
    cptr->qsol = inbufr[count];
    clist = clist->nxt;
  }
  return;
}
//end receive_com1.c

//file receive_com2.c
/*
+ file receive_com2 -  Receive pa.
     Description :
        Receives pa from inbufr and puts it in the halo cells. Side is
        specified by left.
*/
void receive_com2(int ilev, bool left, float inbufr[])
{
  int count, i;
  struct cell *cptr;
  struct cell_list *clist;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Primitives (pa) and zero dqdt for safety in mg.c */
    for (i = 0; i < nq; i++){
      count++;
      cptr->pa[i] = inbufr[count];
      cptr->dqdt[i] = 0.0;
    }
    clist = clist->nxt;
  }
  return;
}
//end receive_com2.c

//file receive_com3.c
/*
+ file receive_com3 -  Receive qa and dqdt.
     Description :
        Receives  qa and  dqdt from  inbufr and  puts it  in  the halo
        cells. Side is specified by left.
*/
void receive_com3(int ilev, bool left, float inbufr[])
{
  int count, i;
  struct cell *cptr;
  struct cell_list *clist;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    for (i = 0; i < nq; i++){
      /* Conservatives (qa)*/
      count++;
      cptr->qa[i] = inbufr[count];
      /* dqdt */
      count++;
      cptr->dqdt[i] = inbufr[count];
    }
    clist = clist->nxt;
  }
  return;
}
//end receive_com3.c

//file receive_com4.c
/*
+ file receive_com4 -  Receive dpdx.
     Description :
        Receives dpdx from inbufr and  puts it in the halo cells. Side
        is specified by left.
*/
void receive_com4(int ilev, bool left, float inbufr[])
{
  int count, i;
  struct cell *cptr;
  struct cell_list *clist;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    for (i = 0; i < nq; i++){
      /* dpdx */
      count++;
      cptr->diff[i] = inbufr[count];
    }
    clist = clist->nxt;
  }
  return;
}
//end receive_com4.c

//file receive_com5.c
/*
+ file receive_com5 -  Receive qa.
     Description :
        Receives qa from inbufr and puts it in the halo cells. Side is
        specified by left.
*/
void receive_com5(int ilev, bool left, float inbufr[])
{
  int count, i;
  struct cell *cptr;
  struct cell_list *clist;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Conservatives (qa)*/
    for (i = 0; i < nq; i++){
      count++;
      cptr->qa[i] = inbufr[count];
    }
    clist = clist->nxt;
  }
  return;
}
//end receive_com5.c

//file receive_com6.c
/*
+ file receive_com6 -  Receive dqdx.
     Description :
        Receives dqdx from inbufr and  puts it in the halo cells. Side
        is specified by left.
*/
void receive_com6(int ilev, bool left, float inbufr[])
{
  int count, i;
  struct cell *cptr;
  struct cell_list *clist;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Whether or not we have differences */
    count++;
    cptr->gotdiff = inbufr[count] > 0;
    for (i = 0; i < nq; i++){
      /* Difference (dqdx)*/
      count++;
      cptr->diff[i] = inbufr[count];
    }
    clist = clist->nxt;
  }
  return;
}
//end receive_com6.c

//file receive_com7.c
/*
+ file receive_com7 -  Receive coarse/fine interpolant from inbufr.
     Description :
        Receives  the flux coarse/fine  interpolant from  the boundary
        joins in inbufr. Side is specified by left.
*/
void receive_com7(int ilev, bool left, float inbufr[])
{
  int id, i, count;
  struct cell_list *clist;
  struct cell *cptr;
  struct join *jptr;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    for (id = 1; id < nd; id++){
      jptr = cptr->ljoin[id-1];
      if (jptr->bound == cibound){
	/* Coarse/fine boundary join */
	for (i = 0; i < nq; i++){
	  count++;
	  jptr->pa0[i] = inbufr[count];
	  count++;
	  jptr->dpdt[i] = inbufr[count];
	}
	for (i = 0; i < nq*(nd - 1); i++){
	  count++;
	  jptr->dpdx[i] = inbufr[count];
	}
      }
    }
    clist = clist->nxt;
  }
  return;
}
//end receive_com7.c

//file receive_com8.c
/*
+ file receive_com8 -  Receive refined from inbufr.
     Description :
        Receives  1 in inbufr  if the  corresponding boundary  cell is
        refined, 0 if it is not.  Puts 1.0 in outbufr if boundary cell
        is refined, -1.0 if it is not. Side is specified by left.
*/
void receive_com8(int ilev, bool left, float inbufr[])
{
  int count;
  struct cell_list *clist;
  struct cell *cptr;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    count++;
    cptr->refined = inbufr[count] > 0;
    clist = clist->nxt;
  }
  return;
}
//end receive_com8.c

//file receive_com9.c
/*
+ file receive_com9 -  Receive halo join boundary flag from inbufr.
     Description :
        Receives  the outer  boundary  flag in  outbufr.   This is  to
        ensure that the halo cells get the right boundary flags.  Side
        is specified by left.
*/
void receive_com9(int ilev, bool left, float inbufr[])
{
  int count;
  struct cell_list *clist;
  struct cell *cptr;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo;
  count = -1;
  while (clist != NULL){
    cptr = clist->cptr;
    count++;
    if (left)
      cptr->ljoin[nd-1]->bound = (int)inbufr[count];
    else
      cptr->rjoin[nd-1]->bound = (int)inbufr[count];
    clist = clist->nxt;
  }
  return;
}
//end receive_com9.c

//file zero_halo_cor.c
/*
+ zero_halo_cor - Zero halo cell coarse/fine correction.
     Description :
        This  zeros  the  halo  cell join  correction  at  coarse/fine
        boundaries on level ilev;
*/
void zero_halo_cor(int ilev)
{
  int id, i;
  struct cell *cptr;
  struct cell_list *clist;
  struct join *jptr;
  /* Left halo */
  clist = plev[ilev]->lhalo;
  while (clist != NULL){
    cptr = clist->cptr;
    for (id = 1; id <= nd; id++){
      /* Left joins */
      jptr = cptr->ljoin[id-1];
      if ((jptr->bound == cibound) && (jptr->cor != NULL)){
	/* Coarse/fine boundary */
	jptr->done = 0;
	for (i = 0; i < nq; i++)
	  jptr->cor[i] = 0.0;
      }
      /* Right joins */
      jptr = cptr->rjoin[id-1];
      if ((jptr->bound == cibound) && (jptr->cor != NULL)){
	/* Coarse/fine boundary */
	jptr->done = 0;
	for (i = 0; i < nq; i++)
	  jptr->cor[i] = 0.0;
      }
    }
    clist = clist->nxt;
  }
  /* Right halo */
  clist = plev[ilev]->rhalo;
  while (clist != NULL){
    cptr = clist->cptr;
    for (id = 1; id <= nd; id++){
      /* Left joins */
      jptr = cptr->ljoin[id-1];
      if ((jptr->bound == cibound)&& (jptr->cor != NULL)){
	/* Coarse/fine boundary */
	for (i = 0; i < nq; i++)
	  jptr->cor[i] = 0.0;
      }
      /* Right joins */
      jptr = cptr->rjoin[id-1];
      if ((jptr->bound == cibound)&& (jptr->cor != NULL)){
	/* Coarse/fine boundary */
	for (i = 0; i < nq; i++)
	  jptr->cor[i] = 0.0;
      }
    }
    clist = clist->nxt;
  }
  return;
}
//end zero_halo_cor.c

//file ref_proc_bound.c
/*
+ ref_proc_bound - Refine processor boundary cells.
     Description :
        Creates cells on level  ilev that should exist on neighbouring
        processors.    The   lists   are   in   plev[ilev]->lref   and
        plev[ilev]->rref.
*/
void ref_proc_bound(int ilev)
{
  int lncell, rncell, npe, ios, maxcom, nproc;
  float *inbufr, *outbufr, dx, xnd, levpbd[MAXPE];
  struct cell *cptr;
  (void)cptr;
  struct cell_list *clist;
  bool info;
  /* Count cells to be sent to left processor */
  clist = plev[ilev]->lref;
  lncell = 0;
  while (clist != NULL){
    lncell++;
    clist = clist->nxt;
  }
  /* Count cells to be sent to right processor */
  clist = plev[ilev]->rref;
  rncell = 0;
  while (clist != NULL){
    rncell++;
    clist = clist->nxt;
  }
  /* Calculate global maxcom */
  maxcom = (int)global_max((float)max(lncell, rncell));
  /* Create communications arrays */
  inbufr = new float[(nd+nq+2)*maxcom + 1 + SAFETY];
  outbufr =new float[(nd+nq+2)*maxcom + 1 + SAFETY];
  /* Zero ins/outs arrays */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Mesh spacing */
  dx = lev[ilev]->dxbar[nd-1];
  /* Calculate processor boundaries on this level */
  levpbd[0] = xl[nd-1];
  for (npe = 1; npe <  limpe; npe++)
    levpbd[npe] = xl[nd-1] + 
      dx*(int)((pbd[npe] - xl[nd-1] + 0.5*dx)/dx);
  levpbd[limpe] = xr[nd-1];
  /* Find left processor */
  nproc = -1;
  xnd = lev[ilev]->xl[nd-1] - 0.5*dx;
  /* Check wrap-around (periodic) */
  if ((xnd < xl[nd-1]) && (bcl[nd-1] == bcper)){
    xnd = xr[nd-1] - 0.5*lev[ilev]->dxbar[nd-1];
    for (npe = mype + 1; npe < limpe; npe++){
      if ((xnd > levpbd[npe]) && (xnd < levpbd[npe+1])){
	nproc = npe;
	break;
      }
    }
  }
  else{
    /* No wrap-around */
    for (npe = mype - 1; npe >= 0; npe--){
      if ((xnd > levpbd[npe]) && (xnd < levpbd[npe+1])){
	nproc = npe;
	break;
      }
    }
  }
  if (lncell > 0){
    /* Set information to left processor */
    outs[nproc] = (nd + nq + 2)*lncell + 1;
    /* No of such cells */
    outbufr[0] = (float)lncell;
    /* Put cell information in outbufr*/
    put_ref(ilev, plev[ilev]->lref, outbufr);
  }
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Receive boundary cell information from right processor */
  info = false;
  for (npe = 0; npe < limpe; npe++){
    if (ins[npe] > 0){
      info = true;
      break;
    }
  }
  if (info)
  /* Receive information */
    receive_ref(ilev, inbufr);
  /* Zero ins/outs arrays */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Cell nd position */
  xnd =  lev[ilev]->xr[nd-1] + 0.5*dx;
  /* Find right processor */
  nproc = -1;
  /* Check wrap-around (periodic) */
  if ((xnd > xr[nd-1]) && (bcl[nd-1] == bcper)){
    xnd = xl[nd-1] + 0.5*lev[ilev]->dxbar[nd-1];
    for (npe = mype - 1; npe >=0; npe--){
      if ((xnd > levpbd[npe]) && (xnd < levpbd[npe+1])){
	nproc = npe;
	break;
      }
    }
  }
  else{
    /* No wrap-around */
    for (npe = mype + 1; npe < limpe; npe++){
      if ((xnd > levpbd[npe]) && (xnd < levpbd[npe+1])){
	nproc = npe;
	break;
      }
    }
  }
  if (rncell > 0){
    /* Set information to right processor */
    outs[nproc] = (nd + nq + 2)*rncell + 1;
    /* No of such cells */
    outbufr[0] = (float)rncell;
    /* Put cell information in outbufr*/
    put_ref(ilev, plev[ilev]->rref, outbufr);
  }
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Receive boundary cell information from left processor */
  info = false;
  for (npe = 0; npe < limpe; npe++){
    if (ins[npe] > 0){
      info = true;
      break;
    }
  }
  if (info)
  /* Receive information */
    receive_ref(ilev, inbufr);
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  return;
}
//end ref_proc_bound.c

//file put_ref.c
/*
+ put_ref.c - Fill outbufr with refine information.
     Description :
        Fills  outbufr with  information about  cells that  are  to be
        refined on  the neighbouring processor using  the list pointed
        to by  start. It removes these cells and the list.
*/
void put_ref(int ilev, struct cell_list *clist, float outbufr[])
{
  int count, id, i;
  struct cell *cptr;
  struct cell_list *tclist;
  count = 0;
  while (clist != NULL){
    cptr = clist->cptr;
    /* Child cell number */
    count++;
    outbufr[count] = (float)cptr->ichild;
    /* Cell volume */
    count++;
    outbufr[count] = cptr->vol;
    /* Coordinates */
    for (id = 1; id <= nd; id++){
      count++;
      outbufr[count] = cptr->xc[id-1];
    }
    /* Solution */
    for (i = 0; i < nq; i++){
      count++;
      outbufr[count] = cptr->qa[i];
    }
    /* Remove cell */
    rem_cell(ilev, cptr);
    /* Remove cell_list */
    tclist = clist;
    clist = clist->nxt;
    delete tclist;
  }
  return;
}
//end put_ref.c

//file receive_ref.c
/*
+ put_receive.c - Refine cells from neighbour.
     Description :
        Uses  the information from  inbufr to  create cells  that were
        generated by ref_cell on  a neighbour processor, but belong on
        this processor.
*/
void receive_ref(int ilev, float inbufr[])
{
  int icell, ncell, count, id, i;
  struct cell *cptr;
  struct cell_list *start, *clist, *tclist;
  ncell = (int)inbufr[0];
  count = 0;
  start = NULL;
  for (icell = 1; icell <= ncell; icell++){
    cptr = new_cell();
    clist = new cell_list;
    clist->cptr = cptr;
    clist->nxt = start;
    start = clist;
    /* Child cell number */
    count++;
    cptr->ichild = (int)inbufr[count];
    /* Cell volume */
    count++;
    cptr->vol = inbufr[count];
    /* Coordinates */
    for (id = 1; id <= nd; id++){
      count++;
      cptr->xc[id-1] = inbufr[count];
    }
    /* Solution */
    for (i = 0; i < nq; i++){
      count++;
      cptr->qa[i] = inbufr[count];
    }
    /* Connect to parent */
    connect_par(ilev, cptr);
    /* Update cell list */
    if (lev[ilev]->lcell != NULL){
      /* Cells already exist on this level */
      lev[ilev]->lcell->nxt = cptr;
      cptr->prev = lev[ilev]->lcell;
      lev[ilev]->lcell = cptr;
    }
    else{
      /* First cells on this level */
      lev[ilev]->fcell= cptr;
      lev[ilev]->lcell= cptr;
      cptr->prev = NULL;
    }
    cptr->nxt = NULL;
  }
  /* Connect cells and remove list */
  clist = start;
  while (clist != NULL){
    cptr = clist->cptr;
    for (id = 1; id <= nd; id++)
      connect_cell(ilev, id, lev[ilev]->fjoin[id-1], 
		   lev[ilev]->ljoin[id-1], cptr);
    tclist = clist;
    clist = clist->nxt;
    delete tclist;
  }
  return;
}
//end receive_ref.c

//file deref_proc_bound.c
/*
+ deref_proc_bound - Derefine processor boundary cells.
     Description :
        Deletes boundary  cells on  level ilev whose  parent is  to be
        derefined. This  is to  deal with fine/coarse  boundaries that
        coincide with  a processor  boundary.  Looks for  parent cells
        with qsol=deref. These might be halo cells, so com_halo has to
        set  qsol  before  this   is  called.  If  the  cells  contain
        particles,  then they  are added  to  the list  pointed to  by
        start.
*/
void deref_proc_bound(int ilev, struct part_list *&start)
{
  struct cell *cptr;
  struct cell_list *clist;
  struct part *pptr, *npptr;
  struct part_list *plist;
  /* Remove left boundary cells */
  clist = plev[ilev]->lbound;
  while (clist != NULL){
    cptr = clist->cptr;
    if (cptr->par != NULL){
      if (cptr->par->qsol == deref){
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
	/* Remove cell */
	rem_cell(ilev, cptr);
      }
    }
    clist = clist->nxt;
  }
  /* Remove right boundary cells */
  clist = plev[ilev]->rbound;
  while (clist != NULL){
    cptr = clist->cptr;
    if (cptr->par != NULL){
      if (cptr->par->qsol == deref){
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
	/* Remove cell */
	rem_cell(ilev, cptr);
      }
    }
    clist = clist->nxt;
  }
  return;
}
//end deref_proc_bound.c

//file com_cor.c
/*
+ file com_cor -  Communicate flux correction.
     Description :
        Communicates the flux correction,  cor, from the halo joins to
        the boundary cell joins.
*/
void com_cor(int ilev)
{
  int count, npe, ios;
  float *inbufr, *outbufr;
  struct cell *cptr;
  (void)cptr;
  struct cell_list *clist;
  (void)clist;
  /* Create communications arrays */
  inbufr = new float[(nq + 1)*plev[ilev]->maxcom + 1 + SAFETY];
  outbufr = new float[(nq + 1)*plev[ilev]->maxcom + 1 + SAFETY];
  /* Zero ins/outs arrays */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  if ((plev[ilev]->lsproc >= 0) && (lev[ilev+1]->fcell != NULL)){
    /* Put information for left processor in outbufr */
    if (lev[ilev]->fcell == NULL)
      /* No cells on this level */
      outbufr[0] = -1.0;
    else
      /* Cells on this level */
      outbufr[0] = 1.0;
    count = 0;
    put_cor(ilev, true, count, outbufr);
    outs[plev[ilev]->lsproc] = count + 1;
  }
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Receive boundary cell information from right processor */
  count = -1;
  for (npe = 0; npe < limpe; npe++){
    if (ins[npe] > 0)
      receive_cor(ilev, false, count, inbufr);
  }
  /* Zero ins/outs arrays */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Zero ins/outs arrays */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  if ((plev[ilev]->rsproc >= 0) && (lev[ilev+1]->fcell != NULL)){
    /* Put information for right processor in outbufr */
    if (lev[ilev]->fcell == NULL)
      /* No cells on this level */
      outbufr[0] = -1.0;
    else
      /* Cells on this level */
      outbufr[0] = 1.0;
    count = 0;
    put_cor(ilev, false, count, outbufr);
    outs[plev[ilev]->rsproc] = count + 1;
  }
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Receive boundary cell information from left processor */
  count = -1;
  for (npe = 0; npe < limpe; npe++){
    if (ins[npe] > 0)
      receive_cor(ilev, true, count, inbufr);
  }
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  return;
}
//end com_cor.c

//file put_cor.c
/*
+ file put_cor -  Put flux correction in outbufr.
     Description :
        Puts the flux  correction etc from the halo  joins in outbufr.
        Side is specified by left.
*/
void put_cor(int ilev, bool left, int &count, float outbufr[])
{
  int id, i;
  struct cell_list *clist;
  struct cell *cptr;
  struct join *jptr;
  if (left)
    /* Left halo */
    clist = plev[ilev]->lhalo;
  else
    /* Right halo */
    clist = plev[ilev]->rhalo; 
  while (clist != NULL){
    cptr = clist->cptr;
    if (left)
      /* Left halo - do right nd joins */
      jptr = cptr->rjoin[nd-1];
    else
      /* Right halo - do left nd joins */
      jptr = cptr->ljoin[nd-1];
    if (jptr->bound == cibound){
      /* Coarse/fine boundary join */
      count++;
      outbufr[count] = (float)jptr->done;
      for (i = 0; i < nq; i++){
	count++;
	outbufr[count] = jptr->cor[i];
      }
    }
    /* Do joins in nd plane */
    for (id = 1; id < nd; id++){
      jptr = cptr->ljoin[id-1];
      if (jptr->bound == cibound){
	/* Coarse/fine boundary join */
	for (i = 0; i < nq; i++){
	  count++;
	  outbufr[count] = jptr->cor[i];
	}
      }
    }
    clist = clist->nxt;
  }
  return;
}
//end put_cor.c

//file receive_cor.c
/*
+ file receive_cor -  Receive cor.
     Description :
        Receives the flux correction, cor in inbufr and puts it in the
        boundary cell joins. Side is specified by left.
*/
void receive_cor(int ilev, bool left, int &count, float inbufr[])
{
  int id, i, done;
  bool cells;
  struct cell *cptr;
  struct cell_list *clist;
  struct join *jptr;
  if (left)
    /* Left boundary */
    clist = plev[ilev]->lbound;
  else
    /* Right boundary */
    clist = plev[ilev]->rbound;
  /* Whether sending processor has cells on this level */
  count++;
  cells = (inbufr[count] > 0);
  while (clist != NULL){
    cptr = clist->cptr;
    if (left)
      /* Left boundary - do left nd joins */
      jptr = cptr->ljoin[nd-1];
    else
      /* Right boundary - do right nd joins */
      jptr = cptr->rjoin[nd-1];
    if (jptr->bound == cibound){
      /* Coarse/fine boundary join */
      count++;
      done = (int)inbufr[count];
      if (done > jptr->done){
	if (cells){
	  jptr->done = done;
	  for (i = 0; i < nq; i++){
	    count++;
	    jptr->cor[i] = inbufr[count];
	  }
	}
	else{
	  jptr->done += done;
	  for (i = 0; i < nq; i++){
	    count++;
	    jptr->cor[i] += inbufr[count];
	  }
	}
      }
      else{
	for (i = 0; i < nq; i++)
	  count++;
      }
    }
    /* Do joins in nd plane */
    for (id = 1; id < nd; id++){
      jptr = cptr->ljoin[id-1];
      if (jptr->bound == cibound){
	/* Coarse/fine boundary join */
	for (i = 0; i < nq; i++){
	  count++;
	  jptr->cor[i] += inbufr[count];
	}
      }
    }
    clist = clist->nxt;
  }
  return;
}
//end receive_cor.c

//file re_part.c
/*
+ re_part - Repartition solution.
     Description :
        This  uses cal_part to  get the  new processor  boundaries and
        then uses these to repartition the grid between processors.
*/
void re_part()
{
  int  id, i, icell, ilev, ncell, npe, maxcom, ios, count, nouts,
    mncell[MAXLEV][MAXPE], nvalues[MAXPE], *touts;
  float *inbufr, *outbufr, ocpu, ncpu;
  struct cell *cptr, *tcptr=0;
  struct cell_list *clist, *tclist,*start[MAXLEV][MAXPE];
  struct part *pptr;
  struct part_list *pmov, *plist;
  /* Set timer */
  parc_clock(&ocpu);
  /* Clear halo and boundary cells */
  for (ilev = 0; ilev < nlevs; ilev++)
    clear_halo(ilev);
  /* Get work array */
  work();
  /* Calculate new processor boundaries */
  cal_part();
  /* Zero cell counts, set start to NULL */
  for (ilev = 0; ilev < nlevs; ilev++){
    for (npe = 0; npe < limpe; npe++){
      mncell[ilev][npe] = 0;
      start[ilev][npe] = NULL;
    }
  }
  /* Find cells to be moved */
  for (ilev = 0; ilev < nlevs; ilev++){
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      /* Find destination processor */
      for (npe = 0; npe < limpe; npe++){
	if ((cptr->xc[nd-1] >= pbd[npe])
	  && (cptr->xc[nd-1] < pbd[npe+1]))
	  break;
      }
      if (npe != mype){
	/* Cell is to be moved - increment mncell */
	mncell[ilev][npe]++;
	/* Add to removal list */
	clist = new cell_list;
	clist->nxt = start[ilev][npe];
	start[ilev][npe] = clist;
	clist->cptr = cptr;
      }
      cptr = cptr->nxt;
    }
  }
  /* No of cells to be moved from this processor */
  ncell = 0;
  for (npe = 0; npe < limpe; npe++){
    for (ilev = 0; ilev < nlevs; ilev++)
      ncell += mncell[ilev][npe];
  }
  /* Put outs for this processor in touts */
  touts = new int[limpe*limpe];
  maxcom = 0;
  for (npe = 0; npe < limpe; npe++){
    nvalues[npe] = limpe;
    touts[npe] = (limpe - 1)*nlevs;
    for (ilev = 0; ilev < nlevs; ilev++)
      touts[npe] += (nq + 3)*mncell[ilev][npe];
    maxcom += touts[npe];
  }
  /* Global maxcom for outs */
  maxcom = (int)global_max((float)maxcom);
  /* Merge touts array */
  merge_array(touts, nvalues);
  /* Master finds maximum maxcom (outs/ins) */
  if (mype == 0){
    for (npe = 0; npe < limpe; npe++){
      count = 0;
      for (i = npe; i < limpe*limpe; i += limpe)
	count += touts[i];
      maxcom = max(maxcom, count);
    }
  }
  delete [] touts;
  broadcast(&maxcom);
  /* Zero ins/outs */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Set up communications arrays */
  inbufr = new float[maxcom + SAFETY];
  outbufr = new float[maxcom + SAFETY];
  /* Fill outbuf */
  count = -1;
  for (npe = 0; npe < limpe; npe++){
    nouts = count;
    if (npe != mype){
      for (ilev = 0; ilev < nlevs; ilev++){
	/* Number of cells moved on this level */
	count++;
	outbufr[count] = (float)mncell[ilev][npe];
	/* Cell data on this level */
	clist = start[ilev][npe];
	while (clist != NULL){
	  cptr = clist->cptr;
	  /* Cell coordinates */
	  for (id = 1; id <= nd; id++){
	    count++;
	    outbufr[count] = cptr->xc[id-1];
	  }
	  /* Solution */
	  for (i = 0; i < nq; i++){
	    count++;
	    outbufr[count] = cptr->qa[i];
	  }
	  clist = clist->nxt;
	}
      }
      outs[npe] = count - nouts;
    }
  }
  /* Remove cells, removal list and make particle list if necessary */
  pmov = NULL;
  for (ilev = 0; ilev < nlevs; ilev ++){
    for (npe = 0; npe < limpe; npe++){
      clist = start[ilev][npe];
      while (clist != NULL){
	cptr = clist->cptr;
	if ((particle) && (cptr->pptr != NULL)){
	  /* Add particles to list */
	  pptr = cptr->pptr;
	  while (pptr != NULL){
	    plist = new part_list;
	    plist->pptr = pptr;
	    plist->pproc = npe;
	    plist->nxt = pmov;
	    pmov = plist;
	    pptr = pptr->nxt;
	  }
	}
	rem_cell(ilev, cptr);
	tclist = clist;
	clist = clist->nxt;
	delete tclist;
      }
    }
  }
  /* Communicate cell information */
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Receive cell information */
  count = -1;
  for (npe = 0; npe < limpe; npe++){
    /* Check for information from processor npe */
    if (npe != mype){
      /* For each level */
      for (ilev = 0; ilev < nlevs; ilev++){
	count++;
	ncell= (int)inbufr[count];
	if (ncell != 0){
	  /* Create and fill cells on this level */
	  cptr = new_cell();
	  /* Check existing list */
	  if (lev[ilev]->lcell == NULL){
	    /* First cells on this level */
	    lev[ilev]->fcell = cptr;
	    lev[ilev]->lcell = cptr;
	  }
	  else
	    /* Add to end of cell list */
	    lev[ilev]->lcell->nxt = cptr;
	  /* Fill new cells */
	  for (icell = 1 ; icell <= ncell; icell++){
	    /* Clean cell for safety */
	    clean_cell(cptr);
	    for (id = 1; id <= nd; id++){
	      /* Coordinates */
	      count++;
	      cptr->xc[id-1] = inbufr[count];
	    }
	    /* Cell volume */
	    cptr->vol = cell_vol(ilev, cptr->xc);
	    for (i = 0; i < nq; i++){
	      /* Solution */
	      count++;
	      cptr->qa[i] = inbufr[count];
	    }
	    /* Cell is filled */
	    cptr->empty = false;
	    cptr->nxt = new_cell();
	    tcptr = cptr;
	    cptr = cptr->nxt;
	  }
	  /* Terminate list */
	  del_cell(tcptr->nxt);
	  tcptr->nxt = NULL;
	  lev[ilev]->lcell = tcptr;
	}
      }
    }
  }
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  if (particle)
    /* Get particles from other processors */
    pstart = mov_part(false, pmov);
  else
    /* Particle not switched on */
    pstart = NULL;
  /* Reset processor boundaries */
  set_proc_bound();
  parc_clock(&ncpu);
  if ((debug) && (mype == 0)){
    cout.setf(ios::scientific);
    cout << "repart: that took " << setprecision(4)
	 << ncpu - ocpu << " secs\n";
  }
  /* Reset parc buffer size */
  parc_reset();
  /* Defrag */
  de_frag();
  return;
}
//end re_part.c

//file set_proc_bound.c
/*
+ set_proc_bound - Set processor boundaries.
     Description :
        This  sets the processor  boundaries for  all levels.  Must be
        called after repartitioning and when creating new grid levels.
*/
void set_proc_bound()
{
  int ilev;
  float dx;
  for (ilev = 0; ilev < nlevs; ilev++){
    /* Mesh spacing in direction nd */
    dx = lev[ilev]->dxbar[nd-1];
    /* Left boundary */
    lev[ilev]->xl[nd-1] = xl[nd-1] + 
      dx*(int)((pbd[mype] - xl[nd-1] + 0.5*dx)/dx);
    /* Right boundary */
    lev[ilev]->xr[nd-1] = xl[nd-1] + 
      dx*(int)((pbd[mype+1] - xl[nd-1] + 0.5*dx)/dx);
    /* Recalculate max no of cells in nd direction on each level */
    lev[ilev]->nc[nd-1] = (int)((lev[ilev]->xr[nd-1] 
      - lev[ilev]->xl[nd-1] + 0.1*dx)/dx);
  }
  return;
}
//end set_proc_bound.c

//file mov_part.c
/*
+ mov_part - Move/copy particles between processors.
     Description :
        If all=false, this moves the  particles in the list pointed to
        by pmov to the appropriate  processor. In that case pproc must
        be  set in  the particle  list element.  If all=true,  then it
        ignores  pmov,  creates  a  list  of  all  particles  on  this
        processor  and copies them  to the  master processor  list. It
        returns  the start  of  the  new list  of  particles for  this
        processor (NULL  if no  particles moved onto  this processor).
        Note that  it does not  assign particles to cells.   Called in
        re_part, re_grid, refine, plotp and step at each level.
*/
struct part_list *mov_part(bool all, struct part_list *pmov)
{
  int  id, i, ipart, npart, npe, ilev, maxcom, ios, count, nouts,
    nvalues[MAXPE], mnpart[MAXPE], *touts;
  float *inbufr, *outbufr;
  struct cell *cptr;
  (void)cptr;
  struct part *pptr;
  struct part_list *plist, *tplist, *start[MAXPE];
  /* Zero particle counts, set start to NULL */
  for (npe = 0; npe < limpe; npe++){
    mnpart[npe] = 0;
    start[npe] = NULL;
  }
  if (all)
    /* Ignore pmov -- make list of all particles on this processor */
    mk_partlist(plist, tplist); 
  else
    /* Use pmov */
    plist = pmov;
  /* Particle counts and destination processor */
  while (plist != NULL){
    pptr = plist->pptr;
    if (all)
      /* Copy all to master */
      npe = 0;
    else
      /* Else use pproc */
      npe = plist->pproc;
    if (npe < 0){
      cout << "bad destination processor for part at ";
      for (id = 1; id <= nd; id++)
	cout << " " << pptr->xp[id-1];
      cout << " from proc " << mype << " to " << npe << '\n';
    }
    /* Counts */
    mnpart[npe]++;
    /* Add to destination processor list */
    tplist = new part_list;
    tplist->nxt = start[npe];
    start[npe] = tplist;
    tplist->pptr = pptr;
    /* Delete original particle list */
    tplist = plist;
    plist = plist->nxt;
    delete tplist;
  }
  /* Put outs for this processor in touts */
  touts = new int[limpe*limpe];
  maxcom = 0;
  for (npe = 0; npe < limpe; npe++){
    nvalues[npe] = limpe;
    touts[npe] = limpe - 1 + (2*nd + 4 + MAXPROP)*mnpart[npe];
    maxcom += touts[npe];
  }
  /* Global maxcom for outs */
  maxcom = (int)global_max((float)maxcom);
  /* Merge touts array */
  merge_array(touts, nvalues);
  /* Master finds maximum maxcom (outs/ins) */
  if (mype == 0){
    for (npe = 0; npe < limpe; npe++){
      count = 0;
      for (i = npe; i < limpe*limpe; i += limpe)
	count += touts[i];
      maxcom = max(maxcom, count);
    }
  }
  broadcast(&maxcom);
  /* Zero ins/outs */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Set up communications arrays */
  inbufr = new float[maxcom + SAFETY];
  outbufr = new float[maxcom + SAFETY];
  /* Fill outbuf */
  count = -1;
  for (npe = 0; npe < limpe; npe++){
    nouts = count;
    if (npe != mype){
      /* Number of particles moved to this processor */
      count++;
      outbufr[count] = (float)mnpart[npe];
      /* Particle data */
      plist = start[npe];
      while (plist != NULL){
	pptr = plist->pptr;
	/* Particle coordinates and velocities */
	for (id = 1; id <= nd; id++){
	  count++;
	  outbufr[count] = pptr->xp[id-1];
	  count++;
	  outbufr[count] = pptr->vp[id-1];
	}
	/* Particle time, mass and thermal energy */
	count++;
	outbufr[count] = pptr->t;
	count++;
	outbufr[count] = pptr->mass;
	count++;
	outbufr[count] = pptr->en;
	/* Other particle properties */
	for (i = 0; i < MAXPROP; i++){
	  count++;
	  outbufr[count] = pptr->props[i];
	}
	if (all)
	  /* Particle level */
	  ilev = pptr->ilev;
	else{
	  /* Not all -- delete particle */
	  delete pptr;
	  /* Just set particle level to -1 */
	  ilev = -1;
	}
	/* Write particle level */ 
	count++;
	outbufr[count] = (float)ilev;
	/* Delete list element */
	tplist = plist;
	plist = plist->nxt;
	delete tplist;
      }
      outs[npe] = count - nouts;
    }
  }
  /* Communicate particle information */
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Receive particle information */
  count = -1;
  /* Initialise particle list */
  pmov = NULL;
  for (npe = 0; npe < limpe; npe++){
    /* Check for information from processor npe */
    if (npe != mype){
      /* Number of particles */
      count++;
      npart = (int)inbufr[count];
      if (npart != 0){
	for (ipart = 1 ; ipart <= npart; ipart++){
	  /* Create particles and add to list */
	  pptr = new part;
	  /* Set particle pointers to NULL (safety) */
	  pptr->prev = NULL;
	  pptr->nxt = NULL;
	  pptr->cptr = NULL;
	  plist = new part_list;
	  plist->pptr = pptr;
	  plist->nxt = pmov; 
	  pmov = plist;
	  /* Particle coordinates and velocities */
	  for (id = 1; id <= nd; id++){
	    count++;
	    pptr->xp[id-1] = inbufr[count]; 
	    count++;
	    pptr->vp[id-1] = inbufr[count]; 
	  }
	  /* Particle time, mass and thermal energy */
	  count++;
	  pptr->t = inbufr[count]; 
	  count++;
	  pptr->mass = inbufr[count]; 
	  count++;
	  pptr->en = inbufr[count]; 
	  /* Other particle properties */
	  for (i = 0; i < MAXPROP; i++){
	    count++;
	    pptr->props[i] = inbufr[count];
	  }
	  /* Particle level */
	  count++;
	  pptr->ilev = (int)inbufr[count]; 
	}
      }
    }
  }
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  delete [] touts;
  return(pmov);
}
//end mov_part.c

//file get_part_proc.c
/*
+ get_part_proc - Get processor for particles
     Description :
        Determines  the  processors  for  the particles  in  the  list
        pointed  to  by start  i.e.  sets  the  ->proc field  in  this
        list. It also sets up the communications buffers.
*/
void get_part_proc(struct part_list *start)
{
  int maxcom, ilev, npe, minproc, maxproc, id, count, ocount, ios, 
  npart, i;
  float *inbufr, *outbufr, xp[MAXD];
  struct cell *cptr;
  struct part *pptr;
  struct part_list *plist; 
  /* Count particles and get range of possible destination processors */
  minproc = MAXPE;
  maxproc = -MAXPE;
  npart = 0;
  plist = start;
  while (plist != NULL){
    pptr = plist->pptr;
    npart++;
    /* Find possible destination processors */
    for (ilev = 0; ilev < nlevs; ilev++){
      for (npe = 0; npe < limpe; npe++){
	if ((levpbd[ilev][npe] <= pptr->xp[nd-1]) &&
	    (pptr->xp[nd-1] < levpbd[ilev][npe+1])){
	  minproc = min(npe, minproc);
	  maxproc = max(npe, maxproc);
	}
      }
    }
    plist = plist->nxt;
  }
  /* Calculate buffer sizes */
  maxcom = (maxproc - minproc + 1);
  maxcom = (2*nd + 4 + MAXPROP)*npart*maxcom;
  maxcom = (int)global_max(float(maxcom));
  /* Set up communications arrays */
  inbufr = new float[maxcom + SAFETY];
  outbufr = new float[maxcom + SAFETY];
  for (npe = 0; npe < limpe; npe++){
    /* Zero ins/outs (safety) */
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Assemble information */
  count = -1;
  for (npe = minproc; npe <= maxproc; npe++){
    /* For each particle */
    plist = start;
    while (plist != NULL){
      pptr = plist->pptr;
      /* Set bad processor as default */
      plist->pproc = -1;
      if (npe != mype){
	/* Put particle coordinates in outbufr */
	for (id = 1; id <= nd; id++){
	  outs[npe]++;
	  count++;
	  outbufr[count] = pptr->xp[id-1];
	}
      }
      plist = plist->nxt;
    }
  }
  /* Exchange data */
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Set outs */
  for (npe = 0; npe < limpe; npe++)
    outs[npe] = ins[npe]/nd;
  count = -1;
  ocount = -1;
  for (npe = 0; npe < limpe; npe++){
    /* Extract particle coordinates */
    for (i = 0; i < outs[npe]; i++){ 
      for (id = 1; id <= nd; id++){
	count++;
	xp[id-1] = inbufr[count];
      }
      ocount++;
      /* Check if particle is on mype */
      cptr = find_fcell(ilev, xp);
      if ((cptr != NULL) && (!cptr->halo))
	/* Particle is on mype */
	outbufr[ocount] = 1.0;
      else
	/* Particle is not on mype */
	outbufr[ocount] = -1.0;
    }
  }
  /* Exchange data */
  ios = 0;
  all2all(&limpe, outs, outbufr, ins, inbufr, &ios);
  /* Extract data */
  count = -1;
  for (npe = 0; npe < limpe; npe++){
    npart = ins[npe];
    plist = start;
    for (i = 1; i <= npart; i++){
      count++;
      if ((plist->pproc < 1) && (inbufr[count] > 0.0))
	plist->pproc = npe;
      plist = plist->nxt;
    }
  }
  /* Check list */
  plist = start;
  while (plist != NULL){
    pptr = plist->pptr;
    if (plist->pproc < 0){
      cout << "getproc part at ";
      for (id = 1; id <= nd; id++)
	cout << " " << pptr->xp[id-1];
      cout << " From proc " << mype << " to " << plist->pproc << '\n';
    }
    plist = plist->nxt;
  }
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  return;
}
//end get_part_proc.c

//file broadcasti0.c
/*
+ broadcasti0 - Broadcast an integer.
     Description :
        The  master processor  sends  an integer  value  to the  other
        processors.  Converts to float and uses parc_broadcastr.
*/
void broadcast(int *value)
{
  int ndim=1, dims[1], ios;
  float buf[1];
  dims[0] = 1;
  buf[0] = (float)*value;
  ios = 0;
  parc_broadcastr(buf, &ndim, dims, &ios);
  *value = (int)buf[0];
  return;
}
//end broadcasti0.c

//file broadcasti1.c
/*
+ broadcasti1 - Broadcast an integer array.
     Description :
        The master processor sends  nvalues integers stored in the buf
        array  to the  other processors.  Converts to  float  and uses
        parc_broadcastr.
*/
void broadcast(int buf[], int nvalues)
{
  int i, ndim=1, dims[1], ios;
  float *rbuf;
  dims[0] = nvalues;
  ios = 0;
  rbuf = new float[nvalues];
  for (i = 0; i < nvalues; i++)
    rbuf[i] = (float)buf[i];
  parc_broadcastr(rbuf, &ndim, dims, &ios);
  for (i = 0; i < nvalues; i++)
    buf[i] = (int)rbuf[i];
  delete [] rbuf;
  return;
}
//end broadcasti1.c

//file broadcastr0.c
/*
+ broadcastr0 - Broadcast a float.
     Description :
        The  master  processor  sends  an  real  value  to  the  other
        processors.  Implemented using parc_broadcastr.
*/
void broadcast(float *value)
{
  int ndim=1, dims[1], ios;
  float buf[1];
  dims[0] = 1;
  buf[0] = *value;
  ios = 0;
  parc_broadcastr(buf, &ndim, dims, &ios);
  *value = buf[0];
  return;
}
//end broadcastr0.c

//file broadcastr1.c
/*
+ broadcastr1 - Broadcast a real array.
     Description :
        The master processor sends nvalues  stored in the buf array to
        the other processors.  Implemented using broadcastr.
*/
void broadcast(float buf[], int nvalues)
{
  int ndim=1, dims[1], ios;
  dims[0] = nvalues;
  ios = 0;
  parc_broadcastr(buf, &ndim, dims, &ios);
  return;
}
//end broadcastr1.c

//file broadcastb0.c
/*
+ broadcastb0 - Broadcast a boolean.
     Description :
        The  master  processor  sends  a boolean  value  to  the  other
        processors.  Implemented using broadcasti0.
*/
void broadcast(bool *value)
{
  float temp;
  if (*value)
    temp = 1.0;
  else
    temp = -1.0;
  broadcast(&temp);
  *value = temp > 0;
  return;
}
//end broadcastb0.c

//file broadcastc0.c
/*
+ broadcastc0 - Broadcast a string.
     Description :
        The  master processor sends  the string  value[] to  the other
        processors.  Implemented using broadcastc0.
*/
void broadcast(char buf[])
{
  int ios;
  ios = 0;
  parc_broadcast0c(buf, &ios);
  return;
}
//end broadcastc0.c

//file merge_arrayi.c
/*
+ merge_arrayi - Merge an integer array.
     Description :
        The master  processor receives nvalues integers  stored in the
        buf array from  each of the other processors.   It then merges
        these into the buf array  in processor order, taking the first
        nvalues from its own buf  array. For the other processors, the
        buf array is unchanged. Implemented by converting to float and
        using merge_arrayr.
*/
void merge_array(int buf[], int nvalues[])
{
  int i, ntot, npe;
  float *rbuf;
  /* Total number of values */
  ntot = 0;
  for (npe = 0; npe < limpe; npe++)
    ntot += nvalues[npe];
  /* Convert to float */
  rbuf = new float[ntot];
  for (i = 0; i < nvalues[mype]; i++)
    rbuf[i] = (float)buf[i];
  /* Call merge_arrayr */
  merge_array(rbuf, nvalues);
  /* Master puts values in buffer and converts back to int */
  if (mype == 0){
    for (i = 0; i < ntot; i++)
      buf[i] = (int)rbuf[i];
  }
  delete [] rbuf;
  return;
}
//end merge_arrayi.c

//file merge_arrayr.c
/*
+ merge_arrayr - Merge a real array.
     Description :
        The  master processor receives  nvalues[npe] values  stored in
        the buf array  from processor npe.  It then  merges these into
        the buf array in processor order, taking the nvalues[0] values
        from its  own buf  array.  For the  other processors,  the buf
        array is unchanged. Implemented using all2all.
*/
void merge_array(float buf[], int nvalues[])
{
  int i, ntot, maxcom, nval0, npe, ios;
  (void)maxcom;
  float *inbufr;
  /* Set ins/outs */
  for (npe = 0; npe < limpe; npe ++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Create output/input buffers */
  ntot = 0;
  for (npe = 1; npe < limpe; npe++)
    ntot += nvalues[npe];
  maxcom = ntot + nvalues[0];
  inbufr = new float[ntot+SAFETY];
  if (mype != 0)
    outs[0] = nvalues[mype];
  /* Exchange data */
  ios = 0;
  all2all(&limpe, outs, buf, ins, inbufr, &ios);
  if (mype == 0){
    /* Put values in buffer */
    nval0 = nvalues[0];
    for (i = 0; i < ntot; i++)
      buf[i+nval0] = inbufr[i];
  }
  /* Delete communications array */
  delete [] inbufr;
  return;
}
//end merge_arrayr.c

//file all2all.c
/*
+ all2all - All to all real data exchange.
     Description :
        Exchanges data between all processor. npe is the
        number  of  processors (limpe),  outs[npe]  is  the number  of
        values  to be  sent by  processor npe,  outbuf[] contains  the
        values in  processor order, ins[npe]  is the number  of values
        received from  processor npe, inbuf[] are  the received values
        in in processor order, ios is a fail integer.
*/
void all2all(int *npe, int outs[], float outbuf[], int ins[],
  float inbuf[], int *ios)
{
  float step_time;
  /* Increment wall time */
  parc_clock(&step_time);
  wall_time += step_time - ostep_time;
  parc_xdatar(npe, outs, outbuf, ins, inbuf, ios);
  /* Set time */
  parc_clock(&ostep_time);
  return;
}
//end all2all.c

//file global_max.c
/*
+ global_max - Calculate global maximum.
     Description :
        Given value from a processor,  it returns the maximum of value
        over all processors.
*/
float global_max(float value)
{
  int npe, nvalues[MAXPE];
  float buf[MAXPE];
  /* Calculate global maxcom */
  buf[0] = value;
  for (npe = 0; npe < limpe; npe++)
    nvalues[npe] = 1;
  merge_array(buf, nvalues);
  if (mype == 0){
    /* Master calculates global maximum */
    for (npe = 1; npe < limpe; npe++)
      value = max(buf[npe], value);
    /* Broadcast global maximum */
    broadcast(&value);
  }
  else{
    /* Syblings receive global maximum */
    broadcast(&value);
  }
  return(value);
}
// end global_max.c

//file global_sum.c
/*
+ global_sum - Calculate global sum.
     Description :
        Given  value from  a processor,  it returns  the sum  of these
        values over all processors.
*/
float global_sum(float value)
{
  int npe, nvalues[MAXPE];
  float buf[MAXPE], sum=0.0;
  buf[0] = value;
  for (npe = 0; npe < limpe; npe++)
    nvalues[npe] = 1;
  merge_array(buf, nvalues);
  if (mype == 0){
    /* Master calculates sum */
    sum = 0.0;
    for (npe = 0; npe < limpe; npe++)
      sum += buf[npe];
    /* Broadcast global sum */
    broadcast(&sum);
  }
  else{
    /* Syblings receive global sum */
    broadcast(&sum);
  }
  return(sum);
}
// end global_sum.c
