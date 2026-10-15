#undef MAIN
#include <iostream>
#include <iomanip>
#include <math.h>
#include <new>
#include <cstdlib>
#include <cstring>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "user.h"
#include "funcs.h"
#define mcom1 1
#define mcom2 2
#define mcom3 3
#define MINCR 500
struct mcell_list {
  struct mcell *mcptr;
  struct mcell_list *nxt;
};
struct mcell *fcell[MAXLEV];
struct mcell_list  *lhalo[MAXLEV], *rhalo[MAXLEV], *lbound[MAXLEV],
  *rbound[MAXLEV], *fcbound[MAXLEV], *ghost[MAXLEV];    
int comsize, mugtop, nebd[MAXLEV];
long long int work;
float *minbufr, *moutbufr;
int mode;

extern void usrpq(float qa[], float pa[]);
extern void usrqp(float pa[], float qa[]);
extern void usrfxq(struct cell *cptr);
extern void usrbd(int id, int bound, join *jptr, float pl[],
  float pr[]);
extern void usrsaf(float dt, struct cell *cptyr, float so[]);
extern void usrnu(float pa[], float nu[]);
extern void usrref(cell* cptr, float qref[]);
extern void usrpref(mcell* mcptr, float pref[]);
extern float cell_vol(int ilev, float xc[]);
extern int dpdx_index(int jid, int id);
extern int it2_to_i(int it[]);
extern float cell_vol(int ilev, float xc[]);
extern float global_max(float value);
extern float global_sum(float value);
extern void broadcast(int buf[], int nvalues);
extern void broadcast(float buf[], int nvalues);
extern void merge_array(int buf[], int nvalues[]);
void merge_array(float buf[], int nvalues[]);
void do_mug(int tmode, int flev, float dt);
void cal_grav();
void create_mug();
struct mcell *new_mcell();
void clean_mcell(struct mcell *mcptr);
void clean_mug();
void clear_mug();
void connect_mcell(int ilev, struct mcell *mcptr);
void mk_mclist(struct cell_list *clist, struct mcell_list
  *&start);
void cal_matrix_cof(int ilev, float dt);
void init_imp_mug(int flev, float dt);
void init_grav_mug();
void cal_bound_pot();
float sum_pot(float xc[]);
void relax(int ilev, int nstep);
void map_up_mcell(struct mcell *pmcptr);
void map_down_mcell(int ilev, struct mcell *pmcptr);
bool cal_mdiffl(int id, struct mcell *mcptr, float diff[]);
bool cal_mdiffr(int id, struct mcell *mcptr, float diff[]);
void set_mug_source(int ilev);
void correct_mug(int ilev);
float maxpot();
float norm(int ilev);
float trunc(int ilev);
void mcom_halo(int mcom, int ilev);
void put_mcom(int mcom, int ilev, bool left);
void put_mcom1(int ilev, int &count,
  struct mcell_list *mclist);
void put_mcom2(int ilev, bool left, int &count,
  struct mcell_list *mclist);
void put_mcom3(int ilev, int &count,
  struct mcell_list *mclist);
void receive_mcom(int mcom, int ilev, bool left);
void receive_mcom1(int ilev, bool left);
void receive_mcom2(int ilev, bool left);
void receive_mcom3(int ilev, bool left);
//file do_mug.c
/*
  +  do_mug -  do implicit diffusion or gravity  with multigrid.
  Description :     
  Uses  a multigrid  to  solve an  implicit  diffusion equation  or
  poisson  equation  on level  flev  (mode=0  - gravity,  mode=1  -
  diffusion). Uses Gauss-Seidel for relaxation and levels 0 .. flev
  for the multigrid.   The solution is in  special mini-cells (type
  mcell).  It assumes that the mcell hierarchy has been created.
*/
void do_mug(int tmode, int flev, float dt)
{
  //char key[SZNAM];  
  int id, ncm, ilev, i, nmug;
  const int clocks = CLOCKS_PER_SEC;
  float tnorm, pot, ttrunc, twork, ocpu;
  struct cell *cptr;
  struct mcell *mcptr;
  bool done, up;
  mode = tmode;
  if (nomug)
    /* Create multigrid */
    create_mug();
  /* Initialise */
  ttrunc = 0.0;
  pot = 0.0;
  ocpu = (float)clock()/(float)clocks;
  initdcpu += mabs((float)clock()/(float)clocks) - mabs(ocpu);
  work = 0;
  nmug = 1;
  if (mode == 0){
    /* Gravity */
    init_grav_mug();
    if (nlevs == 1){
      /* Single grid is special */
      ncm = 0;
      for (id = 1; id <= nd; id++)
	ncm = max(ncm, nc[id-1]);
      /* Just relax */
      relax(0, 2*ncm);
      tnorm = norm(0);
      if (mugprt){
	twork = global_max(work);
	if (mype == 0){
	  cout.setf(ios::scientific); 
	  cout << "lev " << 0 << " residual " << setprecision(4) <<tnorm
	       << " work " << (float)twork <<'\n';
	}
      }
    }
    else{
      /* Initialise mug cycle */
      ilev = alev;
      done = false;
      up = true;
      while (!done){
	/* Relax */
	relax(ilev, nrel);
	if (mugprt){
	  tnorm = norm(ilev);
	  twork = (float)global_max((float)work);
	  if (mype == 0){
	    cout.setf(ios::scientific);
	    cout << setprecision(4);
	    cout << "lev " << ilev << " residual " << tnorm <<
	      " work " << twork << '\n';
	  }
	}
	if (up){
	  if (ilev > 0){
	    ilev--;
	    /* Set multigrid source terms */
	    set_mug_source(ilev);
	  }
	  else{
	    /* Correct */
	    correct_mug(ilev);
	    ilev++;
	    up = false;
	  }
	}
	else{
	  if (ilev != alev){
	    /* Correct */
	    correct_mug(ilev);
	    ilev++;
	  }
	  else{
	    /* Check for convergence */
	    pot = maxpot();
	    if ((norm(alev) < max(normrat*pot, trunc(alev))) ||
		(nmug > maxmug))
	      done = true;
	    else{
	      nmug++;
	      up = true;
	      ilev--;

	      /* Set multigrid source terms */
	      set_mug_source(ilev);
	    }
	  }
	}
      }
    }
    if (mugprt){
      tnorm = norm(alev);
      if (nlevs > 1)
	ttrunc = trunc(alev);
      if (mype == 0){
	if (nmug < maxmug)
	  cout << "mug converged nmug " << nmug << " ";
	else
	  cout << "mug not converged nmug " << nmug << " ";;	  
	if (nlevs == 1)
	  cout << "lev " << nlevs - 1 << " residual " << tnorm
	       << '\n';
	else
	  cout << "lev " << alev << " residual " << tnorm
	       << " trunc " << ttrunc << " maxpot " << pot << '\n';
      }
    }
    /* Put solution in full cell qa[nq-1] */
    for (ilev = alev; ilev >= 0; ilev--){
      if (limpe > 1){
	if (ilev != alev)
	  /* Communicate ilev+1 halo qa (mcom1) */
	  mcom_halo(mcom1, ilev + 1);
      }
      mcptr = fcell[ilev];
      while (mcptr != NULL){
	if (ilev != alev)
	  map_up_mcell(mcptr);
	cptr = mcptr->cptr;
	cptr->qa[iqpot] = mcptr->qa[iqpot];
	mcptr = mcptr->nxt;
      }
    }
  }
  else {
    /* Diffusion */
    init_imp_mug(flev, dt);    
    if (flev == 0){
      /* Coarsest grid is special */
      ncm = 0;
      for (id = 1; id <= nd; id++)
	ncm = max(ncm, nc[id-1]);
      /* Just relax */
      relax(0, 2*ncm);
      tnorm = norm(0);
      if (mugprt){
	twork = global_max(work);
	if (mype == 0){
	  cout.setf(ios::scientific); 
	  cout << "lev " << 0 << " residual " << setprecision(4) <<tnorm
	       << " work " << (float)twork <<'\n';
	}
      }
    }
    else{
      /* Initialise mug cycle */
      ilev = flev;
      done = false;
      up = true;
      while (!done){
	/* Relax */
	relax(ilev, nrel);
	if (mugprt){
	  tnorm = norm(ilev);
	  twork = (float)global_max((float)work);
	  if (mype == 0){
	    cout.setf(ios::scientific);
	    cout << setprecision(4);
	    cout << "lev " << ilev << " residual " << tnorm <<
	      " work " << twork << '\n';
	  }
	}
	if (up){
	  if (ilev > 0){
	    ilev--;
	    /* Set multigrid source terms */
	    set_mug_source(ilev);
	  }
	  else{
	    /* Correct */
	    correct_mug(ilev);
	    ilev++;
	    up = false;
	  }
	}
	else{
	  if (ilev != flev){
	    /* Correct */
	    correct_mug(ilev);
	    ilev++;
	  }
	  else{
	    /* Check for convergence */
	    ttrunc = trunc(flev);
	    tnorm = norm(flev);
	    if ((tnorm <  normrat*ttrunc) || (nmug > maxmug))
	      done = true;
	    else{
	      nmug++;
	      up = true;
	      ilev--;
	      /* Set multigrid source terms */
	      set_mug_source(ilev);
	    }
	  }
	}
      }
    }
    if (mugprt){
      if (mype == 0)
	//      tnorm = norm(flev);
      if (nlevs > 1)
	//	ttrunc = trunc(flev);
      if (mype == 0){
	if (nmug < maxmug)
	  cout << "mug converged nmug " << nmug << " ";
	else
	  cout << "mug not converged nmug " << nmug << " ";
	if (nlevs == 1)
	  cout << "lev " << nlevs - 1 << " residual " << tnorm
	       << '\n';
	else
	  cout << "lev " << flev << " residual " << tnorm
	       << " trunc " << ttrunc << '\n';
      }
    }
    /* Put solution change in flev full cell */
    mcptr = fcell[flev];
    while (mcptr != NULL){
      cptr = mcptr->cptr;
      /* Convert to conserved variables */
      usrqp(mcptr->qa, mcptr->qa);
      for (i = 0; i < nq; i++)
	cptr->dqdt[i] += mcptr->qa[i] - cptr->qa[i];		// TODO only do for cells that aren't solid
      mcptr = mcptr->nxt;
    }
  }
  return;
}
//end do_mug.c


//file cal_grav.c
/*
+  cal_grav -  calculate gravity with multigrid.
     Description :
     Uses  a multigrid  to solve  the Poisson  equation on  the finest
     level. Uses Gauss-Seidel for relaxation  and levels 0 .. alev for
     the  multigrid.   The solution  is  in  special mini-cells  (type
     mcell). It assumes that the mcell hierarchy has been created.
*/
void cal_grav()
{
  int id, ncm, ilev;
  const int clocks = CLOCKS_PER_SEC;
  float tnorm, pot, ttrunc, twork, ocpu;
  struct cell *cptr;
  struct mcell *mcptr;
  bool done, up;
  /* Initialise */
  ttrunc = 0.0;
  pot = 0.0;
  ocpu = (float)clock()/(float)clocks;
  init_grav_mug();
  initdcpu += mabs((float)clock()/(float)clocks) - mabs(ocpu);
  work = 0;
  if (nlevs == 1){
    /* Single grid is special */
    ncm = 0;
    for (id = 1; id <= nd; id++)
      ncm = max(ncm, nc[id-1]);
    /* Just relax */
    relax(0, 2*ncm);
    tnorm = norm(0);
    if (mugprt){
      twork = global_max(work);
      if (mype == 0){
	cout.setf(ios::scientific); 
	cout << "lev " << 0 << " residual " << setprecision(4) <<tnorm
	     << " work " << (float)twork <<'\n';
      }
    }
  }
  else{
    /* Initialise mug cycle */
    ilev = alev;
    done = false;
    up = true;
    while (!done){
      /* Relax */
      relax(ilev, nrel);
      if (mugprt){
	tnorm = norm(ilev);
	twork = (float)global_max((float)work);
	if (mype == 0){
	  cout.setf(ios::scientific);
	  cout << setprecision(4);
	  cout << "lev " << ilev << " residual " << tnorm <<
	    " work " << twork << '\n';
	}
      }
      if (up){
	if (ilev > 0){
	  ilev--;
	  /* Set multigrid source terms */
	  set_mug_source(ilev);
	}
	else{
	  /* Correct */
	  correct_mug(ilev);
	  ilev++;
	  up = false;
	}
      }
      else{
	if (ilev != alev){
	  /* Correct */
	  correct_mug(ilev);
	  ilev++;
	}
	else{
	  /* Check for convergence */
	  pot = maxpot();
	  if (norm(alev) <= max(normrat*pot, trunc(alev)))
	    done = true;
	  else{
	    up = true;
	    ilev--;
	    /* Set multigrid source terms */
	    set_mug_source(ilev);
	  }
	}
      }
    }
  }
  if (mugprt){
    if (mype == 0)
      cout << "mug converged ";
    tnorm = norm(alev);
    if (nlevs > 1)
      ttrunc = trunc(alev);
    if (mype == 0){
      if (nlevs == 1)
	cout << "lev " << nlevs - 1 << " residual " << tnorm
	     << '\n';
      else
	cout << "lev " << alev << " residual " << tnorm
	     << " trunc " << ttrunc << " maxpot " << pot << '\n';
    }
  }
  /* Put solution in full cell qa[nq-1] */
  for (ilev = alev; ilev >= 0; ilev--){
    if (limpe > 1){
      if (ilev != alev)
	/* Communicate ilev+1 halo qa (mcom1) */
	mcom_halo(mcom1, ilev + 1);
    }
    mcptr = fcell[ilev];
    while (mcptr != NULL){
      if (ilev != alev)
	map_up_mcell(mcptr);
      cptr = mcptr->cptr;
      cptr->qa[nq-1] = mcptr->qa[iqpot];
      mcptr = mcptr->nxt;
    }
  }
  return;
}
//end cal_grav.c






//file create_mug.c
/*
+  create_mug -  create grid for multigrid solution.
     Description :
     Creates  and  connects  the  mini-cells for  the  multigrid.   It
     assumes  that any previous  multigrid has  been cleared.  It also
     sets up the communications array, minbufr and moutbufr.  */
void create_mug()
{
  int ilev, icell, id, nbound;
  //float pa[MAXQ], so[MAXQ];
  struct cell *cptr;
  struct cell_list *clist;
  //struct join *jptr;
  struct mcell *mcptr;
  struct mcell_list *mclist, *start;
  /* Clean multigrid (safety) */
  clean_mug();
  /* Determine comsize (size of communications arrays) */
  comsize = 0;
  for (ilev = 0; ilev < nlevs; ilev++){
    /* Count left boundary cells */
    nbound = 0;
    clist = plev[ilev]->lbound;
    while (clist != NULL){
      nbound++;
      clist = clist->nxt;
    }
    comsize = max(comsize, nbound*plev[ilev]->nlproc);
    nbound = 0;
    /* Count right boundary cells */
    clist = plev[ilev]->rbound;
    while (clist != NULL){
      nbound++;
      clist = clist->nxt;
    }
    comsize = max(comsize, nbound*plev[ilev]->nrproc);
  }
  /* Create input/output buffers */
  comsize = nq*(int)global_max((float)comsize);
  minbufr = new float[comsize+SAFETY];
  moutbufr = new float[comsize+SAFETY];
  /* Create/update cell array */
  mugtop = celltop + MINCR;
  if (mcarr == NULL){
    /* Create */
    mcarr = new mcell *[mugtop];
    mcarrtop = mugtop;
    /* Create cells */
    for (icell = 0; icell < mcarrtop ; icell++)
      mcarr[icell] = new mcell;
  }
  else if (mcarrtop < mugtop){
    /* Increase size */
    mcarr = (struct mcell **)realloc((void *)mcarr, 
	      mugtop*sizeof(mcell *));
    /* Create additional cells */
    for (icell = mcarrtop; icell < mugtop ; icell++)
      mcarr[icell] = new mcell;
    mcarrtop = mugtop;
  }
  /* Clean minicells */
  for (icell = 0; icell < mcarrtop; icell++)
    clean_mcell(mcarr[icell]);
  /* Create mini-cell lists on each level */
  for (ilev = 0; ilev < nlevs; ilev++){
    /* Interior cells */
    fcell[ilev] = NULL;
    /* Run through full cells in reverse order */
    cptr = lev[ilev]->lcell;
    while (cptr != NULL){
      /* Add mini-cell to list */
      mcptr = new_mcell();
      mcptr->nxt = fcell[ilev];
      fcell[ilev] = mcptr;
      /* Connect mini-cell to full cell */
      cptr->mcptr = mcptr;
      mcptr->cptr = cptr;
	// TODO set solid flag

      /* Not ghost cell */
      mcptr->ghost = false;
      /* Volume of mini-cell */
      mcptr->vol = cptr->vol;
      /* Coordinates of mini-cell */
      for (id = 1; id <= nd; id++)
	mcptr->xc[id-1] = cptr->xc[id-1];
      cptr = cptr->prev;
    }
    /* Create left boundary cell list */
    start = NULL;
    clist = plev[ilev]->lbound;
    while (clist != NULL){
      mcptr = clist->cptr->mcptr;
      mclist = new mcell_list;
      mclist->mcptr = mcptr;
      mclist->nxt = start;
      start = mclist;
      clist = clist->nxt;
    }
    lbound[ilev] = start;
    /* Create right boundary cell list */
    start = NULL;
    clist = plev[ilev]->rbound;
    while (clist != NULL){
      mcptr = clist->cptr->mcptr;
      mclist = new mcell_list;
      mclist->mcptr = mcptr;
      mclist->nxt = start;
      start = mclist;
      clist = clist->nxt;
    }
    rbound[ilev] = start;
    /* Create left halo cells and list */
    mk_mclist(plev[ilev]->lhalo, lhalo[ilev]);
    /* Create right halo cells and list */
    mk_mclist(plev[ilev]->rhalo, rhalo[ilev]);
    /* Connect interior cells */
    mcptr = fcell[ilev];
    while (mcptr != NULL){
      connect_mcell(ilev, mcptr);
      mcptr = mcptr->nxt;
    }
    /* Connect left halo cells */
    mclist = lhalo[ilev];
    while (mclist != NULL){
      mcptr = mclist->mcptr;
      mcptr->ghost = false;
      connect_mcell(ilev, mcptr);      
      mclist = mclist->nxt;
    }
    /* Connect right halo cells */
    mclist = rhalo[ilev];
    while (mclist != NULL){
      mcptr = mclist->mcptr;
      mcptr->ghost = false;
      connect_mcell(ilev, mcptr);
      mclist = mclist->nxt;
    }
  }
  nomug = false;
  return;
}
//end create_mug.c

//file new_mcell.c
/*
+  new_mcell - create a new mcell.
     Description :
     Grabs memory  for a  new cell from  the mcarr. Increases  size of
     mcarr if necessary. Returns a pointer to the new cell.
*/
struct mcell *new_mcell()
{
  int icell;
  /* Check size of mcarr */
  if (mcarrtop < mugtop + MINCR){
    /* Resize carr array */
    mcarr = (struct mcell **)realloc((void *)mcarr, 
      (mcarrtop + MINCR)*sizeof(mcell *));
    for (icell = mcarrtop; icell < mcarrtop + MINCR; icell++){
      mcarr[icell] = new mcell;
      clean_mcell(mcarr[icell]);
    }
    mcarrtop += MINCR;
  }
  /*Grab memory from mcarr */
  mugtop++;
  return(mcarr[mugtop]);
}
//end new_mcell.c

//file clean_mcell.c
/*
+  clean_mcell - clean mcell
     Description :
     This cleans a  cell so that it  is in the same state  as after it
     has just been allocated.
*/
void clean_mcell(struct mcell *mcptr)
{
  int i, id, icell;
  /* Clear pointers */
  mcptr->nxt = NULL;
  for (icell = 0; icell < nrcell; icell++)
    mcptr->child[icell] = NULL;
  mcptr->par = NULL;
  for (id = 1; id <= nd; id++){
    mcptr->lcell[id-1] = NULL;
    mcptr->rcell[id-1] = NULL;
  }
  mcptr->cptr = NULL;
  for(i = 0; i < nq; i++){  
    mcptr->qa[i] = 0.0;
    mcptr->qt[i] = 0.0;
    mcptr->mso[i] = 0.0;
    mcptr->pso[i] = 0.0;
  }
  return;
}
//end clean_mcell.c

//file clean_mug.c
/*
+  clean_mug -  clean mug.
     Description :
     Removes  the mug  halo, boundary,  fcbound, bound  lists  and the
     communications buffers, minbufr and  moutbufr.  It does not clear
     mcarr or the minicells.  This is done by clear_mug.
*/
void clean_mug()
{
  int ilev, icell;
  struct mcell_list *mclist, *tmclist;
  if (nomug){
    /* Set lists to NULL etc (safety) and return if no multigrid */
    for (ilev = 0; ilev <= alev; ilev++){
      lhalo[ilev] = NULL;
      rhalo[ilev] = NULL;
      lbound[ilev] = NULL;
      rbound[ilev] = NULL;
      fcbound[ilev] = NULL;
      ghost[ilev] = NULL;
      nebd[ilev] = 0;
    }
    return;
  }
  /* Remove halo and boundary lists */
  for (ilev = 0; ilev <= alev; ilev++){
    /* Left halo list */
    mclist = lhalo[ilev];
    while (mclist != NULL){
      tmclist = mclist;
      mclist = mclist->nxt;
      delete tmclist;
    }
    lhalo[ilev] = NULL;
    /* Right halo list */
    mclist = rhalo[ilev];
    while (mclist != NULL){
      tmclist = mclist;
      mclist = mclist->nxt;
      delete tmclist;
    }
    rhalo[ilev] = NULL;
    /* Left boundary list */
    mclist = lbound[ilev];
    while (mclist != NULL){
      tmclist = mclist;
      mclist = mclist->nxt;
      delete tmclist;
    }
    lbound[ilev] = NULL;
    /* Right boundary list */
    mclist = rbound[ilev];
    while (mclist != NULL){
      tmclist = mclist;
      mclist = mclist->nxt;
      delete tmclist;
    }
    rbound[ilev] = NULL;
    /* Clear fcbound lists */
    mclist = fcbound[ilev];
    while (mclist != NULL){
      tmclist = mclist;
      mclist = mclist->nxt;
      delete tmclist;
    }
    fcbound[ilev] = NULL;
    /* Clear bound lists */
    mclist = ghost[ilev];
    while (mclist != NULL){
      tmclist = mclist;
      mclist = mclist->nxt;
      delete tmclist;
    }
    ghost[ilev] = NULL;
    nebd[ilev] = 0;
  }
  /* Clean up minicells */
  if (mcarr != NULL){
    for (icell = 0; icell < mcarrtop; icell++)
      clean_mcell(mcarr[icell]);
  }
  /* Set cell minicell pointers to NULL */
  for (icell = 1; icell <= celltop; icell++)
    carr[icell]->mcptr = NULL;
  /* Delete communications buffers */
  free(minbufr);
  free(moutbufr);
  minbufr = NULL;
  moutbufr = NULL;
  /* No multigrid */
  nomug = true;
  return;
}
//end clean_mug.c

//file clear_mug.c
/*
+  clear_mug -  clear mug.
     Description :
     Removes the  multigrid completely i.e. cleans it  and then clears
     the minicells and the mcarr  array. To be used when switching off
     implicit and clearing a model.
*/
void clear_mug()
{
  int icell;
  /* Clean multigrid (safety) */
  clean_mug();
  /* Clear mcarr array */
  if (mcarr != NULL){
    for (icell = 0; icell < mcarrtop; icell++)
      delete mcarr[icell];
    delete [] mcarr;
    mcarr = NULL;
  }
  /* Zero mugtop and mcarrtop (safety) */
  mugtop = 0;
  mcarrtop = 0;
}
//end clear_mug.c

//file connect_mcell.c
/*
+  connect_mcell -  connect a mini-cell.
     Description :
     Connects the mini-cell mcptr to  its neighbours and parent. Cells
     are created at  coarse/fine boundaries and their  parents are put
     in  the  list pointed  to  by  fcbound[ilev-1].  Ghost  cells  at
     external  boundaries   are  put  in   the  list  pointed   to  by
     ghost[ilev].
*/
void connect_mcell(int ilev, struct mcell *mcptr)
{
  bool doit;
  int id, icell, tid, it[MAXD], ichild;
  float xc[MAXD];
  struct cell *cptr, *left, *right, *pcptr, *ncptr;
  struct mcell *tmcptr, *nmcptr;
  struct mcell_list *mclist;
  /* Corresponding full cell */
  cptr = mcptr->cptr;
  /* Connect to Parent */
  pcptr = cptr->par;
  if (pcptr != NULL){
    mcptr->par = pcptr->mcptr;
    mcptr->par->child[cptr->ichild] = mcptr;
  }
  /* Connect to neighbours */
  for (id = 1; id <= nd; id++){
    /* Left neighbour */
    left = cptr->ljoin[id-1]->left;
    /* Initially NULL (safety) */
    mcptr->lcell[id-1] = NULL;
    if (left != NULL)
      /* Interior join */
      mcptr->lcell[id-1] = left->mcptr;
    else if (cptr->ljoin[id-1]->bound == lebound){
      /* Left exterior boundary - check type */
      if ((bcl[id-1] != bcref) && (!cptr->halo)){
	/* Not symmetric and not halo - create ghost cell */ 
	tmcptr = new_mcell();
	tmcptr->ghost = true;
	/* Cell coordinates and volume */
	for (tid = 1; tid <= nd; tid++){
	  if (tid != id)
	    tmcptr->xc[tid-1] = mcptr->xc[tid-1];
	  else
	    tmcptr->xc[tid-1] = mcptr->xc[tid-1]
	      - lev[ilev]->dxbar[tid-1];
	}
	tmcptr->vol = cell_vol(ilev, tmcptr->xc);
	/* Connect to ilev neighbour */
	mcptr->lcell[id-1] = tmcptr;
	tmcptr->rcell[id-1] = mcptr;
	/* Add to ghost list */
	mclist = new mcell_list;
	mclist->nxt = ghost[ilev];
	ghost[ilev] = mclist;
	mclist->mcptr = tmcptr;
      }
    }
    else if (ilev > 0){
      /* Check for coarse-fine boundary join - get ilev-1 neighbour */
      ncptr = cptr->par->ljoin[id-1]->left;
      if (ncptr != NULL){
	if ((!ncptr->refined) && (pcptr->refined)){
	  /* Coarse-fine boundary join - create if does not exist */
	  nmcptr = ncptr->mcptr;
	  /* Cell coordinates */
	  for (tid = 1; tid <= nd; tid++){
	    if (tid != id)
	      xc[tid-1] = mcptr->xc[tid-1];
	    else
	      xc[tid-1] = mcptr->xc[tid-1]
		- lev[ilev]->dxbar[tid-1];
	  }
	  if (xc[id-1] < xl[id-1])
	    /* Periodic id coordinate */
	    xc[id-1] += xr[id-1] - xl[id-1];
	  /* Find child index */
	  for (tid = 1; tid <= nd; tid++){
	    if (xc[tid-1] < nmcptr->xc[tid-1])
	      it[tid-1] = 1;
	    else
	      it[tid-1] = 2;
	  }
	  ichild = it2_to_i(it);
	  tmcptr = nmcptr->child[ichild];
	  if (tmcptr == NULL){
	    /* Check if parent is in boundary list */
	    doit = true;
	    for (icell = 0; icell < nrcell; icell++){
	      if (nmcptr->child[icell] != NULL){
		doit = false;
		break;
	      }
	    }
	    if (doit){
	      /*Add parent to boundary list */
	      mclist = new mcell_list;
	      mclist->nxt = fcbound[ilev-1];
	      fcbound[ilev-1] = mclist;
	      mclist->mcptr = nmcptr;
	    }
	    /* Create cell */
	    tmcptr = new_mcell();
	    /* Connect to parent */
	    nmcptr->child[ichild] = tmcptr;
	    for (tid = 1; tid <= nd; tid++)
	      tmcptr->xc[tid-1] = xc[tid-1];
	    tmcptr->vol = cell_vol(ilev, tmcptr->xc);
	  }
	  /* Connect to ilev neighbour */
	  mcptr->lcell[id-1] = tmcptr;
	  tmcptr->rcell[id-1] = mcptr;
	}
      }
    }
    /* Right neighbour */
    right = cptr->rjoin[id-1]->right;
    /* Initially NULL (safety) */
    mcptr->rcell[id-1] = NULL;
    if (right != NULL)
      /* Interior join */
      mcptr->rcell[id-1] = right->mcptr;
    else if (cptr->rjoin[id-1]->bound == rebound){
      /* Right exterior boundary - check type */
      if ((bcr[id-1] != bcref) && (!cptr->halo)){
	/* Not symmetric and not halo - create ghost cell */ 
	tmcptr = new_mcell();
	tmcptr->ghost = true;	
	/* Cell coordinates and volume */
	for (tid = 1; tid <= nd; tid++){
	  if (tid != id)
	    tmcptr->xc[tid-1] = mcptr->xc[tid-1];
	  else
	    tmcptr->xc[tid-1] = mcptr->xc[tid-1]
	      + lev[ilev]->dxbar[tid-1];
	}
	tmcptr->vol = cell_vol(ilev, tmcptr->xc);
	/* Connect to ilev neighbour */
	mcptr->rcell[id-1] = tmcptr;
	tmcptr->lcell[id-1] = mcptr;
	/* Add to ghost list */
	mclist = new mcell_list;
	mclist->nxt = ghost[ilev];
	ghost[ilev] = mclist;
	mclist->mcptr = tmcptr;
      }
    }
    else if (ilev > 0){
      /* Check for coarse-fine boundary join - get ilev-1 neighbour */
      ncptr = cptr->par->rjoin[id-1]->right;
      if (ncptr != NULL){
	if ((!ncptr->refined) && (pcptr->refined)){
	  /* Coarse-fine boundary join - create if does not exist */
	  nmcptr = ncptr->mcptr;
	  /* Cell coordinates */
	  for (tid = 1; tid <= nd; tid++){
	    if (tid != id)
	      xc[tid-1] = mcptr->xc[tid-1];
	    else
	      xc[tid-1] = mcptr->xc[tid-1]
		+ lev[ilev]->dxbar[tid-1];
	  }
	  if (xc[id-1] > xr[id-1])
	  /* Periodic id coordinate */
	    xc[id-1] += xl[id-1] - xr[id-1];
	  /* Find child index */
	  for (tid = 1; tid <= nd; tid++){
	    if (xc[tid-1] < nmcptr->xc[tid-1])
	      it[tid-1] = 1;
	    else
	      it[tid-1] = 2;
	  }
	  ichild = it2_to_i(it);
	  tmcptr = nmcptr->child[ichild];
	  if (tmcptr == NULL){
	    /* Check if parent is in boundary list */
	    doit = true;
	    for (icell = 0; icell < nrcell; icell++){
	      if (nmcptr->child[icell] != NULL){
		doit = false;
		break;
	      }
	    }
	    if (doit){
	      /*Add parent to boundary list */
	      mclist = new mcell_list;
	      mclist->nxt = fcbound[ilev-1];
	      fcbound[ilev-1] = mclist;
	      mclist->mcptr = nmcptr;
	    }
	    /* Create cell */
	    tmcptr = new_mcell();
	    /* Connect to parent */
	    nmcptr->child[ichild] = tmcptr;
	    for (tid = 1; tid <= nd; tid++)
	      tmcptr->xc[tid-1] = xc[tid-1];
	    tmcptr->vol = cell_vol(ilev, tmcptr->xc);
	  }
	  /* Connect to ilev neighbour */
	  mcptr->rcell[id-1] = tmcptr;
	  tmcptr->lcell[id-1] = mcptr;
	}
      }
    }
  }
  return;
}
//end connect_mcell.c

//file mk_mclist.c
/*
+  mk_mclist -  make mcell list.
     Description :
     Creates  and initialises  a halo  mcell list  on level  ilev. The
     start of the halo/boundary cell list is in clist and the start of
     the mcell  list is put in  start.  It assumes that  the full cell
     halo list exists.
*/
void mk_mclist(struct cell_list *clist, struct mcell_list
  *&start)
{
  int id;
  struct cell *cptr;
  struct mcell *mcptr;
  struct mcell_list *mclist;
  start = NULL;
  while (clist != NULL){
    cptr = clist->cptr;
    mcptr = new_mcell();
    mclist = new mcell_list;
    mclist->mcptr = mcptr;
    mclist->nxt = start;
    start = mclist;
    /* Connect mini-cell to full cell */
    cptr->mcptr = mcptr;
    mcptr->cptr = cptr;
    /* Coordinates of mini-cell */
    for (id = 1; id <= nd; id++)
      mcptr->xc[id-1] = cptr->xc[id-1];
    /* Volume of mini-cell */
    mcptr->vol = cptr->vol;
    clist = clist->nxt;
  }
}
//end mk_mclist.c

//file init_imp_mug.c
/*
+  init_mug -  initialise solution for implicit multigrid.
     Description :
     Initialise  the   mcells  with   the  solution  and   the  matrix
     coefficients.   Cells on  flev are  initialised with  the current
     solution, those on coarser levels with a time interpolated value.
     For levels coarser  than flev it also uses  the time interpolated
     qa to  reset the pa  in the full cells.  This assumes that  pa is
     correct for  the full cells on  flev.
*/
void init_imp_mug(int flev, float dt)
{
  int ilev, i, id;
  float tfrac, qa[MAXQ], pa[MAXQ], pl[MAXQ], pr[MAXQ];
  struct cell *cptr;
  struct mcell *mcptr, *tmcptr;
  struct mcell_list *mclist;
  /* Fill internal cells on levels levc-1 to flev */
  for (ilev = levc - 1; ilev <= flev; ilev++){
    /* Interpolation time */
    if ((ilev == flev) || (order == 1))
      tfrac = 0.0;
    else
      tfrac = (lev[flev]->nt + dt*lev[levc-1]->ndt/dt0 
	       - lev[ilev]->nt)/lev[ilev]->ndt;
    mcptr = fcell[ilev];
    while (mcptr != NULL){
      cptr = mcptr->cptr;
      /* Calculate solution */
      for (i = 0; i < nq; i++)
	qa[i] = cptr->qa[i] + tfrac*cptr->dqdt[i];
      if (ilev != flev)
	/* Calculate primitives and store in full cell */
	usrpq(qa, cptr->pa);
      else{
	/* Zero multigrid source terms for flev */
	for (i = 0; i < nq; i++)
	  mcptr->mso[i] = 0.0;
      }
      /* Calculate pa */
      usrpq(qa, pa);
      for (i = 0; i < nq; i++){
	/* Fill mini-cell qa with pa (step calculates pa) */	
	mcptr->qa[i] = pa[i];
	/* Set physical source term (old pa) */
	mcptr->pso[i] = pa[i];
      }
      mcptr = mcptr->nxt;
    }
  }
  /* Fill internal cells on levels levc-2 to 0 */
  for (ilev = levc - 2; ilev >= 0; ilev--){
    if (limpe > 1)
      /* Communicate ilev+1 halo qa (mcom1)  */
      mcom_halo(mcom1, ilev+1);
    mcptr = fcell[ilev];
    while (mcptr != NULL){
      /* Map up solution */
      map_up_mcell(mcptr);
      /* Calculate primitives and store in full cell */
      usrpq(mcptr->qa, mcptr->cptr->pa);
      mcptr = mcptr->nxt;
    }
  } 
  /* Fill external ghost cells */
  for (ilev = 0; ilev <= flev; ilev++){
    mclist = ghost[ilev];
    while (mclist != NULL){
      mcptr = mclist->mcptr;
      for (id = 1; id <= nd; id++){
	if (mcptr->xc[id-1] < xl[id-1]){
	  /* Left exterior ghost cell */
	  tmcptr = mcptr->rcell[id-1];
	  //usrpq(tmcptr->qa, pr);	// BM
	  usrbd(id, lebound, tmcptr->cptr->ljoin[id-1], mcptr->qa, pr);
	  break;
	}
	else if (mcptr->xc[id-1] > xr[id-1]){
	  /* Right exterior ghost cell */
	  tmcptr = mcptr->lcell[id-1];
	  //usrpq(tmcptr->qa, pl);	// BM
	  usrbd(id, rebound, tmcptr->cptr->rjoin[id-1], pl, mcptr->qa);
	  break;
	}
      }
      mclist = mclist->nxt;
    }
  }
  /* Map down internal boundary cells */
  for (ilev = 1; ilev < flev; ilev++){
    if (limpe > 1){
      /* Communicate ilev halo qa (mcom1)  */
      mcom_halo(mcom1, ilev);
      /* Communicate ilev halo dqdx (mcom2)*/
      mcom_halo(mcom2, ilev);
    }
    mclist = fcbound[ilev];
    while (mclist != NULL){
      mcptr = mclist->mcptr;
      map_down_mcell(ilev, mcptr);
      mclist = mclist->nxt;
    }
  }
  /* Calculate  matrix coefficients */
  for (ilev= 0; ilev <= flev; ilev++){
    cal_matrix_cof(ilev, dt);
  }
  return;
}  
//end init_imp_mug.c

//file init_grav_mug.c
/*
+  init_grav mug -  initialise solution for gravity multigrid.
     Description :
     Initialise  the mcells  with  the solution.   Cells  on alev  are
     initialised with  the current  solution, those on  coarser levels
     with a  time interpolated  value. Note that  it assumes  that the
     offset for the density is 0 and nq-1 for the potential.
*/
void init_grav_mug()
{
  int ilev;
  float tfrac, rho, pmass;
  struct cell *cptr;
  struct mcell *mcptr;
  struct mcell_list *mclist;
  struct part *pptr; 
  /* First fill of internal cells on levels levc-1 to alev */
  for (ilev = levc - 1; ilev <= alev; ilev++){ 
    mcptr = fcell[ilev];
    /* Interpolation time */
    tfrac = (float)(lev[alev]->nt - lev[ilev]->nt)
      /max(1,lev[ilev]->ndt); 
    while (mcptr != NULL){
      cptr = mcptr->cptr;
	/* Zero multigrid source (safety) */
      mcptr->mso[iqpot] = 0.0;
      if (mabs(tfrac) > tiny)
	/* Do density time interpolation */
	rho = cptr->qa[iqd] + tfrac*cptr->dqdt[iqd];
      else
	/* Use cell density */
	rho = cptr->qa[iqd];
      /* Put density in qa[iqpot] (temporary) */
      mcptr->qa[iqpot] = rho;
      if (particle){
	/* Add particle mass */
	pptr = mcptr->cptr->pptr;
	pmass = 0.0;
	while (pptr != NULL){
	  pmass += pptr->mass;
	  pptr = pptr->nxt;
	}
	mcptr->qa[iqpot] += pmass/mcptr->vol;
      }
      /* Put solution in qt (temporary) */
      mcptr->qt[iqpot] = cptr->qa[iqpot];
      mcptr = mcptr->nxt;
    }
  }
  /* Fill pso for internal cells on levels levc-2 to 0 with density */
  for (ilev = levc - 2; ilev >= 0; ilev--){
    if (limpe > 1)
      /* Communicate ilev+1 halo qa (mcom1)  */
      mcom_halo(mcom1, ilev+1);
    mcptr = fcell[ilev];
    while (mcptr != NULL){
      /* Map up solution */
      map_up_mcell(mcptr);
      mcptr->pso[iqpot] = mcptr->qa[iqpot];
      mcptr = mcptr->nxt;
    }
  }
  /* Final fill of internal cells on levels levc-1 to alev */
  for (ilev = levc - 1; ilev <= alev; ilev++){ 
    mcptr = fcell[ilev];
    while (mcptr != NULL){
      /* Set Poisson source */
      mcptr->pso[iqpot] = - 4.0*gkpi*mcptr->qa[iqpot]*mcptr->vol;
      /* Set solution */
      mcptr->qa[iqpot] = mcptr->qt[iqpot];
      mcptr = mcptr->nxt;
    }
  }
  /* Final fill of internal cells on levels levc-2 to 0 */
  for (ilev = levc - 2; ilev >= 0; ilev--){
    if (limpe > 1)
      /* Communicate ilev+1 halo qa (mcom1)  */
      mcom_halo(mcom1, ilev+1);
    mcptr = fcell[ilev];
    while (mcptr != NULL){
      /* Set Poisson source */
      mcptr->pso[iqpot] = - 4.0*gkpi*mcptr->pso[iqpot]*mcptr->vol;
      /* Map up solution */
      map_up_mcell(mcptr);
      mcptr = mcptr->nxt;
    }
  }
  /* Set boundary potential */
  cal_bound_pot();
  /* Map down coarse/fine boundary cells */
  for (ilev = 1; ilev < alev; ilev++){
    mclist = fcbound[ilev];
    if (limpe > 1){
      /* Communicate ilev halo qa (mcom1)  */
      mcom_halo(mcom1, ilev);
      /* Communicate ilev halo dqdx (mcom2)*/
      mcom_halo(mcom2, ilev);
    }
    while (mclist != NULL){
      mcptr = mclist->mcptr;
      map_down_mcell(ilev, mcptr);
      mclist = mclist->nxt;
    }
  }
  return;
}
// end init_grav_mug.c

//file cal_bound_pot.c
/*
+  cal_bound_pot -  calculate boundary potential.
     Description :
     If boundpot=true, then this  calculates the boundary potential on
     all levels  for boundaries  that are  not periodic  or symmetric.
     This is done by summing the  potential over levb for all external
     boundary cells  using the bound[ilev] lists.   If boundpot=false,
     then it  just sets the  potential to  zero in the  external ghost
     cells.
*/
void cal_bound_pot()
{
  int ilev, ncell[MAXPE], nvalues[MAXPE], ipe, ntot, nmax, id, icell,
    count, ios, noff, shift, maxcom;
  float xc[MAXD];
  struct mcell_list *mclist;
  struct mcell *mcptr;
  if (boundpot){
    if (limpe == 1){
      /* Single processor - just sum */
      for (ilev = 0; ilev <= alev; ilev++){
	mclist = ghost[ilev];
	while (mclist != NULL){
	  mcptr = mclist->mcptr;
	  mcptr->qa[iqpot] = sum_pot(mcptr->xc);
	  mclist = mclist->nxt;
	}
      }
    }
    else{
      /* Multiple processors */
      for (ilev = 0; ilev <= alev; ilev++){
	/* Get no of external boundary cells from other processors */
	ncell[0] = nebd[ilev];
	for (ipe = 0; ipe < limpe; ipe++)
	  nvalues[ipe] = 1;
	merge_array(ncell, nvalues);
	/* Broadcast no of external boundary cells on each processor */
	broadcast(ncell, limpe);
	/* Calculate global total and max */
	ntot = 0;
	nmax = -1;
	for (ipe = 0; ipe < limpe; ipe++){
	  ntot += ncell[ipe];
	  nmax = max(nmax, ncell[ipe]);
	}
	if (ntot > 0){
	  /* Check communication buffers */
	  maxcom = max(limpe*nmax,ntot*nd);
	  if (maxcom > comsize){
	    comsize = maxcom;
	    minbufr = (float *)realloc((void *)minbufr, 
	      (size_t)((comsize+1000)*sizeof(float)));
	    moutbufr = (float *)realloc((void *)moutbufr, 
	      (size_t)((comsize+1000)*sizeof(float)));
	  }
	  /* Put external boundary cell coords in moutbufr */
	  mclist = ghost[ilev];
	  count = -1;
	  while (mclist != NULL){
	    mcptr = mclist->mcptr;
	    for (id = 1; id <= nd; id++){
	      count++;
	      moutbufr[count] = mcptr->xc[id-1];
	    }
	    mclist = mclist->nxt;
	  }
	  for (ipe = 0; ipe < limpe; ipe++)
	    nvalues[ipe] = nd*ncell[ipe];
	  /* Merge coordinates */
	  merge_array(moutbufr, nvalues);
	  /* Broadcast coordinates */
	  broadcast(moutbufr, nd*ntot);
	  /* Contribution to potential from each processor */
	  count = -1;
	  for (icell = 0; icell < ntot; icell++){
	    for (id = 1; id <= nd; id++){
	      count++;
	      xc[id-1] = moutbufr[count];
	    }
	    moutbufr[icell] = sum_pot(xc);
	  }
	  /* Extract local values */
	  noff = 0;
	  for (ipe = 0; ipe < mype; ipe++)
	    noff += ncell[ipe];
	  count = -1;
	  mclist = ghost[ilev];
	  while (mclist != NULL){
	    mcptr = mclist->mcptr;
	    count++;
	    mcptr->qa[iqpot] = moutbufr[noff + count];
	    mclist = mclist->nxt;
	  }
	  /* Delete local values from moutbufr */
	  shift = ncell[mype]; 
	  for (icell = noff; icell < ntot - shift; icell++)
	    moutbufr[icell] = moutbufr[icell + shift];
	  /* Send contributions to each processor */
	  for (ipe = 0; ipe < limpe; ipe++){
	    ins[ipe] = 0;
	    outs[ipe] = ncell[ipe];
	  }
	  outs[mype] = 0;
	  ios = 0;
	  parc_xdatar(&limpe, outs, moutbufr, ins, minbufr, &ios);
	  /* Add potentials from other processors */
	  mclist = ghost[ilev];
	  count = -1;
	  while (mclist != NULL){
	    mcptr = mclist->mcptr;
	    count++;
	    for (ipe = 0; ipe < limpe - 1; ipe++)
	      mcptr->qa[iqpot] += minbufr[ipe*shift + count];
	    mclist = mclist->nxt;
	  }
	}
      }
    }
  }
  else{
    /* Just set them all to zero */
    for (ilev = 0; ilev <= alev; ilev++){
      mclist = ghost[ilev];
      while (mclist != NULL){
	mclist->mcptr->qa[iqpot] = 0.0;
	mclist = mclist->nxt;
      }
    }
  }
  return;
}
//end cal_bound_pot.c

//file sum_pot.c
/*
+  sum_pot -  sum boundary potential.
     Description :
     Sums over the levb cells to calculate the potential at xc.
*/
float sum_pot(float xc[])
{
  int id;
  float pot, dist2, tv;
  struct mcell *mcptr;
  /* Calculate boundary potential at xc by summing over levb */
  pot = 0.0;
  /* Sum over levb cells */
  mcptr = fcell[levb];
  while (mcptr != NULL){
    /* Calculate square of distance */
    dist2 = 0;
    for (id = 1; id <= nd; id++){
      tv = mcptr->xc[id-1] - xc[id-1];
      dist2 += tv*tv;
    }
    /* Increment potential */
    if (nd == 3)
      /* 3 dimensions */
      pot -= mcptr->cptr->qa[0]*mcptr->vol/sqrt(dist2);
    else if ((nd == 1) && (strcmp(geom, "RTP") == 0))
      /* Spherical symmetry */
      pot -= mcptr->cptr->qa[0]*mcptr->vol/fabs(xc[0]);      
    else
      /*  2 dimensions cartesian */
      pot += mcptr->cptr->qa[0]*mcptr->vol*logf(dist2);
    mcptr = mcptr->nxt;
  }
  return(pot);
}
//end sum_pot.c

//file cal_matrix_cof.c
/*
+  cal_matrix_cof -  Calculate matrix coefficients.
     Description :
     Calculate the implicit matrix  coefficients for the mini-cells on
     level ilev and and store them in the mcofp and mcofm arrays.
*/
void cal_matrix_cof(int ilev, float dt)
{
  int id, i, icell;
  float vol, areal, arear, nu[MAXQ], nul[MAXQ], nur[MAXQ], 
    dcoff[MAXQ], tv;
  struct cell *cptr;
  struct join *ljoin, *rjoin;
  struct mcell *mcptr, *pmcptr, *left, *right;
  struct mcell_list *mclist;
  /* Cell diffusion coeffts (store in mcptr->diff temporarily) */
  /* Internal cells */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    /* Diffusion coefft (note mcptr->qa is primitive variables) */    
    usrnu(mcptr->qa, nu);		// TODO skip inside solid body
    for (i = 0; i < nq; i++)
      mcptr->diff[i] = nu[i];
    mcptr = mcptr->nxt;
  }
  /* Internal boundary cells */
  if (ilev > 1){
    mclist = fcbound[ilev-1];
    while (mclist != NULL){
      pmcptr = mclist->mcptr;
      for (icell = 0; icell < nrcell; icell++){
	mcptr = pmcptr->child[icell];
	if (mcptr != NULL){
	  /* Diffusion coefft (note mcptr->qa is primitive variables) */
    	  usrnu(mcptr->qa, nu);
	  for (i = 0; i < nq; i++)
	    mcptr->diff[i] = nu[i];
	}
      }
      mclist = mclist->nxt;
    }
  }
  /* External boundary cells */
  mclist = ghost[ilev];
  // Temp Fix 28/01/25
  mclist = NULL;
  while (mclist != NULL){
    mcptr = mclist->mcptr;
    usrnu(mcptr->qa, nu);	// TODO skip inside solid body
    /* Diffusion coefft (note mcptr->qa is primitive variables) */
    for (i = 0; i < nq; i++)
      mcptr->diff[i] = nu[i];
    mclist = mclist->nxt;
  }
  if (limpe > 1)
    /* Communicate ilev halo nu (mcom3)  */
    mcom_halo(mcom3, ilev);			// TODO communicate state of cell being solid
  /* Calculate matrix coefficients */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    /* Corresponding full cell and volume*/
    cptr = mcptr->cptr;
    vol = cptr->vol;
    /* Initialise diagonal coefficient */
    for (i = 0; i < nq; i++)
      dcoff[i] = 1.0;
    /* For each direction */
    for (id = 1; id <= nd; id++){
      /* tv (= dt/vol*dx) */
      tv = dt/(vol*lev[ilev]->dxbar[id-1]);
      /* Left/right joins and mcells */
      ljoin = cptr->ljoin[id-1];
      rjoin = cptr->rjoin[id-1];
      left = mcptr->lcell[id-1];
      right = mcptr->rcell[id-1];
      /* Left/right areas */
      areal = ljoin->area;
      arear = rjoin->area;
      /* Left diffusion coefficients */
      if (left != NULL){
	// TODO check for left being a solid


	/* Internal */
	for (i = 0; i < nq; i++)
	  nul[i] = left->diff[i];
      }
      else{
	/* Left external symmetry boundary */
	for (i = 0; i < nq; i++)
	  nul[i] = mcptr->diff[i];
      }
      /* Right diffusion coefficients */
      if (right != NULL){
	// TODO check for right being a solid



	/* Internal */
	for (i = 0; i < nq; i++)
	  nur[i] = right->diff[i];
      }
      else{
	/* Right external symmetry boundary */
	for (i = 0; i < nq; i++)
	  nur[i] = mcptr->diff[i];
      }
      /* Interface diffusion coefficients (harmonic mean) */
      for (i = 0; i < nq; i++){
	nul[i] = 2.0/(1.0/max(mcptr->diff[i], tiny)
	  + 1.0/max(nul[i], tiny));
	nur[i] = 2.0/(1.0/max(mcptr->diff[i], tiny)
	  + 1.0/max(nur[i], tiny));
      }
      for (i = 0; i < nq; i++){
	/* Increment central coefft. */
	dcoff[i] += tv*(areal*nul[i] + arear*nur[i]);
	/* Off-diagonal matrix coeffts. */
	mcptr->mcofm[id-1][i] = tv*areal*nul[i];
	mcptr->mcofp[id-1][i] = tv*arear*nur[i];
      }
    }
    for (i = 0; i < nq; i++){
      /* Invert diagonal coefft.*/
      dcoff[i] = 1.0/dcoff[i];
      /* Multiply physical source term by inverse diagonal coefft. */
      mcptr->pso[i] *= dcoff[i];        
      /* Multiply off-diagonal coeffts. by inverse diagonal coefft. */
      for (id = 1; id <= nd; id++){
	mcptr->mcofm[id-1][i] *= dcoff[i];
	mcptr->mcofp[id-1][i] *= dcoff[i];
      }
    }
    mcptr = mcptr->nxt;
  }
  return;
}
//end cal_matrix_cof.c

//file set_mug_source.c
/*
+  set_mug_source -  set multigrid source terms.
     Description : 
     Calculates  the  multigrid  source  terms to  be  used  on  level
     ilev. It  assumes that  the residuals on  level ilev+1  have been
     set. mode=0 - gravity, mode=1 diffusion.

*/
void set_mug_source(int ilev)
{
  int i, il, ih;
  float tv;
  struct mcell *mcptr;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else{
    /* Diffusion */
    il = 0;
    ih = nq;
  }    
  if (limpe > 1)
    /* Communicate ilev+1 halo qa (mcom1) */
    mcom_halo(mcom1, ilev + 1);
  /* Zero multigrid source terms on ilev, map up qa from ilev+1 */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    for(i = il; i < ih; i++)
      mcptr->mso[i] = 0.0;
    map_up_mcell(mcptr);
    mcptr = mcptr->nxt;
  }
  /* Relax for one step on ilev */
  relax(ilev, 1);
  /* Put -qt in source, old qa (qa-qt) in qt, zero qa on ilev  */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    if (mcptr->cptr->refined){
      /* Set source term if refined */
      for (i = il; i < ih; i++)
	mcptr->mso[i] = - mcptr->qt[i];
    }
    /* Reset qt and qa */
    for (i = il; i < ih; i++){
	mcptr->qt[i] = mcptr->qa[i] - mcptr->qt[i];
	mcptr->qa[i] = 0;
    }
    mcptr = mcptr->nxt;
  }
  /* Exchange qa and qt on ilev+1 */
  mcptr = fcell[ilev+1];
  while (mcptr != NULL){
    for (i = il; i < ih; i++){
      tv = mcptr->qa[i];
      mcptr->qa[i] = mcptr->qt[i];
      mcptr->qt[i] = tv;
    }
    mcptr = mcptr->nxt;
  }
  if (limpe > 1)
    /* Communicate ilev+1 halo qa (mcom1) */
    mcom_halo(mcom1, ilev + 1);
  /* Map up qa from ilev+1  */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    map_up_mcell(mcptr);
    mcptr = mcptr->nxt;
  }
  /* Add qa (residual) to source, reset qa on ilev */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    for (i = il; i < ih; i++){
      mcptr->mso[i] += mcptr->qa[i];
      mcptr->qa[i] = mcptr->qt[i];
    }
    mcptr = mcptr->nxt;
  }
  /* Copy qt to qa on ilev+1 */
  mcptr = fcell[ilev+1];
  while (mcptr != NULL){
    for (i = il; i < ih; i++)
      mcptr->qa[i] = mcptr->qt[i];
    mcptr = mcptr->nxt;
  }
  return;
}
//end set_mug_source.c

//file correct_mug.c
/*
+  correct_mug -  correct from ilev to ilev+1.
     Description :
     Corrects the solution on ilev+1  after a number of relaxations on
     ilev. This  is done  by mapping down  the difference  between the
     solutions on the two grids and using this to correct the solution
     on ilev+1. mode=0 - gravity, mode=1 diffusion.
*/
void correct_mug(int ilev)
{
  int i, il, ih;
  struct mcell *mcptr;
  struct mcell_list *mclist;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else{
    /* Diffusion */
    il = 0;
    ih = nq;
  }  
  /* Copy qa to qt on ilev */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    for (i = il; i < ih; i++)
      mcptr->qt[i] = mcptr->qa[i];
    mcptr = mcptr->nxt;
  }
  if ((mode == 0) && boundpot){
    /* Copy ilev external boundary cell qa to qt, zero qa */
    mclist = ghost[ilev];
    while (mclist != NULL){
      mcptr = mclist->mcptr;
      mcptr->qt[iqpot] = mcptr->qa[iqpot];
      mcptr->qa[iqpot] = 0.0;
      mclist = mclist->nxt;
    }
  }
  if (limpe > 1)
    /* Communicate ilev+1 halo qa (mcom1) */
    mcom_halo(mcom1, ilev + 1);
  /* Map up qa from ilev+1 */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    map_up_mcell(mcptr);
    mcptr = mcptr->nxt;
  }
  /* Put qt-qa in qa on ilev */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    for (i = il; i < ih; i++)
      mcptr->qa[i] = mcptr->qt[i] - mcptr->qa[i];
    mcptr = mcptr->nxt;
  }
  /* Copy qa to qt on ilev+1 */
  mcptr = fcell[ilev+1];
  while (mcptr != NULL){
    for (i = il; i < ih; i++)
      mcptr->qt[i] = mcptr->qa[i];
    mcptr = mcptr->nxt;
  }
  if (limpe > 1){
    /* Communicate ilev halo qa (mcom1) */
    mcom_halo(mcom1, ilev);
    /* Communicate ilev halo dqdx (mcom2) */
    mcom_halo(mcom2, ilev);
  }
  /* Map down qa to ilev+1 */
  if (limpe > 1){
    /* Left halo cells */
    mclist = lhalo[ilev];
    while (mclist != NULL){
      map_down_mcell(ilev, mclist->mcptr);
      mclist = mclist->nxt;
    }
    /* Right halo cells */
    mclist = rhalo[ilev];
    while (mclist != NULL){
      map_down_mcell(ilev, mclist->mcptr);
      mclist = mclist->nxt;
    }
  }
  /* Interior cells */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    map_down_mcell(ilev, mcptr);
    mcptr = mcptr->nxt;
  }
  /* Correct ilev+1 */
  mcptr = fcell[ilev+1];
  while (mcptr != NULL){
    /* Put qa+qt in qa */
    for (i = il; i < ih; i++)
      mcptr->qa[i] += mcptr->qt[i];
    mcptr = mcptr->nxt;
  }
  /* Restore solution on ilev */
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    for (i = il; i < ih; i++)
      mcptr->qa[i] = mcptr->qt[i];
    mcptr = mcptr->nxt;
  }
  if (limpe > 1){
    /* Communicate ilev halo qa (mcom1) */
    mcom_halo(mcom1, ilev); 
    /* Communicate ilev halo dqdx (mcom2) */
    mcom_halo(mcom2, ilev);
  }
  /* Map down coarse/fine boundary cells on ilev */
  mclist = fcbound[ilev];
  while (mclist != NULL){
    mcptr = mclist->mcptr;
    map_down_mcell(ilev, mcptr);
    mclist = mclist->nxt;
  }
  return;
}
//end correct_mug.c

//file maxpot.c
/*
+  maxpot - maximum potential for multigrid control.
     Description :
        Returns  the maximum of the potential on all grids.
*/
float maxpot()
{
  int ilev;
  float pot;
  struct mcell *mcptr;
  pot = -big;
  for (ilev = 0; ilev <= alev; ilev++){
    mcptr = fcell[ilev];
    while (mcptr != NULL){
      pot = max(mabs(mcptr->qa[iqpot]), pot);
      mcptr = mcptr->nxt;
    }
  }
  pot = global_max(pot);
  return(pot);
}
//end maxpot.c

//file relax.c
/*
+  relax -  perform multigrid relaxations.
     Description :
     Perform nstep relaxations on grid  ilev. mode=0 - gravity, mode=1
     - diffusion.
*/
void relax(int ilev, int nstep)
{
  int istep, id, i;
  bool gotl, gotr;
  float dq, cofl, cofr, fact, dx, qn[MAXQ], ql[MAXQ], qr[MAXQ],
    diffl[MAXQ], diffr[MAXQ], tv;
  //struct cell *cptr;
  //struct join *jptr;
  struct mcell *mcptr, *left, *right;
  for (istep =1; istep <= nstep; istep++){
    if (limpe > 1)
      /* Communicate ilev halo qa (mcom1) */
      mcom_halo(mcom1, ilev);
    mcptr = fcell[ilev];
    if (mode == 0){
      /* Gravity */
      while (mcptr != NULL){
	/* Initialise */
	dq = 0.0;
	fact = 0.0;
	work++;
	/* For each direction */ 
	for (id = 1; id <= nd; id++){
	  /* Left/right mcells */
	  left = mcptr->lcell[id-1];
	  right = mcptr->rcell[id-1];
	  if (left != NULL)
	    /* Left interior cell */
	    ql[iqpot] = left->qa[iqpot];
	  else
	    /* Left external symmetry boundary */
	    ql[iqpot] = mcptr->qa[iqpot];
	  if (right != NULL)
	    /* Right interior cell */
	    qr[iqpot] = right->qa[iqpot];
	  else
	    /* Right external symmetry boundary */
	    qr[iqpot] = mcptr->qa[iqpot];
	  /* Calculate increment */
	  dx = lev[ilev]->dxbar[id-1];
	  cofl = mcptr->cptr->ljoin[id-1]->area/dx;
	  cofr = mcptr->cptr->rjoin[id-1]->area/dx;
	  dq += cofl*ql[iqpot] + cofr*qr[iqpot];
	  fact += cofl + cofr;
	}
	/* Add physical source */
	dq = (dq + mcptr->pso[iqpot])/fact;
	/* Increment solution and residual */
	mcptr->qt[iqpot] = dq - mcptr->qa[iqpot];
	mcptr->qa[iqpot] = dq;
	mcptr = mcptr->nxt;
      }
      /* Add multigrid source */
      mcptr = fcell[ilev];
      while (mcptr != NULL){
	dq = mcptr->mso[iqpot];
	mcptr->qt[iqpot] += dq;
	mcptr->qa[iqpot] += dq;
	mcptr = mcptr->nxt;
      }
    }
    else {
      /* Diffusion */
      while (mcptr != NULL){
	/* Initialise */
	for (i = 0; i < nq; i++)
	  qn[i] = 0.0;
	work++;					// TODO cells that are solid should do no work
	/* For each direction */ 
	for (id = 1; id <= nd; id++){
	  /* Left/right mcells */
	  left = mcptr->lcell[id-1];
	  right = mcptr->rcell[id-1];
	  if ((left != NULL) && !left->ghost){
		// TODO special treatment for solid body?  We may already cover this with cal_matrix_cof function.... 
	    /* Left interior cell */
	    for (i = 0; i < nq; i++)
	      ql[i] = left->qa[i];
	  }
	  else{
	    /* Left external boundary */
	    gotl = cal_mdiffl(id, mcptr, diffl);	    
	    for (i = 0; i < nq; i++)
	      ql[i] = mcptr->qa[i] - diffl[i];
	  }
	  if ((right != NULL) && !right->ghost){
	    /* Right interior cell */
	    for (i = 0; i < nq; i++)
	      qr[i] = right->qa[i];
	  }
	  else{
		// TODO special treatment for solidbody
	    /* Right external boundary */
	    gotr = cal_mdiffr(id, mcptr, diffr);	    	    
	    for (i = 0; i < nq; i++)
	      qr[i] = mcptr->qa[i] + diffr[i];
	  }
	  /* Increment new solution */
	  for (i = 0; i < nq; i++)
	    qn[i] += mcptr->mcofm[id-1][i]*ql[i] + 
	      mcptr->mcofp[id-1][i]*qr[i];
	}
	/* Add physical source, set solution and residual */
	for (i = 0; i < nq; i++){
	  qn[i] += mcptr->pso[i];
	  mcptr->qt[i] = qn[i] - mcptr->qa[i];
	  mcptr->qa[i] = qn[i];
	}
	mcptr = mcptr->nxt;
      }
      /* Add multigrid source */
      mcptr = fcell[ilev];
      while (mcptr != NULL){
	for (i = 0; i < nq; i++){
	  tv = mcptr->mso[i];
	  mcptr->qt[i] += tv;
	  mcptr->qa[i] += tv;
	}
	mcptr = mcptr->nxt;
      }
    }
  }
  /* Explicitly register a "use" to quiet the compiler warning BM */
  (void)gotl;
  (void)gotr;
  return;
}
//end relax.c

//file map_up_mcell.c
/*
+  map_up_mcell - map up from child mcells.
     Description :
     Maps  up from  the child  cells given  the pointer  pcptr to  the
     parent cell. mode=0 - gravity, mode=1 diffusion.
*/
void map_up_mcell(struct mcell *pmcptr)
{
  int icell, i, il, ih;
  float divvolp, qa[MAXQ];
  struct mcell *mcptr;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else{
    /* Diffusion */
    il = 0;
    ih = nq;
  }    
  /* Only do refined cells */
  if (!pmcptr->cptr->refined)
    return;
  /* Map up solution */
  divvolp = 1.0/pmcptr->vol;
  for (i = il; i < ih; i++)
    qa[i] = 0.0;
  for (icell = 0; icell < nrcell; icell++){
    mcptr = pmcptr->child[icell];
    for (i = il; i < ih; i++)
      qa[i] += mcptr->qa[i]*mcptr->vol*divvolp;			// TODO be mindful of cells that are solid vs not solid
  }
  for (i = il; i < ih; i++)
    pmcptr->qa[i] = qa[i];
  return;
}
//end map_up_mcell.c

//file map_down_mcell.c
/*
+  map_down_mcell - maps down from parent mcells.
     Description :
     Maps down from the parent cell given by the pointer pmcptr to the
     parent cell.  Note that it assumes that a NULL left/right pointer
     either means a halo or  external symmetry boundary cell. mode=0 -
     gravity, mode=1 diffusion.
*/
void map_down_mcell(int ilev, struct mcell *pmcptr)
{
  int id, i, il, ih, rind, icell;
  bool gotl, gotr;
  float pdivdx, pxc[MAXD], cxc[MAXD], dqdx[MAXD][MAXQ], diffl[MAXQ], 
    diffr[MAXQ], qc[MAXQ], vc, vi, vo, tv;
  struct mcell *mcptr, *left, *right;
  struct cell *pcptr;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else{
    /* Diffusion */
    il = 0;
    ih = nq;
  }    
  /* Corresponding full cell */
  pcptr = pmcptr->cptr;
  /* Get gradients */
  for (id = 1; id <= nd; id++){
    /* Coarse cell coordinate */
    pxc[id-1] = pcptr->xc[id-1];
    /* Coarse level inverse mesh spacing */
    pdivdx = 1.0/lev[ilev]->dxbar[id-1];
    /* Left difference */
    left = pmcptr->lcell[id-1];
    if (mode == 0){
      /* Gravity */
      if (left != NULL){
	/* Left cell exists */
	diffl[iqpot] = pmcptr->qa[iqpot] - left->qa[iqpot];
      }
      else if ((id == nd) && (pcptr->halo))
	/* Left halo cell with difference - use cell diff */
	diffl[iqpot] = pmcptr->diff[iqpot];
      else{
	/* Left external boundary - check for symmetry */
	if (bcl[id-1] == bcref)
	  /* Symmetric */
	  diffl[iqpot] = 0.0;
      }
    }
    else
      /* Diffusion */
      gotl = cal_mdiffl(id, pmcptr, diffl);
    /* Right difference */
    right = pmcptr->rcell[id-1];
    if (mode == 0){
      /* Gravity */
      if (right != NULL)
	/* Right cell exists */
	diffr[iqpot] = right->qa[iqpot] - pmcptr->qa[iqpot];
      else if ((id == nd) && (pcptr->halo))
	/* Right halo cell with difference - use cell diff */
	diffr[iqpot] = pmcptr->diff[iqpot];
      else{
	/* Right external boundary - check for symmetry */
	if (bcr[id-1] == bcref)
	  /* Symmetric */
	  diffr[iqpot] = 0.0;
      }
    }
    else
      /* Diffusion */
      gotr = cal_mdiffr(id, pmcptr, diffr);
    for (i = il; i < ih; i++)
      dqdx[id-1][i] = 0.5*(diffl[i] + diffr[i])*pdivdx;
  }
  /* Set central qa 	// TODO be mindful of cells that are solid vs not solid */
  for (i = il; i < ih; i++)
    qc[i] = pmcptr->qa[i];
  if (strcmp(geom, "RZP") == 0){
    /* Cylindrical polars - correct for cell volume */
    vc = pmcptr->vol;
    /* Inner child cell volume */
    if (pmcptr->child[0] != NULL)
      vi = pmcptr->child[0]->vol;
    else{
      for (id = 1; id <= nd; id++)
	cxc[id-1] = pcptr->xc[id-1] - 0.5*lev[ilev+1]->dxbar[id-1];
      vi = cell_vol(ilev+1, cxc);
    }
    /* Outer child cell volume */
    if (strcmp(coord[0], "R") == 0){
      rind = 0;
      if (pmcptr->child[1] != NULL)
	vo = pmcptr->child[1]->vol;
      else{
	for (id = 1; id <= nd; id++)
	  cxc[id-1] = pcptr->xc[id-1] + 0.5*lev[ilev+1]->dxbar[id-1];
	vo = cell_vol(ilev+1, cxc);
      }
    }
    else{
      rind = 1;
      if (pmcptr->child[2] != NULL)
	vo = pmcptr->child[2]->vol;
      else{
	for (id = 1; id <= nd; id++)
	  cxc[id-1] = pcptr->xc[id-1] + 0.5*lev[ilev+1]->dxbar[id-1];
	vo = cell_vol(ilev+1, cxc);
      }
    }
    tv = lev[ilev+1]->dxbar[rind]*(vo - vi)/vc;
    for (i = il; i < ih; i++)
      qc[i] -= tv*dqdx[rind][i];
  }
  /* Map down solution */
  for (icell = 0; icell < nrcell; icell++){
    mcptr = pmcptr->child[icell];
    if (mcptr != NULL){
      /* Central solution */
      for (i = il; i < ih; i++)
	mcptr->qa[i] = qc[i];
      for (id = 1; id <= nd; id++){
	for (i = il; i < ih; i++)
	  mcptr->qa[i] += (mcptr->xc[id-1] - pxc[id-1])*dqdx[id-1][i];
      }
    }
  }
  /* Explicitly register a "use" to quiet the compiler warning BM */
  (void)gotl;
  (void)gotr;  
  return;
}
//end map_down_mcell.c

//file cal_mdiffl.c
/*
+ cal_mdiffl - calculates left mcell gradient
     Description :
        Calculates the left mcell gradient in in direction id, for the
        mcell pointed to by mcptr.   Returns true if the ordinary left
        cell  exists on  this level.   The difference  is returned  in
        diff). Note qa in mcell is pa. Only used for diffusion.
*/
bool cal_mdiffl(int id, struct mcell *mcptr, float diff[])
{
  int i;
  float ql[MAXQ], qr[MAXQ];
  bool got;
  struct cell *cptr;  
  struct cell *left;
  struct join *ljoin;  
  cptr = mcptr->cptr;
  ljoin = cptr->ljoin[id-1];
  left = ljoin->left;
  if (left == NULL){
    if (ljoin->bound == lebound){
      /* Left external boundary */
      got = true;
      for (i = 0; i < nq; i++)
	qr[i] = mcptr->qa[i];
      usrbd(id, lebound, ljoin, ql, qr);
      for (i = 0; i < nq; i++)
	diff[i] = qr[i] - ql[i];
    }
    else if ((id == nd) && (cptr->halo) && (cptr->gotdiff)){
      /* Left halo cell with difference - use cell diff */
      got = true;
      for (i = 0; i < nq; i++)
	diff[i] = mcptr->diff[i];
    }
    else{
      got = false;
      for (i = 0; i < nq; i++)
	diff[i] = 0.0;
    }
  }
  else{
    struct mcell *mleft;
    got = true;
    mleft = left->mcptr;
    for (i = 0; i < nq; i++)
      diff[i] = mcptr->qa[i] - mleft->qa[i];
  }
  return(got);
}
//end cal_mdiffl.c

//file cal_mdiffr.c
/*
+ cal_mdiffr - calculates right mcell gradient
     Description :
        Calculates the  right mcell gradient  in in direction  id, for
        the mcell pointed  to by mcptr.  Returns true  if the ordinary
        right cell exists  on this level.  The  difference is returned
        in diff). Note qa in mcell is pa. Only used for diffusion.
*/
bool cal_mdiffr(int id, struct mcell *mcptr, float diff[])
{
  int i;
  float ql[MAXQ], qr[MAXQ];
  bool got;
  struct cell *cptr;
  struct cell *right;
  struct join *rjoin;
  cptr = mcptr->cptr;  
  rjoin = cptr->rjoin[id-1];
  right = rjoin->right;
  if (right == NULL){
    if (rjoin->bound == rebound){
      /* Right external boundary */
      got = true;
      for (i = 0; i < nq; i++)
	ql[i] = mcptr->qa[i];
      usrbd(id, rebound, rjoin, ql, qr);
      for (i = 0; i < nq; i++)
	diff[i] = qr[i] - ql[i];
    }
    else if ((id == nd) && (cptr->halo) && (cptr->gotdiff)){
      /* Right halo cell with difference - use cell diff */
      got = true;
      for (i = 0; i < nq; i++)
	diff[i] = mcptr->diff[i];
    }
    else{
      got = false;
      for (i = 0; i < nq; i++)
	diff[i] = 0.0;
    }
  }
  else{
    struct mcell *mright;
    got = true;
    mright = right->mcptr;
    for (i = 0; i < nq; i++)
      diff[i] = mright->qa[i] - mcptr->qa[i];
  }
  return(got);
}
//end cal_mdiffr.c

//file norm.c
/*
+  norm - norm for multigrid control.
     Description :
     Returns the norm of the  residual on ilev (per variable).  mode=0
     - gravity, mode=1 diffusion.
  TODO account for internal geom
*/
float norm(int ilev)
{
  int i, il, ih;
  float res[MAXQ], vol, tvol, pref[MAXQ], tv;
  struct mcell *mcptr;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else{
    /* Diffusion */
    il = 0;
    ih = nq;
  }
  for (i = il; i < ih; i++)  
    res[i] = 0.0;
  tvol = 0.0;
  mcptr = fcell[ilev];
  while (mcptr != NULL){
    usrpref(mcptr, pref);
    vol = mcptr->vol;
    tvol += vol;
    for (i = il; i < ih; i++)
      res[i] += vol*mabs(mcptr->qt[i])/pref[i];
    mcptr = mcptr->nxt;
  }
  tv = -1.0;
  for (i = il; i < ih; i++)
    tv = max(res[i], tv);
  tv = global_sum(tv);
  tvol = global_sum(tvol);
  tvol = max(tiny, tvol);
  return(tv/tvol);  
}
//end norm.c

//file trunc.c
/*
+  trunc - truncation error  for multigrid control.
     Description :
     Returns an  estimate of  the truncation error  on level  ilev. It
     assumes that  the multigrid  source term has  been set.  mode=0 -
     gravity, mode=1 diffusion.
  TODO account for internal geom
*/
float trunc(int ilev)
{
  int i, il, ih;
  float mso[MAXQ], vol, tvol, pref[MAXQ], tv;
  struct mcell *mcptr;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else{
    /* Diffusion */
    il = 0;
    ih = nq;
  }
  for (i = il; i < ih; i++)  
    mso[i] = 0.0;  
  tvol = 0.0;
  mcptr = fcell[ilev-1];
  while (mcptr != NULL){
    usrpref(mcptr, pref);
    vol = mcptr->vol;
    tvol += vol;
    for (i = il; i < ih; i++)
      mso[i] += vol*mabs(mcptr->mso[i])/pref[i];
    mcptr = mcptr->nxt;
  }
  tv = -1.0;
  for (i = 0; i < nq; i++)
    tv = max(mso[i], tv);   
  /* Return 0.25 of average */
  tv = global_sum(tv);
  tvol = global_sum(tvol);
  tvol = max(tiny, tvol);
  return(0.25*tv/tvol);
}
//end trunc.c

//file mcom_halo.c
/*
+  mcom_halo - communicate solution to halo cells.
     Description :     
     This fills the halo cells on level ilev with information from the
     boundary cells on  the neighbouring processor. The  nature of the
     information  is determined  by mcom  as follows:  mcom1 qa  mcom2
     dqdx. mode=0 - gravity, mode=1 diffusion.
*/
void mcom_halo(int mcom, int ilev)
{
  int ios;
  /* Return for limpe=1 (safety) */
  if (limpe == 1)
    return;
  /* Assemble information from left boundary cells */
  put_mcom(mcom, ilev, true); 
  /* Send information to left processor */
  ios = 0;
  parc_xdatar(&limpe, outs, moutbufr, ins, minbufr, &ios);
  /* Put received data in right halo cells */
  receive_mcom(mcom, ilev, false);
  /* Assemble information from right boundary cells */
  put_mcom(mcom, ilev, false);
  /* Send information to right processor */
  ios = 0;
  parc_xdatar(&limpe, outs, moutbufr, ins, minbufr, &ios);
  /* Put received data in left halo cells */
  receive_mcom(mcom, ilev, true);
  return;
}
//end mcom_halo.c

//file put_mcom.c
/*
+  put_mcom -  put mcom values in outbufr.
     Description :
     Puts the values determined by  mcom in outbufr. mode=0 - gravity,
     mode=1 diffusion.
*/
void put_mcom(int mcom, int ilev, bool left)
{
  int npe, nprocs, count;
  struct mcell_list *start;
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
    start = lbound[ilev];
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
    start = rbound[ilev];
    /* Start of processor list */
    plist = plev[ilev]->rproc;
    splist = plist;
    /* No of right processors */
    nprocs = plev[ilev]->nrproc;
  }
  /* Set count */
  count = -1;
  /* Put data in outbufr as appropriate */
  if (mcom == mcom1){
    while (plist != NULL){
      put_mcom1(ilev, count, start);
      plist = plist->nxt;
    }
  }
  else if (mcom == mcom2){
    while (plist != NULL){
      put_mcom2(ilev, left, count, start);
      plist = plist->nxt;
    }
  }
  else if ((mcom == mcom3) && (mode == 1)){
    while (plist != NULL){
      put_mcom3(ilev, count, start);
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
//end put_mcom.c

//file put_mcom1.c
/*
+  put_mcom1 -  put qa in moutbufr.
     Description :
     Puts qa from  the boundary cells in moutbufr.   Side is specified
     by left. mode=0 - gravity, mode=1 diffusion.
*/
void put_mcom1(int ilev, int &count, struct mcell_list *mclist)
{
  int i, il, ih;
  struct mcell *mcptr;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else {
    /* Diffusion */
    il = 0;
    ih = nq;
  }
  while (mclist != NULL){
    mcptr = mclist->mcptr;
    /* Solution */
    for (i = il; i < ih; i++){
      count++;
      moutbufr[count] = mcptr->qa[i];
    }
    mclist = mclist->nxt;
  }
  return;
}
//end put_mcom1.c

//file put_mcom2.c
/*
+  put_mcom2 -  put dqdx in moutbufr.
     Description :
     Puts dqdx from the boundary cells in moutbufr.  Side is specified
     by left. mode=0 - gravity, mode=1 diffusion.     
*/
void put_mcom2(int ilev, bool left, int &count, 
  struct mcell_list *mclist)
{
  int i, il, ih;
  float diff[MAXQ];
  struct mcell *mcptr, *mleft, *mright;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else{
    /* Diffusion */
    il = 0;
    ih = nq;
  }  
  while (mclist != NULL){
    mcptr = mclist->mcptr;
    /* Calculate differences */
    if (left){
      /* Left processor boundary - right difference */
      mright = mcptr->rcell[nd-1];
      if (mright != NULL){
	/* Right cell exists */
	for (i = il; i < ih; i++)
	  diff[i] = mright->qa[i] - mcptr->qa[i];
      }
      else{
	/* Right external boundary - check for symmetry */
	if (bcr[nd-1] == bcref)
	  /* Symmetric */
	for (i = il; i < ih; i++)	  
	  diff[i] = 0.0;
	else
	  /* Zero boundary condition */
	for (i = il; i < ih; i++)	  
	  diff[i] = - 2.0*mcptr->qa[i];	
      }
    }
    else{
      /* Right processor boundary - left difference */
      mleft = mcptr->lcell[nd-1];
      if (mleft != NULL)
	/* Left cell exists */
	for (i = il; i < ih; i++)
	  diff[i] = mcptr->qa[i] - mleft->qa[i];
      else{
	/* Left external boundary - check for symmetry */
	if (bcl[nd-1] == bcref)
	  /* Symmetric */
	  for (i = il; i < ih; i++)	  
	    diff[i] = 0.0;
	else
	  /* Zero boundary condition */
	  for (i = il; i < ih; i++)	  
	    diff[i] = 2.0*mcptr->qa[i];
      }
    }
    /* Difference (dqdx) */
    for (i = il; i < ih; i++){
      count++;
      moutbufr[count] = diff[i];
    }
    mclist = mclist->nxt;
  }
  return;
}
//end put_mcom2.c

//file put_mcom3.c
/*
+  put_mcom3 -  put effective diffusion coefficient in moutbufr.
     Description :     
     Puts the effective diffusion  coefficient from the boundary cells
     in moutbufr. Side is specified by left. Only used for diffusion. 
*/
void put_mcom3(int ilev, int &count,
  struct mcell_list *mclist)
{
  int i;
  float nu[MAXQ];
  struct mcell *mcptr;
  while (mclist != NULL){
    mcptr = mclist->mcptr;
    /* Diffusion coefft (note mcptr->qa is primitive variables) */
    usrnu(mcptr->qa, nu);
    for (i = 0; i < nq; i++){
      count++;
      moutbufr[count] = nu[i];
    }
    mclist = mclist->nxt;
  }
  return;
}
//end put_mcom3.c

//file receive_mcom.c
/*
+  receive_mcom -  receive mcom values from minbufr.
     Description :
     Receives the values  determined by mcom in minbufr  and puts them
     in the halo or boundary cells.  Side is specified by left. mode=0
     - gravity, mode=1 diffusion.
*/
void receive_mcom(int mcom, int ilev, bool left)
{
  /* Return if no data */
  if (left && (plev[ilev]->lsproc < 0))
    return;
  if (!left && (plev[ilev]->rsproc < 0))
    return;
  if (mcom == mcom1)
    receive_mcom1(ilev, left);
  else if (mcom == mcom2)
    receive_mcom2(ilev, left);
  else if ((mcom == mcom3) && (mode == 1))
    receive_mcom3(ilev, left);
  return;
}
//end receive_mcom.c

//file receive_mcom1.c
/*
+  receive_mcom1 -  receive qa.
     Description :
     Receives qa from  minbufr and puts it in the  halo cells. Side is
     specified by left. mode=0 - gravity, mode=1 diffusion.
*/
void receive_mcom1(int ilev, bool left)
{
  int count, i, il, ih;
  struct mcell *mcptr;
  struct mcell_list *mclist;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else {
    /* Diffusion */
    il = 0;
    ih = nq;
  }  
  if (left)
    /* Left halo */
    mclist = lhalo[ilev];
  else
    /* Right halo */
    mclist = rhalo[ilev];
  count = -1;
  while (mclist != NULL){
    mcptr = mclist->mcptr;
    /* Solution */
    for (i = il; i < ih; i++){
      count++;
      mcptr->qa[i] = minbufr[count];
    }
    mclist = mclist->nxt;
  }
  return;
}
//end receive_mcom1.c

//file receive_mcom2.c
/*
+  receive_mcom2 -  receive dqdx.
     Description :
     Receives dqdx from minbufr and puts it in the halo cells. Side is
     specified by left. mode=0 - gravity, mode=1 diffusion.     
*/
void receive_mcom2(int ilev, bool left)
{
  int count, i, il, ih;
  struct mcell *mcptr;
  struct mcell_list *mclist;
  if (mode == 0){
    /* Gravity */
    il = iqpot;
    ih = iqpot + 1;
  }
  else {
    /* Diffusion */
    il = 0;
    ih = nq;
  }    
  if (left)
    /* Left halo */
    mclist = lhalo[ilev];
  else
    /* Right halo */
    mclist = rhalo[ilev];
  count = -1;
  while (mclist != NULL){
    mcptr = mclist->mcptr;
    /* Difference (dqdx) */
    for (i = il; i < ih; i++){
      count++;
      mcptr->diff[i] = minbufr[count];
    }
    mclist = mclist->nxt;
  }
  return;
}
//end receive_mcom2.c

//file receive_mcom3.c
/*
+  receive_mcom3 -  receive mu.
     Description :
     Receives the  diffusion coefficient from  minbufr and puts  it in
     diff  in the  halo cells.  Side is  specified by  left.
*/
void receive_mcom3(int ilev, bool left)
{
  int count, i;
  struct mcell *mcptr;
  struct mcell_list *mclist;
  if (left)
    /* Left halo */
    mclist = lhalo[ilev];
  else
    /* Right halo */
    mclist = rhalo[ilev];
  count = -1;
  while (mclist != NULL){
    mcptr = mclist->mcptr;
    /* Diffusion coefficient */
    for (i = 0; i < nq; i++){
      count++;
      mcptr->diff[i] = minbufr[count];
    }
    mclist = mclist->nxt;
  }
  return;
}
//end receive_mcom3.c



// TODO need to communicate whether cell is solid or not...
