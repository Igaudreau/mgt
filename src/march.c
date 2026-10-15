#undef MAIN
#include <iostream>
#include<iomanip>
#include <cstring>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "funcs.h"
extern void merge(int ilev);
extern void correct_flux(int ilev);
extern void re_grid(int ilev);
extern void plot(char *key);
extern void usrpq(float qa[], float pa[]);
//extern void usrmon();
extern void usrfxq(struct cell * cptr);
extern void usrbd(int id, int bound, join *jptr, float pl[], 
  float pr[]);
extern void step(int ilev);
extern void cal_grav();
extern void clean_mug();
extern void de_frag();
extern int re_part();
extern void merge_array(int buf[], int nvalues[]);
extern void merge_array(float buf[], int nvalues[]);
extern void com_halo(int com, int ilev, float tfrac);
extern float global_max(float value);
extern void broadcast(float buf[], int nvalues);
extern int dpdx_index(int jid, int id);
extern void report();
extern bool writem(char filename[]);
void march(int nsteps);
void integrate(int ilev);
void adv_sol(int ilev);
void set_join_int(int ilev, int jid, join *jptr);
void cal_diffl(int ilev, int id, float tfrac, struct cell *cptr, 
  float diffl[]);
void cal_diffr(int ilev, int id, float tfrac, struct cell *cptr, 
  float diffr[]);

//file march.c
/*
+  march -  march solution for nsteps.
      Description  :
         This calls integrate to march  the solution if target is less
         than the current  time. It plots and prints  under control of
         admin and debug.
*/
void march(int nsteps)
{
  char key[SZNAM], tail[SZNAM];
  int istep, ilev, nvals[MAXPE], npe, procupdat[MAXPE], onupdat;
  float dwall, twork, nwall, owall, minwall, maxwall, tv, pcpu[MAXPE],
    procwork[MAXPE], procworkwt[MAXPE];
  bool quit = false;
  if ((target > 0.0) && (target <= t)){
    if (mype == 0)
      cout << "target reached" << '\n';
    return;
  }
  /* Initialise timing */
  nupdat = 0;
  onupdat = 0;
  dwall = 0.0;
  wall_time = 0.0;
  parc_clock(&ostep_time);
  if (limpe > 1 ){
    /* Set nvals for merge_array */
    for (npe = 0; npe < limpe; npe++)
      nvals[npe] = 1;
  }
  for (istep = 1; istep <= nsteps; istep++){
    /* Check target */
    if (t < target){
      /* Initialise timers for timestep */
      pcpu[0] = - wall_time;
      parc_clock(&owall);
      onupdat = nupdat;
      /* Integrate */
      //      mdq = -1.0;
      for (ilev = 0; ilev < nlevs; ilev++)
        lcnv[ilev] = cnv;
      integrate(levc-1);
      /* Merge inactive levels */
      for (ilev = levc - 2; ilev >= 0; ilev--)
	merge(ilev);
      if (limpe > 1){
	pcpu[0] = max(pcpu[0] + wall_time, 1.0e-5);
	merge_array(pcpu, nvals);
	broadcast(pcpu, limpe);
	if (timpart){
	  /* Update workwt */
	  tv = 0.0;
	  for (npe = 0; npe < limpe; npe++)
	    tv += max(pcpu[npe], 1.0e-5);
	  /* Sam's method - doesn't seem to balance nicely for BM
	  workwt *= sqrt(max(pcpu[mype], 1.0e-5)*limpe/tv); */
          /* Update workwt using old method - BM*/
	  workwt *= sqrt(pcpu[mype]*limpe/tv);
	  workwt = max(0.1, min(10.0, workwt));
	}
	else
	  /* Set to 1.0 */
	  workwt = 1.0;
      }
      else
	/* Set to 1.0 */
	workwt = 1.0;
      t = t + dt0;
      if ((limpe > 1) && repart)
	re_part();
      else if ((nlevs > 2) && defrag)
	de_frag();
      /* Dump */
      if ((dump) && (mype ==0)){
	cout.setf(ios::scientific);
	cout <<  "new t " << setprecision(4) << t
	     << " dt " << dt0 << " after step " << setw(4) 
	     << istep << '\n';
      }
      /* Admin */
      if (admin){
	strcpy(key, "P");
	plot(key);
      }
      /* Monitor */
      //usrmon();
      /* Increment times */
      parc_clock(&nwall);
      dwall += nwall - owall;
      /* Check cpu time */
      quit = cputime - nwall < 2.0*(nwall - owall);
      if (quit)
	break;
    }
    else
      break;
  }
  /* Timing */
  dwall = max(dwall, 1.0e-5);
  if (limpe > 1){
    /* Get work on all processors */
    procwork[0] = (float)nupdat;
    merge_array(procwork, nvals);
    twork = 0.0;
    if (mype ==0){
      for (npe = 0; npe < limpe; npe++)
	twork += procwork[npe];
    }
    /* Calculate  timings */
    minwall = big;
    maxwall = -1.0;
    for (npe = 0; npe < limpe; npe++){
      minwall = min(minwall, pcpu[npe]);
      maxwall = max(maxwall, pcpu[npe]);
    }
    dwall = global_max(dwall);
  }
  else{
    minwall = dwall;
    maxwall = dwall;
    twork = (float)nupdat;
    pcpu[0] = dwall;
  }
  /* prof - verbose timings */
  if (prof){
    procupdat[0] = nupdat - onupdat;
    if (limpe > 1){
      procworkwt[0] = workwt;
      merge_array(procupdat, nvals);
      merge_array(procworkwt, nvals);
    }
    if (mype == 0){
      cout << "profile for last step" << '\n';
      cout << "   max/min cpu time ratio " << maxwall/max(minwall,1.0e-5)
        << '\n';
      cout << "   Processor      Wall time        Updates      Speed " <<
        "      Workwt" << '\n';
      for (npe = 0; npe < limpe; npe++){
	cout.setf(ios::scientific);
	cout << setprecision(2);
	cout << setw(4) << "   " << npe << "            " << pcpu[npe] <<
          "     " <<setw(10) << procupdat[npe] << "    " << 
          procupdat[npe]/pcpu[npe] << "      " << procworkwt[npe] <<'\n';
      }
    }
  }
  if (mype == 0){
    cout.setf(ios::scientific);
    if (mugi){
    float mlcnv = -1.0;
    int milev = -1;
      for (ilev = 0; ilev < nlevs; ilev++){
	if (lcnv[ilev] > mlcnv){
	  mlcnv = lcnv[ilev];
	  milev = ilev;
	}
      }
      cout << "max fractional change " << mdq << " max cnv " <<
	mlcnv << " ilev " << milev << '\n';
    }
    tnupdat += (long long int)twork;
    tcpu += dwall;
    cout.setf(ios::scientific);
    cout << "wall clock speed " << twork/max(1.0e-5,dwall) << '\n';
    cout << "that took " << dwall << " secs time now " << t << '\n';
  }
  if (quit){
    /* Write dump file */
    strcpy(key, "dump");
    sprintf(tail, ".%e", t);
    strcat(key, tail);
    writem(key);
    report();
  }
  return;
}
//end march.c

//file integrate.c
/*
+  integrate -  integrate solution.
      Description  :
         This uses  step to integrate  the solution on level  ilev. It
         calls itself to do subgrids if they exist and handles merging
         and regridding.
*/
void integrate(int ilev)
{
  const int clocks = CLOCKS_PER_SEC;
  float ocpu;
  char key[SZNAM];
  /* Advance level by one step */
  step(ilev);
  /* Plot if admin */
  if (admin){
    strcpy(key, "P");
    plot(key);
  }
  if (ilev == nlevs-1){
    /* Finest level */
    if (grav || mugi){
      ocpu = (float)clock()/(float)clocks;
      mdcpu += mabs((float)clock()/(float)clocks) - mabs(ocpu);
    }
    return;
  }
  else{
    /* Check sublevels */
    while (lev[ilev+1]->nt < lev[ilev]->nt)
      integrate(ilev+1);
    if (correct)
      correct_flux(ilev);
    if (regrid)
      re_grid(ilev);
    merge(ilev);
    if (mugi || grav)
      /* Clean multigrid */
      clean_mug();
  }
  return;
}
//end integrate.c

//file adv_sol.c
/*
+  adv_sol -  advance solution.
      Description  :
         This advances the solution (adds cptr->dqdt to cptr->qa).  It
         also sets  up the  boundary interpolant at  internal boundary
         joins.
*/
void adv_sol(int ilev)
{
  int i, id;
  struct cell *cptr;
  struct join *jptr;
  struct cell_list *clist;
  if (limpe > 1)
    /* Communicate halo qa and dqdt (com3) */
    com_halo(com3, ilev, 0.0);
  if (ilev < nlevs){
    /* Set up internal boundary interpolant for internal joins */
    for (id = 1; id <= nd; id++){
      jptr = lev[ilev]->fjoin[id-1];
      while (jptr != NULL){
	if (jptr->bound == cibound)
	  set_join_int(ilev, id, jptr);
	jptr = jptr->nxt;
      }
    }
    if (limpe > 1){
      /* Communicate coarse/fine interpolant for halo joins */
      com_halo(com7, ilev, 0.0);
      if (lev[ilev]->fcell == NULL){
	/* Do halo nd joins */
	clist = plev[ilev]->lhalo;
	while (clist != NULL){
	  jptr = clist->cptr->rjoin[nd-1];
	  if (jptr->bound == cibound)
	    set_join_int(ilev, nd, jptr);
	  clist = clist->nxt;
	}
      }
    }
  }
  /* Advance solution */
  cptr = lev[ilev]->fcell;
  mdq = -1;
  while (cptr != NULL){
    nupdat++;
    for (i = 0; i < nq; i++){
      mdq = max(mdq,
        mabs(cptr->dqdt[i]/max(1.0, cptr->qa[i])));      
      cptr->qa[i] += cptr->dqdt[i];
    }
    /* Fix solution (qa) */
    usrfxq(cptr);
    cptr = cptr->nxt;
  }
  /* Set level viscous Courant number */
  mdq = global_max(mdq);
  lcnv[ilev] = max(0.9, min(mcnv, sqrt(fdq/max(tiny, mdq))*cnv));  
  return;
}
//end adv_sol.c

//file set_join_int.c
/*
+ set_join_int - setup join interpolant.
     Description :
        This sets up the interpolant at join jptr, which is the coarse
        join at  a coarse/fine boundary.  The interpolant  is put into
        the pa0, dpdt and dpdx arrays. It assumes that the solution at
        the old time  level is in qa and that the  solution at the new
        time level  is qa +  dqdt. It checks  for halo cells  and uses
        their gradients.
*/
void set_join_int(int ilev, int jid, join *jptr)
{
  int i, id, il;
  bool gotl, gotr, ggotl, ggotr;
  float divdx[MAXD], pl0[MAXQ], pr0[MAXQ], qal[MAXQ],qar[MAXQ], 
    pl[MAXQ], pr[MAXQ], gradl[MAXQ], gradr[MAXQ], ggradl[MAXQ], 
    ggradr[MAXQ];
  struct cell *left, *right;
  struct join *ljoin, *rjoin;
    /* Left/right pa in direction perpendicular to join */
  left = jptr->left;
  right = jptr->right;
  for (i = 0; i < nq; i++){
    qal[i] = left->qa[i];
    qar[i] = right->qa[i];
  }
  usrpq(qal, pl0);
  usrpq(qar, pr0);
  /* Central state at old time */
  for (i = 0; i < nq; i++)
    jptr->pa0[i] = 0.5*(pl0[i] + pr0[i]);
  if (nd > 1){
    /* Gradients at old time */
    for (id = 1; id <= nd; id++){
      divdx[id-1] = 0.5/lev[ilev]->dxbar[id-1];
      if (id != jid){
	/* Do left of internal boundary join */
	ljoin = left->ljoin[id-1];
	if (ljoin->left != NULL){
	  gotl = true;
	  /* Internal join -- get left gradient on left */
	  usrpq(ljoin->left->qa, pl);
	  for (i = 0; i < nq; i++)
	    gradl[i] = pl0[i] - pl[i];
	}
	else{
	  if (ljoin->bound == lebound){
	    /* Left external boundary join */
	    gotl = true;
	    usrbd(id, lebound, ljoin, pl, pl0);
	    for (i = 0; i < nq; i++)
	      gradl[i] = pl0[i] - pl[i];
	  }
	  else
	    /* Left internal boundary join */
	    gotl = false;
	}
	rjoin = left->rjoin[id-1];
	if (rjoin->right != NULL){
	  gotr = true;
	  /* Internal join -- get right gradient on left */
	  usrpq(rjoin->right->qa, pr);
	  for (i = 0; i < nq; i++)
	    gradr[i] = pr[i] - pl0[i];
	}
	else{
	  if (rjoin->bound == rebound){
	    /* Right external boundary join */
	    gotr = true;
	    usrbd(id, rebound, rjoin, pl0, pr);
	    for (i = 0; i < nq; i++)
	      gradr[i] = pr[i] - pl0[i];
	  }
	  else
	    /* Right internal boundary join */
	    gotr = false;
	}
	/* Final gradient on left */
	if (gotl && gotr){
	  ggotl = true;
	  ggotr = true;
	  for (i = 0; i < nq; i++)
	    ggradl[i] = av(gradl[i], gradr[i]);
	}
	else if (gotl){
	  ggotl = true;
	  ggotr = false;
	  for (i = 0; i < nq; i++)
	    ggradl[i] = gradl[i];
	}
	else if (gotr){
	  ggotr = true;
	  ggotl = false;
	  for (i = 0; i < nq; i++)
	    ggradl[i] = gradr[i];
	}
	else{
	  ggotl = false;
	  ggotr = false;
	}
	/* Do right of internal boundary join */
	ljoin = right->ljoin[id-1];
	if (ljoin->left != NULL){
	  /* Internal join -- get left gradient on right*/
	  gotl = true;
	  usrpq(ljoin->left->qa, pl);
	  for (i = 0; i < nq; i++)
	    gradl[i] = pr0[i] - pl[i];
	}
	else{
	  if (ljoin->bound == lebound){
	    /* Left external boundary join */
	    gotl = true;
	    usrbd(id, lebound, ljoin, pl, pr0);
	    for (i = 0; i < nq; i++)
	      gradl[i] = pr0[i] - pl[i];
	  }
	  else
	    /* Left internal boundary join */
	    gotl = false;
	}
	rjoin = right->rjoin[id-1];
	if (rjoin->right != NULL){
	  /* Internal join -- get right gradient on right */
	  gotr = true;
	  usrpq(rjoin->right->qa, pr);
	  for (i = 0; i < nq; i++)
	    gradr[i] = pr[i] - pr0[i];
	}
	else{
	  if (rjoin->bound == rebound){
	    /* Right external boundary join */
	    gotr = true;
	    usrbd(id, rebound, rjoin, pr0, pr);
	    for (i = 0; i < nq; i++)
	      gradr[i] = pr[i] - pr0[i];
	  }
	  else
	    /* Right internal boundary join */
	    gotr = false;
	}
	/* Final gradient on right */
	if (gotl && gotr){
	  ggotl = true;
	  ggotr = true;
	  for (i = 0; i < nq; i++)
	    ggradr[i] = av(gradl[i], gradr[i]);
	}
	else if (gotl){
	  ggotl = true;
	  ggotr = false;
	  for (i = 0; i < nq; i++)
	    ggradr[i] = gradl[i];
	}
	else if (gotr){
	  ggotl = false;
	  ggotr = true;
	  for (i = 0; i < nq; i++)
	    ggradr[i] = gradr[i];
	}
	else{
	  ggotl = false;
	  ggotr = false;
	}
	/* Put gradient in dpdx */
	il = dpdx_index(jid, id);
	if (ggotl && ggotr){
	  for (i = 0; i < nq; i++)
	    jptr->dpdx[i+il] = (ggradl[i] + ggradr[i])*divdx[id-1];
	}
	else if (ggotl){
	  for (i = 0; i < nq; i++)
	    jptr->dpdx[i+il] = 2.0*ggradl[i]*divdx[id-1];
	}
	else if (ggotr){
	  for (i = 0; i < nq; i++)
	    jptr->dpdx[i+il] = 2.0*ggradr[i]*divdx[id-1];
	}
	else{
	  for (i = 0; i < nq; i++)
	    jptr->dpdx[i+il] = 0.0;
	}
      }
    }
  }
  /* Central state time derivative */
  for (i = 0; i < nq; i++){
    qal[i] = left->qa[i] + left->dqdt[i];
    qar[i] = right->qa[i] + right->dqdt[i];
  }
  usrpq(qal, pl0);
  usrpq(qar, pr0);
  for (i = 0; i < nq; i++)
    jptr->dpdt[i] = 0.5*(pl0[i] + pr0[i]) - jptr->pa0[i];
  return;
}
//end set_join_int.c

//file cal_diffl.c
/*
+  cal_diffl - calculate left difference for primitives.
      Description  :
         Calculates left difference of primitive variables in the cell
         given  by  cptr  for  level  ilev.   Result  is  returned  in
         diffl[MAXQ].
*/
void cal_diffl(int ilev, int id, float tfrac, struct cell *cptr, 
  float diffl[])
{
  int i, il, tid;
  float pl[MAXQ], pr[MAXQ];
  struct join *ljoin, *pjptr;
  struct cell *left;
  /* Zero for safety */
  for (i = 0; i < nq; i++)
    diffl[i] = 0.0;
  ljoin = cptr->ljoin[id-1];
  left = ljoin->left;
  if (left != NULL){
    /* Left cell exists */
    for (i = 0; i < nq; i++)
      diffl[i] = cptr->pa[i] - left->pa[i];
  }
  else{
    /* Left needs boundary treatment */
    if (ljoin->bound == lebound){
      /* Left external boundary */
      for (i = 0; i < nq; i++)
	pr[i] = cptr->pa[i];
      usrbd(id, lebound, ljoin, pl, pr);
      for (i = 0; i < nq; i++)
	diffl[i] = pr[i] - pl[i];
    }
    else if (ljoin->bound == lfibound){
      /* Left internal boundary */
      pjptr = ljoin->par;
      for (i = 0; i < nq; i++)
	pl[i] = pjptr->pa0[i] + tfrac*pjptr->dpdt[i];
      for (tid = 1; tid <= nd; tid++){
	if (id != tid){
	  il = dpdx_index(id, tid);
	  for (i = 0; i < nq; i++)
	    pl[i] += (ljoin->xj[tid-1] - pjptr->xj[tid-1])*
	      pjptr->dpdx[i+il];
	}
      }
      for (i = 0; i < nq; i++)
	diffl[i] = 2.0*(cptr->pa[i] - pl[i]);
    }
  }
  return;
}
//end cal_diffl.c

//file cal_diffr.c
/*
+  cal_diffr -  calculate right difference for primitives.
      Description  :
         Calculates  right difference  of primitive  variables  in the
         cell given  by cptr  for level ilev.   Result is  returned in
         diffr[MAXQ].
*/
void cal_diffr(int ilev, int id, float tfrac, struct cell *cptr, 
  float diffr[])
{
  int i, il, tid;
  float pl[MAXQ], pr[MAXQ];
  struct join *rjoin, *pjptr;
  struct cell *right;
  /* Zero for safety */
  for (i = 0; i < nq; i++)
    diffr[i] = 0.0;
  rjoin = cptr->rjoin[id-1];
  right = rjoin->right;
  if (right != NULL){
    /* Right cell exists */
    for (i = 0; i < nq; i++)
      diffr[i] = right->pa[i] - cptr->pa[i];
  }
  else{
    /* Right needs special treatment */
    if (rjoin->bound == rebound){
      /* Right external boundary */
      for (i = 0; i < nq; i++)
	pl[i] = cptr->pa[i];
      usrbd(id, rebound, rjoin, pl, pr);
      for (i = 0; i < nq; i++)
	diffr[i] = pr[i] - pl[i];
    }
    else if (rjoin->bound == rfibound){
      /* Right internal boundary */
      pjptr = rjoin->par;
      for (i = 0; i < nq; i++)
	pr[i] = pjptr->pa0[i] + tfrac*pjptr->dpdt[i];
      for (tid = 1; tid <= nd; tid++){
	if (id != tid){
	  il = dpdx_index(id, tid);
	  for (i = 0; i < nq; i++)
	    pr[i] += (rjoin->xj[tid-1] - pjptr->xj[tid-1])*
	      pjptr->dpdx[i+il];
	}
      }
      for (i = 0; i < nq; i++)
	diffr[i] = 2.0*(pr[i] - cptr->pa[i]);
    }
  }
  return;
}
//end cal_diffr.c
