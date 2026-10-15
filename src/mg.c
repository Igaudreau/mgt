#undef MAIN
#undef USER
#define mcom1 1
#include <iostream>
#include<iomanip>
#include <sys/times.h>
#include<math.h>
#include <cstring>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "user.h"
#include "funcs.h"

extern void com_halo(int com, int ilev, float tfrac);
extern void zero_halo_cor(int ilev);
extern int dpdx_index(int jid, int id);
extern void adv_sol(int ilev);
extern void cal_diffl(int ilev, int id, float tfrac, struct cell *cptr, 
  float diffl[]);
extern void cal_diffr(int ilev, int id, float tfrac, struct cell *cptr, 
  float diffr[]);
extern struct cell *find_fcell(int &ilev, float xc[]);
extern struct part_list *mov_part(bool all, struct part_list *pmov);
extern void get_part_proc(struct part_list *start);
extern struct part_list *assign_part(struct part_list *start);
extern void rm_part(struct part *pptr);
extern void usrbd(int id, int bound, join *jptr, float pl[], 
float pr[]);
extern void usrsaf(int norder, float dt, cell *cptr, float s[]);
extern void usrpq(float qa[], float pa[]);
extern void usrqp(float pa[], float qa[]);
extern float usrss(float pa[]);
extern float usrte(float pa[]);
extern bool usrps(cell *cptr, float &s);
extern void usrfxp(float pa[]);
extern void usrmu(float pa[], float mu[]);
extern void usrpartpop(int ilev, float tcur);
extern void usrprate(struct part *pptr, float dt, float pa[], 
  float dvp[], float &den, float &dm);
extern float usrdmdt(float pmass, float pa[]);
extern float usrYE(float pa[]);
extern float usrAbar(float pa[]);
extern float usrein(float temp, float rho, float Abar, float Ye);
extern float usrpg2te(float tempi, float rho, float pg, float Abar,
  float Ye);
extern void usrreact(float pa[], float ydot[]);
extern float usrshockpos(bool verbose, int dir, float xn, float qa[],
  float pa[]);
extern float global_max(float value);
extern float global_sum(float value);
extern void mcom_halo(int mcom, int ilev);
extern void create_mug();
extern void init_mug();
extern void cal_grav();
extern void usrSootFoil(cell *cptr);
void step(int ilev);
void hflux(int dir, float pvl[], float pvr[], float f[]);
void rieman(float dl, float ul, float pl, float cl, float dr, float ur,
  float pr, float cr, float &d, float &u, float &p);
float wave(float p, float pi);
void vflux(int id, int ilev, struct join *jptr, float tfrac, float f[]);
void vel_diff(int id, struct join *jptr, float tfrac, float vdiff[]);
void apply_grav(float dt);
void evol_part(int ilev, float tcur, float dt, float tfrac, 
  struct cell *cptr, struct part_list *&pmov);

//file step.c
/*
+  step -  time step  solution.
      Description  : 
         This  determines  the  stable  timestep  and  integrates  the
         solution over  that timestep.   This done by  determining the
         fluxes  at  each  of  the  joins, incrementing  dqdt  in  the
         left/right cells  and then  using dqdt to  update qa  in each
         cell. It does 1, 2, 3 dimensions.
*/
void step(int ilev)
{
  bool first;
  (void)first;
  int i, il, id, tid, ndtc=0, ndt, ntf=0, ntc=0, maxjaccall, maxreactcall,
    nreactcall, maxnsub, ncell;
  float rate, hrate, vrate, grate, ss, divdx=0, divrho, tv, dt, f[MAXQ], 
    dth, pl[MAXQ], pr[MAXQ], diffl[MAXQ], diffr[MAXQ], tfrac=0, area, 
    dtar, mu[MAXQ], so[MAXQ], rho, tcur=0, srate, ydot[MAXSPECIES+1],
    spos;
  (void)rho;
  struct cell *cptr, *left, *right;
  struct join *jptr, *pjptr;
  struct cell_list *clist;
  struct part *pptr;
  struct part_list *pmov;
  (void)pptr;
  /* Check for atmosphere */
  spos = -big;
  if (strncmp(cas, "ATMOS", 5) == 0){
    /* Find shockpos */
    spos = usrshockpos(false, 1, 0, pl, pr);
    if (debug)
      cout << "ilev " << ilev << " spos " << spos << '\n';
  }
  /* Store current level in curlev */
  curlev = ilev;
  /* Zero gsl ode counters and set max nsub */
  maxnsub = -1;
  maxjaccall = -1;
  maxreactcall = -1;
  nreactcall = 0;
  ncell = 0;
  /* Calculate primitives and rate etc.*/
  hrate = -1.0;
  vrate = -1.0;
  grate = -1.0;
  srate = -1.0;
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    /* Count cells */
    ncell++;
    /* Primitives */
    usrpq(cptr->qa, cptr->pa);
    /* Zero dqdt */
    for (i = 0; i < nq; i++)
      cptr->dqdt[i] = 0.0;
    ss = usrss(cptr->pa);
    for (id = 1; id <= nd; id++){
      divdx = cptr->rjoin[id-1]->area/cptr->vol;
      hrate = max(hrate, (mabs(cptr->pa[iqu0+id]) + ss)*divdx);
    }
    if (react){
      usrreact(cptr->pa, ydot);
      /* Isotope rate */      
      for (i = 0; i < niso; i++)
	srate = max(srate, fabs(Ai[i]*ydot[i]));
      /* Energy rate */
      srate = max(srate, fabs(ydot[niso]/cptr->pa[iqte]));
    }
    if (viscous){
      divrho = 1.0/cptr->qa[iqd];
      tv = -1.0;
      usrmu(cptr->pa, mu);
      for (i = 0; i <= iqal0 + scad; i++)
	tv = max(tv, mu[i]);
      tv *= divrho;
      for (i = iqal0 + scad + 1; i <= iqal0 + nal; i++)
	tv = max(tv, mu[i]);
      vrate = max(vrate, 4.0*tv*divdx*divdx);
      /* Zero velocity differences */
      for (id = 0; id < nd; id++){
	for (i = 0; i < nd; i++)
	  cptr->dudx[i][id] = 0.0;
      }
    }
    if (grav)
      grate = max(grate, sqrt(gconst*cptr->qa[iqd]));
    cptr = cptr->nxt;
  }
  /* Divide rates by cn, cnv, cng and calculate global rates */
  hrate /= cn;
  vrate /= cnv;
  grate /= cng;
  srate /= cns;
  rate = max(hrate, vrate);
  rate = max(hrate, srate);  
  rate = max(rate, grate);
  if (limpe > 1){
    rate = global_max(rate);
    hrate = global_max(hrate);
    vrate = global_max(vrate);
    srate = global_max(srate);    
    grate = global_max(grate);
  }
  /* Check if level exists */
  if (rate < 0.0){
    /* No cells on this level */
    lev[ilev]->ndt = lev[ilev-1]->ndt;
    lev[ilev]->nt += lev[ilev-1]->ndt;
    return;
  }
  /* Timestep */
  if (ilev == levc - 1){
    /* Base timestep -- set base timestep and all times */
    dt0 = min(target - t, cnf/rate);
    dt = dt0;
    ndt = ipow(MAXLEV);
    lev[levc-1]->ndt = ndt; 
    for (i = 0; i < nlevs; i++)
      lev[i]->nt = 0;
  }
  else{
    /* Timestep for fine level */
    ndtc = lev[ilev-1]->ndt;
    ndt = ndtc;
    dt = (float)(dt0*ndt/lev[levc-1]->ndt);
    while(dt > 1.0/rate){
      dt *=0.5;
      ndt = ndt >> 1;
    }
    /* Coarse and fine level times etc */
    ntc = lev[ilev-1]->nt - ndtc;
    ntf = lev[ilev]->nt;
    if (ntc + ndtc < ntf + ndt){
      ndt = ntc + ndtc - ntf;
      dt = (float)(dt0*ndt/lev[levc-1]->ndt);
    }
    lev[ilev]->ndt = ndt;
    /* Interpolation time */
    tfrac = (float)(ntf - ntc)/ndtc;
  }
  if (particle){
    /* Calculate current time for lagrangian particles */
    tcur = t + dt0*(float)lev[ilev]->nt/(float)lev[levc-1]->ndt;
    /* Set start of particle move lists to NULL */
    pmov = NULL;
    /* Update particle population */
    usrpartpop(ilev, tcur);
  }
  /* First order time step */
  if (order == 2)
    dth = 0.5*dt;
  else
    dth = dt;
  if (limpe > 1){
  /* Zero flux correction for halo joins */
    zero_halo_cor(ilev);
    /* Communicate halo pa (com2) */
    com_halo(com2, ilev, tfrac);
  }
  if ((ilev == alev) && grav){
    /* Do gravity */
    if (nomug){
      create_mug();
      cal_grav();
    }
    apply_grav(0.5*dt);
  }
  if (!hydrostatic){
    /* For each direction */
    for (id = 1; id <= nd; id++){
      /*Calculate fluxes and add to dqdt */
      jptr = lev[ilev]->fjoin[id-1];
      while (jptr != NULL){
	if (jptr->bound < lebound){
	  /* No special boundary treatment */
	  left = jptr->left;
	  right = jptr->right;
	  for (i = 0; i < nq; i++){
	    pr[i] = right->pa[i];
	    pl[i] = left->pa[i];
	  }
	  hflux(id, pl, pr, f);
	  area = jptr->area;
	  if (viscous)
	    vflux(id, ilev, jptr, tfrac, f);
	  for (i = 0; i < nq; i++){
	    left->dqdt[i] -=  f[i]*area;
	    right->dqdt[i] +=  f[i]*area;
	  }
	  if (jptr->bound == cibound){
	    if (order == 1){
	      /* Coarse internal boundary -- put flux in cor */
	      dtar = dt*area;
	      /* Set count for flux corrections (parallel) */
	      jptr->done = 1;
	      for (i = 0; i < nq; i++)
		jptr->cor[i] = -f[i]*dtar;
	    }
	  }
	}
	else{
	  /* Needs boundary treatment */
	  if (jptr->bound == lebound){
	    /* Left external boundary */
	    right = jptr->right;
	    for (i = 0; i < nq; i++)
	      pr[i] = right->pa[i];
	    usrbd(id, lebound, jptr, pl, pr);
	    hflux(id, pl, pr, f);
	    area = jptr->area;
	    if (viscous)
	      vflux(id, ilev, jptr, tfrac, f);
	    for (i = 0; i < nq; i++)
	      right->dqdt[i] +=  f[i]*area;
	  }
	  else if (jptr->bound == lfibound){
	    /* Left internal boundary */
	    right = jptr->right;
	    pjptr = jptr->par;
	    for (i = 0; i < nq; i++){
	      pr[i] = right->pa[i];
	      pl[i] = pjptr->pa0[i]  + tfrac*pjptr->dpdt[i];
	    }
	    for (tid = 1; tid <= nd; tid++){
	      if (id != tid){
		il = dpdx_index(id, tid);
		for (i = 0; i < nq; i++)
		  pl[i] += (jptr->xj[tid-1] - pjptr->xj[tid-1])*
		    pjptr->dpdx[i+il];
	      }
	    }
	    for (i = 0; i < nq; i++)
	      pl[i] = 2.0*pl[i] - pr[i];
	    hflux(id, pl, pr, f);
	    area = jptr->area;
	    if (viscous)
	      vflux(id, ilev, jptr, tfrac, f);
	    for (i = 0; i < nq; i++)
	      right->dqdt[i] +=  f[i]*area;
	    if (order == 1){
	      /* Fine internal boundary -- put flux in cor */
	      dtar = dt*area;
	      /* Set count for flux corrections (parallel) */
	      pjptr->done++;
	      for (i = 0; i < nq; i++)
		pjptr->cor[i] += f[i]*dtar;
	    }
	  }
	  if (jptr->bound == rebound){
	    /* Right external boundary join */
	    left = jptr->left;
	    for (i = 0; i < nq; i++)
	      pl[i] = left->pa[i];
	    usrbd(id, rebound, jptr, pl, pr);
	    hflux(id, pl, pr, f);
	    area = jptr->area;
	    if (viscous)
	      vflux(id, ilev, jptr, tfrac, f);
	    for (i = 0; i < nq; i++)
	      left->dqdt[i] -=  f[i]*area;
	  }
	  else if (jptr->bound == rfibound){
	    /* Right internal boundary */
	    left = jptr->left;
	    pjptr = jptr->par;
	    for (i = 0; i < nq; i++){
	      pl[i] = left->pa[i];
	      pr[i] = pjptr->pa0[i] + tfrac*pjptr->dpdt[i];
	    }
	    for (tid = 1; tid <= nd; tid++){
	      if (id != tid){
		il = dpdx_index(id, tid);
		for (i = 0; i < nq; i++)
		  pr[i] += (jptr->xj[tid-1] - pjptr->xj[tid-1])*
		    pjptr->dpdx[i+il];
	      }
	    }
	    for (i = 0; i < nq; i++)
	      pr[i] = 2.0*pr[i] - pl[i];
	    hflux(id, pl, pr, f);
	    area = jptr->area;
	    if (viscous)
	      vflux(id, ilev, jptr, tfrac, f);
	    for (i = 0; i < nq; i++)
	      left->dqdt[i] -=  f[i]*area;
	    if (order == 1){
	      /* Fine internal boundary -- put flux in cor */
	      dtar = area*dt;
	      /* Set count for flux corrections (parallel) */
	      pjptr->done++;
	      for (i = 0; i < nq; i++)
		pjptr->cor[i] += f[i]*dtar;
	    }
	  }
	}
	jptr = jptr->nxt;
      }
    }
  }
  /* Add sources, put solution change in  dqdt */
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    /* Calculate sources */
    jaccall = 0;
    reactcall = 0;    
    usrsaf(1, dth, cptr, so);
    if (order == 1){
      maxnsub = max(maxnsub, nsub);
      maxjaccall = max(maxjaccall, jaccall);
      maxreactcall = max(maxreactcall, reactcall);
      nreactcall += reactcall;
    }    
    tv = dth/cptr->vol;
    if (spos + par[9] < 0.0){
      for (i = 0; i < nq; i++)
	cptr->dqdt[i] = tv*cptr->dqdt[i] + dth*so[i];
    }
    else{
      for (i = 0; i < nq; i++)
	cptr->dqdt[i] = 0.0;
    }
    /* Special step control for STEADY model */
    if (strncmp(cas, "STEADY", 6) == 0){
      /* Check for near sonic point */
      float oss, omach, nmach, dmach, pa[MAXQ], qa[MAXQ];
      oss = usrss(cptr->pa);
      if (cptr->pa[iqu0+1] > 0.9*oss){      
	for (i = 0; i < nq; i++)
	  qa[i] = cptr->qa[i] + cptr->dqdt[i];
	omach = cptr->pa[iqu0+1]/oss;
	usrpq(qa, pa);	
	nmach = pa[iqu0+1]/usrss(pa);
	dmach  = nmach - omach;
	if (dmach > 0.5*(1.0 - nmach)){
	  tv = 0.001*(1.0 - nmach)/dmach;
	  for (i = 0; i < nq; i++)
	    cptr->dqdt[i] *= tv;
	  dt0 *= tv;
	  cn *=tv;
	  if (omach + dmach*tv > 0.9999)
	    target = t;
	}
      }
    }
    if ((particle) && (order == 1)){
      /* Evolve lagrangian particles */
      if (cptr->pptr != NULL)
	evol_part(ilev, tcur, dt, tfrac, cptr, pmov);
    }    
    cptr = cptr->nxt;
  }
  if (order == 1){
    /* Advance solution to full time and return */
    adv_sol(ilev);
    if (particle){
      /* Assign local particles */
      pmov = assign_part(pmov);
      if (limpe > 1){
	/* Set destination processors for non-local particles */
	get_part_proc(pmov);
	/* Move particles and assign to cells */
	pmov = mov_part(false, pmov);
	if (assign_part(pmov) != NULL)
	  cout << "mg pmov not NULL on proc " << mype << '\n';
      }
    }
    lev[ilev]->nt += ndt;
    if ((ilev == alev) && grav){
      /* Do gravity */
      cal_grav();
      apply_grav(0.5*dt);
    }
    if (debug){
      if (limpe > 1){
	maxjaccall = global_max(maxjaccall);
	maxreactcall = global_max(maxreactcall);	
	maxnsub = global_max(maxnsub);
	nreactcall = global_sum(nreactcall);
	ncell = global_sum(ncell);	
      }
      if (mype == 0){      
	cout.setf(ios::scientific);
	cout << "level " << ilev << setprecision(4) << " dt "
	     << dt << " hyp limit " << 1.0/hrate
	     << " so limit " << 1.0/srate;            
	if (viscous)
	  cout << " vis limit " << setprecision(4) << 
	    1.0/max(vrate, 1.0e-10);
	if (grav)
	  cout << " grav limit " << setprecision(4) << 
	    1.0/max(grate, 1.0e-10);
	if (ilev == alev)
	  cout << " time frac " << 
	    (float)lev[ilev]->nt/(float)lev[levc-1]->ndt;
	cout << '\n';
	cout << "maxjacall " << maxjaccall <<
          " maxreactcall " << maxreactcall <<" maxnsub " << maxnsub <<
	  " nreactcall " << nreactcall << " avreactcall "<<
	  nreactcall/ncell << '\n';
      }
    }
    return;
  }
  /* Second order Step */
  /* Second order interpolation time */
  if (ilev > levc)
    tfrac = (float)(ntf + ndt/2 - ntc)/ndtc;
  /* Zero dqdt, dudx */
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    /* Primitives */
    for (i = 0; i < nq; i++)
      cptr->dqdt[i] += cptr->qa[i];
    usrpq(cptr->dqdt, cptr->pa);
    /* Zero dqdt */
    for (i = 0; i < nq; i++)
      cptr->dqdt[i] = 0.0;
    /* Zero velocity differences */
    for (id = 0; id < nd; id++){
      for (i = 0; i < nd; i++)
	cptr->dudx[i][id] = 0.0;
    }
    cptr = cptr->nxt;
  }
  if (limpe > 1){
    /* Communicate halo pa, dpdx (com2, com4) */
    com_halo(com2, ilev, tfrac);
    com_halo(com4, ilev, tfrac);
  }
  if (!hydrostatic){
    /* For each direction */
    for (id = 1; id <= nd; id++){
      /*Calculate average differences for local cells */
      cptr = lev[ilev]->fcell;
      while (cptr != NULL){
	/* Calculate left/right differences */
	cal_diffl(ilev, id, tfrac, cptr, diffl);
	cal_diffr(ilev, id, tfrac, cptr, diffr);
	/* Average difference */
	for (i = 0; i < nq; i++)
	  cptr->diff[i] = av(diffr[i], diffl[i]);
	cptr = cptr->nxt;
      }
      if ((limpe > 1) && (id == nd)){
	if ((mype != 0) || (bcl[nd-1] == bcper)){
	  /*Calculate average differences for left halo cells */
	  clist = plev[ilev]->lhalo;
	  while (clist != NULL){
	    cptr = clist->cptr;
	    cal_diffr(ilev, id, tfrac, cptr, diffr);
	    for (i = 0; i < nq; i++)
	      cptr->diff[i] = av(diffr[i], cptr->diff[i]);
	    clist = clist->nxt;
	  }
	}
	if ((mype != limpe-1) || (bcl[nd-1] == bcper)){
	  /*Calculate average differences for right halo cells */
	  clist = plev[ilev]->rhalo;
	  while (clist != NULL){
	    cptr = clist->cptr;
	    cal_diffl(ilev, id, tfrac, cptr, diffl);
	    for (i = 0; i < nq; i++)
	      cptr->diff[i] = av(diffl[i], cptr->diff[i]);
	    clist = clist->nxt;
	  }
	}
      }
      /* Calculate fluxes and increment dqdt */
      jptr = lev[ilev]->fjoin[id-1];
      divdx = 1.0/lev[ilev]->dxbar[id-1];
      while (jptr != NULL){
	if (jptr->bound < lebound){
	  /* No special boundary treatment */
	  left = jptr->left;
	  right = jptr->right;
	  for (i = 0; i < nq; i++){
	    pl[i] = left->pa[i] + 0.5*left->diff[i];
	    pr[i] = right->pa[i] - 0.5*right->diff[i];
	  }
	  hflux(id, pl, pr, f);
	  area = jptr->area;
	  if (viscous)
	    vflux(id, ilev, jptr, tfrac, f);
	  for (i = 0; i <  nq; i++){
	    left->dqdt[i] -=  f[i]*area;
	    right->dqdt[i] +=  f[i]*area;
	  }
	  if (jptr->bound == cibound){
	    /* Coarse internal boundary -- put flux in cor */
	    dtar = dt*area;
	    /* Set count for flux corrections (parallel) */
	    jptr->done = 1;
	    for (i = 0; i < nq; i++)
	      jptr->cor[i] = -f[i]*dtar;
	  }
	}
	else{
	  /* Needs boundary treatment */
	  if (jptr->bound == lebound) {
	    /* Left external boundary*/
	    right = jptr->right;
	    for (i = 0; i < nq; i++)
	      pr[i] = right->pa[i] - 0.5*right->diff[i];
	    usrbd(id, lebound, jptr, pl, pr);
	    if (bcl[id-1] == bcfix)
	      for (i = 0; i < nq; i++)
		pl[i] += 0.5*right->diff[i];
	    hflux(id, pl, pr, f);
	    area = jptr->area;
	    if (viscous)
	      vflux(id, ilev, jptr, tfrac, f);
	    for (i = 0; i < nq; i++)
	      right->dqdt[i] +=  f[i]*area;
	  }
	  else if (jptr->bound == lfibound){
	    /* Left internal boundary */
	    right = jptr->right;
	    pjptr = jptr->par;
	    for (i = 0; i < nq; i++)
	      pl[i] = pjptr->pa0[i] + tfrac*pjptr->dpdt[i];
	    for (tid = 1; tid <= nd; tid++){
	      if (id != tid){
		il = dpdx_index(id, tid);
		for (i = 0; i < nq; i++)
		  pl[i] += (jptr->xj[tid-1] - pjptr->xj[tid-1])*
		    pjptr->dpdx[i+il];
	      }
	    }
	    for (i = 0; i < nq; i++)
	      pr[i] = right->pa[i] - 0.5*right->diff[i];
	    hflux(id, pl, pr, f);
	    area = jptr->area;
	    if (viscous)
	      vflux(id, ilev, jptr, tfrac, f);
	    for (i = 0; i < nq; i++)
	      right->dqdt[i] +=  f[i]*area;
	    /* Fine internal boundary -- put flux in cor */
	    dtar = dt*area;
	    /* Set count for flux corrections (parallel) */
	    pjptr->done++;
	    for (i = 0; i < nq; i++)
	      pjptr->cor[i] += f[i]*dtar;
	  }
	  else if (jptr->bound == rebound){
	    /* Right external boundary */
	    left = jptr->left;
	    for (i = 0; i < nq; i++)
	      pl[i] = left->pa[i] + 0.5*left->diff[i];
	    usrbd(id, rebound, jptr, pl, pr);
	    if (bcr[id-1] == bcfix)
	      for (i = 0; i < nq; i++)
		pr[i] -= 0.5*left->diff[i];
	    hflux(id, pl, pr, f);
	    area = jptr->area;
	    if (viscous)
	      vflux(id, ilev, jptr, tfrac, f);
	    for (i = 0; i < nq; i++)
	      left->dqdt[i] -=  f[i]*area;
	  }
	  else if (jptr->bound == rfibound){
	    /* Right internal boundary */
	    left = jptr->left;
	    pjptr = jptr->par;
	    for (i = 0; i < nq; i++)
	      pr[i] = pjptr->pa0[i] + tfrac*pjptr->dpdt[i];
	    for (tid = 1; tid <= nd; tid++){
	      if (id != tid){
		il = dpdx_index(id, tid);
		for (i = 0; i < nq; i++)
		  pr[i] += (jptr->xj[tid-1] - pjptr->xj[tid-1])*
		    pjptr->dpdx[i+il];
	      }
	    }
	    for (i = 0; i < nq; i++)
	      pl[i] = left->pa[i] + 0.5*left->diff[i];
	    hflux(id, pl, pr, f);
	    area = jptr->area;
	    if (viscous)
	      vflux(id, ilev, jptr, tfrac, f);
	    for (i = 0; i < nq; i++)
	      left->dqdt[i] -=  f[i]*area;
	    /* Fine internal boundary -- put flux in cor */
	    dtar = dt*area;
	    /* Set count for flux corrections (parallel) */
	    pjptr->done++;
	    for (i = 0; i < nq; i++)
	      pjptr->cor[i] += f[i]*dtar;
	  }
	}
	jptr = jptr->nxt;
      }
    }
  }
  /* Add sources, put solution change in  dqdt */
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    /* Calculate sources */
    jaccall = 0;
    reactcall = 0;        
    usrsaf(2, dt, cptr, so);
    maxnsub = max(maxnsub, nsub);
    maxjaccall = max(maxjaccall, jaccall);
    maxreactcall = max(maxreactcall, reactcall);
    nreactcall += reactcall;    
    tv = dt/cptr->vol;
    if (spos + par[9] < 0.0){
      for (i = 0; i < nq; i++)
	cptr->dqdt[i] = tv*cptr->dqdt[i] + dth*so[i];
    }
    else{
      for (i = 0; i < nq; i++)
	cptr->dqdt[i] = 0.0;
    }
    if (particle){
      /* Evolve lagrangian particles */
      if (cptr->pptr != NULL)
	evol_part(ilev, tcur, dt, tfrac, cptr, pmov);
    }
    cptr = cptr->nxt;
  }
  /* Advance solution to full time */
  adv_sol(ilev);
  if (particle){
    /* Assign local particles */
    pmov = assign_part(pmov);
    if (limpe > 1){
      /* Set destination processors for non-local particles */
      get_part_proc(pmov);
      /* Move particles and assign to cells */
      pmov = mov_part(false, pmov);
      if (assign_part(pmov) != NULL)
	cout << "mg pmov not NULL on proc " << mype << '\n';
    }
  }
  lev[ilev]->nt += ndt;
  if ((ilev == alev) && grav){
    /* Do gravity */
    cal_grav();
    apply_grav(0.5*dt);
  }
  if (debug){
    if (limpe > 1){    
      maxjaccall = global_max(maxjaccall);
      maxreactcall = global_max(maxreactcall);      
      maxnsub = global_max(maxnsub);
      nreactcall = global_sum(nreactcall);
      ncell = global_sum(ncell);	      
    }
    if (mype == 0){      
      cout.setf(ios::scientific);
      cout << "level " << ilev << setprecision(4) << " dt "
	   << dt << " hyp limit " << 1.0/hrate
	   << " so limit " << 1.0/srate;            
      if (viscous)
	cout << " vis limit " << setprecision(4) << 
	  1.0/max(vrate, 1.0e-10);
      if (grav)
	cout << " grav limit " << setprecision(4) << 
	  1.0/max(grate, 1.0e-10);
      if (ilev == alev)
	cout << " time frac " << 
	  (float)lev[ilev]->nt/(float)lev[levc-1]->ndt;
      cout << '\n';
      cout << "maxjaccall " << maxjaccall <<
	" maxreactcall " << maxreactcall <<" maxnsub " << maxnsub <<
	  " nreactcall " << nreactcall << " avreactcall "<<
	  nreactcall/ncell << '\n';      
    }
  }
  if (sootfoil) {
    cptr = lev[ilev]->fcell;
    while (cptr != NULL) {
      usrSootFoil(cptr);
      cptr = cptr->nxt;
    }
  }
  return;
}
//end step.c

//file hflux.c
/*
+ hflux - Fluxes from Riemann problem.
     Description :
        This routine calculates the fluxes at an interface with pvl on
        left, pvr on  right.  Uses the linear solver  if god=true, hll
        if hll=true, Kurganov-Tadmor if kt=true.
*/
void hflux(int dir, float pvl[], float pvr[], float f[])
{
  int i, id, iqu;
  float tel, ter, cel, cer, cl, cr, pl, pr, dl, dr, p, u, d, te, ud, ul,
    ur, ev, pi[MAXQ], Abar, Ye, vel[MAXQ], fl[MAXQ], fr[MAXQ], qvl[MAXQ],
    qvr[MAXQ], sl, sr, s, mflux, cvis, vis, tv;
  /* Zero all fluxes (safety) */
  for (i = 0; i < nq; i++){
    f[i] = 0;
    fl[i] = 0;
    fr[i] = 0;
  }
  /* Calculate pressure,  density, velocity */
  iqu = iqu0 + dir;
  if (fxp){
    usrfxp(pvr);
    usrfxp(pvl);
  }
  else{
    pvl[iqd] = max(dfloor,pvl[iqd]);
    pvr[iqd] = max(dfloor,pvr[iqd]);
    pvl[iqe] = max(pfloor,pvl[iqe]);
    pvr[iqe] = max(pfloor,pvr[iqe]);
  }
  /* Left state */
  ul = pvl[iqu];
  dl = pvl[iqd];
  pl = pvl[iqe];
  tel = pvl[iqte];
  cel = usrss(pvl);
  cl = dl*cel;
  /* Right state */
  ur = pvr[iqu];
  dr = pvr[iqd];
  pr = pvr[iqe];
  ter = pvr[iqte]; 
  cer = usrss(pvr);
  cr = dr*cer;
  if (god){
    /* Linear Riemann solver */
    p = (cr*pl + cl*pr - cr*cl*(ur - ul))/(cr + cl);
    u = (cr*ur + cl*ul - (pr - pl))/(cr + cl);
    if (u >= 0.0){
      if ((ul - cel) > 0.0){
	d = dl;
	u = ul;
	p = pl;
      }
      else
	d = cl*dl/(cl + (u - ul)*dl);
    }
    else{
      if ((ur + cer) <= 0.0){
	d = dr;
	u = ur;
	p = pr;
      }
      else
	d = cr*dr/(cr - (u - ur)*dr);
    }
    /* Transverse velocities and scalars in resolved state */
    if (u > 0.0){
      for (id = 1; id <= nd; id++)
	vel[id-1] = pvl[iqu0+id];
      for (i = iqal0 + 1; i <= iqal0 + scad; i++)
	pi[i] = pvl[i];
    }
    else{
      for (id = 1; id <= nd; id++)
	vel[id-1] = pvr[iqu0+id];
      for (i = iqal0 + 1; i <= iqal0 + scad; i++)
	pi[i] = pvr[i];
    }
    vel[iqu-1] = u;
    p = max(pfloor, p);
    /* Mass flux and energy flux */
    mflux = d*u;
    f[iqd] = mflux;
    f[iqe] = p;
    /* Add velocity bits into momentum fluxes and energy */
    for (id = 1; id <= nd; id++){
      i = iqu0 + id;
      ud = vel[id-1];
      f[i] = mflux*ud;
      f[iqe] = f[iqe] + 0.5*d*ud*ud;
    }
    /* Supernova  EOS */      
    Abar = usrAbar(pi);
    Ye = usrYE(pi);      
    te = usrpg2te(0.5*(tel + ter), d, p, Abar, Ye);
    ev = usrein(te, d, Abar, Ye);
    /* Add this into energy flux, add pressure into normal momentum */
    f[iqe] = u*(f[iqe] + ev);
    f[iqu] = f[iqu] + p;
    /* k-epsilon */
    if (keps){
      /* Fluxes of turbulence variables */
      if (u > 0){
	for (i = iqk; i <= iqeps; i++)
	  f[i] = mflux*pvl[i];
	tv = 2.0*d*pvl[iqk]/3.0;
	f[iqu] += tv;
	f[iqe] += tv*u;
      }
      else{
	for (i = iqk; i <= iqeps; i++)
	  f[i] = mflux*pvr[i];
	tv = 2.0*d*pvr[iqk]/3.0;
	f[iqu] += tv;
	f[iqe] += tv*u;
      }
    }
    /* Scalars and internal energy */
    if (u > 0){
      for (i = iqal0 + 1; i <= iqal0 + scad; i++)
	f[i] = mflux*pvl[i];
      if (inten)
	/* Advect internal energy */
	f[iqie] = mflux*pvl[iqie];
    }
    else{
      for (i = iqal0 + 1; i <= iqal0 + scad; i++)
	f[i] = mflux*pvr[i];
      if (inten)
	/* Advect internal energy */
	f[iqie] = mflux*pvr[iqie];
    }
  }
  else{
    /* HLL or KT - Calculate left/right conservatives */
    usrqp(pvl, qvl);
    usrqp(pvr, qvr);
    /* Calculate left/right wavespeeds */
    cel = usrss(pvl);
    cer = usrss(pvr);
    sl = pvl[iqu] - cel;
    sr = pvr[iqu] + cer;
    /* Calculate left flux */
    u = pvl[iqu];
    p = pvl[iqe];
    d = pvl[iqd];
    mflux = d*u;
    /* Mass flux */
    fl[iqd] = mflux;
    for (id = 1; id <= nd; id++)
      vel[id-1] = pvl[iqu0+id];
    vel[iqu-1] = u;
    /* Momentum flux */
    fl[iqu] = p;
    for (id = 1; id <= nd; id++){
      i = iqu0 + id;
      fl[i] += mflux*vel[id-1];
    }
    /* Energy flux */
    fl[iqe] = u*(qvl[iqe] + p);
    if (keps){
      /* Fluxes of turbulence variables */
      for (i = iqk; i <= iqeps; i++)
	fl[i] = mflux*pvl[i];
      tv = 2.0*d*pvl[iqk]/3.0;
      fl[iqu] += tv;
      fl[iqe] += tv*u;
    }
    /* Scalars and internal energy */
    for (i = iqal0 + 1; i <= iqal0 + scad; i++)
      fl[i] = mflux*pvl[i];
    if (inten)
      /* Advect internal energy */
      fl[iqie] = mflux*pvl[iqie];
    /* Calculate right flux */
    u = pvr[iqu];
    p = pvr[iqe];
    d = pvr[iqd];
    mflux = d*u;
    /* Mass flux */
    fr[iqd] = mflux;
    for (id = 1; id <= nd; id++)
      vel[id-1] = pvr[iqu0+id];
    vel[iqu-1] = u;
    /* Momentum flux */
    fr[iqu] = p;
    for (id = 1; id <= nd; id++){
      i = iqu0 + id;
      fr[i] += mflux*vel[id-1];
    }
    /* Energy flux */
    fr[iqe] = u*(qvr[iqe] + p);
    if (keps){
      /* Fluxes of turbulence variables */
      for (i = iqk; i <= iqeps; i++)
	fr[i] = mflux*pvr[i];
      tv = 2.0*d*pvr[iqk]/3.0;
      fr[iqu] += tv;
      fr[iqe] += tv*u;
    }
    /* Scalars and possibly internal energy */
    for (i = iqal0 + 1; i <= iqal0 + scad; i++)
      fr[i] = mflux*pvr[i];
    if (inten)
      /* Internal energy */
      fr[iqie] = mflux*pvr[iqie];
    if (hll){
      /* HLL fluxes */
      if (sl > 0.0){
	/* Take left flux */
	for ( i = 0; i <= iqal0 + scad; i++)
	  f[i] = fl[i];
	if (inten)
	  /* Internal energy flux */
	  f[iqie] = fl[iqie];
      }
      else if (sr < 0.0){
	/* Take right flux */
	for ( i = 0; i <= iqal0 + scad; i++)
	  f[i] = fr[i];
	if (inten)
	  /* Internal energy flux */
	  f[iqie] = fr[iqie];
      }
      else{
	/* Take HLL flux */
	for ( i = 0; i <= iqal0 + scad; i++)
	  f[i] = (sr*fl[i] - sl*fr[i] + sr*sl*(qvr[i] - qvl[i]))
	    /(sr - sl);
	if (inten)
	  /* Internal energy flux */
	  f[iqie] = (sr*fl[iqie] - sl*fr[iqie] + 
	    sr*sl*(qvr[iqie] - qvl[iqie]))/(sr - sl);
      }
    }
    else if (kt){
      /* Kurganov-Tadmor fluxes */
      s = max(mabs(sl), mabs(sr));
      for ( i = 0; i <= iqal0+ scad; i++)
	f[i] = 0.5*(fr[i] + fl[i] - s*(qvr[i] - qvl[i]));
      if (inten)
	/* Internal energy flux */
	f[iqie] = 0.5*(fr[iqie] + fl[iqie] - 
		       s*(qvr[iqie] - qvl[iqie]));
    }
  }
  /* Artificial dissipation */
  cvis = 2.0/(1.0/(pvl[iqd]*cel) + 1.0/(pvr[iqd]*cer));
  vis = avisp*cvis;
  for (id = 1; id <= nd; id++){
    i = iqu0+ id;
    tv = vis*(pvl[i] - pvr[i]);
    f[i] += tv;
    f[iqe] +=  0.5*tv*(pvl[i] + pvr[i]);
  }
  f[iqe] += avise*cvis*(usrte(pvl) - usrte(pvr));
  return;
}
//end hflux.c

//file vflux.c
/*
+ vflux - Viscous fluxes.
     Description :
        This  routine calculates  the viscous  fluxes at  an interface
        pointed to by  jptr. These are added to  the hyperbolic fluxes
        in f[]. It is invoked if viscous is true.

*/
void vflux(int id, int ilev, struct join *jptr, float tfrac, float f[])
{
  int i, tid, il;
  float divdx[MAXD], vcof, pl[MAXQ], pr[MAXQ], pi[MAXQ], 
    dv[MAXD][MAXD], vdiffl[MAXD], vdiffr[MAXD], divv, fv[MAXD], 
    mu[MAXQ], mul[MAXQ], mur[MAXQ];
  struct cell *left, *right;
  struct join *pjptr;
  /* Inverse mesh spacings */
  for (tid = 1; tid <= nd; tid++)
    divdx[tid-1] = 1.0/lev[ilev]->dxbar[tid-1];
  /* Left,  right and interface states */
  if (jptr->bound < lebound){
    /* No special boundary treatment */
    left = jptr->left;
    right = jptr->right;
    for (i = 0; i < nq; i++){
      pl[i] = left->pa[i];
      pr[i] = right->pa[i];
      pi[i] = 0.5*(pl[i] + pr[i]);
    }
    /* Velocity gradients for left/right cells */
    for (i = 0; i < nd; i++){
      left->dudx[i][id-1] += pi[iqu0+i+1]*divdx[id-1];
      right->dudx[i][id-1] -= pi[iqu0+i+1]*divdx[id-1];
    }
    /* Velocity gradients at interface */
    for (tid = 1; tid <= nd; tid++){
      /* Normal velocity gradients */
      dv[tid-1][id-1] = (pr[iqu0+tid] - pl[iqu0+tid])*divdx[id-1];
      /* Transverse velocity gradients */
      if (id != tid){
	vel_diff(tid, left->ljoin[tid-1], tfrac, vdiffl);
	vel_diff(tid, left->rjoin[tid-1], tfrac, vdiffr);
	for (i = 1; i <= nd; i++)
	  dv[i-1][tid-1] = vdiffl[i-1] + vdiffr[i-1];
	vel_diff(tid, right->ljoin[tid-1], tfrac, vdiffl);
	vel_diff(tid, right->rjoin[tid-1], tfrac, vdiffr);
	for (i = 1; i <= nd; i++){
	  dv[i-1][tid-1] += vdiffl[i-1] + vdiffr[i-1];
	  dv[i-1][tid-1] *= 0.25*divdx[tid-1];
	}
      }
    }
  }
  else{
    /* Needs boundary treatment */
    if (jptr->bound == lebound){
      /* Left external boundary */
      right = jptr->right;
      for (i = 0; i < nq; i++)
	pr[i] = right->pa[i];
      usrbd(id, lebound, jptr, pl, pr);
      /* Interface state */
      for (i = 0; i < nq; i++)
	pi[i] = 0.5*(pl[i] + pr[i]);
      /* Velocity differences for right cell */
      for (i = 0; i < nd; i++)
	right->dudx[i][id-1] -= pi[iqu0+i+1]*divdx[id-1];
      /* Velocity gradients at interface */
      for (tid = 1; tid <= nd; tid++){
	/* Normal velocity gradients */
	dv[tid-1][id-1] = (pr[iqu0+tid] - pl[iqu0+tid])*divdx[id-1];
	/* Transverse velocity gradients */
	if (id != tid){
	  if (bcl[id-1] == bcwall){
	    for (i = 1; i <= nd; i++)
	      dv[i-1][tid-1] = 0.0;
	  }
	  else{
	    vel_diff(tid, right->ljoin[tid-1], tfrac, vdiffl);
	    vel_diff(tid, right->rjoin[tid-1], tfrac, vdiffr);
	    for (i = 1; i <= nd; i++)
	      dv[i-1][tid-1] = 0.5*(vdiffl[i-1] + vdiffr[i-1])
		*divdx[tid-1];
	  }
	}
      }
    }
    else if (jptr->bound == lfibound){
      /* Left internal boundary */
      right = jptr->right;
      pjptr = jptr->par;
      for (i = 0; i < nq; i++){
	pr[i] = right->pa[i];
	pi[i] = pjptr->pa0[i]  + tfrac*pjptr->dpdt[i];
      }
      for (tid = 1; tid <= nd; tid++){
	if (id != tid){
	  il = dpdx_index(id, tid);
	  for (i = 0; i < nq; i++)
	    pi[i] += (jptr->xj[tid-1] - pjptr->xj[tid-1])*
	      pjptr->dpdx[i+il];
	}
      }
      /* Notional left state */
      for (i = 0; i < nq; i++)
	pl[i] = 2.0*pi[i] - pr[i]; 
      /* Velocity differences for right cell */
      for (i = 0; i < nd; i++)
	right->dudx[i][id-1] -= pi[iqu0+i+1]*divdx[id-1];
      /* Velocity gradients at interface */
      for (tid = 1; tid <= nd; tid++){
	/* Normal velocity gradients */
	dv[tid-1][id-1] = 2.0*(pr[iqu0+tid] - pi[iqu0+tid])*divdx[id-1];
	/* Transverse velocity gradients */
	if (id != tid){
	  il = dpdx_index(id, tid);
	  for (i = 1; i <= nd; i++)
	    dv[i-1][tid-1] = pjptr->dpdx[iqu0+i+il];
	}
      }
    }
    else if (jptr->bound == rebound){
      /* Right external boundary */
      left = jptr->left;
      for (i = 0; i < nq; i++)
	pl[i] = left->pa[i];
      usrbd(id, rebound, jptr, pl, pr);
      /* Interface state */
      for (i = 0; i < nq; i++)
	pi[i] = 0.5*(pl[i] + pr[i]);
      /* Velocity differences for left cell */
      for (i = 0; i < nd; i++)
	left->dudx[i][id-1] += pi[iqu0+i+1]*divdx[id-1];
      /* Velocity gradients at interface */
      for (tid = 1; tid <= nd; tid++){
	/* Normal velocity gradients */
	dv[tid-1][id-1] = (pr[iqu0+tid] - pl[iqu0+tid])*divdx[id-1];
	/* Transverse velocity gradients */
	if (id != tid){
	  if (bcr[id-1] == bcwall){
	    for (i = 1; i <= nd; i++)
	      dv[i-1][tid-1] = 0.0;
	  }
	  else{
	    vel_diff(tid, left->ljoin[tid-1], tfrac, vdiffl);
	    vel_diff(tid, left->rjoin[tid-1], tfrac, vdiffr);
	    for (i = 1; i <= nd; i++)
	      dv[i-1][tid-1] = 0.5*(vdiffl[i-1] + vdiffr[i-1])
		*divdx[tid-1];
	  }
	}
      }
    }
    else if (jptr->bound == rfibound){
      /* Right internal boundary */
      left = jptr->left;
      pjptr = jptr->par;
      for (i = 0; i < nq; i++){
	pl[i] = left->pa[i];
	pi[i] = pjptr->pa0[i] + tfrac*pjptr->dpdt[i];
      }
      for (tid = 1; tid <= nd; tid++){
	if (id != tid){
	  il = dpdx_index(id, tid);
	  for (i = 0; i < nq; i++)
	    pi[i] += (jptr->xj[tid-1] - pjptr->xj[tid-1])*
	      pjptr->dpdx[i+il];
	}
      }
      /* Notional right state */
      for (i = 0; i < nq; i++)
	pr[i] = 2.0*pi[i] - pl[i]; 
      /* Velocity differences for left cell */
      for (i = 0; i < nd; i++)
	left->dudx[i][id-1] += pi[iqu0+i+1]*divdx[id-1];
      /* Velocity gradients at interface */
      for (tid = 1; tid <= nd; tid++){
	/* Normal velocity gradients */
	dv[tid-1][id-1] = 2.0*(pi[iqu0+tid] - pl[iqu0+tid])*divdx[id-1];
	/* Transverse velocity gradients */
	if (id != tid){
	  il = dpdx_index(id, tid);
	  for (i = 1; i <= nd; i++)
	    dv[i-1][tid-1] = pjptr->dpdx[iqu0+i+il];
	}
      }
    }
  }
  /* Calculate fluxes */
  /* Diffusion coefficients (Harmonic mean) */
  usrmu(pl, mul);
  usrmu(pr, mur);
  for (i = 0; i < nq; i++){
    if ((mabs(mul[i]) > tiny) && (mabs(mur[i]) > tiny))
      mu[i] = 2.0/(1.0/mul[i] + 1.0/mur[i]);
    else
      mu[i] = 0.0;
  }
  /* Divergence of v */
  divv = 0.0;
  for (i = 0; i < nd; i++)
    divv += dv[i][i];
  if (strcmp(geom, "RZP") == 0){
    if ((strcmp(coord[0], "R") == 0) && (id == 1))
      divv += pi[iqu0+1] /max(jptr->xj[0], 0.5/divdx[0]);
    else if ((strcmp(coord[1], "R") == 0) && (id == 2))
      divv += pi[iqu0+2] /max(jptr->xj[1], 0.5/divdx[1]);
  }
  /* Momentum fluxes */
  vcof = mu[iqu0+1];
  for (i = 0; i < nd; i++)
    fv[i] = - vcof*(dv[i][id-1] + dv[id-1][i]);
  fv[id-1] += 2.0*vcof*divv/3.0;
  for (i = 1; i <= nd; i++){
    f[iqu0+i] += fv[i-1];
    /* Energy flux */
    f[iqe] += pi[iqu0+i]*fv[i-1];
  }
  /* Thermal conduction */
  f[iqe] -= mu[iqe]*(usrte(pr) - usrte(pl))*divdx[id-1];
  /* Diffusion of remaining variables */
  for (i = iqe + 1; i < nq; i++)
    f[i] -= mu[i]*(pr[i] - pl[i])*divdx[id-1];
  return;
}
//end vflux.c

//file vel_diff.c
/*
+ vel_diff - Velocity differences at interface.
     Description :
        This  routine  calculates   the  velocity  difference  at  the
        interface (join) pointed to by jptr. The result is returned in
        vdiff[MAXQ]) with  grad[0] =  x velocity difference  etc. This
        can be turned  into a gradient by dividing  by the appropriate
        dx for the grid and direction.
*/

void vel_diff(int id, struct join *jptr, float tfrac, float vdiff[])
{
  int i, tid, il;
  float pl[MAXQ], pr[MAXQ];
  struct cell *left, *right;
  struct join *pjptr;
  if (jptr->bound < lebound){
    /* No special boundary treatment */
    for (tid = 1; tid <= nd; tid++)
      vdiff[tid-1] = jptr->right->pa[iqu0+tid] 
	- jptr->left->pa[iqu0+tid];
  }
  else{
    /* Needs boundary treatment */
    if (jptr->bound == lebound){
      /* Left external boundary */
      right = jptr->right;
      for (i = 0; i < nq; i++)
	pr[i] = right->pa[i];
      usrbd(id, lebound, jptr, pl, pr);
      for (tid = 1; tid <= nd; tid++)
	vdiff[tid-1] = pr[iqu0+tid] - pl[iqu0+tid];
    }
    else if (jptr->bound == lfibound){
      /* Left internal boundary */
      right = jptr->right;
      pjptr = jptr->par;
      for (tid = 1; tid <= nd; tid++)
	vdiff[tid-1] = right->pa[iqu0+tid] - pjptr->pa0[iqu0+tid]
		      - tfrac*pjptr->dpdt[iqu0+tid];
      for (tid = 1; tid <= nd; tid++){
	if (id != tid){
	  il = dpdx_index(id, tid);
	  for (i = 1; i <= nd; i++)
	    vdiff[i-1] -= (jptr->xj[tid-1] - pjptr->xj[tid-1])*
	      pjptr->dpdx[iqu0+i+il];
	}
      }
      for (tid = 1; tid <= nd; tid++)
	vdiff[tid-1] *= 2.0;
    }
    else if (jptr->bound == rebound){
      /* Right external boundary */
      left = jptr->left;
      for (i = 0; i < nq; i++)
	pl[i] = left->pa[i];
      usrbd(id, rebound, jptr, pl, pr);
      for (tid = 1; tid <= nd; tid++)
	vdiff[tid-1] = pr[iqu0+tid] - pl[iqu0+tid];
    }
    else if (jptr->bound == rfibound){
      /* Right internal boundary */
      left = jptr->left;
      pjptr = jptr->par;
      for (tid = 1; tid <= nd; tid++)
	vdiff[tid-1] = pjptr->pa0[iqu0+tid]
          + tfrac*pjptr->dpdt[iqu0+tid] - left->pa[iqu0+tid];
      for (tid = 1; tid <= nd; tid++){
	if (id != tid){
	  il = dpdx_index(id, tid);
	  for (i = 1; i <= nd; i++)
	    vdiff[i-1] += (jptr->xj[tid-1] - pjptr->xj[tid-1])*
	      pjptr->dpdx[iqu0+i+il];
	}
      }
      for (tid = 1; tid <= nd; tid++)
	vdiff[tid-1] *= 2.0;
    }
  }
  return;
}
//end vel_diff.c

//file apply_grav.c
/*
+  apply_grav -  apply gravitational force.
     Description :
     Applies the  gravitational force to all levels  using the current
     potential. It assumes that  the mcell hierarchy has been created.
*/
void apply_grav(float dt)
{
  int ilev, id;
  float tv, ql, qr, dmom, quo, qun;
  struct cell *cptr;
  struct mcell *mcptr, *left, *right;
  /* Set coefficient in Poisson equation */
  tv = 0.5*gconst;
  for (ilev = levc; ilev <= alev; ilev++){
    if (limpe > 1)
      /* Communicate ilev halo qa (com1) */
      mcom_halo(mcom1, ilev);
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      if (!cptr->refined){
	/* Only do unrefined cells */
	mcptr = cptr->mcptr;
	for (id = 1; id <= nd; id++){
	  /* Left potential */
	  left = mcptr->lcell[id-1];
	  if (left != NULL)
	    /* left cell exists */
	    ql = left->qa[iqpot];
	  else
	    /* Left external symmetry boundary */
	    ql = mcptr->qa[iqpot];
	  /* Right potential */
	  right = mcptr->rcell[id-1];
	  if (right != NULL)
	    /* Right interior cell */
	    qr = right->qa[iqpot];
	  else
	    /* Right external symmetry boundary*/
	    qr = mcptr->qa[iqpot];
	  dmom = - cptr->qa[iqd]*tv*(qr - ql)*dt/
	    (lev[ilev]->dxbar[id-1]);
	  quo = cptr->qa[iqu0+id];
	  qun = quo + dmom;
	  cptr->qa[iqu0+id] = qun;
	  cptr->qa[iqe] += 0.5*(qun*qun - quo*quo)/cptr->qa[iqd];
	}
      }
      cptr = cptr->nxt;
    }
  }
  return;
}
//end apply_grav.c

//file evol_part.c
/*
+ evol_part - Evolve lagrangian particles.
     Description :
        This evolves the  lagrangian particles in the  cell pointed to
        by  cptr.   Particles  that   need  to   be  moved   to  other
        cells/processors are put in the list pointed to by pmov.
*/
void evol_part(int ilev, float tcur, float dt, float tfrac, 
  struct cell *cptr, struct part_list *&pmov)
{
  int id, i, tilev, pproc, npe;
  (void)tilev, (void)pproc, (void)npe;
  float diffl[MAXQ], diffr[MAXQ], dpdx[MAXQ][MAXD], xc[MAXD], dx, 
    pc[MAXQ], pa[MAXQ], xp0[MAXD], vp0[MAXD], dvp[MAXD], xp[MAXD], 
    vp[MAXD], en0, en, m0, m, den, dm, dtp, dtph;
  bool outcell, outdomain, above;
  (void)above;
  struct cell *ncptr, *tcptr;
  (void)ncptr, (void)tcptr;
  struct part *pptr, *npptr;
  struct part_list *plist;
  for (id = 1; id <= nd; id++){
    /* Calculate fluid primitive derivatives */
    cal_diffl(ilev, id, tfrac, cptr, diffl);
    cal_diffr(ilev, id, tfrac, cptr, diffr);
    dx = lev[ilev]->dxbar[id-1];
    for (i = 0; i < nq; i++)
      dpdx[i][id-1] = av(diffl[i], diffr[i])/dx;
    /* Cell centre coordinates */
    xc[id-1] = cptr->xc[id-1];
  }
  /* Central fluid state */
  for (i = 0; i < nq; i++)
    pc[i] = cptr->pa[i];
  /* Run through particles */
  pptr = cptr->pptr;
  while (pptr != NULL){
    pptr = cptr->pptr;
    /* Next particle */
    npptr = pptr->nxt;
    /* Set particle timestep */
    dtp = max(0.0, tcur + dt - pptr->t);
    /* Half timestep */
    dtph = 0.5*dtp;
    /* Update particle time*/
    pptr->t += dtp;
    /* Initial position and velocity */
    for (i = 0; i < nd; i++){
      xp0[i] = pptr->xp[i];
      vp0[i] = pptr->vp[i];
    }
    /* Initial Thermal energy */
    en0 = pptr->en;
    /* Initial mass */
    m0 = pptr->mass;
    /* Fluid state at initial position */
    for (i = 0; i < nq; i++){
      pa[i] = pc[i];
      for (id = 1; id <= nd; id++)
	pa[i] += (xp0[id-1] - xc[id-1])*dpdx[i][id-1];
    }
    /* Half step */
    usrprate(pptr, dtph, pa, dvp, den, dm);
    for (id = 1; id <= nd; id++){
      xp[id-1] = xp0[id-1] + dtph*vp0[id-1];
      //vp[id-1] = vp0[id-1] + dtph*dvp[id-1];
      vp[id-1] = pa[iqu0+id];// Use fluid velocity
    }
    m = m0 + dtph*dm;
    en = en0 + dtph*den;
    /* Fluid state at half time position */
    for (i = 0; i < nq; i++){
      pa[i] = pc[i];
      for (id = 1; id <= nd; id++)
	pa[i] += (xp[id-1] - xc[id-1])*dpdx[i][id-1];
    }
    /* Set half step particle variables for usrprate */
    for (id = 1; id <= nd; id++){
      pptr->xp[id-1] = xp[id-1];
      vp[id-1] = pa[iqu0+id];      
    }
    pptr->mass = m;
    pptr->en = en;
    /* Full step */
    usrprate(pptr, dtp, pa, dvp, den, dm);
    for (id = 1; id <= nd; id++){
      xp[id-1] = xp0[id-1] + dtp*vp[id-1];
      //vp[id-1] = vp0[id-1] + dtp*dvp[id-1];
      vp[id-1] = pa[iqu0+id];// Use fluid velocity     
    }
    m = m0 + dtp*dm;
    en = en0 + dtp*den;
    /* Move particle off nd cell boundary if necessary */
    if (mabs(xp[nd-1] - cptr->ljoin[nd-1]->xj[nd-1]) < 
	1.0e-3*lev[ilev]->dxbar[nd-1]){
      xp[nd-1] = cptr->ljoin[nd-1]->xj[nd-1] - 
        1.0e-3*lev[ilev]->dxbar[nd-1];
    }
    if (mabs(xp[nd-1] - cptr->rjoin[nd-1]->xj[nd-1]) < 
	1.0e-3*lev[ilev]->dxbar[nd-1]){
      xp[nd-1] = cptr->rjoin[nd-1]->xj[nd-1] + 
        1.0e-3*lev[ilev]->dxbar[nd-1];
    }
    /* Check for symmetry boundaries */
    for (id = 1; id <= nd; id++){
      if ((xp[id-1] <= xl[id-1]) && (bcl[id-1] == bcref)){
	/* Left symmetry boundary -- reset velocity and position */
	xp[id-1] = 2.0*xl[id-1] - xp[id-1]; 
	vp[id-1] = -vp[id-1];
      }
      if ((xp[id-1] >= xr[id-1]) && (bcr[id-1] == bcref)){
	/* Right symmetry boundary -- reset velocity and position */
	xp[id-1] = 2.0*xr[id-1] - xp[id-1]; 
	vp[id-1] = -vp[id-1];
      }
    }
    /* Update particle variables */
    for (id = 1; id <= nd; id++){
      pptr->xp[id-1] = xp[id-1];
      pptr->vp[id-1] = vp[id-1];
    }
    pptr->mass = m;
    pptr->en = en;
    /* Check if particle is still in cell */
    outcell = false;
    for (id = 1; id <= nd; id++){
      if (xp[id-1] < cptr->ljoin[id-1]->xj[id-1]){
	outcell = true;
	break;
      }
      if (xp[id-1] > cptr->rjoin[id-1]->xj[id-1]){
	outcell = true;
	break;
      }
    }
    if (outcell){
      /* Remove from cell */
      rm_part(pptr);
      /* Check if in domain */
      outdomain = false;
      for (id = 1; id <= nd; id++){
	if (xp[id-1] <= xl[id-1]){
	  outdomain = true;
	  break;
	}
	if (xp[id-1] >= xr[id-1]){
	  outdomain = true;
	  break;
	}
      }
      if (outdomain){
	/* Outside domain -- delete */
	cout << "deleting particle at ";
	for (id = 1; id <= nd; id++)
	  cout << " " << pptr->xp[id-1];
	cout << '\n';
	delete pptr;
      }
      else{
	/* Add to processor move list */
	plist = new part_list;
	plist->pptr = pptr;
	plist->pproc = -1;
	plist->nxt = pmov; 
	pmov = plist;
	cout.setf(ios::scientific);
	cout << std::setprecision(10);
      }
    }
    pptr = npptr;
  }
  return;
}
//end evol_part.c
