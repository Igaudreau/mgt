#undef MAIN
#define PLOT
#include <iostream>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <new>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "funcs.h"
#include "plot.h"
extern bool mvcmd(char key[]);
extern bool get_cmd(char key[]);
extern bool get_n(char key[]);
extern bool get_r(float &value);
extern bool get_i(int &value);
extern void report();
extern bool cal_diffql(int id, struct cell *cptr, float diff[]);
extern bool cal_diffqr(int id, struct cell *cptr, float diff[]);
extern bool usrps(struct cell *cptr, float &s);
extern bool usrpartps(char pitem[], struct part *pptr, float &s);
extern float usrvel(int id, struct cell *cptr);
extern float usrbfield(int id, struct cell *cptr);
extern float usrgfield(int id, struct cell *cptr);
extern void usrga(int id, struct cell *cptr, struct cell *left,
  struct cell *right);
extern void i_to_it(int ndim, int nt[], int icell, int it[]);
extern void readxq(char []);
extern int it_to_i(int ndim, int nt[], int it[]);
extern struct part_list *mov_part(bool all, struct part_list *pmov);
extern void mk_partlist(struct part_list *&start, 
  struct part_list *&end);
extern void merge_array(int buf[], int nvalues[]);
extern void merge_array(float buf[], int nvalues[]);
extern void broadcast(bool *value);
extern void broadcast(int buf[], int nvalues);
extern float global_sum(float value);
extern float global_max(float value);
extern struct cell *find_fcell(int &ilev, float xc[]);
void com_halo(int com, int ilev, float tfrac);
void plact();
void pldac();
bool plcmd(char key[]);
void plot(char key[]);
void plotg(char name[], char key[], int *ios);
void plotp(char name[], char key[], int *ios);
void levint(int ilev, char key[]);
void fill(struct cell *cptr, bool vvect, bool bvect, int ncell, 
  int nt[], float divdx[]);
void wr_graph();

//file plact.c
/*
+ plact - activate plotting module (graphics)
    Description :
       Call graphics activation routine and set defaults.
*/
void plact()
{
  int id;
  /* Activate implementation package */
  if (mype == 0)
    model_prepare();
  /* Default plotting item (1st conserved variable) */
  strcpy(item, "QA1");
  /* Kind of value (solution itself) */
  strcpy(kind, "Q");
  /* Set coordinate ranges */
  for (id = 0; id < MAXD; id++){
    xlg[id] = -big;
    xrg[id] = big;
  }
  /* No masks */
  for (id = 1; id <= nd;id++){
    mcoord[id-1] = false;
    mval[id-1] = 0.0;
  }
  /* Do not write out graphics file */
  wrgraph = false;
  return;
}
//end plact.c

//file pldac.c
/*
+ pldac - deactivate plotting module (pl version)
    Description :
       call graphics  deactivation routine (this may  close plot files
       etc.)
*/
void pldac()
{
  /* Just deactivate implementation package */
  model_flush();
  return;
}
//end pldac.c

//file plcmd.c
/*
+ plcmd - plotting command interface
    Description :
       Implements  rad-specific graphics commands  in a  muvi context.
       We want  to give the program ademo  versions of g, p,  v, s and
       plot first  class status because  a user will normally  want to
       look  at  models  produced  by  this  program.   therefore,  we
       intercept these commands here, and implement them properly.  we
       also intercept  the surrogate commands,  gg, pp, vv,  ss, pplot
       and  pass them onto  the muvi  command interpreter  by changing
       key!!!  this is just about ok, but relies on rather special use
       of the mvcmd function.
*/
bool plcmd(char key[])
{
  bool found, fail;
  int oplmax, tlev, id;
  char filename[SZNAM];
  found = true;
  /* Intercept muvi commands */
  if ((strcmp(key, "GG") == 0) ||
     (strcmp(key, "PP") == 0) ||
     (strcmp(key, "SS") == 0) ||
     (strcmp(key, "VV") == 0) ||
     (strcmp(key, "PPLOT") == 0)){
    found = false;
  }
  /* Plots - show current plotting parameters */
  else if (strcmp(key, "PLOTS") == 0){
    cout << "rad/plot parameters:-" << '\n';
    cout << "item " << item << " kind " << kind << '\n';
  }
  /* Plot item - set plotting item (pg, rho, etc.) */
  else if (strcmp(key, "PLOT") == 0){
    if (!get_cmd(key)){
      cout << "Error reading plot item" << '\n';
      report();
    }
    else
      strcpy(item, key);
  }
  /* Kind kind - set kind of plot (q=solution, r=residual, etc.) */
  else if (strcmp(key, "KIND") == 0){
    if (!get_cmd(key)){
      cout << "Error reading plot kind" << '\n';
      report();
    }
    else
      strcpy(kind, key);
  }
  /* Readxq - read .xq file and plot */
  else if (strcmp(key, "READXQ") == 0){
    if (!get_n(key)){
      cout << "Error : reading filename" << '\n';
      report();
    }
    else{
      strcpy(filename, key);
      readxq(filename);
    }
  }
  /* Export - export plot to file (*.exp) */
  else if (strcmp(key, "EXPORT") == 0){
    if (!get_n(key)){
      cout << "Error : reading filename" << '\n';
      report();
    }
    else{
      strcpy(expname, key);
      strcat(expname, ".xq");
      fail = false;
      /* Master checks if file exists */
      if (mype == 0){
	if ((expfile=fopen(expname, "r"))){
	  cout << "Error: file exists" << '\n';
	  fclose(expfile);
	  fail = true;
	}
      }
      /* All processors need to know about failure */ 
      broadcast(&fail);
      if (fail)
	report();
      else{
	/* Master checks if file can be opened */
	if (mype == 0){
	  fail = false;
	  if (!(expfile = fopen(expname, "w"))){
	    cout << "Error: unable to open graphics file" <<  '\n';
	    fail = true;
	  }
	}
	/* All processors need to know about failure */ 
	broadcast(&fail);
	if (fail)
	  report();
	else
	  wrgraph = true;
      }
    }
  }
  /* XR - x coordinate range */
  else if (strcmp(key, "XR") == 0){
    int i;
    for (i = 0; i < nd; i++){
      if (strcmp(coord[i], "X") == 0)
	break;
    }
    if (!get_r(xlg[i])){
      if (mype == 0)      
	cout << "Error : reading left limit" << '\n';
      report();
    }
    else if (!get_r(xrg[i])){
      if (mype == 0)
	cout << "Error : reading right limit" << '\n';
      report();
    }
  }
  /* NOXR - default x coordinate range */
  else if (strcmp(key, "NOXR") == 0){
    int i;
    for (i = 0; i < nd; i++){
      if (strcmp(coord[i], "X") == 0)
	break;
    }
    xlg[i] = -big;
    xrg[i] = big;
  }
  /* YR - y coordinate range */
  else if (strcmp(key, "YR") == 0){
    int i;
    for (i = 0; i < nd; i++){
      if (strcmp(coord[i], "Y") == 0)
	break;
    }
    if (!get_r(xlg[i])){
      if (mype == 0)      
	cout << "Error : reading left limit" << '\n';
      report();
    }
    else if (!get_r(xrg[i])){
      if (mype == 0)            
	cout << "Error : reading right limit" << '\n';
      report();
    }
  }
  /* NOYR - default x coordinate range */
  else if (strcmp(key, "NOYR") == 0){
    int i;
    for (i = 0; i < nd; i++){
      if (strcmp(coord[i], "Y") == 0)
	break;
    }
    xlg[i] = -big;
    xrg[i] = big;
  }
  /* ZR - z coordinate range */
  else if (strcmp(key, "ZR") == 0){
    int i;
    for (i = 0; i < nd; i++){
      if (strcmp(coord[i], "Z") == 0)
	break;
    }
    if (!get_r(xlg[i])){
      if (mype == 0)      
	cout << "Error : reading left limit" << '\n';
      report();
    }
    else if (!get_r(xrg[i])){
      if (mype == 0)            
	cout << "Error : reading right limit" << '\n';
      report();
    }
  }
  /* NOZR - default x coordinate range */
  else if (strcmp(key, "NOZR") == 0){
    int i;
    for (i = 0; i < nd; i++){
      if (strcmp(coord[i], "Z") == 0)
	break;
    }
    xlg[i] = -big;
    xrg[i] = big;
  }
  /* RR - r coordinate range */
  else if (strcmp(key, "RR") == 0){
    if ((strcmp(geom,"RZP") != 0) && (strcmp(geom,"RTP") != 0)){
      if (mype == 0)
	cout << "Error : not RZP or RTP" << '\n';
    }
    else{
      int i;
      for (i = 0; i < nd; i++){
	if (strcmp(coord[i], "R") == 0)
	  break;
      }
      if (!get_r(xlg[i])){
	if (mype == 0)      	
	  cout << "Error : reading left limit" << '\n';
	report();
      }
      else if (!get_r(xrg[i])){
	if (mype == 0)      	
	  cout << "Error : reading right limit" << '\n';
	report();
      }
    }
  }
  /* NORR - default r coordinate range */
  else if (strcmp(key, "NORR") == 0){
    if ((strcmp(geom,"RZP") != 0) && (strcmp(geom,"RTP") != 0)){
      if (mype == 0)
	cout << "Error : not RZP or RTP" << '\n';
    }
    else{
      int i;
      for (i = 0; i < nd; i++){
	if (strcmp(coord[i], "R") == 0)
	  break;
      }
      xlg[i] = -big;
      xrg[i] = big;
    }
  }
  else if (strcmp(key, "MASK") == 0){
    if (!get_cmd(cname)){
      if (mype == 0)            
	printf("Error:  reading coordinate name.");
      report();
      return(true);
    }
    else{
      found = false;
      for (id = 1; id <= nd; id++){
	if (strcmp(coord[id-1], cname) == 0){
	  mcoord[id-1] = true;
	  mid = id;
	  found = true;
	  break;
	}
      }
      if (!found){
	if (mype == 0)      	
	  cout << "Error:  unknown coordinate." << '\n';
	report();
	return(true);
      }
      else if (!get_r(mval[mid-1])){
	if (mype == 0)      	
	  cout << "Error:  reading mask coordinate value." << '\n';
	report();
	return(true);
      }
    }
  }
  else if (strcmp(key, "MASKI") == 0){
    if (!get_i(id)){
      if (mype == 0)            
	printf("Error:  reading coord id.");
      report();
      return(true);
    }
    else if (!get_r(mval[id-1])){
      if (mype == 0)            
	cout << "Error:  reading mask coordinate value." << '\n';
      report();
      return(true);
    }
    else
      mcoord[id-1] = true;
  }
  else if (strcmp(key, "UNMASK") == 0){
    if (!get_cmd(cname)){
      if (mype == 0)            
	printf("Error:  reading coordinate name.");
      report();
      return(true);
    }
    else{
      if (strcmp(cname, "*") == 0){
	found = true;
	for (id = 1; id <= MAXD; id++)
	  mcoord[id-1] = false;
      }
      else{
	found = false;
	for (id = 1; id <= MAXD; id++){
	  if (strcmp(coord[id-1], cname) == 0){
	    mcoord[id-1] = false;
	    found = true;
	    break;
	  }
	}
      }
      if (!found){
	if (mype == 0)      	
	  cout << "Error:  unknown coordinate." << '\n';
	report();
	return(true);
      }
    }
  }
  /* P, G - plot solution, grid */
  else if ((strcmp(key, "P") == 0) ||
          (strcmp(key, "G") == 0) ||
          (strcmp(key, "V") == 0) ||
          (strcmp(key, "B") == 0) ||
	   (strcmp(key, "GV") == 0)){
    if (nom){
      if (mype == 0)
	cout <<  "no model available." << '\n';
    }
    else if (strcmp(key, "G") != 0)
      plot(key);
    else{
      oplmax = plmax;
      for (tlev = oplmax; tlev >= plmin; tlev--){
	if (alev >= tlev){
	  plmax = tlev;
	  plot(key);
	}
      }
      plmax = oplmax;
    }
  }
  /* Otherwise, not found (yet) */
  else
    found = false;
  /* If not found, try for commands in the muvi interpreter */
  if ((!found) && (mvcmd(key)))
    found = true;
  return(found);
}
//end plcmd.c

//file plot.c
/*
+ plot - plot solution, grid or velocities
    Description :
       This routine  puts a parody of  the current rad  model into the
       muvi model database.   the parody is based on  key which can be
       one of gg, pp, vv  (grid,scalar,vector); any other key leads to
       the entire model being transferred.
*/
void plot(char key[])
{
  int static model_count(0);
  int ios;
  char name[SZNAM];
  /*Invent a temporary model name */
  sprintf(name, "rad#%d", 10000000+model_count);
  model_count += 1;
  ios = 0;
  if (strncmp(item,"PPART",5) == 0)
    /* Plot particles */
    plotp(name, key, &ios);
  else
    /* Plot solution, grid or velocities */
    plotg(name, key, &ios);
  /* If that worked, visualize this */
  if ((mype == 0) && (!wrgraph)){
    if (ios == 0){
      vis_begin();
      if (strcmp(key, "G") == 0)
	model_seta((char*)"grid", &ios);
      else if (strcmp(key, "P") == 0)
	model_seta((char*)"scalar", &ios);
      else if ((strcmp(key, "V") == 0) || (strcmp(key, "B") == 0)
	       || (strcmp(key, "GV") == 0))  
	model_seta((char*)"vector", &ios);
      else
	std::cout << "Error : illegal plot option" << '\n';
      model_vis(&ios);
      vis_end();
    }
    /* Remove the model (or die trying) */
    ios = 0;
    model_remove(name, &ios);
  }
  wrgraph = false;
  return;
}
//end plot.c

//file plotg.c
/*
+ plotg - plot solution, grid or velocities
     Description :
        This plots the current model using muvi.
*/
void plotg(char name[], char key[], int *ios)
{
  char co[SZCO], un[SZCO], physic[SZNAM], vl[SZCQ], vu[SZCQ],
    vvlab[SZNAM];
  bool doit;
  float qqb[MAXQB], qmin, qmax;
  int tlev, i, j, k, id, it[MAXD], ssb[MAXSSB], nqb, vmlt[3][3], 
    nvalues[MAXPE], buf[2*MAXPE], npe, count, tid;
  struct cell *cptr;
  /* Set ios flag */
  *ios = 0;
  /* Check for masks */
  ndm = nd;
  for (id = 1; id <= nd; id++){
    if (mcoord[id-1] == true){
      /* Move mask off cell edge (temporary copy) */
      mvalt[id-1] = xl[id-1] + lev[plmax]->dxbar[id-1]*(int)((mval[id-1]
        - xl[id-1])/lev[plmax]->dxbar[id-1])
        + 0.01*lev[plmax]->dxbar[id-1];
      ndm--;
    }
  }
  /* G or P -- scalar field only */
  if ((strcmp(key, "G") == 0) ||
      (strcmp(key, "P") == 0)){
    nnq = 1;
    if (mype == 0){
      strcpy(vvlab, item);
      model_set0s((char*)"plot", vvlab, ios);
      /* Composite strings for variable labels */
      strcpy(vl, ";");strcat(vl, item);strcat(vl, ";");
    }
  }
  /* V, B, GV -- vector field) */
  else if  ((strcmp(key, "V") == 0) || (strcmp(key, "B") == 0)
	    || (strcmp(key, "GV") == 0)){
    nnq = ndm;
    if (mype == 0){
      strcpy(vl, ";");
      for (id = 1; id <= nd; id++){
	if (!mcoord[id-1]){
	  strcat(vl, "U");
	  strcat(vl, coord[id-1]);
	  strcat(vl, ";");
	}
      }
    }
  }
  else{
    cout << "Error:  unexpected key." << '\n';
    *ios = 1;
  }
  if (*ios == 1)
    return;
  /* Geometry */
  /* Construct composite strings for coordinates, units */
  if (mype == 0){
    strcpy(co, ";");
    strcpy(un, ";");
    for (id = 1; id <= nd; id++){  
      if (!mcoord[id-1]){
	strcat(co, coord[id-1]);
	strcat(un, ";");
	strcat(co, ";");
      }
    }
    for (i = 0; i <= nnq; i++)
      vu[i] = ';';
    vu[nnq+1] = '\0';
    /* Wash out mlt */
    for (i = 0; i < 3*nnq; i++)
      vmlt[0][i] = 0;
    /* Physics label */
    strcpy(physic, "gas");
  }
  /* Finest grid mesh spacing */
  for (id = 1; id <= nd; id++)
      dxt[id-1] = lev[plmax]->dxbar[id-1];
  /* Set local plotting domain and mesh spacing (possibly masked) */
  tid = 1;
  for (id = 1; id <= nd; id++){
    if (!mcoord[id-1]){
      xlgt[tid-1] = max(xlg[id-1], xl[id-1]);
      xlgt[tid-1] = xl[id-1] + dxt[id-1]*(int)((xlgt[tid-1] - xl[id-1])
        /dxt[id-1]);
      xrgt[tid-1] = min(xrg[id-1], xr[id-1]);
      xrgt[tid-1] = xr[id-1] - dxt[id-1]*(int)((xr[id-1] - xrgt[tid-1])
        /dxt[id-1]);
      dx[tid-1] = dxt[id-1];
      tid++;
    }
  }
  /* Calculate sizes etc*/
  tid = 1;
  for (id = 1; id <= nd; id++){
    /* Local values - only unmasked coordinates */
    if (!mcoord[id-1]){
      lxl[tid-1] = max(lev[plmax]->xl[id-1], xlgt[tid-1]);
      nl[tid-1] = max(0, (int)((min(xrgt[tid-1], lev[plmax]->xr[id-1])
        - lxl[tid-1] + 0.1*dx[tid-1])/dx[tid-1]));
      tid++;
    }
  }
  nsa = 1;
  for (id = 1; id <= ndm; id++)
    nsa *= nl[id-1];
  /* Determine sizes of graphics arrays */
  if (limpe > 1){
    /* Check mask direction */
    if (mcoord[nd-1]){
      /* nd coordinate is masked - check if on processor */
      if ((mvalt[nd-1] < lev[plmax]->xl[nd-1])
	  || (mvalt[nd-1] > lev[plmax]->xr[nd-1]))
	/* Not on processor - no contribution */
	nsa = 0;
    }
    mergnsa = (int)global_sum((float)nsa);
  }
  else
    mergnsa = nsa;
  nqa = nsa*nnq;
  mergnqa = mergnsa*nnq;
  /* Make graphics arrays */
  if (mype == 0){
    ssa = new int[2*mergnqa+SAFETY];
    qqa = new float[2*mergnqa+SAFETY];
  }
  else{
    ssa = new int[nqa+SAFETY];
    qqa = new float[nqa+SAFETY];
  }
  /* Fill graphics arrays with grid data */
  for (j = 0; j < nsa; j++)
    ssa[j] = bblank;
  for (j = 0; j < nqa; j++)
    qqa[j] = 0.0;         
  if ((strncmp(item, "GA", 2) ==0) || (strcmp(key, "GV") == 0)){
    /* Gravitational force is special - put in pa[iqu0+id] */
    for (tlev = plmin; tlev <= plmax; tlev++){
      if (limpe > 1)
	/* Communicate tlev halo qa (com5)  */
	com_halo(com5, tlev, 0.0);
      cptr = lev[tlev]->fcell;
      while (cptr != NULL){
	for (id = 1; id <= nd; id++)
	  usrga(id, cptr, cptr->ljoin[id-1]->left, 
		cptr->rjoin[id-1]->right);
	cptr = cptr->nxt;
      }
    }
  }
  if (strcmp(key, "G") == 0){
    cptr = lev[plmax]->fcell;
    while (cptr != NULL){
      doit = true;
      /* Check mask */
      for (id = 1; id <= nd; id++){
	if (mcoord[id-1]){
	  if ((cptr->xc[id-1] < mvalt[id-1] - 0.5*dxt[id-1])
	    || (cptr->xc[id-1] >= mvalt[id-1] + 0.5*dxt[id-1]))
	    doit = false;
	}
      }
      tid = 1;
      for (id = 1; id <= nd; id++){
	if (!mcoord[id-1]){
	  it[tid-1] = (int)((cptr->xc[id-1] - lxl[tid-1])/dx[tid-1])
	    + 1;
	  if ((it[tid-1]< 1) || (it[tid-1] > nl[tid-1]))
	    doit = false;
	  tid++;
	}
      }
      if (doit){
	i = it_to_i(ndm, nl, it) - 1;
	ssa[i] = ffluid;
	cptr->mat = ffluid;		// TODO new version has this removed - why?
      }
      cptr = cptr->nxt;
    }
  }
  else if (strcmp(key, "P") == 0){
    /* Check for level interpolation */
    if (plmin != plmax){
      /* Interpolate from coarser levels */
      for (tlev = plmin; tlev < plmax; tlev++)
	levint(tlev, key);
    }
    else{
      if (limpe > 1)
	/* Communicate ilev halo qa (com5)  */
	com_halo(com5, plmax, 0.0);
    }
    cptr = lev[plmax]->fcell;
    while (cptr != NULL){
      doit = true;
      /* Check mask */
      for (id = 1; id <= nd; id++){
	if (mcoord[id-1]){
	  if ((cptr->xc[id-1] < mvalt[id-1] - 0.5*dxt[id-1])
	    || (cptr->xc[id-1] >= mvalt[id-1] + 0.5*dxt[id-1]))
	    doit = false;
	}
      }
      tid = 1;
      for (id = 1; id <= nd; id++){
	if (!mcoord[id-1]){
	  it[tid-1] = (int)((cptr->xc[id-1] - lxl[tid-1])/dx[tid-1])
	    + 1;
	  if ((it[tid-1]< 1) || (it[tid-1] > nl[tid-1]))
	    doit = false;
	  tid++;
	}
      }
      if (doit){
	i = it_to_i(ndm, nl, it) - 1;
	if (usrps(cptr, qqa[i]))
	  ssa[i] = cptr->mat;		// TODO new version has = ffluid; why?
      }
      cptr = cptr->nxt;
    }
  }
  else if ((strcmp(key, "V") == 0) || (strcmp(key, "B") == 0)
	   || (strcmp(key, "GV") == 0)){
    /* Check for level interpolation */
    if (plmin != plmax)
      /* Interpolate from coarser levels */
      for (tlev = plmin; tlev <= plmax; tlev++)
	levint(tlev, key);
    cptr = lev[plmax]->fcell;
    while (cptr != NULL){
      doit = true;
      /* Check mask */
      for (id = 1; id <= nd; id++){
	if (mcoord[id-1]){
	  if ((cptr->xc[id-1] < mvalt[id-1] - 0.5*dxt[id-1])
	    || (cptr->xc[id-1] >= mvalt[id-1] + 0.5*dxt[id-1]))
	    doit = false;
	}
      }
      tid = 1;
      for (id = 1; id <= nd; id++){
	if (!mcoord[id-1]){
	  it[tid-1] = (int)((cptr->xc[id-1] - lxl[tid-1])/dx[tid-1])
	    + 1;
	  if ((it[tid-1]< 1) || (it[tid-1] > nl[tid-1]))
	    doit = false;
	  tid++;
	}
      }
      if (doit){
	i = it_to_i(ndm, nl, it) - 1;
	ssa[i] = ffluid;
	k = ndm*i - 1;
	if (strcmp(key, "V") == 0){
	  for (id = 1; id <= nd; id++){
	    if (!mcoord[id-1]){
	      k = k + 1;
	      qqa[k] = usrvel(id, cptr);
	    }
	  }
	}
	else if (strcmp(key, "B") == 0){
	  for (id = 1; id <= nd; id++){
	    if (!mcoord[id-1]){
	      k = k + 1;
	      qqa[k] = usrbfield(id, cptr);
	    }
	  }
	}
	else if (strcmp(key, "GV") == 0){
	  for (id = 1; id <= nd; id++){
	    if (!mcoord[id-1]){
	      k = k + 1;
	      qqa[k] = usrgfield(id, cptr);
	    }
	  }
	}
      }
      cptr = cptr->nxt;
    }
  }
  /* Total domain size */
  if (limpe > 1){
    if (!mcoord[nd-1])
      nl[ndm-1] = (int)global_sum((float)nl[ndm-1]);
  }
  if (mype == 0){
    /* Make master coordinate arrays */
    nxad = 0;
    for (id = 1; id <= ndm; id++)
      nxad += nl[id-1];
    xad = new float[nxad];
    xed = new float[nxad+ndm];
    /* Fill master coordinate arrays */
    nxad = - 1;
    nxed = - 1;
    for (id = 1; id <= ndm; id++){
      ncd[0][id-1] = nl[id-1];
      /* Left boundary */
      nxed++;
      xed[nxed] = xlgt[id-1];
      for (k = 1; k <= nl[id-1]; k++){
	nxad++;
	xad[nxad] = (k - 0.5)*dx[id-1] + xlgt[id-1];
	nxed++;
	xed[nxed] = k*dx[id-1] + xlgt[id-1];
      }
    }
  }
  if (limpe > 1){
    /* Get sizes from other processors */
    buf[0] = nsa;
    buf[1] = nqa;
    for (npe = 0; npe < limpe; npe++)
      nvalues[npe] = 2;
    merge_array(buf, nvalues);
    broadcast(buf, 2*limpe);
    /* Merge ssa */
    for (npe = 0; npe < limpe; npe++)
      nvalues[npe] = buf[2*npe];
    merge_array(ssa, nvalues);
    /* Merge qqa */
    for (npe = 0; npe < limpe; npe++)
      nvalues[npe] = buf[2*npe+1];
    merge_array(qqa, nvalues);
  }
  if (mype == 0){
    /* Put model */
    nsa = mergnsa;
    nqa = mergnqa;
    nqb = 0;
    count = 0;
    count = 0;
    for (j = 0; j < mergnsa; j++){
      if (ssa[j] == ffluid)
	count += 1;
    }
    qmin = big;
    qmax = - big;
    if ((strcmp(key, "P") == 0) ||
	(strcmp(key, "V") == 0) || 
        (strcmp(key, "B") == 0)){
      for (j = 0; j < mergnsa; j++){
	if (ssa[j] == ffluid){
	  qmin = min(qmin,qqa[j]);
	  qmax = max(qmax,qqa[j]);
	}
      }
    }
    if ((count == 0) && (!wrgraph)){
      cout << "Error:  nothing to plot!" << '\n';
      *ios = 1;
    }
    else if ((fabs(qmin-qmax) <= 1.0e-34) && (!wrgraph)){
      cout << "note:  scalar field is degenerate with  " <<
	qmin << '\n';
      *ios = 1;
    }
    else if (wrgraph){
      wr_graph();
      fclose(expfile);
    }
    else{
      model_put_ugrid(name, title, geom, t, ndm, co, un, ncd, xed, xad,
        physic, nnq, vl, vu, vmlt, nqa, ssa, qqa,  nqb, ssb, qqb, ios);
    }
  }
  /* Delete graphics  arrays */
  delete [] ssa;
  delete [] xed;
  delete [] xad;
  delete [] qqa;

  return;
} 
//end plotg.c

//file plotp.c
/*
+ plotp - plot particles
     Description :
        This plots the particles using muvi.
*/
void plotp(char name[], char key[], int *ios)
{
  char co[SZCO], un[SZCO], physic[SZNAM], vl[SZCQ], vu[SZCQ],
    vvlab[SZNAM], pitem[SZNAM];
  //bool doit;
  float *qpt, *xpa, qmin, qmax;
  int i, j, ndm, id, nnq, vmlt[3][3], *spt, npart;
  struct part *pptr;
  struct part_list *plist, *tplist, *start, *end, *ostart;
  ostart = NULL;
  xpa = NULL;
  qpt = NULL;
  spt = NULL;
  qmin = big;
  qmax = -big;
  i = -1; j = -1;
  npart = 0;
  /* Check for masks */
  ndm = nd;
  for (id = 1; id <= nd; id++){
    if (mcoord[id-1] == true)
      ndm--;
  }
  /* Get particles from other processors */
  ostart = mov_part(true, NULL);
  if (mype == 0){
    /* Add particles on master */
    mk_partlist(start, end);
    if (start == NULL)
      start = ostart;
    else
      end->nxt = ostart;
    /* Count particles */
    npart = 0;
    plist = start;
    while (plist != NULL){
      npart ++;
      plist = plist->nxt;
    }
    if (npart == 0){
      cout << "note:  no particles" << '\n';
      return;
    }
    /* Create particle plot arrays */
    xpa = new float[ndm*npart+SAFETY];
    qpt = new float[npart+SAFETY];
    spt = new int[npart+SAFETY];
    /* Initialize the arrays BM 
    for (i = 0; i < ndm*npart+SAFETY; i++)
      xpa[i] = 0.0;
    for (i = 0; i < npart+SAFETY; i++) {
      qpt[i] = 0.0;
      spt[i] = 0;
    }*/
    /* Particle plot item */
    //for (i = 5; i <= strlen(item); i++)
    int item_len = static_cast<int>(strlen(item)); 
    for (i = 5; i <= item_len; i++)
      pitem[i-5] = item[i];
    /* Fill particle plot array */
    qmin = big;
    qmax = -big;
    plist = start;
    i = -1;
    j = -1;
    while (plist != NULL){
      pptr = plist->pptr;
      for (id = 0; id < nd; id++){
	if (mcoord[id] != true){
	  i++;
	  xpa[i] = pptr->xp[id];
	}
      }
      j++;
      spt[j] = ffluid;
      if (!usrpartps(pitem, pptr, qpt[j]))
	qpt[j] = (float)(j + 1); 
      qmin = min(qmin,qpt[j]);
      qmax = max(qmax,qpt[j]);
      plist = plist->nxt;
    }
  }
  /* Delete other processor particles and list */
  plist = ostart;
  while (plist != NULL){
    /* Delete particle */
    delete plist->pptr;
    /* Delete list element */
    tplist = plist;
    plist = plist->nxt;
    delete tplist;
  }
  /* Only master plots */
  if (mype == 0){
    if (end != NULL){
      /* Delete master particle list */
      end->nxt = NULL;
      plist = start;
      while (plist != NULL){
	/* Delete list element */
	tplist = plist;
	plist = plist->nxt;
	delete tplist;
      }
    }
    cout << npart << " particles " << '\n';
    if (qmax - qmin < tiny){
      /* Degenerate scalar -- add virtual particle */
      for (id = 1; id <= nd; id++)
	xpa[i+id] = xr[id-1] + lev[0]->dxbar[id-1];
      spt[j+1] = ffluid;
      qpt[j+1] = 1.1*qpt[0];
      npart++;
    }
    nnq = 1;
    /* Set ios flag */
    *ios = 0;
    model_seta((char*)"scalar", ios);
    strcpy(vvlab, item);
    model_set0s((char*)"plot", vvlab, ios);
    /* Composite strings for variable labels */
    strcpy(vl, ";");strcat(vl, item);strcat(vl, ";");
    /* Construct composite strings for coordinates, units */
    strcpy(co, ";");
    strcpy(un, ";");
    for (id = 1; id <= nd; id++){  
      if (!mcoord[id-1]){
	strcat(co, coord[id-1]);
	strcat(un, ";");
	strcat(co, ";");
      }
    }
    strcpy(vu,"; ;");
    /* Wash out vmlt */
    vmlt[0][0] = 0;
    vmlt[0][1] = 1;
    vmlt[0][2] = 0;
    /* Physics label */
    strcpy(physic, "gas");
    model_put_scatter(name, title, t, geom, ndm, co, un, physic,
      nnq, vl, vu, vmlt,  npart, xpa, spt, qpt, ios);
    /* Delete plotting arrays */
    delete [] xpa;
    delete [] qpt;
    delete [] spt;
    /* Reset to NULL for memory safety BM */
    xpa = NULL;
    qpt = NULL;
    spt = NULL;
  }
  return;
} 
//end plotp.c

//file levint.c
/*
+ levint - Interpolates solution from lev ilev for plotting.
     Description :
        This interpolates  the solution from  level ilev to  plmax and
        puts the results in the qqa array. The actual interpolation is
        done by fill.  It assumes that the level plmax plotting arrays
        have already been made.
*/
void levint(int ilev, char key[])
{
  int id, ncell, nt[MAXD];
  float divdx[MAXD]; 
  bool bvect, vvect;
  struct cell *cptr;
  struct cell_list *clist;
  /* No of fine cells, mesh spacing, size of finest grid  etc */
  ncell = ipow(nd*(plmax - ilev));
  for (id = 1; id <= nd; id++){
    divdx[id-1] = 1.0/lev[ilev]->dxbar[id-1];
    nt[id-1] = ipow(plmax - ilev);
  }
  /* Velocity vector */
  vvect = (strcmp(key, "V") == 0);
    /* Magnetic field vector */
  bvect = (strcmp(key, "B") == 0);
  if (limpe > 1){
    /* Communicate ilev halo qa dqdx (com5, com6) */
    com_halo(com5, ilev, 0.0);
    com_halo(com6, ilev, 0.0);
    /* Run through left halo cells */
    clist = plev[ilev]->lhalo;
    while (clist != NULL){
      cptr = clist->cptr;
      /* Check for refinement */
      if (!cptr->refined)
	/* Fill qqa */
	fill(cptr, vvect, bvect, ncell, nt, divdx);
      clist = clist->nxt;
    }
    /* Run through right halo cells */
    clist = plev[ilev]->rhalo;
    while (clist != NULL){
      cptr = clist->cptr;
      /* Check for refinement */
      if (!cptr->refined)
	/* Fill qqa */
	fill(cptr, vvect, bvect, ncell, nt, divdx);
      clist = clist->nxt;
    }
  }
  /* Run through internal cells */
  cptr = lev[ilev]->fcell;
  while (cptr != NULL){
    /* Check for refinement */
    if (!cptr->refined)
      /* Fill qqa */
      fill(cptr, vvect, bvect, ncell, nt, divdx);
    cptr = cptr->nxt;
  }
  return;
}
//end levint.c

//file fill.c
/*
+ fill - Fills qqa array.
     Description :
        This  fills  the  qqa  array  with the  solution  obtained  by
        interpolating from the parent  cell pointed to by cptr. Called
        by levint.   It assumes that  the level plmax  plotting arrays
        have already been made.
*/
void fill(struct cell *cptr, bool vvect, bool bvect, int ncell, 
  int nt[], float divdx[])
{
  bool doit, gotl, gotr;
  int icell, i, id, it[MAXD], j, k, tid;
  float  xlb[MAXD], xc[MAXD], pxc[MAXD], qa0[MAXQ], qat[MAXQ], 
    diffl[MAXQ], diffr[MAXQ], grad[MAXQ][MAXD];
  for (id = 1; id <= nd; id++){
    /* Parent cell position */
    pxc[id-1] = cptr->xc[id-1];
    /* Parent cell left boundary */
    xlb[id-1] = pxc[id-1] - 0.5/divdx[id-1];
    /* Calculate gradients */
    gotl = cal_diffql(id, cptr, diffl);
    gotr = cal_diffqr(id, cptr, diffr);
    if (gotl && gotr){
      for (i = 0; i < nq; i++)
	grad[i][id-1] = av(diffl[i], diffr[i])*divdx[id-1];
    }
    else if (gotr){
      for (i = 0; i < nq; i++)
	grad[i][id-1] = diffr[i]*divdx[id-1];
    }
    else{
      for (i = 0; i < nq; i++)
	grad[i][id-1] = diffl[i]*divdx[id-1];
    }
  }
  /* Run through cells on level plmax */
  for (icell = 1; icell <= ncell; icell++){
    /* Subcell positions */
    i_to_it(nd, nt, icell, it);
    for (id = 1; id <= nd; id++)
      xc[id-1] = xlb[id-1] + (it[id-1] - 0.5)*dxt[id-1];
    /* xad index */
    doit = true;
    /* Check mask */
    for (id = 1; id <= nd; id++){
      if (mcoord[id-1]){
	if ((xc[id-1] < mvalt[id-1] - 0.5*dxt[id-1])
	  || (xc[id-1] >= mvalt[id-1] + 0.5*dxt[id-1]))
	  doit = false;
      }
    }
    tid = 1;
    for (id = 1; id <= nd; id++){
      if (!mcoord[id-1]){
	it[tid-1] = (int)((xc[id-1] - lxl[tid-1])/dx[tid-1])
	  + 1;
	if ((it[tid-1]< 1) || (it[tid-1] > nl[tid-1]))
	  doit = false;
	tid++;
      }
    }
    if (doit){
      j = it_to_i(ndm, nl, it) - 1;
      /* Check for valid j */
      if ((j >= 0) && (j < nsa)){
	ssa[j] = ffluid;
	/* Central state */
	for (i = 0; i < nq; i++)
	  qa0[i] = cptr->qa[i];
	/* Interpolate using grad */
	for (i = 0; i < nq; i++){
	  for (id = 1; id <= nd; id++)
	    qa0[i] += grad[i][id-1]*(xc[id-1] - pxc[id-1]);
	}
	/* Copy of cptr->qa */
	for (i = 0; i < nq; i++){
	  qat[i] = cptr->qa[i];
	  cptr->qa[i] = qa0[i];
	}
	if (vvect){
	  /* Velocity vector */
	  k = ndm*j - 1;
	  for (id = 1; id <= nd; id++){
	    if (!mcoord[id-1]){
	      k = k + 1;
	      qqa[k] = usrvel(id, cptr);
	    }
	  }
	}
	else if (bvect){
	  /* Magnetic field vector */
	  k = ndm*j - 1;
	  for (id = 1; id <= nd; id++){
	    if (!mcoord[id-1]){
	      k = k + 1;
	      qqa[k] = usrbfield(id, cptr);
	    }
	  }
	}
	else{
	  /* Scalar */
	  if (usrps(cptr, qqa[j]))
	    ssa[j] = cptr->mat;
	  else
	    ssa[j] = bblank;
	}
	/* Restore cptr->qa */
	for (i = 0; i < nq; i++)
	  cptr->qa[i] = qat[i];
      }
    }
  }
  return;
}
//end fill.c

//file wr_graph.c
/*
+ wr_graph - Write graphics arrays to file.
     Description :
        This  writes   the  graphics  array   as  text  to   the  file
        expname.exp.  The first line contains the size in the form nx,
        ny, nz.  Subsequent lines contain  x coord, y coord,  z coord,
        value for each cell.
*/
void wr_graph()
{
  int id, i, it[MAXD];
  float xc[MAXD];
  /* Write Dimensions */
  fprintf(expfile, "%d ", ndm);
  /* Write sizes */
  for (id = 1; id <= ndm; id++)
    fprintf(expfile, "%d ", nl[id-1]);
  fprintf(expfile, "\n");
    /* Use "export" as title */
  fprintf(expfile, "%s \n", "export");
  /* Write coordinate names */
  for (id = 1; id <= nd; id++)
    if (!mcoord[id-1])
      fprintf(expfile, "%s ", coord[id-1]);
  /* Write plot item */
  fprintf(expfile, "%s \n ", item);
  /* Run through graphics array */
  for (i = 0; i < nqa; i++){
    i_to_it(ndm, nl, i+1, it);
    /* Coordinates */
    for (id = 1; id <= ndm; id++){
      xc[id-1] = xlgt[id-1] + (it[id-1] - 0.5)*dx[id-1];
      fprintf(expfile, "%.6E ", xc[id-1]);
    }
    /* Value */
      fprintf(expfile, "%.6E \n", qqa[i]);
  }
  return;
}
//end wr_graph.c
