#undef MAIN
#undef PLOT
#include <iostream>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "const.h"
#include "struct.h"
#include "para.h"
#include "funcs.h"
#include "plot.h"
struct datap {
  float xc[MAXD], value;
  struct datap *nxt;
};
char coordt[MAXD][SZITEM];
extern bool open_dat(const char tail[], char name[], bool techo);
extern bool get_n(char name[]);
extern bool get_cmd(char name[]);
extern bool get_r(float &value);
extern bool get_i(int &value);
extern void report();
void readxq(char file[]);
void redeye(char name[], char file[], int *ios);
using namespace std;
//file readxq.c
/*
+ readxq - read and plot tabulated data
     Description :
        This routine reads a 1-dimensional array of scalars from named
        text file, and plots them.
*/
void readxq(char file[])
{
  int static model_count(0);
  int ios;
  char name[SZNAM];
  /* invent a temporary model name */
  sprintf(name, "rad#%d", 10000000+model_count);
  model_count += 1;
  ios = 0;
  /* read data etc and load model */
  redeye(name, file, &ios);
  /* Master does it */
  if (mype == 0){
    if (ios == 0){
      vis_begin();
      model_seta((char*)"scalar", &ios);
      model_vis(&ios);
      vis_end();
    }
    /* remove the model */
    ios = 0;
    model_remove(name, &ios);
  }
  return;
}
//end readxq.c

//file redeye.c
/* 
+ redeye - plot solution, grid or velocities
      Description :
*/
void redeye(char name[], char file[], int *ios)
{
  char co[SZCO], un[SZCO], physic[SZNAM], vl[SZCQ], vu[SZCQ],
    title[SZNAM], geom[SZNAM], key[SZNAM];
  float qqb[MAXQB], tt, dxt[MAXD], xl[MAXD], xr[MAXD];
  int ndt, ga_size, nl[MAXD], i, id, tid, k, ssb[MAXSSB], nqb, 
    vmlt[3][3];
  bool doit;
  struct datap *start, *cdata, *tcdata;
  /* Master opens file */
  if (mype == 0){
    if (!open_dat(".xq", file, false)){
      cout << "Error: cannot open file" << '\n';
      return;
    }
  }
  if (!get_i(ndt)){
    if (mype == 0)
      cout << "Error: reading dimensions" << '\n';
    report();
    return;
  }
  for (id = 1; id <= ndt; id++){
    if (!get_i(nl[id-1])){
      if (mype == 0)
	cout << "Error: reading array dimension " << id << '\n';
      report();
      return;
    }
  }
  if (!get_n(title)){
    if (mype == 0)
      cout << "Error: reading title" << '\n';
    report();
    return;
  }
  for (id = 1; id <= ndt; id++){
    if (!get_cmd(coordt[id-1])){
      if (mype == 0)
	cout << "Error: reading coordinate name " << id << '\n';
      report();
      return;
    }
  }
  if (!get_cmd(key)){
    if (mype == 0)
      cout << "Error: reading variable name" << '\n';
    report();
    return;
  }
  /* Get data points */
  ga_size = 0;
  start = new datap;
  cdata = start;
  for (id = 1; id <= ndt; id++){
    xl[id-1] = big;
    xr[id-1] = -big;
  }
  while (get_r(cdata->xc[0])){
    for (id = 1; id < ndt; id++)
      get_r(cdata->xc[id]);
    for (id = 1; id <= ndt; id++){
      xl[id-1] = min(xl[id-1], cdata->xc[id-1]);
      xr[id-1] = max(xr[id-1], cdata->xc[id-1]);
    }
    get_r(cdata->value);
    ga_size++;
    cdata->nxt = new datap;
    tcdata = cdata;
    cdata = cdata->nxt;
  }
  /* Terminate list */
  if (ga_size > 1){
    delete tcdata->nxt;
    tcdata->nxt = NULL;
  }
  else{
    delete start;
    if (mype == 0)
      cout << "Error: no data" << '\n';
    return;
  }
  /* Only master does anything  */
  if (mype == 0){
    /* One dimensional is special */
    if (ndt == 1)
      nl[0] = ga_size;
    /* Create coordinate arrays */
    nxad = 0;
    for (id = 1; id <= ndt; id++)
      nxad += nl[id-1];
    xad = new float[nxad];
    xed = new float[nxad+ndt];
    /* Create graphics arrays */
    ssa = new int[ga_size];
    qqa = new float[ga_size];
    /* Calculate mesh spacing, set domain boundaries to cell edges */
    for (id = 1; id <= ndt; id++){
      dxt[id-1] = (xr[id-1] - xl[id-1])/(float)(nl[id-1] - 1);
      xl[id-1] -= 0.5*dxt[id-1];
      xr[id-1] += 0.5*dxt[id-1];
    }
    /* Check for masks */
    ndm = ndt;
    for (id = 1; id <= ndt; id++){
      if (mcoord[id-1] == true){
	/* Move mask off cell edge (temporary copy) */
	mvalt[id-1] = xl[id-1] + dxt[id-1]*(int)((mval[id-1]
	  - xl[id-1])/dxt[id-1]) + 0.01*dxt[id-1];
	ndm--;
      }
    }
    /* Set local plotting domain and mesh spacing (possibly masked) */
    tid = 1;
    for (id = 1; id <= ndt; id++){
      if (!mcoord[id-1]){
	xlgt[tid-1] = max(xlg[id-1], xl[id-1]);
	xlgt[tid-1] = xl[id-1] + 
          dxt[id-1]*(int)((xlgt[tid-1] - xl[id-1])/dxt[id-1]);
	  xrgt[tid-1] = min(xrg[id-1], xr[id-1]);
	xrgt[tid-1] = min(xrg[id-1], xr[id-1]);
	xrgt[tid-1] = xr[id-1] - 
          dxt[id-1]*(int)((xr[id-1] - xrgt[tid-1])
	  /dxt[id-1]);
	dx[tid-1] = dxt[id-1];
	nl[tid-1] = max(0, (int)((xrgt[tid-1] - xlgt[tid-1]
	  + 0.1*dxt[id-1])/dxt[id-1]));
	tid++;
      }
    }
    /* In one dimension count points (possible unequal mesh spacing) */
    if (ndt == 1){
      cdata = start;
      nl[0] = 0;
      while (cdata != NULL){
	if ((cdata->xc[0] >= xlgt[0]) && 
	    (cdata->xc[0] <= xrgt[0]))
	  nl[0]++;
	cdata = cdata->nxt;
      }
    }
    /* Fill coordinate arrays */
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
    /* Fill graphics arrays */
    cdata = start;
    i = -1;
    while (cdata != NULL){
      doit = true;
      tid = 1;
      for (id = 1; id <= ndt; id++){
	if ((cdata->xc[id-1] < xlgt[tid-1]) || 
	    (cdata->xc[id-1] > xrgt[tid-1]))
	  doit = false;
	if (mcoord[id-1]){
	  if ((cdata->xc[id-1] < mvalt[id-1] - 0.5*dxt[id-1])
	      || (cdata->xc[id-1] >= mvalt[id-1] + 0.5*dxt[id-1]))
	    doit = false;
	}
	else
	  tid++;
      }
      if (doit){
	i++;
	/* One dimensional uses x coord */
	if (ndt == 1)
	  xad[i] = cdata->xc[0];
	qqa[i] = cdata->value;
	ssa[i] = ffluid;
      }
      cdata = cdata->nxt;
    }
    /* Plot */
    nsa = i + 1;
    nnq = 1;
    nqa = nnq*nsa;
    nqb = 0;
    tt = 0.0;
    /* Reset xed for one dimensional */
    if (ndt == 1){
      xed[0] = xlgt[0];
      for (i = 0; i < nsa - 1; i++)
	xed[i+1] = xad[i] + 0.5*dx[0];
      xed[nsa] = xrgt[0];
    }
    /* Geometry */
    strcpy(geom, "XYZ");
    /* Construct composite strings for coordinates, units */
    strcpy(co, ";");
    strcpy(un, ";");
    for (id = 1; id <= ndm; id++){  
      strcat(co, coordt[id-1]);
      strcat(un, ";");
      strcat(co, ";");
    }
    for (i = 0; i <= nnq; i++)
      vu[i] = ';';
    vu[nnq+1] = '\0';
    /* Plot item */
    model_set0s((char*)"plot", key, ios);
    /* Composite strings for variable labels */
    strcpy(vl, ";");strcat(vl, key);strcat(vl, ";");
    /* Wash out mlt */
    for (i = 0; i < 3*nnq; i++)
      vmlt[0][i] = 0;
    /* Physics label */
    strcpy(physic, " ");
    model_put_ugrid(name, title, geom, tt, ndm, co, un, ncd, xed, xad, 
      physic, nnq, vl, vu, vmlt, nqa, ssa, qqa,  nqb, ssb, qqb, ios);
    /* Delete graphics arrays */
    delete [] ssa;
    delete [] xed;
    delete [] xad;
    delete [] qqa;
  }
  /* Delete data list */
  cdata = start;
  while (cdata != NULL){
    tcdata = cdata;
    cdata = cdata->nxt;
    delete tcdata;
  }
  return;
}
//end redeye.c

