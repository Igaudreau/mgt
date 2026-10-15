/* modelv.c:  model visualisation module 

     This module performs the abstractions from a model dataset,
     and fills a grid/scatter "vis" structure, which the vis module
     can plot.  Its external routuines all impose the Error Status
    scheme.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "const.h"
#include "model.h"
#include "vis.h"
#define MAXMASK MAXD
#define MAXVECTOR 6

static char scalar_name[SZNAM],
  mask_name[MAXMASK][SZNAM],vector_coord[MAXVECTOR][SZNAM],
  vector_scalar[MAXVECTOR][SZNAM];
static bool noscalar,mask_used[MAXMASK],
  mask_isindex[MAXMASK],mask_time,mask_time_isindex,
  nos,nov,nog,pick_point,vector_used[MAXVECTOR],
  nolog_scalar;
static int mask_index[MAXMASK],mask_time_index,pick_index;
static float mask_value[MAXMASK],mask_time_value;

static bool nomodelv=TRUE;
static char modelv_type[SZNAM];
static grid_model_struct *modelv_gp;
static scatter_model_struct *modelv_sp;

static bool first_time=TRUE;


/* EXPLICIT DECLARATIONS ADDED BY BM */
int model_reset();
int model_show();
void model_set0s(char name[], char value[], int *ios);
int model_set0i(char name[], int value, int *ios);
void model_set0r(char name[], float value, int *ios);
int model_set0l(char name[], bool value, int *ios);
int model_seta(char name[], int *ios);
int model_vis_either(scatter_model_struct *sp, scatter_vis_struct *svp, grid_vis_struct *gvp, int *ios);
int model_vis_crunch(scatter_model_struct *sp, grid_vis_struct *gvp, int *ios);
int model_vis_grid(grid_model_struct *gp, grid_vis_struct *gvp, int *ios);
int model_vis_scatter(scatter_model_struct *sp, scatter_vis_struct *svp, int *ios);
int model_selectvg(grid_model_struct *gp, int *ios);
int model_selectvs(scatter_model_struct *sp, int *ios);
int model_deselect();

extern bool csame(char a[], char b[]);
extern void ccopy(char a[], int offset, int n, char b[]);
extern bool cabrev(char a[], char b[], int *offset);
extern void ccat(char a[], int lenb, char b[]);
extern int ds_method(char met[], int *ios);
extern void user_setup_model(char geom[], int nd, char co[], char un[], 
  char phys[], int nq, char vl[], char vu[], int vmlt[][3], int *ios);
extern void user_check_scalar(char name[], int *ios);
extern void user_get_scalar(float q[], char name[], float *value);
extern bool geo_met(char geom[], int nd);
extern void geo_x(int lumd, float *z, float *x);
extern void vis_grid(grid_vis_struct *gvp, int *ios);
extern void vis_scatter(scatter_vis_struct *svp, int *ios);


/* ADMIN */

int model_reset() {
   int i;

   if (first_time) {
      for (i=0; i < MAXMASK; i++) {
         mask_used[i] = FALSE;
      }
      mask_time = TRUE;
      mask_time_isindex = TRUE;
      mask_time_index = 0;
      pick_point = FALSE;
      for (i=0; i < MAXVECTOR; i++) {
         vector_used[i] = FALSE;
      }
      noscalar = TRUE;
      nolog_scalar = TRUE;
      first_time = FALSE;
   }
   nos=TRUE;
   nov=TRUE;
   nog=TRUE;
   return 0;
}

int model_show() {
   int i,count;

   printf("Model parameters:-\n");
   if (noscalar) {
      printf("  No scalar selected\n");
   } else {
      printf("  Scalar '%s' selected\n", scalar_name);
   }
      if (!nolog_scalar) {
         printf("  Log(scalar)\n");
      }
   printf("  Vector pairs [ ");
   for (i=0; i < MAXVECTOR; i++) {
      if (vector_used[i]) {
         printf("%s=%s ", vector_coord[i], vector_scalar[i]);
      }
   }
   printf("]\n");
   count = 0;
   if (mask_time) {
      printf("  Mask time by ");
      if (mask_time_isindex) {
         printf("frame index %d\n", mask_time_index);
      } else {
         printf("nearest time %f\n", mask_time_value);
      }
   } else {
     printf("  No time mask\n");
   }
   for (i=0; i < MAXMASK; i++) {
      if (mask_used[i]) {
         count = count + 1;
         printf("  Mask coordinate '%s' ",mask_name[i]);
         if (mask_isindex[i]) {
            printf("index=%d\n", mask_index[i]);
         } else {
            printf("value=%g\n", mask_value[i]);
         }
      }
   }
   if (count < 1) {
      printf("  No coordinate masks\n");
   }
   if (pick_point) {
      printf("  No pick point\n");
   } else {
      printf("  Pick point ");
      printf("index=%d\n", pick_index);
   }
   return 0;
}

/* SET ATTRIBUTE ROUTINES */

void model_set0s(char name[], char value[], int *ios)
{
   int i,offset,ii;
   bool found;

   if (*ios) {
   } else if (csame("history",name)) {
      if (nomodelv) {
         printf("error:  no model selected\n");
         *ios = ERR;
      } else if (csame("grid",modelv_type)) {
         if (modelv_gp->nhist >= MAXHIST) {
            printf("error:  too many history reecords\n");
            *ios = ERR;
         } else {
            ccopy(value, 0, SZLIN, modelv_gp->hist[modelv_gp->nhist]);
            modelv_gp->nhist = modelv_gp->nhist + 1;
         }
      } else if (csame("scatter",modelv_type)) {
         if (modelv_sp->nhist >= MAXHIST) {
            printf("error:  too many history reecords\n");
            *ios = ERR;
         } else {
            ccopy(value, 0, SZLIN, modelv_sp->hist[modelv_sp->nhist]);
            modelv_sp->nhist = modelv_sp->nhist + 1;
         }
      }
   } else if (csame("plot",name)) {
      ccopy(value, 0, SZNAM, scalar_name);
      noscalar = FALSE;
   } else if (cabrev("vplot",name,&offset)) {
      found = FALSE;
      for (i=0; i < MAXVECTOR; i++) {
         if (vector_used[i]) {
            if (csame(&(name[offset]),vector_coord[i])) {
               ii = i;
               found = TRUE;
               break;
            }
         }
      }
      if (!found) {
         for (i=0; i < MAXVECTOR; i++) {
            if (!vector_used[i]) {
               ii = i;
               found = TRUE;
               break;
            }
         }
      }
      if (found) {
         ccopy(&(name[offset]), 0, SZNAM, vector_coord[ii]);
         ccopy(value, 0, SZNAM, vector_scalar[ii]);
         vector_used[ii] = TRUE;
      } else {
         printf("error:  too many coordinate vectors\n");
         *ios = ERR;
      }
   } else {
      printf("error:  unknown string attribute '%s' (model_set0s)\n",
        name);
      *ios = ERR;
   }
}

int model_set0i(name, value, ios)
   char name[];
   int *ios;
   int value;
{
   int i,offset,ii;
   bool found;

   if (*ios) {
   } else if (csame("maskti",name)) {
      if (value < 0) {
         printf("error:  illegal mask time index\n");
         *ios = ERR;
      } else {
         mask_time = TRUE;
         mask_time_isindex = TRUE;
         mask_time_index = value;
      }
   } else if (csame("picki",name)) {
      if (value < 0) {
         printf("error:  illegal pick index\n");
         *ios = ERR;
      } else {
         pick_point = TRUE;
         pick_index = value;
      }
   } else if (cabrev("maski",name,&offset)) {
      found = FALSE;
      for (i=0; i < MAXMASK; i++) {
         if (mask_used[i]) {
            if (csame(&(name[offset]),mask_name[i])) {
               ii = i;
               found = TRUE;
               break;
            }
         }
      }
      if (!found) {
         for (i=0; i < MAXMASK; i++) {
            if (!mask_used[i]) {
               ii = i;
               found = TRUE;
               break;
            }
         }
      }
      if (found) {
         ccopy(&(name[offset]), 0, SZNAM, mask_name[ii]);
         mask_isindex[ii] = TRUE;
         mask_index[ii] = value;
         mask_used[ii] = TRUE;
      } else {
         printf("error:  too many coordinate masks\n");
         *ios = ERR;
      }
   } else {
      printf("error:  unknown integer attribute '%s' (model_set0i)\n",
        name);
   }
   return 0;
}

void model_set0r(char name[], float value, int *ios)
{
   int i,offset,ii;
   bool found;
   if (*ios) {
   } else if (csame("maskt",name)) {
      mask_time = TRUE;
      mask_time_isindex = FALSE;
      mask_time_value = value;
   } else if (cabrev("mask",name,&offset)) {
      found = FALSE;
      for (i=0; i < MAXMASK; i++) {
         if (mask_used[i]) {
            if (csame(&(name[offset]),mask_name[i])) {
               ii = i;
               found = TRUE;
               break;
            }
         }
      }
      if (!found) {
         for (i=0; i < MAXMASK; i++) {
            if (!mask_used[i]) {
               ii = i;
               found = TRUE;
               break;
            }
         }
      }
      if (found) {
         ccopy(&(name[offset]), 0, SZNAM, mask_name[ii]);
         mask_isindex[ii] = FALSE;
         mask_value[ii] = value;
         mask_used[ii] = TRUE;
      } else {
         printf("error:  too many coordinate masks\n");
         *ios = ERR;
      }
   } else {
      printf("error:  unknown real attribute '%s' (model_set0r)\n",
        name);
   }
}

int model_set0l(name, value, ios)
   char name[];
   int *ios;
   bool value;
{

   if (*ios) {
   } else {
      printf("error:  unknown logical attribute '%s' (model_set0l)\n",
        name);
   }
   return 0;
}

int model_seta(char name[], int *ios)
{
   int offset,i;
   bool found;

   if (*ios) {
   } else if (csame("binary_data",name)) {
      ds_method("binary", ios);
   } else if (csame("text_data",name)) {
      ds_method("text", ios);
   } else if (csame("unpick",name)) {
      pick_point = FALSE;
   } else if (csame("unmaskt",name)) {
      mask_time = FALSE;
   } else if (cabrev("unmask",name,&offset)) {
      found = FALSE;
      for (i=0; i < MAXMASK; i++) {
         if (mask_used[i]) {
            if (csame(&(name[offset]), mask_name[i])) {
               mask_used[i] = FALSE;
               found = TRUE;
               break;
            } else if (csame(&(name[offset]), "*")) {
               mask_used[i] = FALSE;
               found = TRUE;
            }
         }
      }
      if (!found) {
         printf("error:  cannot unmask coordinate '%s'\n", 
           &(name[offset]));
         *ios = ERR;
      }
   } else if (csame("scalar",name)) {
      nos = FALSE;
   } else if (csame("noscalar",name)) {
      nos = TRUE;
   } else if (csame("vector",name)) {
      nov = FALSE;
   } else if (csame("novector",name)) {
      nov = TRUE;
   } else if (csame("grid",name)) {
      nog = FALSE;
   } else if (csame("nogrid",name)) {
      nog = TRUE;
   } else if (csame("log",name)) {
      nolog_scalar = FALSE;
   } else if (csame("nolog",name)) {
      nolog_scalar = TRUE;
   } else {
      printf("error:  unknown assertion '%s' (model_seta)\n",
        name);
   }
   return 0;
}

/* MODEL ABSTRACTION */

int model_vis_either(sp, svp, gvp, ios)
   scatter_model_struct *sp;
   scatter_vis_struct *svp;
   grid_vis_struct *gvp;
   int *ios;
{
/* pick means we convert to 1-dimensional time grid */
   if (*ios) {
   } else if (pick_point) {
      if (mask_time) {
         printf("warning:  ignore mask_time scatter in pick mode!\n");
      }
      model_vis_crunch(sp, gvp, ios);
   } else {
      model_vis_scatter(sp, svp, ios);
   }
   return 0;
}

int model_vis_crunch(sp, gvp, ios)
   scatter_model_struct *sp;
   grid_vis_struct *gvp;
   int *ios;
{
   int id,i,k,scalar_index,ioscalar,ipt;
   scatter_frame_struct *fp;
   char co[SZCO],un[SZCO],vl[SZCQ],vu[SZCQ];
   float q;
   scalar_index = -1;		// BM

   gvp->crunched = TRUE;
/* construct composite coord string */
   ccopy(":", 0, SZCO, co);
   for (id=0; id < sp->nd; id++) {
      ccat(sp->coords[id], SZCO, co);
      ccat(":", SZCO, co);
   }
/* construct composite units string */
   ccopy(":", 0, SZCO, un);
   for (id=0; id < sp->nd; id++) {
      ccat(sp->units[id], SZCO, un);
      ccat(":", SZCO, un);
   }
/* construct composite vlabels string */
   ccopy(":", 0, SZCQ, vl);
   for (i=0; i < sp->nq; i++) {
      ccat(sp->vlabels[i], SZCQ, vl);
      ccat(":", SZCQ, vl);
   }
/* construct composite vunits string */
   ccopy(":", 0, SZCQ, vu);
   for (i=0; i < sp->nq; i++) {
      ccat(sp->vunits[i], SZCQ, vu);
      ccat(":", SZCQ, vu);
   }
/* setup user in preparation for inquisition */
   user_setup_model(sp->geom, sp->nd, co, un, sp->physics, sp->nq,
     vl, vu, sp->vmlt, ios);
/* ok to proceed? */
   if (*ios) {
   } else if (!pick_point) {
      printf("error:  not pick_point mode (model_vis_crunch)\n");
      *ios = ERR;
/* pick point */
   } else {
/*   we ignored grid as it is unstructured */
      gvp->nogrid = TRUE;
/*   pick by index (make it safe) */
      if (pick_index < 0) {
         pick_index = 0;
      } else if (pick_index > (sp->npt-1)) {
         pick_index = sp->npt - 1;
      }
/*   title */
      ccopy(sp->title, 0, SZLIN, gvp->title);
/*   geom */
      ccopy(sp->geom, 0, SZNAM, gvp->geom);
/*   number of dimensions (time only) */
      gvp->nd = 1;
#ifdef CURV
      gvp->curv = FALSE;
#endif
/*   time coordinate labels */
      ccopy("t", 0, SZNAM, gvp->coords[0]);
      ccopy("time", 0, SZNAM, gvp->clabel[0]);
      ccopy(" ", 0, SZNAM, gvp->units[0]);
      gvp->nc[0] = sp->nt;
      gvp->nodal[0] = TRUE;
#ifdef CURV
      gvp->kad[gvp->ndo] = 0;
#endif
/*   number of coordinate dimensions */
      gvp->ndc = 0;
/*   number of frozen coordinates */
      gvp->ndo = sp->nd + sp->ndo;
/*   space coordinate labels */
      for (id=0; id < gvp->ndo; id++) {
         ccopy(sp->coords[id], 0, SZNAM, gvp->coords[id+1]);
         ccopy(sp->clabel[id], 0, SZNAM, gvp->clabel[id+1]);
         ccopy(sp->units[id], 0, SZNAM, gvp->units[id+1]);
#ifdef CURV
         gvp->kad[id] = id+1;
#endif
      }
/*   scalar */
      gvp->noscalar = noscalar || nos;
      if (!gvp->noscalar) {
/*      search for explicitly named variable */
         if (nolog_scalar) {
            ccopy(scalar_name, 0, SZNAM, gvp->qlabel);
         } else {
            //sprintf(gvp->qlabel, "log(%s)", scalar_name);
	    snprintf(gvp->qlabel, sizeof(gvp->qlabel), "log(%.25s)", scalar_name);		// BM for safety
         }
         scalar_index = -1;
         for (i=0; i < sp->nq; i++) {
            if (csame(scalar_name, sp->vlabels[i])) {
               scalar_index = i;
               ccopy(sp->vunits[i], 0, SZNAM, gvp->qunits);
               break;
            }
         }
/*      check for user-defined scalar */
         if (scalar_index < 0) {
            user_check_scalar(scalar_name, &ioscalar);
            if (ioscalar != OK) {
               printf("warning:  scalar=%s undefined\n", scalar_name);
            }
         }
      }
/*   no vector */
      gvp->novector = TRUE;
/*   xe is missing (because nodal) */
      gvp->xe[0] = NIL;
/*   xa memory+store */
      if (!*ios) {
	gvp->xa[0] = (float*)calloc(gvp->nc[0], sizeof(float));
         *ios = (gvp->xa[0] == NULL);
         if (*ios) {
            printf("error:  no room for xa\n");
         }
      }
/*   sa memory */
      if (!*ios) {
         gvp->sa = (int*)calloc(gvp->nc[0], sizeof(int));
         *ios = (gvp->sa == NULL);
         if (*ios) {
            printf("error:  not enough memory for vis/sa\n");
         }
      }
/*   qa memory */
      if (*ios) {
      } else if (!gvp->noscalar) {
         gvp->qa = (float*)calloc(gvp->nc[0], sizeof(float));
         *ios = (gvp->qa == NULL);
         if (*ios) {
            printf("error:  not enough memory for vis/qa\n");
         }
      }
/*   va memory void (because boring) */
      gvp->va = NULL;
/*   grab time coordinates and data */
      if (!*ios) {
         gvp->xmin[0] = BIG;
         gvp->xmax[0] = -BIG;
         gvp->qmin = BIG;
         gvp->qmax = -BIG;
/*
         kpt = pick_index*sp->nd;
*/
         ipt = pick_index*sp->nq;
         for (i=0; i < sp->nt; i++) {
            fp = &(sp->frame[i]);
            gvp->xa[0][i] = fp->t;
            if (fp->t < gvp->xmin[0]) {
               gvp->xmin[0] = fp->t;
            }
            if (fp->t > gvp->xmax[0]) {
               gvp->xmax[0] = fp->t;
            }
            if (!gvp->noscalar) {
               gvp->sa[i] = fp->spt[pick_index];
               if (gvp->sa[i] == FLUID) {
                  if (scalar_index < 0) {
                     user_get_scalar(&(fp->qpt[ipt]),
                       scalar_name, &q);
                     gvp->qa[i] = q;
                  } else {
                     k = sp->nq*pick_index + scalar_index;
                     gvp->qa[i] = fp->qpt[k];
                  }
                  if (nolog_scalar) {
                  } else if (gvp->qa[i] > 0.0) {
                     gvp->qa[i] = log10(gvp->qa[i]);
                  } else {
#ifdef BUG
                     gvp->sa[i] == 0;
#else
                     gvp->sa[i] = 0;
#endif
                  }
                  if (gvp->qa[i] < gvp->qmin) {
                     gvp->qmin = gvp->qa[i];
                  }
                  if (gvp->qa[i] > gvp->qmax) {
                     gvp->qmax = gvp->qa[i];
                  }
               }
            }
         }
/*      determine whether we have ranges */
#ifdef CURV
#ifdef DEGEN_BUG
         gvp->noqr = gvp->qmax <= gvp->qmin;
#else
         if (gvp->qmax == gvp->qmin) {
#ifdef DEBUG
            fprintf(stderr, "model_vis_crunch:  degenerate\n");
#endif
            if (gvp->qmin == 0.0) {
               gvp->qmin = -1.0;
               gvp->qmax = 1.0;
#ifdef DEBUG
               fprintf(stderr, "..degenerate zero\n");
#endif
            } else {
#ifdef DEBUG
               fprintf(stderr, "..degenerate value=%f\n", gvp->qmin);
#endif
               gvp->qmin = gvp->qmax*0.9;
               gvp->qmax = gvp->qmax*1.1;
            }
         } else if (gvp->qmax < gvp->qmin) {
#ifdef DEBUG
            fprintf(stderr, "model_vis_crunch:  no data\n");
#endif
            gvp->qmin = -1.0;
            gvp->qmax = +1.0;
         }
         gvp->noqr = FALSE;
#endif
#else
         gvp->noqr = gvp->qmax > gvp->qmin;
#endif
         gvp->noxr[0] = gvp->xmax[0] > gvp->xmin[0];
      }
   }
/* only submit to vis if error-free */
   if (!*ios) {
      vis_grid(gvp, ios);
   }
/* revoke scalar/vector requests */
   nog = TRUE;
   nos = TRUE;
   nov = TRUE;
   return 0;
}

int model_vis_grid(gp, gvp, ios)
   grid_model_struct *gp;
   grid_vis_struct *gvp;
   int *ios;
{
   int id,ido,jdo,i,mask[MAXD],ix,nx,iy,ny,iz,nz,k,ki,kv,j,jv,jd,it,ns,
     nv,scalar_index,vector_index[MAXD],coord_index[MAXD],
     jdv,idv,indexv,ioscalar;
   bool masked[MAXD];
   float dx,dxmin,dtmin,v,q,dt;
   grid_frame_struct *fp;
   char vector_name[MAXD][SZNAM],co[SZCO],un[SZCO],vl[SZCQ],vu[SZCQ];
  scalar_index = -1;	// BM

#ifdef CURV
   if (gp->curv) {
      if (geo_met(gp->geom, (gp->nd+gp->ndo))) {
      } else {
         printf("Warning:  trouble in model_vis_grid() with geo_met\n");
      }
   }
   gvp->curv = gp->curv;
   ccopy(gp->geom, 0, SZNAM, gvp->geom);
#endif
   gvp->crunched = FALSE;
/* construct composite coord string */
   ccopy(":", 0, SZCO, co);
   for (id=0; id < gp->nd; id++) {
      ccat(gp->coords[id], SZCO, co);
      ccat(":", SZCO, co);
   }
/* construct composite units string */
   ccopy(":", 0, SZCO, un);
   for (id=0; id < gp->nd; id++) {
      ccat(gp->units[id], SZCO, un);
      ccat(":", SZCO, un);
   }
/* construct composite vlabels string */
   ccopy(":", 0, SZCQ, vl);
   for (i=0; i < gp->nq; i++) {
      ccat(gp->vlabels[i], SZCQ, vl);
      ccat(":", SZCQ, vl);
   }
/* construct composite vunits string */
   ccopy(":", 0, SZCQ, vu);
   for (i=0; i < gp->nq; i++) {
      ccat(gp->vunits[i], SZCQ, vu);
      ccat(":", SZCQ, vu);
   }
/* setup user in preparation for inquisition */
   user_setup_model(gp->geom, gp->nd, co, un, gp->physics, gp->nq,
     vl, vu, gp->vmlt, ios);
/* title */
   ccopy(gp->title, 0, SZLIN, gvp->title);
/* grid? */
   gvp->nogrid = nog;
/* scalar */
   gvp->noscalar = noscalar||nos;
   if (!gvp->noscalar) {
      if (nolog_scalar) {
         ccopy(scalar_name, 0, SZNAM, gvp->qlabel);
      } else {
         //sprintf(gvp->qlabel, "log(%s)", scalar_name);
	 snprintf(gvp->qlabel, sizeof(gvp->qlabel), "log(%.25s)", scalar_name);		// BM for safety
      }
/*   search for explicit name */
      scalar_index = -1;
      for (i=0; i < gp->nq; i++) {
         if (csame(scalar_name, gp->vlabels[i])) {
            scalar_index = i;
            ccopy(gp->vunits[i], 0, SZNAM, gvp->qunits);
            break;
         }
      }
/*   check for user-defined scalar */
      if (scalar_index < 0) {
         user_check_scalar(scalar_name, &ioscalar);
         if (ioscalar != OK) {
            printf("warning:  user scalar=%s undefined\n", scalar_name);
         }
      }
   }
/* whether vector is there */
   gvp->novector = nov;
/* construct firm masks from spec */
   gvp->nd = 0;
   ns = 1;
   for (id=0; id < gp->nd; id++) {
      masked[id] = FALSE;
      for (i=0; i < MAXMASK; i++) {
         if (mask_used[i]) {
            if (csame(mask_name[i],gp->coords[id])) {
               masked[id] = TRUE;
               if (mask_isindex[i]) {
                  mask[id] = mask_index[i];
               } else {
                  dxmin = BIG;
                  for (j=0; j < gp->nc[0][id]; j++) {
                     dx = mask_value[i]-gp->xa[0][id][j];
                     dx = dx*dx;
                     if (dx < dxmin) {
                        mask[id] = j;
                        dxmin = dx;
                     }
                  }
               }
               break;
            }
         }
      }
      if (!masked[id]) {
         gvp->nd = gvp->nd + 1;
         coord_index[gvp->nd-1] = id;
         ns = ns*gp->nc[0][id];
         mask[id] = -1;
      }
   }
/* number of space dimensions */
   gvp->ndc = gvp->nd;
/* time mask? */
   if (!mask_time) {
      gvp->nd = gvp->nd + 1;
   }
/* work through axes */
   gvp->ndo = 0;
   jd = 0;
   jdv = 0;
   for (id=0; id < gp->nd; id++) {
/*   find vector index for this coordinate */
      indexv = -1;
      for (i=0; i < MAXVECTOR; i++) {
         if (vector_used[i]) {
            if (csame(gp->coords[id],vector_coord[i])) {
               indexv = i;
               break;
            }
         }
      }
/*   coordinate masked?  add onto end of list */
      if (masked[id]) {
         gvp->ndo = gvp->ndo + 1;
         ido = gvp->ndo - 1;
         jdo = gvp->nd + ido;
#ifdef CURV
         gvp->kad[jdo] = gp->kad[id];
#endif
         ccopy(gp->coords[id], 0, SZNAM, gvp->coords[jdo]);
         ccopy(gp->clabel[id], 0, SZNAM, gvp->clabel[jdo]);
         ccopy(gp->units[id], 0, SZNAM, gvp->units[jdo]);
         gvp->xo[ido] = (gp->xa[0][id])[mask[id]];
         idv = gvp->ndc + ido;
         if (indexv < 0) {
            ccopy(" ", 0, SZNAM, vector_name[idv]);
         } else {
            ccopy(vector_scalar[indexv], 0, SZNAM, vector_name[idv]);
         }
         vector_index[idv] = indexv;
         gvp->vector_index[idv] = jdo;
         ccopy(gp->coords[id], 0, SZNAM, gvp->vlabels[idv]);
/*   suitable coordinate */
      } else {
#ifdef CURV
         gvp->kad[jd] = gp->kad[id];
#endif
         ccopy(gp->coords[id], 0, SZNAM, gvp->coords[jd]);
         ccopy(gp->clabel[id], 0, SZNAM, gvp->clabel[jd]);
         ccopy(gp->units[id], 0, SZNAM, gvp->units[jd]);
         if (indexv < 0) {
            ccopy(" ", 0, SZNAM, vector_name[jdv]);
         } else {
            ccopy(vector_scalar[indexv], 0, SZNAM, vector_name[jdv]);
         }
         vector_index[jdv] = indexv;
         gvp->vector_index[jdv] = jd;
         ccopy(gp->coords[id], 0, SZNAM, gvp->vlabels[jdv]);
         gvp->nc[jd] = gp->nc[0][id];
         if (!*ios) {
            gvp->xe[jd] = (float*)calloc((gvp->nc[jd]+1),sizeof(float));
            *ios = (gvp->xe[jd] == NULL);
            if (*ios) {
               printf("error:  no room for xe\n");
               break;
            }
         }
         if (!*ios) {
            gvp->xa[jd] = (float*)calloc(gvp->nc[jd], sizeof(float));
            *ios = (gvp->xa[jd] == NULL);
            if (*ios) {
               printf("error:  no room for xa\n");
               break;
            }
         }
         if (!*ios) {
            for (i=0; i < gvp->nc[jd]+1; i++) {
               gvp->xe[jd][i] = gp->xe[0][id][i];
            }
            for (i=0; i < gvp->nc[jd]; i++) {
               gvp->xa[jd][i] = gp->xa[0][id][i];
            }
            gvp->xmin[jd] = gvp->xe[jd][0];
            gvp->xmax[jd] = gvp->xe[jd][gvp->nc[jd]];
            gvp->nodal[jd] = FALSE;
         }
         jd = jd + 1;
         jdv = jdv + 1;
      }
   }
/* copy over frozen space coordinates */
   jdv = gvp->ndc;
   for (id=0; id < gp->ndo; id++) {
      jd = gp->nd + id;
      gvp->ndo = gvp->ndo + 1;
      ido = gvp->ndo - 1;
      jdo = gvp->nd + ido;
#ifdef CURV
      gvp->kad[jdo] = gp->kad[jd];
#endif
/*   find vector index for this coordinate */
      indexv = -1;
      for (i=0; i < MAXVECTOR; i++) {
         if (vector_used[i]) {
            if (csame(gp->coords[jd],vector_coord[i])) {
               indexv = i;
               break;
            }
         }
      }
      ccopy(gp->coords[jd], 0, SZNAM, gvp->coords[jdo]);
      ccopy(gp->clabel[jd], 0, SZNAM, gvp->clabel[jdo]);
      ccopy(gp->units[jd], 0, SZNAM, gvp->units[jdo]);
      gvp->xo[ido] = gp->xo[id];
      if (indexv < 0) {
         ccopy(" ", 0, SZNAM, vector_name[jdv]);
      } else {
         ccopy(vector_scalar[indexv], 0, SZNAM, vector_name[jdv]);
      }
      vector_index[jdv] = indexv;
      gvp->vector_index[jdv] = jdo;
      ccopy(gp->coords[jd], 0, SZNAM, gvp->vlabels[jdv]);
      jdv = jdv + 1;
   }
   gvp->ndv = jdv;
/* reject degenerate case */
   if (*ios) {
   } else if (gvp->nd < 1) {
      printf("error:  masking makes dataset degenerate\n");
      *ios = ERR;
/* time frames */
   } else if (mask_time) {
      gvp->ndo = gvp->ndo + 1;
      ido = gvp->ndo - 1;
      jdo = gvp->nd + ido;
#ifdef CURV
      gvp->kad[jdo] = gp->nd + gp->ndo;
#endif
      ccopy("t", 0, SZNAM, gvp->coords[jdo]);
      ccopy("time", 0, SZNAM, gvp->clabel[jdo]);
      ccopy(" ", 0, SZNAM, gvp->units[jdo]);
      if (!mask_time_isindex) {
         dtmin = BIG;
         mask_time_index = -1;
         for (i=0; i < gp->nt; i++) {
            fp = &(gp->frame[i]);
            dt = fp->t - mask_time_value;
            dt = dt*dt;
            if (dt < dtmin) {
               mask_time_index = i;
               dtmin = dt;
            }
         }
      } else {
         if (mask_time_index < 0) {
            printf("warning:  invalid time_index, default zero\n");
            mask_time_index = 0;
         }
      }
      fp = &(gp->frame[mask_time_index]);
      gvp->xo[ido] = fp->t;
   } else {
      mask_time_index = -1;
      jd = gvp->nd - 1;
#ifdef CURV
      gvp->kad[jd] = gp->nd + gp->ndo;
#endif
      ccopy("t", 0, SZNAM, gvp->coords[jd]);
      ccopy("time", 0, SZNAM, gvp->clabel[jd]);
      ccopy(" ", 0, SZNAM, gvp->units[jd]);
      gvp->nodal[jd] = TRUE;
      gvp->nc[jd] = gp->nt;
      ns = ns*gvp->nc[jd];
/*   memory for xe and xa */
      gvp->xe[jd] = (float*)calloc((gvp->nc[jd]+1), sizeof(float));
      *ios = (gvp->xe[jd] == NULL);
      if (*ios) {
         printf("error:  no room for xe\n");
      }
      if (!*ios) {
         gvp->xa[jd] = (float*)calloc(gvp->nc[jd], sizeof(float));
         *ios = (gvp->xa[jd] == NULL);
         if (*ios) {
            printf("error:  no room for xa\n");
         }
      }
      if (!*ios) {
         for (i=0; i < gp->nt; i++) {
            fp = &(gp->frame[i]);
            gvp->xa[jd][i] = fp->t;
            if (i == 0) {
               gvp->xmin[jd] = fp->t;
            } else if (i == gp->nt-1) {
               gvp->xmax[jd] = fp->t;
            }
         }
         if (gp->nt > 1) {
            dt = 0.5*(gvp->xa[jd][1] - gvp->xa[jd][0]);
         } else {
            dt = 0.0;
         }
         gvp->xe[jd][0] = gvp->xa[jd][0] - dt;
         for (i=1; i < gp->nt; i++) {
            gvp->xe[jd][i] = 0.5*(gvp->xa[jd][i-1] + gvp->xa[jd][i]);
         }
         if (gp->nt > 1) {
            dt = 0.5*(gvp->xa[jd][gp->nt-1] - gvp->xa[jd][gp->nt-2]);
         } else {
            dt = 0.0;
         }
         gvp->xe[jd][gp->nt] = gvp->xa[jd][gp->nt-1] + dt;
         gvp->xmin[jd] = gvp->xe[jd][0];
         gvp->xmax[jd] = gvp->xe[jd][gp->nt];
      }
   }
#ifdef CURV
/*
   printf("model_vis_grid:-\n");
   for (id=0; id < gvp->nd; id++) {
      printf("id=%d,coords=%s,kad=%d\n", id, gvp->coords[id],
        gvp->kad[id]);
   }
   for (ido=0; ido < gvp->ndo; ido++) {
      id = gvp->nd + ido;
      printf("id=%d,coords=%s,kad=%d,xo=%f\n", id, gvp->coords[id],
        gvp->kad[id], gvp->xo[ido]);
   }
*/
#ifdef CURV_MUVI
   if (gvp->curv) {
/*
      printf("model_vis_grid:  gvp->geom=%s,nd=%d,ndc=%d,ndo=%d\n",
       gvp->geom, gvp->nd, gvp->ndc, gvp->ndo);
      for (i=0; i < (gvp->nd+gvp->ndo); i++) {
         printf("kad[%d]=%d\n", i, gvp->kad[i]);
      }
*/
/*    [r,z,phi] with some masking */
      if (csame("MUVI_RZP",gvp->geom)) {
/*       Masked to become [r,phi] */
         if ((gvp->nd==2) &&(gvp->ndc==2) &&
           (gvp->ndo==2) && (gvp->kad[gvp->nd]==1)) {
            gvp->ndo = 0;
            gvp->ndc = gvp->nd;
            ccopy("MUVI_RP", 0, SZNAM, gvp->geom);
            printf("Note:  now using surrogate geometry MUVI_RP\n");
            if (!geo_met(gvp->geom, gvp->nd)) {
               printf("model_vis_grid:  trouble with that\n");
            }
/*       Masked to become [r,z] */
         } else if ((gvp->nd==2) &&(gvp->ndc==2) &&
           (gvp->ndo==2) && (gvp->kad[gvp->nd]==2)) {
            gvp->ndo = 0;
            gvp->ndc = gvp->nd;
            gvp->curv = FALSE;
            ccopy("cartesian", 0, SZNAM, gvp->geom);
            printf("Note:  reverting surrogate cartesian geometry\n");
         }
      }
   }
#endif
#endif
/* sa memory */
   if (!*ios) {
      gvp->sa = (int*)calloc(ns, sizeof(int));
      *ios = (gvp->sa == NULL);
      if (*ios) {
         printf("error:  not enough memory for vis/sa\n");
      }
   }
/* qa memory */
   if (*ios) {
   } else if (!gvp->noscalar) {
      gvp->qa = (float*)calloc(ns, sizeof(float));
      *ios = (gvp->qa == NULL);
      if (*ios) {
         printf("error:  not enough memory for vis/qa\n");
      }
   }
/* va memory */
   if (*ios) {
   } else if (!gvp->novector) {
      nv = ns*gvp->ndv;
      gvp->va = (float*)calloc(nv, sizeof(float));
      *ios = (gvp->va == NULL);
      if (*ios) {
         printf("error:  not enough memory for vis/va\n");
      }
   }
/* invent fictive added coordinate sizes */
   nx = gp->nc[0][0];
   if (gp->nd >= 2) {
      ny = gp->nc[0][1];
      if (gp->nd == 3) {
         nz = gp->nc[0][2];
      } else {
         nz = 1;
         mask[2] = 0;
      }
   } else {
      ny = 1;
      nz = 1;
      mask[2] = 0;
      mask[1] = 0;
   }
/* transfer vector indices */
   if (!gvp->novector) {
      for (id=0; id < gvp->ndv; id++) {
         if (vector_index[id] < 0) {
            ccopy(" ", 0, SZNAM, vector_name[id]);
         } else {
            vector_index[id] = -1;
            for (i=0; i < gp->nq; i++) {
               if (csame(vector_name[id],gp->vlabels[i])) {
                  vector_index[id] = i;
                  break;
               }
            }
/*         check for user-defined scalar */
            if (vector_index[id] < 0) {
               user_check_scalar(vector_name[id], &ioscalar);
               if (ioscalar != OK) {
                  printf("warning:  user vector=%s undefined\n",
                    vector_name[id]);
               }
            }
         }
      }
   }
/* abstract switch,scalar and/or vector from original */
   if (!*ios) {
      j = 0;
      for (it=0; it < gp->nt; it++) {
         fp = &(gp->frame[it]);
         if ((!mask_time) || (mask_time_index == it)) {
            for (iz=0; iz < nz; iz++) {
               if ((!masked[2]) || (mask[2] == iz)) {
                  for (iy=0; iy < ny; iy++) {
                     if ((!masked[1]) || (mask[1] == iy)) {
                        for (ix=0; ix < nx; ix++) {
                           if ((!masked[0]) || (mask[0] == ix)) {
/*                           index in original grid array */
                              k = (iz*ny + iy)*nx + ix;
/*                           switch */
                              gvp->sa[j] = fp->sa[k];
/*                           fluid case only */
                              if (gvp->sa[j] == FLUID) {
/*                              scalar field */
                                 if (!gvp->noscalar) {
                                    if (scalar_index < 0) {
                                       ki = k*gp->nq;
                              user_get_scalar(&(fp->qa[ki]),
                                scalar_name, &q);
                                       gvp->qa[j] = q;
                                    } else {
                                       ki = k*gp->nq + scalar_index;
                                       gvp->qa[j] = fp->qa[ki];
                                    }
                                    if (nolog_scalar) {
                                    } else if (gvp->qa[j] > 0.0) {
                                       gvp->qa[j] = log10(gvp->qa[j]);
                                    } else {
                                       gvp->sa[j] = 0;
                                    }
                                 }
/*                              vector field */
                                 if (!gvp->novector) {
                                    for (jd=0; jd < gvp->ndv; jd++) {
                                       jv = j*gvp->ndv + jd;
                                       if (vector_index[jd] < 0) {
                                          ki = k*gp->nq;
                                          user_get_scalar(&(fp->qa[ki]),
                                            vector_name[jd], &q);
                                          gvp->va[jv] = q;
                                       } else {
                                          i = vector_index[jd];
                                          kv = k*gp->nq + i;
                                          gvp->va[jv] = fp->qa[kv];
                                       }
                                    }
                                 }
                              }
/*                           increment output index */
                              j = j + 1;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
/*   find scalar range? */
      if (!gvp->noscalar) {
         gvp->qmin = BIG;
         gvp->qmax = -BIG;
         for (j=0; j < ns; j++) {
            if (gvp->sa[j] == FLUID) {
               if (gvp->qa[j] < gvp->qmin) {
                  gvp->qmin = gvp->qa[j];
               }
               if (gvp->qa[j] > gvp->qmax) {
                  gvp->qmax = gvp->qa[j];
               }
            }
         }
#ifdef DEGEN_BUG
         gvp->noqr = gvp->qmax <= gvp->qmin;
#else
         if (gvp->qmax == gvp->qmin) {
#ifdef DEBUG
            fprintf(stderr, "model_vis_grid:  degenerate data\n");
#endif
            if (gvp->qmin == 0.0) {
#ifdef DEBUG
               fprintf(stderr, "..degenerate zero\n");
#endif
               gvp->qmin = -1.0;
               gvp->qmax = 1.0;
            } else {
#ifdef DEBUG
               fprintf(stderr, "..degenerate value=%f\n", gvp->qmin);
#endif
               gvp->qmin = gvp->qmax*0.9;
               gvp->qmax = gvp->qmax*1.1;
            }
         } else if (gvp->qmax < gvp->qmin) {
#ifdef DEBUG
            fprintf(stderr, "model_vis_grid:  no data\n");
#endif
            gvp->qmin = -1.0;
            gvp->qmax = +1.0;
         }
         gvp->noqr = FALSE;
#endif
      }
/*   find vector range? */
      if (!gvp->novector) {
         for (id=0; id < gvp->ndv; id++) {
            gvp->vmax[id] = 0.0;
         }
         for (j=0; j < ns; j++) {
            if (gvp->sa[j] == FLUID) {
               for (id=0; id < gvp->ndv; id++) {
                  jv = j*gvp->ndv + id;
                  v = gvp->va[jv];
                  if (v < 0.0) {
                     v = -v;
                  }
                  if (v > gvp->vmax[id]){
                     gvp->vmax[id] = v;
                  }
               }
            }
         }
         for (id=0; id < gvp->ndv; id++) {
            gvp->novmax[id] = gvp->vmax[id] == 0.0;
         }
      }
   }
/* revoke scalar/vector requests */
   nog = TRUE;
   nos = TRUE;
   nov = TRUE;
/* only submit to vis if error free */
   if (!*ios) {
      vis_grid(gvp, ios);
   }
   (void)coord_index; /* Silences the unused-but-set-variable warning BM */
   return 0;
}

int model_vis_scatter(sp, svp, ios)
   scatter_model_struct *sp;
   scatter_vis_struct *svp;
   int *ios;
{
   int id,i,k,j,ns,nv,jd,jv,kv,ido,ioscalar,ki,
     scalar_index,vector_index[MAXD],indexv;
   float dt,dtmin,v,q;
   scatter_frame_struct *fp;
   char vector_name[MAXD][SZNAM],co[SZCO],un[SZCO],vl[SZCQ],vu[SZCQ];
   scalar_index = -1;		// BM

/* construct composite coord string */
   ccopy(":", 0, SZCO, co);
   for (id=0; id < sp->nd; id++) {
      ccat(sp->coords[id], SZCO, co);
      ccat(":", SZCO, co);
   }
/* construct composite units string */
   ccopy(":", 0, SZCO, un);
   for (id=0; id < sp->nd; id++) {
      ccat(sp->units[id], SZCO, un);
      ccat(":", SZCO, un);
   }
/* construct composite vlabels string */
   ccopy(":", 0, SZCQ, vl);
   for (i=0; i < sp->nq; i++) {
      ccat(sp->vlabels[i], SZCQ, vl);
      ccat(":", SZCQ, vl);
   }
/* construct composite vunits string */
   ccopy(":", 0, SZCQ, vu);
   for (i=0; i < sp->nq; i++) {
      ccat(sp->vunits[i], SZCQ, vu);
      ccat(":", SZCQ, vu);
   }
/* setup user in preparation for inquisition */
   user_setup_model(sp->geom, sp->nd, co, un, sp->physics, sp->nq,
     vl, vu, sp->vmlt, ios);
/* title */
   ccopy(sp->title, 0, SZLIN, svp->title);
/* geom */
   ccopy(sp->geom, 0, SZNAM, svp->geom);
/* scalar */
   svp->noscalar = noscalar||nos;
   if (!svp->noscalar) {
/*   search for explicit name */
      if (nolog_scalar) {
         ccopy(scalar_name, 0, SZNAM, svp->qlabel);
      } else {
         //sprintf(svp->qlabel, "log(%s)", scalar_name);
         snprintf(svp->qlabel, sizeof(svp->qlabel), "log(%.25s)", scalar_name);	// For safety BM
      }
      scalar_index = -1;
      for (i=0; i < sp->nq; i++) {
         if (csame(scalar_name, sp->vlabels[i])) {
            scalar_index = i;
            ccopy(sp->vunits[i], 0, SZNAM, svp->qunits);
            break;
         }
      }
/*   check for user-defined scalar */
      if (scalar_index < 0) {
         user_check_scalar(scalar_name, &ioscalar);
         if (ioscalar != OK) {
            printf("warning:  scalar=%s undefined\n", scalar_name);
         }
      }
   }
/* copy over space dimensions */
   svp->nd = sp->nd;
   for (id=0; id < svp->nd; id++) {
      ccopy(sp->coords[id], 0, SZNAM, svp->coords[id]);
      ccopy(sp->clabel[id], 0, SZNAM, svp->clabel[id]);
      ccopy(sp->units[id], 0, SZNAM, svp->units[id]);
   }
   svp->ndc = svp->nd;
#ifdef CURV
   svp->curv = sp->curv;
#endif
/* frozen coordinates */
   svp->ndo = sp->ndo;
   for (id=0; id < svp->ndo; id++) {
      jd = svp->nd + id;
      ccopy(sp->coords[jd], 0, SZNAM, svp->coords[jd]);
      ccopy(sp->clabel[jd], 0, SZNAM, svp->clabel[jd]);
      ccopy(sp->units[jd], 0, SZNAM, svp->units[jd]);
      svp->xo[id] = sp->xo[id];
#ifdef CURV
      svp->kad[id] = id;
#endif
   }
/* vector */
   svp->novector = nov;
   if (!svp->novector) {
      svp->ndv = svp->nd + svp->ndo;
      for (id=0; id < svp->ndv; id++) {
/*      find vector index for this coordinate */
         indexv = -1;
         for (i=0; i < MAXVECTOR; i++) {
            if (vector_used[i]) {
               if (csame(svp->coords[id],vector_coord[i])) {
                  indexv = i;
                  break;
               }
            }
         }
         svp->vector_index[id] = id;
         if (indexv < 0) {
            ccopy(" ", 0, SZNAM, vector_name[id]);
            vector_index[id] = -1;
         } else {
            ccopy(vector_scalar[indexv], 0, SZNAM, vector_name[id]);
/*         search for explicit name */
            vector_index[id] = -1;
            for (i=0; i < sp->nq; i++) {
               if (csame(vector_name[id],sp->vlabels[i])) {
                  vector_index[id] = i;
                  break;
               }
            }
/*         check for user-defined scalar */
            if (vector_index[id] < 0) {
               user_check_scalar(vector_name[id], &ioscalar);
               if (ioscalar != OK) {
                  printf("warning:  user vector=%s undefined\n",
                    vector_name[id]);
               }
            }
         }
         ccopy(sp->coords[id], 0, SZNAM, svp->vlabels[id]);
      }
   }
/* find time index (assured) */
   if (*ios) {
   } else if (!mask_time) {
      printf("error:  scatter should either be pick or mask_time\n");
      *ios = ERR;
   } else {
      if (!mask_time_isindex) {
         dtmin = BIG;
         mask_time_index = -1;
         for (i=0; i < sp->nt; i++) {
            fp = &(sp->frame[i]);
            dt = fp->t - mask_time_value;
            dt = dt*dt;;
            if (dt < dtmin) {
               mask_time_index = i;
               dtmin = dt;
            }
         }
      } else if (mask_time_index < 0) {
         mask_time_index = 0;
      } else if (mask_time_index >= sp->nt) {
         mask_time_index = sp->nt - 1;
      }
      fp = &(sp->frame[mask_time_index]);
      svp->ndo = svp->ndo + 1;
      ido = svp->ndo - 1;
      svp->xo[ido] = fp->t;
   }
/* number of pointss */
   svp->npt = sp->npt;
/* xpt memory */
   if (!*ios) {
      ns = svp->nd*svp->npt;
      svp->xpt = (float*)calloc(ns, sizeof(float));
      *ios = (svp->xpt == NULL);
      if (*ios) {
         printf("error:  not enough memory for vis/xpt\n");
      } else {
         for (id=0; id < svp->nd; id++) {
            svp->xmin[id] = BIG;
            svp->xmax[id] = -BIG;
         }
         k = 0;
         for (j=0; j < sp->npt; j++) {
            for (id=0; id < svp->nd; id++) {
               svp->xpt[k] = sp->xpt[k];
               if (svp->xpt[k] < svp->xmin[id]) {
                  svp->xmin[id] = svp->xpt[k];
               }
               if (svp->xpt[k] > svp->xmax[id]) {
                  svp->xmax[id] = svp->xpt[k];
               }
               k = k + 1;
            }
         }
         for (id=0; id < svp->nd; id++) {
            svp->noxr[id] = svp->xmax[id] <= svp->xmin[id];
         }
      }
   }
/* spt memory */
   if (!*ios) {
      ns = svp->npt;
      svp->spt = (int*)calloc(ns, sizeof(int));
      *ios = (svp->spt == NULL);
      if (*ios) {
         printf("error:  not enough memory for vis/spt\n");
      }
   }
/* qpt memory */
   if (*ios) {
   } else if (!svp->noscalar) {
      ns = svp->npt;
      svp->qpt = (float*)calloc(ns, sizeof(float));
      *ios = (svp->qpt == NULL);
      if (*ios) {
         printf("error:  not enough memory for vis/qpt\n");
      }
   }
/* vpt memory */
   if (*ios) {
   } else if (!svp->novector) {
      nv = svp->ndv*svp->npt;
      svp->vpt = (float*)calloc(nv, sizeof(float));
      *ios = (svp->vpt == NULL);
      if (*ios) {
         printf("error:  not enough memory for vis/vpt\n");
      }
   }
/* abstract switch,scalar and/or vector from original */
   if (!*ios) {
/*   single time frame */
      fp = &(sp->frame[mask_time_index]);
      for (j=0; j < sp->npt; j++) {
         svp->spt[j] = fp->spt[j];
         if (svp->spt[j] == FLUID) {
            if (!svp->noscalar) {
               if (scalar_index < 0) {
                  k = j*sp->nd;
                  ki = j*sp->nq;
                  user_get_scalar(&(fp->qpt[ki]),
                    scalar_name, &q);
                  svp->qpt[j] = q;
               } else {
                  k = sp->nq*j + scalar_index;
                  svp->qpt[j] = fp->qpt[k];
               }
               if (nolog_scalar) {
               } else if (svp->qpt[j] > 0.0) {
                  svp->qpt[j] = log10(svp->qpt[j]);
               } else {
                  svp->spt[j] = 0;
               }
            }
            if (!svp->novector) {
               for (jd=0; jd < svp->ndv; jd++) {
                  jv = j*svp->ndv + jd;
                  if (vector_index[jd] < 0) {
                     k = j*sp->nd;
                     ki = j*sp->nq;
                     user_get_scalar(&(fp->qpt[ki]),
                        vector_name[jd], &q);
                     svp->vpt[jv] = q;
                  } else {
                     kv = j*sp->nq + vector_index[jd];
                     svp->vpt[jv] = fp->qpt[kv];
                  }
               }
            }
         }
      }
/*   find scalar range */
      if (!svp->noscalar) {
         svp->qmin = BIG;
         svp->qmax = -BIG;
         for (j=0; j < svp->npt; j++) {
            if (svp->spt[j] == FLUID) {
               if (svp->qpt[j] < svp->qmin) {
                  svp->qmin = svp->qpt[j];
               }
               if (svp->qpt[j] > svp->qmax) {
                  svp->qmax = svp->qpt[j];
               }
            }
         }
      }
/*   find vector range? */
      if (!svp->novector) {
         for (id=0; id < svp->ndv; id++) {
            svp->vmax[id] = 0.0;
         }
         for (j=0; j < svp->npt; j++) {
            if (svp->spt[j] == FLUID) {
               for (id=0; id < svp->ndv; id++) {
                  jv = j*svp->ndv + id;
                  v = svp->vpt[jv];
                  if (v < 0.0) {
                     v = -v;
                  }
                  if (v > svp->vmax[id]){
                     svp->vmax[id] = v;
                  }
               }
            }
         }
         for (id=0; id < svp->ndv; id++) {
            svp->novmax[id] = svp->vmax[id] == 0.0;
         }
      }
   }
/* revoke requested scalar/vector */
   nog = TRUE;
   nos = TRUE;
   nov = TRUE;
/* submit to vis */
   if (!*ios) {
      vis_scatter(svp, ios);
   }
   return 0;
}

int model_selectvg(gp, ios)
   grid_model_struct *gp;
   int *ios;
{
   ccopy("grid", 0, SZNAM, modelv_type);
   modelv_gp = gp;
   nomodelv = FALSE;
   return 0;
}

int model_selectvs(sp, ios)
   scatter_model_struct *sp;
   int *ios;
{
   ccopy("scatter", 0, SZNAM, modelv_type);
   modelv_sp = sp;
   nomodelv = FALSE;
   return 0;
}

int model_deselect() {
   nomodelv = TRUE;
   return 0;
}
