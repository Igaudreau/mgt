/* model.c:  model module 

     This module stores, manipulates and generates plot data for
     one or more model datasets.  Its external routines all impose
     the Error Status scheme.

     It prepares subservient modules such as "ds","vis", etc. .
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"
#include "model.h"
#define VERSION 2.2
#define VERSION_DATE "1-aug-96"
#define MAXHANDLE 4

typedef struct handle_type {
   char name[SZNAM];
   bool used;
   bool open;
   int unit;
   bool writen;
   bool empty;
   char data_type[SZNAM];
   union {
      grid_model_struct grid_model;
      scatter_model_struct scatter_model;
   } data_union;
} handle_struct;

static handle_struct handle[MAXHANDLE];
static char model_name[SZNAM];
static int model_index;
static bool nomodel=TRUE;
static bool dump=FALSE;

#include "vis.h"
static grid_vis_struct grid_vis;
static scatter_vis_struct scatter_vis;


/* EXPLICIT DECLARATIONS ADDED BY BM */
bool model_geom_same(char geomi[], int ndi, char geom[], int nd, char coords[MAXD][SZNAM]);
void model_prepare();
#undef find_bug
#ifdef find_bug
  void model_find_bug_();
#endif
int model_flush();
int model_print();
int model_debug();
int model_debug_grid(grid_model_struct *gp);
int model_debug_scatter(scatter_model_struct *sp);
int model_create(char name[], char file[], int *ios);
int model_close(char name[], int *ios);
int model_remove(char name[], int *ios);
int model_dump(char name[], char file[], int *ios);
int model_load(char name[], char file[], int *ios);
int model_merge(char name1[], char name2[], char name3[], int *ios);
void model_put_ugrid(char name[], char title[], char geom[], float t, 
  int  nd, char coords[], char units[], 
  int nc[][MAXD], float xe[], float xa[], char physics[], int nq,
  char vlabels[], char vunits[],  int vmlt[][3],
  int topn, int sa[], float qa[], int topj, int sb[], float qb[],
  int *ios);
#ifdef PLOTG_SLICE
  int model_put_grid(char name[], char title[], float t, char geom[], int nd, char coords[], char units[], int limlev,
    int nc[][MAXD], float xe[], float xa[], char physics[], int nq, char vlabels[], char vunits[], int vmlt[][3], int nlev,
    int topn, int subn[], int sa[], float qa[], int topj, int subj[], int sb[], float qb[], 
    int ndo, float xo[], int kad[], int *ios);
#else
  int model_put_grid(char name[], char title[], float t, char geom[], int nd, char coords[], char units[], int limlev,
    int nc[][MAXD], float xe[], float xa[], char physics[], int nq, char vlabels[], char vunits[], int vmlt[][3], int nlev,
    int topn, int subn[], int sa[], float qa[], int topj, int subj[], int sb[], float qb[], int *ios);
#endif
int model_put_scatter(char name[], char title[], float t, char geom[], int nd, char coords[], char units[],
  char physics[], int nq, char vlabel[], char vunits[], int vmlt[][3], int npt, float xpt[], int spt[], float qpt[], int *ios);
int model_vis(int *ios);
int model_select(char name[], int *ios);
int model_find(char name[], int *imodel);
int model_new(int *imodel);
#ifdef PLOTG_SLICE
  int model_put_grid_header(grid_model_struct *gp, char title[], char geom[], int nd, char coords[], char units[],
    int limlev, int nc[][MAXD], float xe[], float xa[], char physics[], int nq, char vlabels[], char vunits[], int vmlt[][3], 
    int ndo,float xo[],int kad[], int *ios);
#else
  int model_put_grid_header(grid_model_struct *gp, char title[], char geom[], int nd, char coords[], char units[],
    int limlev, int nc[][MAXD], float xe[], float xa[], char physics[], int nq, char vlabels[], char vunits[], int vmlt[][3], int *ios);
#endif
int model_zap_grid_header(grid_model_struct *gp, int *ios);
int model_check_grid_header(grid_model_struct *gp, char title[], char geom[], int nd, char coords[], char units[],
  int limlev, int nc[][MAXD], float xe[], float xa[], char physics[], int nq, char vlabels[], char vunits[], int vmlt[][3], int *ios);
int model_put_grid_frame(grid_model_struct *gp, int it, float t, int nlev, int topn, int subn[], int sa[], float qa[],
  int topj, int subj[], int sb[], float qb[], int *ios);
int model_zap_grid_frame(grid_model_struct *gp, int it, int *ios);
int model_dump_grid(int unit, grid_model_struct *gp, int *ios);
int model_dump_grid_header(int unit, grid_model_struct *gp,int *ios);
int model_dump_grid_frame(int unit, grid_model_struct *gp, int it, int *ios);
int model_load_grid(int unit, grid_model_struct *gp, int *ios);
int model_load_grid_header(int unit, grid_model_struct *gp, int *ios);
int model_load_grid_frame(int unit, grid_model_struct *gp, int it, int *ios);
int model_put_scatter_header(scatter_model_struct *sp, char title[], char geom[], int nd, char coords[], char units[],
  char physics[], int nq, char vlabels[], char vunits[], int vmlt[][3], int npt, float xpt[], int *ios);
int model_zap_scatter_header(scatter_model_struct *sp, int *ios);
int model_check_scatter_header(scatter_model_struct *sp, char title[], char geom[], int nd, char coords[], char units[],
  char physics[], int nq, char vlabels[], char vunits[], int vmlt[][3], int npt, float xpt[], int *ios);
int model_put_scatter_frame(scatter_model_struct *sp, int it, float t, int npt, int spt[], float qpt[], int *ios);
int model_zap_scatter_frame(scatter_model_struct *sp, int it, int *ios);
int model_dump_scatter(int unit, scatter_model_struct *sp, int *ios);
int model_dump_scatter_header(int unit, scatter_model_struct *sp, int *ios);
int model_dump_scatter_frame(int unit, scatter_model_struct *sp, int it, int *ios);
int model_load_scatter(int unit, scatter_model_struct *sp, int *ios);
int model_load_scatter_header(int unit, scatter_model_struct *sp, int *ios);
int model_load_scatter_frame(int unit, scatter_model_struct *sp, int it, int *ios);
#ifdef CURV
  int model_setup_geom(char geomi[], int ndi, char coordis[MAXD][SZNAM], char geom[], int *nd, char coords[MAXD][SZNAM], int *ndo, float xo[], 
    int kad[MAXD], bool *curv, int *ios);
#else
  model_setup_geom(char geomi[], int ndi, char geom[], int *nd, char coords[MAXD][SZNAM], int *ndo, float xo[], int *ios);
#endif
#ifdef PLOTG_SLICE
int model_masked_geom(char geomi[], int ndi, char coordis[MAXD][SZNAM], char geom[], int *nd, char coords[MAXD][SZNAM], int *ndo, float xo[], 
  int kad[MAXD], bool *curv, int *ios);
#endif


/* EXTERNAL ROUTINES - ADDED BY BM */
extern void sys_date(int str_len, char str[]);
extern void ccopy(char a[], int offset, int n, char b[]);
extern bool csame(char a[], char b[]);
extern void cscan(int lena, char a[], int offset, int lenb, char b[], int *newoff);
extern bool geo_met(char geom[], int nd);
extern int ds_prepare();
extern int ds_flush();
extern int ds_path(char path[], int *ios);
extern int ds_method(char met[], int *ios);
extern int ds_create(char file[], int *unit, int *ios);
extern int ds_open(char file[], int *unit, int *ios);
extern int ds_close(int unit, int *ios);
extern int ds_index(char name[], int no, int lencomb, char comb[]);
extern int ds_put_begin(int unit, char name[], int *ios);
extern int ds_put_end(int unit, char name[], int *ios);
extern int ds_get_begin(int unit, char name[], int *ios);
extern int ds_get_end(int unit, char name[], int *ios);
extern int ds_put0i(int unit, char name[], int value, int *ios);
extern int ds_put1i(int unit, char name[], int nvalue, int value[], int *ios);
extern int ds_put0r(int unit, char name[], float value, int *ios);
extern int ds_put1r(int unit, char name[], int nvalue, float value[], int *ios);
extern int ds_put0b(int unit, char name[], bool value, int *ios);
extern int ds_put1b(int unit, char name[], int nvalue, bool value[], int *ios);
extern int ds_put0s(int unit, char name[], char value[], int *ios);
extern int ds_put0c(int unit, char name[], char value, int *ios);
extern int ds_put1c(int unit, char name[], int nvalue, char value[], int *ios);
extern bool ds_get(int unit, char name[], char type[], int *ios);
extern int ds_get0i(int unit, char name[], int *value, int *ios);
extern int ds_get1i(int unit, char name[], int nvalue, int value[], int *ios);
extern int ds_get0r(int unit, char name[], float *value, int *ios);
extern int ds_get1r(int unit, char name[], int nvalue, float value[], int *ios);
extern int ds_get0b(int unit, char name[], bool *value, int *ios);
extern int ds_get1b(int unit, char name[], int nvalue, bool value[], int *ios);
extern int ds_get0s(int unit, char name[], char value[], int *ios);
extern int ds_get0c(int unit, char name[], char *value, int *ios);
extern int ds_get1c(int unit, char name[], int nvalue, char value[], int *ios);
extern int merror(char text[]);
extern int ds_filex(char path[], char name[], char ext[], int lenfile, char file[]);
extern int ds_check_unit(int unit, int *ios);
extern int ds_put_header(int unit, char name[], char type[], int nvalue, int *ios);
extern int ds_get_prelude(int unit, char name[], char type[], int nvalue, int *ios);
extern int ds_get_header(int unit, int *ios);
extern int ds_check_header(char name[], char type[], int nvalue, int *ios);
extern int ds_indent(int unit, int indent, int *ios);
extern int ds_read_error(int unit, char proc[], int *ios);
extern int ds_write_error(int unit, char proc[], int *ios);
extern int ds_put1c_binary(int unit, int nvalue, char value[], int *ios);
extern int ds_get1c_binary(int unit, int *nvalue, char value[], int *ios);
extern int ds_put1i_binary(int unit, int nvalue, int value[], int *ios);
extern int ds_get1i_binary(int unit, int nvalue, int value[], int *ios);
extern int ds_put1r_binary(int unit, int nvalue, float value[], int *ios);
extern int ds_get1r_binary(int unit, int nvalue, float value[], int *ios);
extern int ds_put1b_binary(int unit, int nvalue, bool value[], int *ios);
extern int ds_get1b_binary(int unit, int nvalue, bool value[], int *ios);
extern int ds_put0s_binary(int unit, char value[], int *ios);
extern int ds_get0s_binary(int unit, int lenvalue, char value[], int *ios);
extern int csg_prepare();
extern int csg_flush();
extern int csg_show();
extern int csg_debug();
extern int model_reset();
extern int model_show();
extern void model_set0s(char name[], char value[], int *ios);
extern int model_set0i(char name[], int value, int *ios);
extern void model_set0r(char name[], float value, int *ios);
extern int model_set0l(char name[], bool value, int *ios);
extern int model_seta(char name[], int *ios);
extern int model_vis_either(scatter_model_struct *sp, scatter_vis_struct *svp, grid_vis_struct *gvp, int *ios);
extern int model_vis_crunch(scatter_model_struct *sp, grid_vis_struct *gvp, int *ios);
extern int model_vis_grid(grid_model_struct *gp, grid_vis_struct *gvp, int *ios);
extern int model_vis_scatter(scatter_model_struct *sp, scatter_vis_struct *svp, int *ios);
extern int model_selectvg(grid_model_struct *gp, int *ios);
extern int model_selectvs(scatter_model_struct *sp, int *ios);
extern int model_deselect();
extern void vis_prepare();
extern void vis_flush();
extern void vis_begin();
extern void vis_end();






bool model_geom_same(geomi, ndi, geom, nd, coords)
   char geomi[],geom[],coords[MAXD][SZNAM];
   int ndi,nd;
{
   char expect[MAXD][SZNAM];
   bool same;
   int id;

   same = TRUE;
   if (ndi != nd) {
      printf("error:  number of dimensions changed\n");
      same = FALSE;
   } else if (csame("xyz",geomi)) {
      if (!csame("cartesian",geom)) {
         printf("error:  geometry changed\n");
         same = FALSE;
      }
      ccopy("x", 0, SZNAM, expect[0]);
      ccopy("y", 0, SZNAM, expect[1]);
      ccopy("z", 0, SZNAM, expect[2]);
   } else if (csame("rzp",geomi)) {
      if (!csame("cylidrical",geom)) {
         printf("error:  geometry changed\n");
         same = FALSE;
      }
      ccopy("r", 0, SZNAM, expect[0]);
      ccopy("z", 0, SZNAM, expect[1]);
      ccopy("phi", 0, SZNAM, expect[2]);
   } else if (csame("rpz",geomi)) {
      if (!csame("cylindrical",geom)) {
         printf("error:  geometry changed\n");
         same = FALSE;
      }
      ccopy("r", 0, SZNAM, expect[0]);
      ccopy("phi", 0, SZNAM, expect[1]);
      ccopy("z", 0, SZNAM, expect[2]);
   } else if (csame("rtp",geomi)) {
      if (!csame("spherical",geom)) {
         printf("error:  geometry changed\n");
         same = FALSE;
      }
      ccopy("r", 0, SZNAM, expect[0]);
      ccopy("theta", 0, SZNAM, expect[1]);
      ccopy("phi", 0, SZNAM, expect[2]);
#ifdef CURV
   } else if (geo_met(geomi, ndi)) {
      for (id=0; id < MAXD; id++) {
         ccopy(coords[id], 0, SZNAM, expect[id]);
      }
#endif
   } else {
      printf("error:  unknown model geometry '%s'\n", geomi);
      same = FALSE;
   }
   if (same) {
      for (id=0; id < MAXD; id++) {
         if (!csame(expect[id],coords[id])) {
            printf("error:  coordinates changed\n");
            printf("..id=%d,expect=%s,got=%s\n", id, expect[id],
              coords[id]);
            same = FALSE;
         }
      }
   }
   return(same);
}

void model_prepare() {
   int i;

   printf("\n");
   printf(" Module MUVI version %4.2f (%s).\n", VERSION, VERSION_DATE);
   printf(" Copyright (C) 1990-97  Mantis Numerics Ltd.\n");
   printf(" 46 The Calls, Leeds LS2 7EY, UK.\n");
   printf("\n");
   ds_prepare();
   csg_prepare();
   vis_prepare();
   for (i= 0; i < MAXHANDLE; i++) {
      handle[i].used = FALSE;
   }
   nomodel = TRUE;
   model_reset();
   return;
}

#undef find_bug
#ifdef find_bug
void model_find_bug_() {
   printf("model_find_bug:  units=(%d,%d,%d)\n", handle[0].unit,
     handle[1].unit, handle[2].unit);
   printf("model_find_bug:  address=%d\n", (int)&handle[0].unit);
   printf("model_find_bug:  address=%d\n", (int)&handle[1].unit);
   printf("model_find_bug:  address=%d\n", (int)&handle[2].unit);
}
#endif
int model_flush() {
   int i,ios;

   ios = OK;
   for (i=0; i < MAXHANDLE; i++) {
      if (handle[i].used) {
         printf("model_flush:  removing i=%d\n", i);
         if (handle[i].open) {
            printf("model_flush:  closing i=%d\n",i);
            printf("model_flush:  unit=%d\n", handle[i].unit);
            model_close(handle[i].name,&ios);
            if (ios) {
               printf("error:  closing model '%s' (model_flush)\n",
                 handle[i].name);
               ios = OK;
            }
         }
         model_remove(handle[i].name,&ios);
         if (ios) {
            printf("error:  removing model '%s' (model_flush)\n",
              handle[i].name);
            ios = OK;
         }
      }
   }
   vis_flush();
   csg_flush();
   ds_flush();
   return 0;
}

int model_print() {
   grid_model_struct *gp;
   scatter_model_struct *sp;
   int i,count;

   printf("Model Database:-\n");
   count = 0;
   for (i=0; i < MAXHANDLE; i++) {
      if (handle[i].used) {
         count = count + 1;
         if (count == 1) {
            printf("%12s %7s %7s %6s %s\n","name", "status", "type",
              "frames", "title");
         }
         printf("%12s", handle[i].name);
         if (handle[i].open) {
            printf(" %7s","file");
         } else {
            printf(" %7s","memory");
         }
         printf(" %7s", handle[i].data_type);
         if (csame("grid",handle[i].data_type)) {
            gp = &(handle[i].data_union.grid_model);
            printf(" %6d",gp->nt);
            printf(" %s\n", gp->title);
         } else if (csame("scatter",handle[i].data_type)) {
            sp = &(handle[i].data_union.scatter_model);
            printf(" %6d",sp->nt);
            printf(" %s\n", sp->title);
         }
      }
   }
   if (nomodel) {
      printf("  No model selected for plotting\n");
   } else {
      printf("  Model select for plotting is '%s'\n", model_name);
   }
   return 0;
}

int model_debug()
{
   handle_struct *hp;
   grid_model_struct *gp;
   scatter_model_struct *sp;

   if (nomodel) {
      printf("No model selected\n");
   } else {
      hp = &(handle[model_index]);
      if (csame("grid",hp->data_type)) {
         gp = &(hp->data_union.grid_model);
         model_debug_grid(gp);
      } else if (csame("scatter",hp->data_type)) {
         sp = &(hp->data_union.scatter_model);
         model_debug_scatter(sp);
      }
   }
   return 0;
}

int model_debug_grid(gp)
   grid_model_struct *gp;
{
   grid_frame_struct *fp;
   int id,jd,i,lev,it,n;

   printf("Grid Model:-\n");
   printf("  version=%4.2f,date=%s\n", gp->version, gp->date);
   printf("  title='%s'\n", gp->title);
   printf("  geom=%s,nd=%d,ndo=%d\n", gp->geom, gp->nd, gp->ndo);
   if (gp->nhist > 0) {
      printf("  history,records=%d\n", gp->nhist);
      for (i=0; i < gp->nhist; i++) {
         printf("  %d: %s\n", i, gp->hist[i]);
      }
   }
   printf("  coordinates:-\n");
   for (id=0; id < (gp->nd+gp->ndo); id++) {
      printf("    id=%d,coords=%s,clabel='%s',units='%s'",
        id, gp->coords[id], gp->clabel[id], gp->units[id]);
      if (id >= gp->nd) {
         jd = id - gp->nd;
         printf(",xo=%E", gp->xo[jd]);
      }
      printf("\n");
   }
   printf("  physics=%s,nq=%d\n", gp->physics, gp->nq);
   for (i=0; i < gp->nq; i++) {
      printf("    i=%d,vlabel=%s,vunits=%s,M=%d,L=%d,T=%d\n",
        i, gp->vlabels[i], gp->vunits[i], gp->vmlt[i][0],
        gp->vmlt[i][1], gp->vmlt[i][2]);
   }
   printf("  limlev=%d,grids:-\n", gp->limlev);
   for (lev=0; lev < (gp->limlev+1); lev++) {
      printf("    lev=%d\n", lev);
      for (id=0; id < gp->nd; id++) {
         n = gp->nc[lev][id];
         printf("      id=%d,nc=%d\n", id, gp->nc[lev][id]);
         printf("        xe_min=%E\n", gp->xe[lev][id][0]);
         printf("        xa_min=%E\n", gp->xa[lev][id][0]);
         printf("        xa_max=%E\n", gp->xa[lev][id][n-1]);
         printf("        xe_max=%E\n", gp->xe[lev][id][n]);
      }
   }
   printf("  nt=%d,frames:-\n", gp->nt);
   for (it=0; it < gp->nt; it++) {
      fp = &(gp->frame[it]);
      printf("    it=%d,t=%E\n", it, fp->t);
   }
   return 0;
}

int model_debug_scatter(sp)
   scatter_model_struct *sp;
{
   scatter_frame_struct *fp;
   int id,jd,i,it,k;

   printf("Scatter Model:-\n");
   printf("  version=%4.2f,date=%s\n", sp->version, sp->date);
   printf("  title='%s'\n", sp->title);
   if (sp->nhist > 0) {
      printf("  history,records=%d\n", sp->nhist);
      for (i=0; i < sp->nhist; i++) {
         printf("  %d: %s\n", i, sp->hist[i]);
      }
   }
   printf("  geom=%s,nd=%d,ndo=%d\n", sp->geom, sp->nd, sp->ndo);
   printf("  coordinates:-\n");
   for (id=0; id < (sp->nd+sp->ndo); id++) {
      printf("    id=%d,coords=%s,clabel='%s',units='%s'",
        id, sp->coords[id], sp->clabel[id], sp->units[id]);
      if (id >= sp->nd) {
         jd = id - sp->nd;
         printf(",xo=%E", sp->xo[jd]);
      }
      printf("\n");
   }
   printf("  physics=%s,nq=%d\n", sp->physics, sp->nd);
   for (i=0; i < sp->nq; i++) {
      printf("    i=%d,vlabel=%s,vunits=%s,M=%d,L=%d,T=%d\n",
        i, sp->vlabels[i], sp->vunits[i], sp->vmlt[i][0],
        sp->vmlt[i][1], sp->vmlt[i][2]);
   }
   printf("  npt=%d,xpt:-\n", sp->npt);
   for (i=0; i < sp->npt; i++) {
      printf("    xpt=");
      for (id=0; id < sp->nd; id++) {
         k = sp->nd*i + id;
         printf(" %E", sp->xpt[k]);
      }
      printf("\n");
      if ((i >= 9) && (i < (sp->npt-1))) {
         printf("    etc...\n");
         break;
      }
   }
   printf("  nt=%d,frames:-\n", sp->nt);
   for (it=0; it < sp->nt; it++) {
      fp = &(sp->frame[it]);
      printf("    it=%d,t=%E\n", it, fp->t);
   }
   return 0;
}


int model_create(name, file, ios)
   char name[],file[];
   int *ios;
{
   int imodel,unit;
   handle_struct *hp;

   if (!*ios) {
/*
      printf("Creating model name '%s' on file '%s'\n", name, file);
*/
      model_find(name, &imodel);
      if (imodel >= 0) {
         printf("error:  model '%s' exists (model_create)\n", name);
         *ios = ERR;
      } else {
         model_new(&imodel);
         if (imodel < 0) {
            printf("error:  no slot for model '%s' (model_create)\n",
              name);
            *ios = ERR;
         } else {
            ds_create(file, &unit, ios);
            if (*ios) {
               printf("error:  creating file for model '%s' ",name);
               printf("(model_create\n");
            } else {
               hp = &(handle[imodel]);
               ccopy(name, 0, SZNAM, hp->name);
               hp->used = TRUE;
               hp->empty = TRUE;
               hp->open = TRUE;
               hp->unit = unit;
               hp->writen = TRUE;
            }
         }
      }
   }
   return 0;
}

int model_close(name, ios)
   char name[];
   int *ios;
{
   int imodel,input_ios;
   handle_struct *hp;

   input_ios = *ios;
   *ios = OK;
/*
   printf("Closing model name '%s'\n", name);
*/
   model_find(name, &imodel);
   if (imodel < 0) {
      printf("error:  model '%s' not found (model_close)\n",name);
      *ios = ERR;
   } else {
      hp = &(handle[imodel]);
      if (hp->open) {
         if (hp->writen) {
            ds_put_end(hp->unit, hp->data_type, ios);
            if (*ios) {
               printf("error:  writing trailer for model '%s' ",name);
               printf("(model_close)\n");
            }
         }
         ds_close(hp->unit, ios);
         hp->open = FALSE;
      } else {
         printf("error:  model '%s' not open (model_close)\n", name);
         *ios = ERR;
      }
      if (*ios) {
         printf("error:  tourble closing model '%s' (model_close)\n",
           name);
         *ios = OK;
      }
      model_remove(name, ios);
   }
   if (input_ios) {
      *ios = input_ios;
   }
   return 0;
}

int model_remove(char name[], int *ios)
{
   grid_model_struct *gp;
   scatter_model_struct *sp;
   int imodel;
   handle_struct *hp;

   if (!*ios) {
/*
      printf("Removing model name '%s'\n", name);
*/
      model_find(name, &imodel);
      if (imodel < 0) {
         printf("error:  model '%s' not found (model_remove)\n", name);
         *ios = ERR;
      } else {
         hp = &(handle[imodel]);
         if (hp->open) {
            printf("error:  model '%s' not open (model_remove)\n",name);
            *ios = ERR;
         } else {
            if (csame("grid",handle[imodel].data_type)) {
               gp = &(handle[imodel].data_union.grid_model);
               model_zap_grid_frame(gp, 0, ios);
               model_zap_grid_header(gp, ios);
            } else if (csame("scatter",handle[imodel].data_type)) {
               sp = &(handle[imodel].data_union.scatter_model);
               model_zap_scatter_frame(sp, 0, ios);
               model_zap_scatter_header(sp, ios);
            }
            hp->used = FALSE;
            if (nomodel) {
            } else if (csame(name,model_name)) {
               nomodel = TRUE;
               model_deselect();
            }
         }
      }
   }
   return 0;
}

int model_dump(name, file, ios)
   char name[],file[];
   int *ios;
{
   int imodel,unit;
   handle_struct *hp;
   grid_model_struct *gp;
   scatter_model_struct *sp;

   if (!*ios) {
/*
      printf("Dumping model name '%s' to file '%s'\n", name, file);
*/
      model_find(name, &imodel);
      if (imodel < 0) {
         printf("error:  model '%s' not found (model_dump)\n", name);
         *ios = ERR;
      } else {
         ds_create(file, &unit, ios);
         if (*ios) {
            printf("error:  creating file '%s' for model '%s' ",
              file, name);
            printf("(model_create)\n");
         } else {
            hp = &(handle[imodel]);
            ds_put_begin(unit, hp->data_type, ios);
            if (*ios) {
               printf("error:  dumping data type for model '%s' ",name);
               printf("(model_dump)\n");
            } else {
               if (csame("grid",hp->data_type)) {
                  gp = &(hp->data_union.grid_model);
                  model_dump_grid(unit, gp, ios);
                  if (*ios) {
                     printf("error:  dumping header for model ");
                     printf("'%s' (model_dump)\n", name);
                  }
               } else if (csame("scatter", hp->data_type)) {
                  sp = &(hp->data_union.scatter_model);
                  model_dump_scatter(unit, sp, ios);
                  if (*ios) {
                     printf("error:  dumping header for model ");
                     printf("'%s' (model_dump)\n", name);
                  }
               }
               ds_put_end(unit, hp->data_type, ios);
            }
            ds_close(unit, ios);
         }
     }
   }
   return 0;
}

int model_load(name, file, ios)
   char file[],name[];
   int *ios;
{
   int imodel,unit;
   handle_struct *hp;
   grid_model_struct *gp;
   scatter_model_struct *sp;
   char nameb[SZNAM],typeb[SZNAM];

   if (!*ios) {
/*
      printf("Loading model name '%s' from file '%s'\n", name, file);
*/
      model_find(name, &imodel);
      if (imodel >= 0) {
         printf("error:  model '%s' exists (model_load)\n", name);
         *ios = ERR;
      } else {
         model_new(&imodel);
         if (imodel < 0) {
            printf("error:  no slot for model '%s' (model_load)\n",
              name);
            *ios = ERR;
         } else {
            ds_open(file, &unit, ios);
            if (*ios) {
               printf("error:  opening file '%s' for model '%s' ",
                 file, name);
               printf("(model_load)\n");
            } else {
               hp = &(handle[imodel]);
               if (!ds_get(unit, nameb, typeb, ios)) {
                  printf("error:  loading data type for model '%s' ",
                    name);
                  printf("(model_load)\n");
               } else {
                  ccopy(name, 0, SZNAM, hp->name);
                  hp->open = FALSE;
                  hp->empty = TRUE;
                  ccopy(typeb, 0, SZNAM, hp->data_type);
                  if (csame("grid", hp->data_type)) {
                     gp = &(hp->data_union.grid_model);
                     model_load_grid(unit, gp, ios);
                     if (*ios) {
                        printf("error:  loading header for model ");
                        printf("'%s' (model_load)\n", name);
                     }
                  } else if (csame("scatter", hp->data_type)) {
                     sp = &(hp->data_union.scatter_model);
                     model_load_scatter(unit, sp, ios);
                     if (*ios) {
                        printf("error:  loading header for model ");
                        printf("'%s' (model_load)\n", name);
                     }
                  } else {
                     printf("error:  data type '%s' unknown ",
                       hp->data_type);
                     printf("for model '%s' (model_load)\n", name);
                     *ios = ERR;
                  }
                  if (*ios) {
                     printf("error:  loading model '%s' ", name);
                     printf("(model_load)\n");
                     hp->used = FALSE;
                  } else {
                     hp->used = TRUE;
                     model_select(name, ios);
                  }
               }
               ds_close(unit, ios);
            }
         }
      }
   }
   return 0;
}

int model_merge(name1, name2, name3, ios)
   char name1[],name2[],name3[];
   int *ios;
{
   printf("error:  function not implemented (model_merge)\n");
   return 0;
}
/*
void model_put_ugrid(name, title, geom, t, nd, coords, units, 
  nc, xe, xa,  physics, nq, vlabels, vunits,  vmlt,
  topn, sa, qa, topj, sb, qb, ios)
   char name[],title[],geom[],coords[],units[],physics[],
     vlabels[],vunits[];
   int nd,nc[][MAXD],nq,vmlt[][3],topn,topj,*ios,sa[],sb[];
   float t,xe[],xa[],qa[],qb[];
*/
void model_put_ugrid(char name[], char title[], char geom[], float t, 
  int  nd, char coords[], char units[], 
  int nc[][MAXD], float xe[], float xa[], char physics[], int nq,
  char vlabels[], char vunits[],  int vmlt[][3],
  int topn, int sa[], float qa[], int topj, int sb[], float qb[],
  int *ios)
{
   int subn[1],subj[1],limlev,nlev;
#ifdef PLOTG_SLICE
   int ndo=0,kad[MAXD],id;
   float xo[MAXD];
#endif
   /*
   printf("address of geom %p", (void *)&geom[0]);printf("\n");
   printf("address of name %p", (void *)&name[0]);printf("\n");
   printf("b address of title %p", (void *)&title[0]);printf("\n");
   printf("e address of title %p", (void *)&title[SZLIN]);printf("\n");
   printf("address of t %p", (void *)&t);printf("\n");
   printf("name ");printf(name);printf("\n");
   printf("title ");printf(title);printf("\n");
   printf("t %f", t);printf("\n");
   printf("geom ");printf(geom);printf("\n");
   printf("nd %d", nd);;printf("\n");
   printf("coords ");printf(coords);printf("\n");
   printf("units ");printf(units);printf("\n");
   printf("physics ");printf(physics);printf("\n");
   for (id = 0; id < nd; id++){   
     printf("nc %d ", nc[0][id]);
   }
   printf("\n");
   for (id = 0; id <= topn; id++){   
     printf("xe %f ", xe[id]);
   }
   printf("\n");
   for (id = 0; id < topn; id++){   
     printf("xa %f ", xa[id]);
   }
   printf("\n");
   printf("nq %d", nq);;printf("\n");
   printf("vlabels ");printf(vlabels);printf("\n");
   printf("vunits ");printf(vunits);printf("\n");
   printf("vunits ");printf(vunits);
   for (id = 0; id < 3; id++){
     printf("vmlt %d ", vmlt[0][id]);printf("\n");
   }
   printf("topn %d", topn);;printf("\n");
   for (id = 0; id < 11; id++){
     printf("sa %d ", sa[id]);
   }
   printf("\n");
   for (id = 0; id < 11; id++){
     printf("qa %f ", qa[id]);
   }
   printf("\n");
   printf("topj %d", topj);;printf("\n");
   for (id = 0; id < 11; id++){
     printf("sb %d ", sb[id]);
   }
   printf("\n");
   for (id = 0; id < 11; id++){
     printf("qb %f ", qb[id]);
   }
   printf("\n");
   */
   limlev = 0;
   nlev = 0;
#ifdef PLOTG_SLICE
   for (id=0; id < nd; id++) {
      kad[id] = id;
   }
   model_put_grid(name, title, t, geom, nd, coords, units, limlev,
     nc, xe, xa, physics, nq, vlabels, vunits, vmlt, nlev,
     topn, subn, sa, qa, topj, subj, sb, qb, 
     ndo, xo, kad, ios);
#else
   model_put_grid(name, title, t, geom, nd, coords, units, limlev,
     nc, xe, xa, physics, nq, vlabels, vunits, vmlt, nlev,
     topn, subn, sa, qa, topj, subj, sb, qb, ios);
#endif
}

#ifdef PLOTG_SLICE
int model_put_grid(name, title, t, geom, nd, coords, units, limlev,
  nc, xe, xa,  physics, nq, vlabels, vunits, vmlt, nlev,
  topn, subn, sa, qa, topj, subj, sb, qb, 
  ndo, xo, kad, ios)
#else
model_put_grid(name, title, t, geom, nd, coords, units, limlev,
  nc, xe, xa,  physics, nq, vlabels, vunits, vmlt, nlev,
  topn, subn, sa, qa, topj, subj, sb, qb, ios)
#endif
   char name[],title[],geom[],coords[],units[],physics[],
     vlabels[],vunits[];
   int nd,limlev,nc[][MAXD],nq,vmlt[][3],nlev,topn,subn[],topj,subj[],
     *ios,sa[],sb[];
   float t,xe[],xa[],qa[],qb[];
#ifdef PLOTG_SLICE
   int ndo,kad[];
   float xo[];
#endif
{
   int imodel;

   handle_struct *hp;
   grid_model_struct *gp;

   if (!*ios) {
/*
      printf("Put grid model name '%s'\n", name);
*/
      model_find(name, &imodel);
      if (imodel < 0) {
         model_new(&imodel);
         if (imodel < 0) {
            printf("error:  no slot for model '%s' (model_put_grid)\n",
              name);
            *ios = ERR;
         } else {
            hp = &(handle[imodel]);
            ccopy(name, 0, SZNAM, hp->name);
            hp->used = TRUE;
            hp->empty = TRUE;
            hp->open = FALSE;
         }
      } else {
         hp = &(handle[imodel]);
      }
      if (!*ios) {
         gp = &(hp->data_union.grid_model);
         if (hp->empty) {
            ccopy("grid", 0, SZNAM, hp->data_type);
#ifdef PLOTG_SLICE
            model_put_grid_header(gp, title, geom, nd, coords, units,
              limlev, nc, xe, xa, physics, nq, vlabels, vunits, vmlt,
              ndo, xo, kad, ios);
#else
            model_put_grid_header(gp, title, geom, nd, coords, units,
              limlev, nc, xe, xa, physics, nq, vlabels, vunits, vmlt,
              ios);
#endif
            if (*ios) {
               printf("error:  storing header for model '%s' ", name);
               printf("(model_put_grid)\n");
            } else if (hp->open) {
               ds_put_begin(hp->unit, hp->data_type, ios);
               if (*ios) {
                  printf("error:  dumping data type for model '%s' ",
                    name);
                  printf("(model_put_grid)\n");
               } else {
                  model_dump_grid_header(hp->unit, gp, ios);
                  if (*ios) {
                     printf("error:  dumping header for model '%s' ",
                       name);
                     printf("(model_put_grid)\n");
                  }
               }
            }
            hp->empty = FALSE;
         } else {
            model_check_grid_header(gp, title, geom, nd, coords, units,
              limlev, nc, xe, xa, physics, nq, vlabels, vunits, vmlt,
              ios);
            if (*ios) {
               printf("error:  inconsistent header for model '%s' ",
                 name);
               printf("(model_put_grid)\n");
            }
         }
         if (*ios) {
         } else if (gp->nt >= MAXT) {
            printf("error:  too many time frames for model '%s' ",
              name);
            printf("(model_put_grid)\n");
            *ios = ERR;
         } else {
            gp->nt = gp->nt + 1;
            model_put_grid_frame(gp, gp->nt-1, t, nlev, topn, subn, sa,
              qa, topj, subj, sb, qb, ios);
            if (*ios) {
               printf("error:  storing frame %d for model '%s' ", 
                 gp->nt-1, name);
               printf("(model_put_grid)\n");
            } else if (hp->open) {
               model_dump_grid_frame(hp->unit, gp, gp->nt-1, ios);
               if (*ios) {
                  printf("error:  dumping frame %d for model '%s' ",
                    gp->nt-1, name);
                  printf("(model_put_grid)\n");
               }
               model_zap_grid_frame(gp, gp->nt-1, ios);
            } else {
               model_select(name, ios);
            }
         }
      }
   }
   return 0;
}

int model_put_scatter(name, title, t, geom, nd, coords, units,
  physics, nq, vlabel, vunits, vmlt, npt, xpt, spt, qpt, ios)
   char name[],title[],geom[],coords[],units[],
     physics[],vlabel[],vunits[];
   float t,xpt[],qpt[];
   int nd,nq,vmlt[][3],npt,*ios,spt[];
{
   int imodel;
   handle_struct *hp;
   scatter_model_struct *sp;

   if (!*ios) {
/*
      printf("Put scatter model '%s'\n", name);
*/
      model_find(name, &imodel);
      if (imodel < 0) {
         model_new(&imodel);
         if (imodel < 0) {
            printf("error:  no slot for model '%s' ", name);
            printf("(model_put_scatter)\n");
            *ios = ERR;
         } else {
            hp = &(handle[imodel]);
            ccopy(name, 0, SZNAM, hp->name);
            hp->used = TRUE;
            hp->empty = TRUE;
            hp->open = FALSE;
         }
      } else {
         hp = &(handle[imodel]);
      }
      if (!*ios) {
         sp = &(hp->data_union.scatter_model);
         if (hp->empty) {
            ccopy("scatter", 0, SZNAM, hp->data_type);
            model_put_scatter_header(sp, title, geom, nd, coords, units,
              physics, nq, vlabel, vunits, vmlt, npt, xpt, ios);
            if (*ios) {
               printf("error:  storing header for model '%s' ",name);
               printf("(model_put_scatter)\n");
            } else if (hp->open) {
               ds_put_begin(hp->unit, hp->data_type, ios);
               if (*ios) {
                  printf("error:  dumping data type for model '%s' ",
                    name);
                  printf("(model_put_scatter)\n");
               } else {
                  model_dump_scatter_header(hp->unit, sp, ios);
                  if (*ios) {
                     printf("error:  dumping header for model '%s' ",
                       name);
                     printf("(model_put_scatter)\n");
                  }
               }
            }
            hp->empty = FALSE;
         } else {
            model_check_scatter_header(sp, title, geom, nd, coords,
              units, physics, nq, vlabel, vunits, vmlt, npt, xpt, ios);
            if (*ios) {
               printf("error:  inconsistent header for model '%s' ",
                 name);
               printf("(model_put_scatter)\n");
            }
         }
         if (*ios) {
         } else if (sp->nt >= MAXT) {
            printf("error:  too may time frames for model '%s' ",
              name);
            printf("(model_put_scatter)\n");
            *ios = ERR;
         } else {
            sp->nt = sp->nt + 1;
            model_put_scatter_frame(sp, sp->nt-1, t, npt, spt, qpt,
              ios);
            if (*ios) {
               printf("error:  storing frame %d for model '%s' ",
                 sp->nt-1, name);
               printf("(model_put_scatter)\n");
            } else if (hp->open) {
               model_dump_scatter_frame(hp->unit, sp, sp->nt-1, ios);
               if (*ios) {
                  printf("error:  dumping frame %d for model '%s' ",
                    sp->nt-1, name);
                  printf("(model_put_scatter)\n");
               }
               model_zap_scatter_frame(sp, sp->nt-1, ios);
            } else {
               model_select(name, ios);
            }
         }
      }
   }
   return 0;
}

/* PLOTTING ROUTINES */

int model_vis(int *ios)
{
   handle_struct *hp;
   grid_model_struct *gp;
   scatter_model_struct *sp;
   grid_vis_struct *gvp;
   scatter_vis_struct *svp;

   if (*ios) {
   } else {
      vis_begin();
      if (!nomodel) {
         hp = &(handle[model_index]);
         if (csame("grid",hp->data_type)) {
            gp = &(hp->data_union.grid_model);
            gvp = &(grid_vis);
            model_vis_grid(gp, gvp, ios);
         } else if (csame("scatter",hp->data_type)) {
            sp = &(hp->data_union.scatter_model);
            svp = &(scatter_vis);
            gvp = &(grid_vis);
            model_vis_either(sp, svp, gvp, ios);
         }
      }
      vis_end();
   }
   return 0;
}

int model_select(name, ios)
   char name[];
   int *ios;
{
   handle_struct *hp;
   grid_model_struct *gp;
   scatter_model_struct *sp;
   int imodel;

   if (!*ios) {
      model_find(name, &imodel);
      if (imodel < 0) {
         printf("error:  model '%s' not found (model_select)\n", name);
         *ios = ERR;
      } else {
         nomodel = FALSE;
         ccopy(name, 0, SZNAM, model_name);
         model_index = imodel;
         hp = &(handle[model_index]);
         if (csame("grid",hp->data_type)) {
            gp = &(hp->data_union.grid_model);
            model_selectvg(gp, ios);
         } else if (csame("scatter", hp->data_type)) {
            sp = &(hp->data_union.scatter_model);
            model_selectvs(sp, ios);
         } else {
             model_deselect();
         }
         model_reset();
      }
   }
   return 0;
}

/* INTERNAL ROUTINES */

int model_find(name, imodel)
   char name[];
   int *imodel;
{
   int i;

   *imodel = (-1);
   for (i=0; i < MAXHANDLE; i++) {
      if (*imodel < 0) {
         if (handle[i].used) {
            if (csame(name,handle[i].name)) {
               *imodel = i;
            }
         }
      }
   }
   return 0;
}

int model_new(imodel)
   int *imodel;
{
   int i;

   *imodel = (-1);
   for (i=0; i < MAXHANDLE; i++) {
      if (*imodel < 0) {
         if (!handle[i].used) {
            *imodel = i;
         }
      }
   }
   return 0;
}

/* INTERNAL GRID */

#ifdef PLOTG_SLICE
int model_put_grid_header(gp, title, geom, nd, coords, units,
  limlev, nc, xe, xa, physics, nq, vlabels, vunits, vmlt, 
  ndo,xo,kad, ios)
#else
model_put_grid_header(gp, title, geom, nd, coords, units,
  limlev, nc, xe, xa, physics, nq, vlabels, vunits, vmlt, ios)
#endif
   grid_model_struct *gp;
   char title[],geom[],coords[],units[],physics[],vlabels[],vunits[];
   int nd,limlev,nc[][MAXD],nq,vmlt[][3],*ios;
   float xe[],xa[];
#ifdef PLOTG_SLICE
   int ndo,kad[];
   float xo[];
#endif
{
   int offset,id,i,j,ke,ka,lev,k;
#ifdef CURV
   char coordis[MAXD][SZNAM];
   int ido;
#endif

   if (dump) {
      printf("model_put_grid_header:  starts\n");
   }
/* version */
   gp->version = VERSION;
/* date now */
   sys_date(SZDATE, gp->date);
/* put title */
   ccopy(title, 0, SZLIN, gp->title);
   if (dump) {
      printf("..title=%s\n", gp->title);
   }
/* no history yet */
   gp->nhist = 0;
/* put geometry */
/* check geometry, setup degreneracy */
#ifdef CURV
   offset = 0;
   for (id=0; id < nd; id++) {
      if (id >= 0 && id < MAXD) {			// for some reason the compiler wants this check BM
        cscan(SZCO,coords,offset,SZNAM,coordis[id],&offset);
      }
#ifdef PLOTG_SLICE
      gp->kad[id] = kad[id];
#endif
/*
      fprintf(stderr, "model_put_grid_header:  id=%d,co=%s\n",
        id, coordis[id]);
*/
   }
#ifdef PLOTG_SLICE
   gp->ndo = ndo;
   offset=0;
   for (ido=0; ido < ndo; ido++) {
      gp->xo[ido] = xo[ido];
      cscan(SZCO,coords,offset,SZNAM,coordis[nd+ido],&offset);
      gp->kad[nd+ido] = kad[nd+ido];
   }
   model_masked_geom(geom, nd, coordis, gp->geom, &gp->nd, gp->coords,
     &gp->ndo, gp->xo, gp->kad, &gp->curv, ios);
#else
   model_setup_geom(geom, nd, coordis, gp->geom, &gp->nd, gp->coords,
     &gp->ndo, gp->xo, gp->kad, &gp->curv, ios);
#endif
/*
   printf("model_put_grid_header:-\n");
   for (id=0; id < gp->nd; id++) {
      printf("id=%d,coords=%s,kad=%d\n", id, gp->coords[id],
       gp->kad[id]);
   }
   for (ido=0; ido < gp->ndo; ido++) {
      id = gp->nd + ido;
      printf("id=%d,coords=%s,kad=%d,xo=%f\n", id, gp->coords[id],
         gp->kad[id], gp->xo[ido]);
   }
*/
#else
   model_setup_geom(geom, nd, gp->geom, &gp->nd, gp->coords, &gp->ndo,
     gp->xo, ios);
#endif
   if (*ios) {
      printf("error:  checking model geometry\n");
   } else if ((nd < 0) || (nd > MAXD)) {
      merror("grid_header:  illegal number of dimensions");
      *ios = ERR;
   } else {
      offset = 0;
      for (id=0; id < nd; id++) {
         cscan(SZCO,coords,offset,SZNAM,gp->clabel[id],&offset);
         if (dump) {
            printf("..clabel[%d]=%s\n", id, gp->clabel[id]);
         }
      }
      offset = 0;
      for (id=0; id < nd; id++) {
         cscan(SZCO,units,offset,SZNAM,gp->units[id],&offset);
         if (dump) {
            printf("..units[%d]=%s\n", id, gp->units[id]);
         }
      }
      for (id=nd; id < MAXD; id++) {
         ccopy(gp->coords[id], 0, SZNAM, gp->clabel[id]);
         ccopy(gp->units[0], 0, SZNAM, gp->units[id]);
      }
   }
/* put physics */
   ccopy(physics, 0, SZNAM, gp->physics);
   gp->nq = nq;
   if (dump) {
      printf("..physics=%s,nq=%d\n", gp->physics, gp->nq);
   }
   if ((nq < 1) || (nq > MAXQ)) {
      merror("grid_header:  illegal number of variables");
      *ios = ERR;
   } else {
      offset = 0;
      for (i=0; i < nq; i++) {
         cscan(SZCQ,vlabels,offset,SZNAM,gp->vlabels[i],&offset);
         if (dump) {
            printf("..vlabels[%d]=%s\n", i, gp->vlabels[i]);
         }
      }
      offset = 0;
      for (i=0; i < nq; i++) {
         cscan(SZCQ,vunits,offset,SZNAM,gp->vunits[i],&offset);
         if (dump) {
            printf("..vunits[%d]=%s\n", i, gp->vunits[i]);
         }
      }
      for (i=0; i < nq; i++) {
         for (k=0; k < 3; k++) {
            gp->vmlt[i][k] = vmlt[i][k];
         }
      }
   }
/* grid levels and coordinates */
   gp->limlev = limlev;
   if (*ios) {
   } else if ((limlev < 0) || (limlev > MAXLEV-1)) {
      merror("grid_header:  illegal number of refined grids");
      *ios = ERR;
   } else {
      ke = 0;
      ka = 0;
      for (lev=0; lev <=limlev; lev++) {
         if (*ios) {
            break;
         } else {
            for (id = 0; id < nd; id++) {
               gp->nc[lev][id] = nc[lev][id];
               gp->xe[lev][id] = (float*)calloc((nc[lev][id]+1),
                 sizeof(float));
               *ios = (gp->xe[lev][id] == NULL);
               if (*ios) {
                  merror("grid_header:  not enough memory for xe");
                  break;
               } else {
                  for (j=0; j < (nc[lev][id]+1); j++) {
                     (gp->xe[lev][id])[j] = xe[ke];
                     ke = ke + 1;
                  }
               }
               gp->xa[lev][id] = (float*)calloc(nc[lev][id],
                 sizeof(float));
               *ios = (gp->xa[lev][id] == NULL);
               if (*ios) {
                  merror("grid_header:  not enough memory for xa");
                  break;
               } else {
/*
                  printf("..lev=%d,id=%d\n", lev, id);
                  printf("..nc=%d\n", nc[lev][id]);
*/
                  for (j=0; j < nc[lev][id]; j++) {
/*
                     printf("..j=%d,ka=%d,xa=%g\n", j, ka, xa[ka]);
*/
                     (gp->xa[lev][id])[j] = xa[ka];
                     ka = ka + 1;
                  }
               }
            }
         }
      }
   }
/* reset number of time frames */
   gp->nt = 0;
   return 0;
}

int model_zap_grid_header(gp, ios)
   grid_model_struct *gp;
   int *ios;
{
   int id,lev;

   for (lev=0; lev <= gp->limlev; lev++) {
      for (id=0; id < gp->nd; id++) {
         if (gp->xe[lev][id] != NULL) {
            cfree(gp->xe[lev][id], (gp->nc[lev][id]+1), sizeof(float));
            gp->xe[lev][id] = NULL;
         }
         if (gp->xa[lev][id] != NULL) {
            cfree(gp->xa[lev][id], gp->nc[lev][id], sizeof(float));
            gp->xa[lev][id] = NULL;
         }
      }
   }
   return 0;
}

int model_check_grid_header(gp, title, geom, nd, coords, units,
  limlev, nc, xe, xa, physics, nq, vlabels, vunits, vmlt, ios)
   grid_model_struct *gp;
   char title[],geom[],coords[],units[],physics[],vlabels[],vunits[];
   int nd,limlev,nc[][MAXD],nq,vmlt[][3],*ios;
   float xe[],xa[];
{
   char temp[SZNAM];
   int offset,i,j,lev,id,ke,ka;

/* check title, no big deal */
   if (!csame(title,gp->title)) {
      merror("grid_header:  title changed");
   }
/* check geometry (serious) */
   if (!model_geom_same(geom,nd,gp->geom,gp->nd, gp->coords)) {
      merror("grid_header:  geometry changed (fatal)");
      *ios = ERR;
   } else if (nd != gp->nd) {
      merror("grid_header:  number of dimensions changed (fatal)");
      *ios = ERR;
   } else {
      offset = 0;
      for (id=0; id < nd; id++) {
         cscan(SZCO,coords,offset,SZNAM,temp,&offset);
         if (!csame(gp->clabel[id],temp)) {
            merror("grid_header:  coords changed (fatal)");
            *ios = ERR;
            break;
         }
      }
      offset = 0;
      for (id=0; id < nd; id++) {
         cscan(SZCO,units,offset,SZNAM,temp,&offset);
         if (!csame(gp->units[id],temp)) {
            merror("grid_header:  units changed (fatal)");
            *ios = ERR;
            break;
         }
      }
   }
/* check physics (serious) */
   if (!csame(physics,gp->physics)) {
      merror("grid_header:  physics changed (fatal)");
      *ios= ERR;
   } else if (nq != gp->nq) {
      merror("grid_header:  number of variables changed (fatal)");
      *ios = ERR;
   } else {
      offset = 0;
      for (i=0; i < nq; i++) {
         cscan(SZCQ,vlabels,offset,SZNAM,temp,&offset);
         if (!csame(gp->vlabels[i],temp)) {
            merror("grid_header:  vlabels changed (fatal)");
            *ios = ERR;
            break;
         }
      }
      offset = 0;
      for (i=0; i < nq; i++) {
         cscan(SZCQ,vunits,offset,SZNAM,temp,&offset);
         if (!csame(gp->vunits[i],temp)) {
            merror("grid_header:  vunits changed (fatal)");
            *ios = ERR;
            break;
         }
      }
      for (i=0; i < nq; i++) {
         if (vmlt[i][0] != gp->vmlt[i][0]) {
            merror("grid_header:  vmlt changed (fatal)");
            *ios = ERR;
            break;
         } else if (vmlt[i][1] != gp->vmlt[i][1]) {
            merror("grid_header:  vmlt changed (fatal)");
            *ios = ERR;
            break;
         } else if (vmlt[i][2] != gp->vmlt[i][2]) {
            merror("grid_header:  vmlt changed (fatal)");
            *ios = ERR;
            break;
         }
      }
   }
/* grid levels and coordinates */
   if (*ios) {
   } else if (limlev != gp->limlev) {
      merror("grid_header:  limiting number of refined levels changed");
      *ios = ERR;
   } else {
      ke = 0;
      ka = 0;
      for (lev=0; lev <=limlev; lev++) {
         if (*ios) {
            break;
         } else {
            for (id = 0; id < nd; id++) {
               if (nc[lev][id] != (gp->nc[lev][id])) {
                  merror("grid_header:  grid size changed");
                  *ios = ERR;
                  break;
               } else {
                  for (j=0; j < (nc[lev][id]+1); j++) {
                     if ((gp->xe[lev][id])[j] != xe[ke]) {
                        merror("grid_header:  xe changed");
                        *ios = ERR;
                        break;
                     }
                     ke = ke + 1;
                  }
               }
               if (*ios) {
               } else {
                  for (j=0; j < nc[lev][id]; j++) {
                     if ((gp->xa[lev][id])[j] != xa[ka]) {
                        merror("grid_header:  xa changed");
                        *ios = ERR;
                        break;
                     }
                     ka = ka + 1;
                  }
               }
            }
         }
      }
   }
   return 0;
}

int model_put_grid_frame(gp, it, t, nlev, topn, subn, sa, qa,
  topj, subj, sb, qb, ios)
   grid_model_struct *gp;
   int it,nlev,topn,subn[],topj,subj[],*ios,sa[],sb[];
   float t,qa[],qb[];
{
   grid_frame_struct *fp;
   int i,j,k;

/* grab memory pointer before we begin */
   fp = &(gp->frame[it]);
/* actual number of grid levels */
   fp->nlev = nlev;
/* numbers of nodes and joins in frame */
   fp->topn = topn;
   fp->topj = topj;
/* memory (sa) */
   if (!*ios) {
      fp->sa = (int*)calloc(topn, sizeof(int));
      *ios = fp->sa == NULL;
      if (*ios) {
         merror("grid_frame:  not enough memory for sa");
      } else {
         for (j=0; j < topn; j++) {
            (fp->sa)[j] = sa[j];
         }
      }
   }
/* memory (qa) */
   if (!*ios) {
      fp->qa = (float*)calloc(topn*gp->nq, sizeof(float));
      *ios = (fp->qa == NULL);
      if (*ios) {
         merror("grid_frame:  not enough memory for qa");
      } else {
         k = 0;
         for (j=0; j < topn; j++) {
            for (i=0; i < gp->nq; i++) {
               (fp->qa)[k] = qa[k];
               if (dump) {
                  printf("model_put... j=%d,i=%d,k=%d,qa=%f\n",
                    j, i, k, (fp->qa)[k]);
               }
               k = k + 1;
            }
         }
      }
   }
/* memory (subn,only if nlev > 0) */
   if (*ios) {
   } else if (gp->limlev > 0) {
      fp->subn = (int*)calloc(topn, sizeof(int));
      *ios = fp->subn == NULL;
      if (*ios) {
         merror("grid_frame:  not enough memory for subn");
      } else {
         for (j=0; j < topn; j++) {
            (fp->subn)[j] = subn[j];
         }
      }
   } else {
      fp->subn = NULL;
   }
/* memory (sb) */
   if (*ios) {
   } else if (topj < 1) {
      fp->sb = NULL;
   } else {
      fp->sb = (int*)calloc(topj, sizeof(int));
      *ios = (fp->sb == NULL);
      if (*ios) {
         merror("grid_frame:  not enough memory for sb");
      } else {
         for (j=0; j < topj; j++) {
            (fp->sb)[j] = sb[j];
         }
      }
   }
/* memory (qb) */
   if (*ios) {
   } else if (topj < 1) {
     fp->qb = NULL;
   } else {
      fp->qb = (float*)calloc(topj*gp->nq, sizeof(float));
      *ios = (fp->qb == NULL);
      if (*ios) {
         merror("grid_frame:  not enough memory for qb");
      } else {
         k = 0;
         for (j=0; j < topj; j++) {
            for (i=0; i < gp->nq; i++) {
               (fp->qb)[k] = qb[k];
               k = k + 1;
            }
         }
      }
   }
/* memory (subj,only if nlev > 0) */
   if (*ios) {
   } else if (topj < 1) {
   } else if (gp->limlev > 0) {
      fp->subj = (int*)calloc(topj, sizeof(int));
      *ios = (fp->subj == NULL);
      if (*ios) {
         merror("grid_frame:  not enough memory for subj");
      } else {
         for (j=0; j < topj; j++) {
            (fp->subj)[j] = subj[j];
         }
      }
   } else {
      fp->subj = NULL;
   }
/* frame time */
   fp->t = t;
   return 0;
}

int model_zap_grid_frame(gp, it, ios)
   grid_model_struct *gp;
   int it,*ios;
{
   grid_frame_struct *fp;
   int i;

   if (!*ios) {
      for (i = gp->nt-1; i >= it; i--) {
         fp = &(gp->frame[i]);
         if (fp->sa != NULL) {
            cfree(fp->sa, fp->topn, sizeof(int));
            fp->sa = NULL;
         }
         if (fp->qa != NULL) {
            cfree(fp->qa, fp->topn*gp->nq, sizeof(float));
            fp->qa = NULL;
         }
         if (fp->sb != NULL) {
            cfree(fp->sb, fp->topj, sizeof(int));
            fp->sb = NULL;
         }
         if (fp->qb != NULL) {
            cfree(fp->qb, fp->topj*gp->nq, sizeof(float));
            fp->qb = NULL;
         }
         if (fp->subn != NULL) {
            cfree(fp->subn, fp->topn, sizeof(int));
            fp->subn = NULL;
         }
         if (fp->subj != NULL) {
            cfree(fp->subj, fp->topj, sizeof(int));
            fp->subj = NULL;
         }
      }
      gp->nt = it;
   }
   return 0;
}

/* INTERNAL GRID DATA DUMP */

int model_dump_grid(unit, gp, ios)
   int unit,*ios;
   grid_model_struct *gp;
{
   int it;

   model_dump_grid_header(unit, gp, ios);
   for (it=0; it < gp->nt; it++) {
      model_dump_grid_frame(unit, gp, it, ios);
   }
   return 0;
}

int model_dump_grid_header(unit, gp,ios)
   int unit,*ios;
   grid_model_struct *gp;
{
   char temp[SZNAM];
   int id,i,lev;

/* common ground */
   ds_put0r(unit, "version", gp->version, ios);
   ds_put0s(unit, "date", gp->date, ios);
   ds_put0s(unit, "title", gp->title, ios);
   ds_put0i(unit, "nhist", gp->nhist, ios);
   for (i=0; i < gp->nhist; i++) {
      ds_put0s(unit, "HISTORY", gp->hist[i], ios);
      if (*ios) {
         printf("error:  writing model history records\n");
         break;
      }
   }
   ds_put0s(unit, "geom", gp->geom, ios);
   ds_put0i(unit, "nd", gp->nd, ios);
   ds_put0i(unit, "ndo", gp->ndo, ios);
   for (id=0; id < (gp->nd+gp->ndo); id++) {
      ds_index("coords", id, SZNAM, temp);
      ds_put0s(unit, temp, gp->coords[id], ios);
   }
   for (id=0; id < (gp->nd+gp->ndo); id++) {
      ds_index("clabel", id, SZNAM, temp);
      ds_put0s(unit, temp, gp->clabel[id], ios);
   }
   for (id=0; id < (gp->nd+gp->ndo); id++) {
      ds_index("units", id, SZNAM, temp);
      ds_put0s(unit, temp, gp->units[id], ios);
   }
   ds_put1r(unit, "xo", gp->ndo, gp->xo, ios);
   ds_put0s(unit, "physics", gp->physics, ios);
   ds_put0i(unit, "nq", gp->nq, ios);
   for (i=0; i < gp->nq; i++) {
      ds_index("vlabels", i, SZNAM, temp);
      ds_put0s(unit, temp, gp->vlabels[i], ios);
   }
   for (i=0; i < gp->nq; i++) {
      ds_index("vunits", i, SZNAM, temp);
      ds_put0s(unit, temp, gp->vunits[i], ios);
   }
   for (i=0; i < gp->nq; i++) {
      ds_index("vmlt", i, SZNAM, temp);
      ds_put1i(unit, temp, 3, gp->vmlt[i], ios);
   }
/* grid specifics */
   ds_put0i(unit, "limlev", gp->limlev, ios);
   for (lev=0; lev < (gp->limlev+1); lev++) {
      for (id=0; id < gp->nd; id++) {
         ds_put0i(unit, "nc", gp->nc[lev][id], ios);
         ds_put1r(unit, "xe", gp->nc[lev][id]+1, gp->xe[lev][id], ios);
         ds_put1r(unit, "xa", gp->nc[lev][id], gp->xa[lev][id], ios);
      }
   }
   return 0;
}

int model_dump_grid_frame(unit, gp, it, ios)
   int unit,it,*ios;
   grid_model_struct *gp;
{
   grid_frame_struct *fp;

   fp = &(gp->frame[it]);
   ds_put_begin(unit, "grid_frame", ios);
   ds_put0r(unit, "t", fp->t, ios);
   ds_put0i(unit, "nlev", fp->nlev, ios);
   ds_put0i(unit, "topn", fp->topn, ios);
   ds_put1i(unit, "sa", fp->topn, fp->sa, ios);
   ds_put1r(unit, "qa", fp->topn*gp->nq, fp->qa, ios);
   if (gp->limlev > 0) {
      ds_put1i(unit, "subn", fp->topn, fp->subn, ios);
   }
   ds_put0i(unit, "topj", fp->topj, ios);
   ds_put1i(unit, "sb", fp->topj, fp->sb, ios);
   ds_put1r(unit, "qb", fp->topj*gp->nq, fp->qb, ios);
   if (gp->limlev > 0) {
      ds_put1i(unit, "subj", fp->topj, fp->subj, ios);
   }
   ds_put_end(unit, "grid_frame", ios);
   return 0;
}

/* INTERNAL GRID DATA LOAD */

int model_load_grid(unit, gp, ios)
   int unit,*ios;
   grid_model_struct *gp;
{
   char nameb[SZNAM],typeb[SZNAM];

   model_load_grid_header(unit, gp, ios);
   gp->nt = 0;
   while (ds_get(unit, nameb, typeb, ios)) {
      if (!csame("begin", nameb)) {
         merror("model_load_grid:  expecting grid_frame");
         *ios = ERR;
         break;
      } else if (!csame("grid_frame", typeb)) {
         merror("model_load_grid:  expecting grid_frame");
         *ios = ERR;
      } else {
         gp->nt = gp->nt + 1;
         model_load_grid_frame(unit, gp, gp->nt-1, ios);
         if (*ios) {
            printf("model_load_grid:  load_grid_frame failed!\n");
         }
      }
   }
   return 0;
}

int model_load_grid_header(unit, gp,ios)
   int unit,*ios;
   grid_model_struct *gp;
{
   char temp[SZNAM];
   int id,i,lev,j;

/* common ground */
   ds_get0r(unit, "version", &(gp->version), ios);
   ds_get0s(unit, "date", gp->date, ios);
   ds_get0s(unit, "title", gp->title, ios);
   ds_get0i(unit, "nhist", &(gp->nhist), ios);
   for (i=0; i < gp->nhist; i++) {
      if (i < MAXHIST) {
         j = i;
      } else {
         j = MAXHIST - 1;
      }
      ds_get0s(unit, "HISTORY", gp->hist[j], ios);
   }
   if (gp->nhist > MAXHIST) {
      printf("note:  some history records lost\n");
      gp->nhist = MAXHIST;
   }
   ds_get0s(unit, "geom", gp->geom, ios);
   ds_get0i(unit, "nd", &(gp->nd), ios);
   ds_get0i(unit, "ndo", &(gp->ndo), ios);
   if (!*ios) {
      for (id=0; id < (gp->nd+gp->ndo); id++) {
         ds_index("coords", id, SZNAM, temp);
         ds_get0s(unit, temp, gp->coords[id], ios);
      }
      for (id=0; id < (gp->nd+gp->ndo); id++) {
         ds_index("clabel", id, SZNAM, temp);
         ds_get0s(unit, temp, gp->clabel[id], ios);
      }
      for (id=0; id < (gp->nd+gp->ndo); id++) {
         ds_index("units", id, SZNAM, temp);
         ds_get0s(unit, temp, gp->units[id], ios);
      }
   }
   ds_get1r(unit, "xo", gp->ndo, gp->xo, ios);
   ds_get0s(unit, "physics", gp->physics, ios);
   ds_get0i(unit, "nq", &(gp->nq), ios);
   if (!*ios) {
      for (i=0; i < gp->nq; i++) {
         ds_index("vlabels", i, SZNAM, temp);
         ds_get0s(unit, temp, gp->vlabels[i], ios);
      }
      for (i=0; i < gp->nq; i++) {
         ds_index("vunits", i, SZNAM, temp);
         ds_get0s(unit, temp, gp->vunits[i], ios);
      }
      for (i=0; i < gp->nq; i++) {
         ds_index("vmlt", i, SZNAM, temp);
         ds_get1i(unit, temp, 3, gp->vmlt[i], ios);
      }
   }
/* grid specifics */
   ds_get0i(unit, "limlev", &(gp->limlev), ios);
   if (*ios) {
   } else if ((gp->limlev < 0) || (gp->limlev >= MAXLEV)) {
      merror("model_load_grid_header:  illegal limlev");
      *ios = ERR;
   } else {
      for (lev=0; lev < (gp->limlev+1); lev++) {
         for (id=0; id < gp->nd; id++) {
            ds_get0i(unit, "nc", &(gp->nc[lev][id]), ios);
            if (*ios) {
               break;
            } else if (gp->nc[lev][id] < 1) {
               merror("model_load_grid_header:  illegal nc");
               *ios = ERR;
            } else {
               gp->xe[lev][id] = (float*)calloc((gp->nc[lev][id]+1),
                 sizeof(float));
               *ios = (gp->xe[lev][id] == NULL);
               ds_get1r(unit, "xe", gp->nc[lev][id]+1, gp->xe[lev][id],
                 ios);
               gp->xa[lev][id] = (float*)calloc(gp->nc[lev][id],
                 sizeof(float));
               *ios = (gp->xa[lev][id] == NULL);
               ds_get1r(unit, "xa", gp->nc[lev][id], gp->xa[lev][id],
                 ios);
             }
         }
      }
   }
   return 0;
}

int model_load_grid_frame(unit, gp, it, ios)
   int unit,it,*ios;
   grid_model_struct *gp;
{
   grid_frame_struct *fp;

/* grab memory pointer before we begin */
   fp = &(gp->frame[it]);
/* frame time */
   ds_get0r(unit, "t", &(fp->t), ios);
/* refined levels */
   ds_get0i(unit, "nlev", &(fp->nlev), ios);
/* numbers of nodes in frame */
   ds_get0i(unit, "topn", &(fp->topn), ios);
/* memory (sa) */
   if (!*ios) {
      fp->sa = (int*)calloc(fp->topn, sizeof(int));
      *ios = fp->sa == NULL;
      if (*ios) {
         merror("grid_frame:  not enough memory for sa");
      } else {
         ds_get1i(unit, "sa", fp->topn, fp->sa, ios);
      }
   }
/* memory (qa) */
   if (!*ios) {
      fp->qa = (float*)calloc(fp->topn*gp->nq, sizeof(float));
      *ios = (fp->qa == NULL);
      if (*ios) {
         merror("grid_frame:  not enough memory for qa");
      } else {
         ds_get1r(unit, "qa", fp->topn*gp->nq, fp->qa, ios);
      }
   }
/* memory (subn,only if limlev > 0) */
   if (*ios) {
   } else if (gp->limlev > 0) {
      fp->subn = (int*)calloc(fp->topn, sizeof(int));
      *ios = fp->subn == NULL;
      if (*ios) {
         merror("grid_frame:  not enough memory for subn");
      } else {
         ds_get1i(unit, "subn", fp->topn, fp->subn, ios);
      }
   } else {
      fp->subn = NULL;
   }
/* number of joins */
   ds_get0i(unit, "topj", &(fp->topj), ios);
/* memory (sb) */
   if (*ios) {
   } else if (fp->topj < 1) {
      fp->sb = NULL;
      ds_get1i(unit, "sb", fp->topj, fp->sb, ios);
   } else {
      fp->sb = (int*)calloc(fp->topj, sizeof(int));
      *ios = (fp->sb == NULL);
      if (*ios) {
         merror("grid_frame:  not enough memory for sb");
      } else {
         ds_get1i(unit, "sb", fp->topj, fp->sb, ios);
      }
   }
/* memory (qb) */
   if (*ios) {
   } else if (fp->topj < 1) {
      fp->qb = NULL;
      //ds_get1r(unit, "qb", fp->topj*gp->nq, fp->sb, ios);
      ds_get1i(unit, "qb", fp->topj * gp->nq, fp->sb, ios);	// Fix by BM
   } else {
      fp->qb = (float*)calloc(fp->topj*gp->nq, sizeof(float));
      *ios = (fp->qb == NULL);
      if (*ios) {
         merror("grid_frame:  not enough memory for qb");
      } else {
         ds_get1r(unit, "qb", fp->topj*gp->nq, fp->qb, ios);
      }
   }
/* memory (subj,olimly if nlev > 0) */
   if (*ios) {
   } else if (fp->topj < 1) {
      fp->subj = NULL;
      if (gp->limlev > 0) {
         ds_get1i(unit, "subj", fp->topj, fp->sb, ios);
      }
   } else if (gp->limlev > 0) {
      fp->subj = (int*)calloc(fp->topj, sizeof(int));
      *ios = (fp->subj == NULL);
      if (*ios) {
         merror("grid_frame:  not enough memory for subj");
      } else {
         ds_get1i(unit, "subj", fp->topj, fp->subj, ios);
      }
   } else {
      fp->subj = NULL;
   }
/* frame end */
   ds_get_end(unit, "grid_frame", ios);
   return 0;
}

/* INTERNAL SCATTER */

int model_put_scatter_header(sp, title, geom, nd, coords, units,
  physics, nq, vlabels, vunits, vmlt, npt, xpt, ios)
   scatter_model_struct *sp;
   char title[],geom[],coords[],units[],physics[],vlabels[],vunits[];
   int nd,vmlt[][3],npt,*ios, nq;  //MIR
   float xpt[];
{
   int offset,id,i,j,k;
#ifdef CURV
   char coordis[MAXD][SZNAM];
#endif

   if (dump) {
      printf("model_put_scatter_header:  starts\n");
   }
/* version */
   sp->version = VERSION;
/* date */
   sys_date(SZDATE, sp->date);
/* put title */
   ccopy(title, 0, SZLIN, sp->title);
   if (dump) {
      printf("..title=%s\n", sp->title);
   }
/* no history yet */
   sp->nhist = 0;
/* check geometry, assign coordinate names */
#ifdef CURV
   offset = 0;
   for (id=0; id < nd; id++) {
      cscan(SZCO,coords,offset,SZNAM,coordis[id],&offset);
   }
   model_setup_geom(geom, nd, coordis, sp->geom, &sp->nd, sp->coords,
     &sp->ndo, sp->xo, sp->kad, &sp->curv, ios);
#else
   model_setup_geom(geom, nd, sp->geom, &sp->nd, sp->coords, &sp->ndo,
     sp->xo, ios);
#endif
   if (*ios) {
   } else if ((nd < 0) || (nd > MAXD)) {
      merror("scatter_header:  illegal number of dimensions");
      *ios = ERR;
   } else {
      offset = 0;
      for (id=0; id < nd; id++) {
         cscan(SZCO,coords,offset,SZNAM,sp->clabel[id],&offset);
         if (dump) {
            printf("..clabel[%d]=%s\n", id, sp->clabel[id]);
         }
      }
      offset = 0;
      for (id=0; id < nd; id++) {
         cscan(SZCO,units,offset,SZNAM,sp->units[id],&offset);
         if (dump) {
            printf("..units[%d]=%s\n", id, sp->units[id]);
         }
      }
      for (id=nd; id < MAXD; id++) {
         ccopy(sp->coords[id], 0, SZNAM, sp->clabel[id]);
         ccopy(sp->units[0], 0, SZNAM, sp->units[id]);
      }
   }
/* put physics */
   ccopy(physics, 0, SZNAM, sp->physics);
   sp->nq = nq;
   if (dump) {
      printf("..physics=%s,nq=%d\n", sp->physics, sp->nq);
   }
   if ((nq < 1) || (nq > MAXQ)) {
      merror("scatter_header:  illegal number of variables");
      *ios = ERR;
   } else {
      offset = 0;
      for (i=0; i < nq; i++) {
         cscan(SZCQ,vlabels,offset,SZNAM,sp->vlabels[i],&offset);
         if (dump) {
            printf("..vlabels[%d]=%s\n", i, sp->vlabels[i]);
         }
      }
      offset = 0;
      for (i=0; i < nq; i++) {
         cscan(SZCQ,vunits,offset,SZNAM,sp->vunits[i],&offset);
         if (dump) {
            printf("..vunits[%d]=%s\n", i, sp->vunits[i]);
         }
      }
      for (i=0; i < nq; i++) {
         for (k=0; k < 3; k++) {
            sp->vmlt[i][k] = vmlt[i][k];
         }
      }
   }
/* point positions */
   if (*ios) {
   } else if (npt < 1) {
      merror("scatter_header:  illegal number of points");
      *ios = ERR;
   } else {
      sp->npt = npt;
      if (dump) {
         printf("..npt=%d\n", sp->npt);
      }
      sp->xpt = (float*)calloc((npt*nd), sizeof(float));
      *ios = (sp->xpt == NULL);
      if (*ios) {
         merror("scatter_header:  not enough memory for xpt");
      } else {
         for (j=0; j < npt; j++) {
            if (dump) {
               printf("..j=%d,xpt=", j);
            }
            for (id=0; id < nd; id++) {
               k = nd*j + id;
               (sp->xpt)[k] = xpt[k];
               if (dump) {
                  printf("%d %f ", k, (sp->xpt)[k]);
               }
            }
            if (dump) {
               printf("\n");
            }
         }
      }
   }
   sp->nt = 0;
   return 0;
}

int model_zap_scatter_header(sp, ios)
   scatter_model_struct *sp;
   int *ios;
{
   if (sp->xpt != NULL) {
      cfree(sp->xpt, (sp->npt*sp->nd), sizeof(float));
      sp->xpt = NULL;
   }
   return 0;
}

int model_check_scatter_header(sp, title, geom, nd, coords, units,
  physics, nq, vlabels, vunits, vmlt, npt, xpt, ios)
   scatter_model_struct *sp;
   char title[],geom[],coords[],units[],physics[],vlabels[],vunits[];
   int nd,nq,vmlt[][3],npt,*ios;
   float xpt[];
{
   char temp[SZNAM];
   int offset,i,j,k,id;

/* check title, no big deal */
   if (!csame(title,sp->title)) {
      merror("model title changed");
   }
/* check geometry (serious) */
   if (!model_geom_same(geom, nd, sp->geom, sp->nd, sp->coords)) {
      merror("model geometry changed (fatal)");
      *ios = ERR;
   } else if (nd != sp->nd) {
      merror("model number of dimensions changed (fatal)");
      *ios = ERR;
   } else {
      offset = 0;
      for (id=0; id < nd; id++) {
         cscan(SZCO,coords,offset,SZNAM,temp,&offset);
         if (!csame(sp->clabel[id],temp)) {
            merror("scatter_grid:  coords changed (fatal)");
            *ios = ERR;
            break;
         }
      }
      offset = 0;
      for (id=0; id < nd; id++) {
         cscan(SZCO,units,offset,SZNAM,temp,&offset);
         if (!csame(sp->units[id],temp)) {
            merror("scatter_grid:  units changed (fatal)");
            *ios = ERR;
            break;
         }
      }
   }
/* check physics (serious) */
   if (!csame(physics,sp->physics)) {
      merror("model physics changed (fatal)");
      *ios= ERR;
   } else if (nq != sp->nq) {
      merror("model number of variables changed (fatal)");
      *ios = ERR;
   } else {
      offset = 0;
      for (i=0; i < nq; i++) {
         cscan(SZCQ,vlabels,offset,SZNAM,temp,&offset);
         if (!csame(sp->vlabels[i],temp)) {
            merror("scatter_grid:  vlabels changed (fatal)");
            *ios = ERR;
            break;
         }
      }
      offset = 0;
      for (i=0; i < nq; i++) {
         cscan(SZCQ,vunits,offset,SZNAM,temp,&offset);
         if (!csame(sp->vunits[i],temp)) {
            merror("scatter_grid:  vunits changed (fatal)");
            *ios = ERR;
            break;
         }
      }
      for (i=0; i < nq; i++) {
         if (vmlt[i][0] != sp->vmlt[i][0]) {
            merror("grid_header:  vmlt changed (fatal)");
            *ios = ERR;
            break;
         } else if (vmlt[i][1] != sp->vmlt[i][1]) {
            merror("grid_header:  vmlt changed (fatal)");
            *ios = ERR;
            break;
         } else if (vmlt[i][2] != sp->vmlt[i][2]) {
            merror("grid_header:  vmlt changed (fatal)");
            *ios = ERR;
            break;
         }
      }
   }
/* check point coordinates (only if OK so far,serious) */
   if (*ios) {
   } else if (npt != sp->npt) {
      merror("scatter_model:  number of points changed (fatal)");
      *ios = ERR;
   } else {
      for (j=0; j < npt; j++) {
         if (*ios) {
            break;
         } else {
            for (id=0; id < nd; id++) {
               k = nd*j+id;
               if (xpt[k] != (sp->xpt)[k]) {
                  merror("scatter_model:  coordinates changed (fatal)");
                  *ios = ERR;
                  break;
               }
            }
         }
      }
   }
   return 0;
}

int model_put_scatter_frame(sp, it, t, npt, spt, qpt, ios)
   scatter_model_struct *sp;
   float t,qpt[];
   int it,npt,*ios,spt[];
{
   scatter_frame_struct *fp;
   int i,j,k;

   fp = &(sp->frame[it]);
   fp->t = t;
   if (dump) {
      printf("..t=%f\n", fp->t);
   }
   fp->spt = (int*)calloc(sp->npt, sizeof(int));
   *ios = fp->spt == NULL;
   if (*ios) {
      merror("not enough frame memory for spt (fatal)");
   } else {
      fp->qpt = (float*)calloc((sp->nq*sp->npt), sizeof(float));
      *ios = fp->qpt == NULL;
      if (*ios) {
         merror("not enough frame memory for qpt (fatal)");
      }
   }
   if (*ios) {
      merror("failed to create scatter frame");
   } else {
      for (j=0; j < npt; j++) {
         (fp->spt)[j] = spt[j];
         if (dump) {
            printf("..j=%d,spt=%c,qpt=",j,(fp->spt)[j]);
         }
         for (i=0; i < sp->nq; i++) {
            k = (sp->nq)*j+i;
            (fp->qpt)[k] = qpt[k];
            if (dump) {
               printf(" %d %f", k, (fp->qpt)[k]);
            }
         }
         if (dump) {
            printf("\n");
         }
      }
   }
   return 0;
}

int model_zap_scatter_frame(sp, it, ios)
   scatter_model_struct *sp;
   int it,*ios;
{
   scatter_frame_struct *fp;
   int i;

   if (!*ios) {
      for (i = sp->nt-1; i >= it; i--) {
         fp = &(sp->frame[i]);
         if (fp->spt != NULL) {
            cfree(fp->spt, sp->npt, sizeof(int));
            fp->spt = NULL;
         }
         if (fp->qpt != NULL) {
            cfree(fp->qpt, sp->nd*sp->npt, sizeof(float));
            fp->qpt = NULL;
         }
      }
      sp->nt = it;
   }
   return 0;
}

/* INTERNAL SCATTER DATA DUMP */

int model_dump_scatter(unit, sp, ios)
   int unit,*ios;
   scatter_model_struct *sp;
{
   int it;

   model_dump_scatter_header(unit, sp, ios);
   for (it=0; it < sp->nt; it++) {
      model_dump_scatter_frame(unit, sp, it, ios);
   }
   return 0;
}

int model_dump_scatter_header(unit, sp,ios)
   int unit,*ios;
   scatter_model_struct *sp;
{
   char temp[SZNAM];
   int id,i;

/* common ground */
   ds_put0r(unit, "version", sp->version, ios);
   ds_put0s(unit, "date", sp->date, ios);
   ds_put0s(unit, "title", sp->title, ios);
   ds_put0i(unit, "nhist", sp->nhist, ios);
   for (i=0; i < sp->nhist; i++) {
      ds_put0s(unit, "HISTORY", sp->hist[i], ios);
      if (*ios) {
         printf("error:  writing model history records\n");
         break;
      }
   }
   ds_put0s(unit, "geom", sp->geom, ios);
   ds_put0i(unit, "nd", sp->nd, ios);
   ds_put0i(unit, "ndo", sp->ndo, ios);
   for (id=0; id < (sp->nd+sp->ndo); id++) {
      ds_index("coords", id, SZNAM, temp);
      ds_put0s(unit, temp, sp->coords[id], ios);
   }
   for (id=0; id < (sp->nd+sp->ndo); id++) {
      ds_index("clabel", id, SZNAM, temp);
      ds_put0s(unit, temp, sp->clabel[id], ios);
   }
   for (id=0; id < (sp->nd+sp->ndo); id++) {
      ds_index("units", id, SZNAM, temp);
      ds_put0s(unit, temp, sp->units[id], ios);
   }
   ds_put1r(unit, "xo", sp->ndo, sp->xo, ios);
   ds_put0s(unit, "physics", sp->physics, ios);
   ds_put0i(unit, "nq", sp->nq, ios);
   for (i=0; i < sp->nq; i++) {
      ds_index("vlabels", i, SZNAM, temp);
      ds_put0s(unit, temp, sp->vlabels[i], ios);
   }
   for (i=0; i < sp->nq; i++) {
      ds_index("vunits", i, SZNAM, temp);
      ds_put0s(unit, temp, sp->vunits[i], ios);
   }
   for (i=0; i < sp->nq; i++) {
      ds_index("vmlt", i, SZNAM, temp);
      ds_put1i(unit, temp, 3, sp->vmlt[i], ios);
   }
/* scatter specifics */
   ds_put0i(unit, "npt", sp->npt, ios);
   ds_put1r(unit, "xpt", (sp->nd)*(sp->npt), sp->xpt, ios);
   return 0;
}

int model_dump_scatter_frame(unit, sp, it, ios)
   int unit,it,*ios;
   scatter_model_struct *sp;
{
   scatter_frame_struct *fp;

   fp = &(sp->frame[it]);
   ds_put_begin(unit, "scatter_frame", ios);
   ds_put0r(unit, "t", fp->t, ios);
   ds_put1i(unit, "spt", sp->npt, fp->spt, ios);
   ds_put1r(unit, "qpt", (sp->npt)*(sp->nq), fp->qpt, ios);
   ds_put_end(unit, "scatter_frame", ios);
   return 0;
}

/* INTERNAL SCATTER DATA LOAD */

int model_load_scatter(unit, sp, ios)
   int unit,*ios;
   scatter_model_struct *sp;
{
   char nameb[SZNAM],typeb[SZNAM];

   model_load_scatter_header(unit, sp, ios);
   sp->nt = 0;
   while (ds_get(unit, nameb, typeb, ios)) {
      if (!csame("begin", nameb)) {
         merror("model_load_scatter:  expecting scatter_frame");
         *ios = ERR;
         break;
      } else if (!csame("scatter_frame", typeb)) {
         merror("model_load_scatter:  expecting scatter_frame");
         *ios = ERR;
      } else {
         sp->nt = sp->nt + 1;
         model_load_scatter_frame(unit, sp, sp->nt-1, ios);
      }
   }
   return 0;
}

int model_load_scatter_header(unit, sp, ios)
   int unit,*ios;
   scatter_model_struct *sp;
{
   char temp[SZNAM];
   int id,i,j;

/* common ground */
   ds_get0r(unit, "version", &(sp->version), ios);
   ds_get0s(unit, "date", sp->date, ios);
   ds_get0s(unit, "title", sp->title, ios);
   ds_get0i(unit, "nhist", &(sp->nhist), ios);
   for (i=0; i < sp->nhist; i++) {
      if (i < MAXHIST) {
         j = i;
      } else {
         j = MAXHIST - 1;
      }
      ds_get0s(unit, "HISTORY", sp->hist[j], ios);
   }
   if (sp->nhist > MAXHIST) {
      printf("note:  some history records lost\n");
      sp->nhist = MAXHIST;
   }
   ds_get0s(unit, "geom", sp->geom, ios);
   ds_get0i(unit, "nd", &(sp->nd), ios);
   ds_get0i(unit, "ndo", &(sp->ndo), ios);
   for (id=0; id < (sp->nd+sp->ndo); id++) {
      ds_index("coords", id, SZNAM, temp);
      ds_get0s(unit, temp, sp->coords[id], ios);
   }
   for (id=0; id < (sp->nd+sp->ndo); id++) {
      ds_index("clabel", id, SZNAM, temp);
      ds_get0s(unit, temp, sp->clabel[id], ios);
   }
   for (id=0; id < (sp->nd+sp->ndo); id++) {
      ds_index("units", id, SZNAM, temp);
      ds_get0s(unit, temp, sp->units[id], ios);
   }
   ds_get1r(unit, "xo", sp->ndo, sp->xo, ios);
   ds_get0s(unit, "physics", sp->physics, ios);
   ds_get0i(unit, "nq", &(sp->nq), ios);
   for (i=0; i < sp->nq; i++) {
      ds_index("vlabels", i, SZNAM, temp);
      ds_get0s(unit, temp, sp->vlabels[i], ios);
   }
   for (i=0; i < sp->nq; i++) {
      ds_index("vunits", i, SZNAM, temp);
      ds_get0s(unit, temp, sp->vunits[i], ios);
   }
   for (i=0; i < sp->nq; i++) {
      ds_index("vmlt", i, SZNAM, temp);
      ds_get1i(unit, temp, 3, sp->vmlt[i], ios);
   }
/* scatter specifics */
   ds_get0i(unit, "npt", &(sp->npt), ios);
   if (*ios) {
   } else if (sp->npt < 1) {
      merror("too few scatter points");
      *ios = ERR;
   } else {
      sp->xpt = (float*)calloc((sp->npt*sp->nd), sizeof(float));
      *ios = (sp->xpt == NULL);
      if (*ios) {
         merror("scatter_header:  not enough memory for xpt");
      } else {
         ds_get1r(unit, "xpt", (sp->nd)*(sp->npt), sp->xpt, ios);
      }
   }
   return 0;
}

int model_load_scatter_frame(unit, sp, it, ios)
   int unit,it,*ios;
   scatter_model_struct *sp;
{
   scatter_frame_struct *fp;

   fp = &(sp->frame[it]);
   ds_get0r(unit, "t", &(fp->t), ios);
   fp->spt = (int*)calloc(sp->npt, sizeof(int));
   *ios = fp->spt == NULL;
   if (*ios) {
      merror("not enough frame memory for spt (fatal)");
   } else {
      fp->qpt = (float*)calloc((sp->nq*sp->npt), sizeof(float));
      *ios = fp->qpt == NULL;
      if (*ios) {
         merror("not enough frame memory for qpt (fatal)");
      }
   }
   if (!*ios) {
      ds_get1i(unit, "spt", sp->npt, fp->spt, ios);
      ds_get1r(unit, "qpt", (sp->npt)*(sp->nq), fp->qpt, ios);
      ds_get_end(unit, "scatter_frame", ios);
   }
   return 0;
}

#ifdef CURV
int model_setup_geom(geomi, ndi, coordis, geom, nd, coords, ndo, xo, 
  kad, curv, ios)
#else
model_setup_geom(geomi, ndi, geom, nd, coords, ndo, xo, ios)
#endif
   char geomi[],geom[],coords[MAXD][SZNAM];
   int ndi,*nd,*ndo,*ios;
#ifdef CURV
   char coordis[MAXD][SZNAM];
   int kad[MAXD];
   bool *curv;
#endif
   float xo[];
{
   int id;

   if (*ios) {
   } else {
#ifdef CURV
      *curv = FALSE;
      *ndo = 0;
#endif
      *nd = ndi;
      if (csame("xyz",geomi)) {
         ccopy("cartesian", 0, SZNAM, geom);
         ccopy("x", 0, SZNAM, coords[0]);
         ccopy("y", 0, SZNAM, coords[1]);
         ccopy("z", 0, SZNAM, coords[2]);
      } else if (csame("rzp",geomi)) {
         ccopy("cylidrical", 0, SZNAM, geom);
         ccopy("r", 0, SZNAM, coords[0]);
         ccopy("z", 0, SZNAM, coords[1]);
         ccopy("phi", 0, SZNAM, coords[2]);
#ifdef CURV
         if (ndi==3) {
            printf("Note:  using surrogate geometry MUVI_RZP\n");
            ccopy("MUVI_RZP", 0, SZNAM, geom);
            if (geo_met(geom, ndi)) {
               *curv = TRUE;
            } else {
               *ios = 1;
            }
         }
#endif
      } else if (csame("rpz",geomi)) {
         ccopy("cylindrical", 0, SZNAM, geom);
         ccopy("r", 0, SZNAM, coords[0]);
         ccopy("phi", 0, SZNAM, coords[1]);
         ccopy("z", 0, SZNAM, coords[2]);
      } else if (csame("rtp",geomi)) {
         ccopy("spherical", 0, SZNAM, geom);
         ccopy("r", 0, SZNAM, coords[0]);
         ccopy("theta", 0, SZNAM, coords[1]);
         ccopy("phi", 0, SZNAM, coords[2]);
#ifdef CURV
      } else if (geo_met(geomi, ndi)) {
         ccopy(geomi, 0, SZNAM, geom);
         for (id=0; id < ndi; id++) {
            ccopy(coordis[id], 0, SZNAM, coords[id]);
         }
         *curv = TRUE;
#endif
      } else {
         printf("error:  unknown model geometry '%s'\n", geom);
         *ios = 1;
      }
#ifdef CURV
      for (id=0; id < ndi; id++) {
         kad[id] = id;
      }
#endif
   }
   return 0;
}

#ifdef PLOTG_SLICE
int model_masked_geom(geomi, ndi, coordis, geom, nd, coords, ndo, xo, 
  kad, curv, ios)
   char geomi[],geom[],coords[MAXD][SZNAM];
   int ndi,*nd,*ndo,*ios;
   char coordis[MAXD][SZNAM];
   int kad[MAXD];
   bool *curv;
   float xo[];
{
   int id,ndom;

   if (*ios) {
   } else {
      *curv = FALSE;
      *nd = ndi;
      ndom = ndi + (*ndo);
      for (id=0; id < ndom; id++) {
         ccopy(coordis[id], 0, SZNAM, coords[id]);
      }
      if (csame("xyz",geomi)) {
         ccopy("cartesian", 0, SZNAM, geom);
      } else if (csame("rzp",geomi)) {
         ccopy("cylidrical", 0, SZNAM, geom);
/*       Convert [r,z,phi] to a surrogate curvilinear geometry */
         if (ndom==3) {
            printf("Note:  using surrogate geometry MUVI_RZP\n");
            ccopy("MUVI_RZP", 0, SZNAM, geom);
            if (geo_met(geom, ndom)) {
               *curv = TRUE;
            } else {
               *ios = 1;
            }
         }
      } else if (csame("rpz",geomi)) {
         ccopy("cylindrical", 0, SZNAM, geom);
      } else if (csame("rtp",geomi)) {
         ccopy("spherical", 0, SZNAM, geom);
      } else if (geo_met(geomi, ndom)) {
         ccopy(geomi, 0, SZNAM, geom);
         *curv = TRUE;
      } else {
         printf("error:  unknown model geometry '%s'\n", geom);
         *ios = 1;
      }
   }
   return 0;
}
#endif
