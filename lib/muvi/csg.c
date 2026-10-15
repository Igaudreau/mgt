/* csg.c:  csg module 

     This module handles the CSG Solid Model dataset.

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"
#include "csg.h"

#define VERSION 0.00

static csg_struct csg;
static bool nocsg=TRUE,csg_tagged=FALSE,nocsg_at,nocsg_or;
static float csg_x,csg_y,csg_z,csg_azim,csg_altd,csg_twistd;


/* EXPLICIT DECLARATIONS ADDED BY BM */
bool csg_check_shape(char shape[], int npar, float par[]);
int csg_prepare();
int csg_flush();
int csg_show();
int csg_debug();
int csg_open(char type[], int *ios);
int csg_close(int *ios);
int csg_reset();
int csg_dump(char file[], int *ios);
int csg_load(char file[], int *ios);
int csg_begin(char tag[], int *ios);
int csg_end(int *ios);
int csg_at(float x, float y, float z, int *ios);
int csg_orient(float azim, float altd, float twistd, int *ios);
int csg_op(char op[], char shape[], int npar, float par[], int *ios);
int csg_npar(char shape[], int *npar, int *ios);
int csg_ar(float at_x, float at_y, float at_z, float azim, float altd, float twistd, float ar[MAXD4][MAXD4]);
int csg_arinv(float at_x, float at_y, float at_z, float azim, float altd, float twistd, float ar[MAXD4][MAXD4]);
int csg_tran(float at_x, float at_y, float at_z, float a[MAXD4][MAXD4]);
int csg_mult(float a[MAXD4][MAXD4], float b[MAXD4][MAXD4], float c[MAXD4][MAXD4]);
int csg_copy(float a[MAXD4][MAXD4], float b[MAXD4][MAXD4]);
int csg_rot(int id, float angled, float a[MAXD4][MAXD4]);
int csg_vmult(float a[MAXD4][MAXD4], float x[MAXD4], float y[MAXD4]);
int csg_vmult3(float a[MAXD4][MAXD4], float x[MAXD4], float y[MAXD4]);
int csg_vis(int *ios);
int csg_set0s(char name[], char value[], int *ios);
int csg_set0i(char name[], int value, int *ios);
int csg_set0r(char name[], float value, int *ios);
int csg_set0l(char name[], bool value, int *ios);
int csg_seta(char name[], int *ios);
int csg_zerom(float a[MAXD4][MAXD4]);


/* EXTERNAL ROUTINES - ADDED BY BM*/
extern void ccopy(char a[], int offset, int n, char b[]);
extern bool csame(char a[], char b[]);
extern int rem(int a,int b);
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
extern void vis_begin(void);
extern void vis_csg(csg_struct *cp, int *ios);
extern void vis_end();



bool csg_check_shape(shape, npar, par)
   char shape[];
   int npar;
   float par[];
{
   bool check;

   check = TRUE;
   if (csame("sphere",shape)) {
      if (npar != 1) {
         printf("error:  sphere has wrong number of parameters\n");
         check = FALSE;
      } else if (par[0] <= 0.0) {
         printf("error: sphere has invalid radius\n");
         check = FALSE;
      }
   } else if (csame("cylinder",shape)) {
      if (npar != 2) {
         printf("error:  cylinder has wrong number of parameters\n");
         check = FALSE;
      } else if (par[0] <= 0.0) {
         printf("error:  cylinder has invalid raqdius\n");
         check = FALSE;
      } else if (par[1] <= 0.0) {
         printf("error:  cylinder has invalid height\n");
         check = FALSE;
      }
   } else if (csame("block",shape)) {
      if (npar != 3) {
         printf("error:  block has wrong numberof parameters\n");
         check = FALSE;
      } else if (par[0] <= 0.0) {
         printf("error:  block has invalid x-length\n");
         check = FALSE;
      } else if (par[1] <= 0.0) {
         printf("error:  block has invalid y-length\n");
         check = FALSE;
      } else if (par[2] <= 0.0) {
         printf("error:  block has invalid z-length\n");
         check = FALSE;
      }
   } else if (csame("rect",shape)) {
      if (npar != 2) {
         printf("error:  rect has wrong numberof parameters\n");
         check = FALSE;
      } else if (par[0] <= 0.0) {
         printf("error:  rect has invalid x-length\n");
         check = FALSE;
      } else if (par[1] <= 0.0) {
         printf("error:  rect has invalid y-length\n");
         check = FALSE;
      }
   } else if (csame("cone",shape)) {
      if (npar != 2) {
         printf("error:  cone has wrong numberof parameters\n");
         check = FALSE;
      } else if (par[0] <= 0.0) {
         printf("error:  cone has invalid radius\n");
         check = FALSE;
      } else if (par[1] <= 0.0) {
         printf("error:  cone has invalid height\n");
         check = FALSE;
      }
#ifdef STL
   } else if (csame("stl",shape)) {
      if (npar%12 != 0) {
         printf("error:  stl has wrong numberof parameters\n");
         check = FALSE;
      }
#endif
   } else {
      printf("error:  unknown shape=%s\n", shape);
      check = FALSE;
   }
   return(check);
}

int csg_prepare() {
   csg_reset();
   nocsg = TRUE;
   return 0;
}

int csg_flush() {
   csg_reset();
   nocsg = TRUE;
   return 0;
}

int csg_show() {
   printf("CSG Database:-\n");
   if (nocsg) {
      printf("  Empty.\n");
   } else {
      printf("  name=%s,type=%s\n", csg.name, csg.type);
      printf("  title='%s'\n", csg.title);
      printf("  tag_count=%d,object_count=%d\n", csg.ncl, csg.nob);
   }
   return 0;
}

int csg_debug() {
   char ops[3][SZNAM];
   int iob,icl,i,j;

   ccopy("ADD", 0, SZNAM, ops[ADD]);
   ccopy("SUB", 0, SZNAM, ops[SUB]);
   csg_show();
   if (!nocsg) {
      printf("  Contents:-\n");
      for (icl=0; icl < csg.ncl; icl++) {
         printf("    tag=%s\n", csg.cltag[icl]);
         for (i=0; i < csg.clnob[icl]; i++) {
            iob = csg.clob0[icl] + i;
            printf("      op=%s,shape=%s,par=",ops[csg.obop[iob]],
              csg.obnam[iob]);
            for (j=0; j < csg.obndat[iob]; j++) {
               printf(" %.3E", csg.obdat[iob][j]);
            }
            printf("\n");
            printf("        at=[%.3E,%.3E,%.3E]\n",
              csg.at[iob][0],csg.at[iob][1],csg.at[iob][2]);
            printf("        azim=%.1f,altd=%.1f,twist=%.1f\n",
              csg.azim[iob],csg.altd[iob],csg.twistd[iob]);
         }
      }
   }
   return 0;
}

int csg_open(type, ios)
   char type[];
   int *ios;
{
   csg_reset();
   ccopy(type, 0, SZNAM, csg.type);
   return 0;
}

int csg_close(ios)
   int *ios;
{
   return 0;
}

int csg_reset() {
   int id;

   csg_x = 0.0;
   csg_y = 0.0;
   csg_z = 0.0;
   nocsg_at = TRUE;
   csg_azim = 0.0;
   csg_altd = 0.0;
   csg_twistd = 0.0;
   nocsg_or = TRUE;
   csg_tagged = FALSE;
   ccopy("anon", 0, SZNAM, csg.name);
   ccopy(" ", 0, SZNAM, csg.type);
   ccopy("Scene", 0, SZLIN, csg.title);
   for (id=0; id < MAXD; id++) {
      csg.noxr[id] = TRUE;
      csg.noxmax[id] = TRUE;
      csg.noxmin[id] = TRUE;
   }
   csg.ncl = 0;
   csg.nob = 0;
   csg.ndat = 0;
   nocsg = FALSE;
   return 0;
}

int csg_dump(file, ios)
   char file[];
   int *ios;
{
   int unit,id,icl,iob;

   if (!*ios) {
      ds_create(file, &unit, ios);
      if (*ios) {
         printf("error:  creating file '%s' for CSG ",
           file);
         printf("(csg_dump)\n");
      } else {
         ds_put_begin(unit, "csg", ios);
         if (*ios) {
            printf("error:  dumping data type for CSG file '%s' ",file);
            printf("(csg_dump)\n");
         } else {
            ds_put0r(unit, "version", VERSION, ios);
            ds_put0s(unit, "type", csg.type, ios);
            ds_put0s(unit, "title", csg.title, ios);
            for (id=0; id < MAXD; id++) {
               ds_put0b(unit, "noxr", csg.noxr[id], ios);
               if (!csg.noxr[id]) {
                  ds_put0r(unit, "xmin", csg.xmin[id], ios);
                  ds_put0r(unit, "xmax", csg.xmax[id], ios);
               }
            }
            ds_put0i(unit, "ncl", csg.ncl, ios);
            if (csg.ncl > 0) {
               for (icl=0; icl < csg.ncl; icl++) {
                  ds_put0s(unit, "cltag", csg.cltag[icl], ios);
                  ds_put0i(unit, "clob0", csg.clob0[icl], ios);
                  ds_put0i(unit, "clnob", csg.clnob[icl], ios);
               }
            }
            ds_put0i(unit, "nob", csg.nob, ios);
            if (csg.nob > 0) {
               for (iob=0; iob < csg.nob; iob++) {
                  ds_put0s(unit, "obnam", csg.obnam[iob], ios);
                  ds_put0i(unit, "obop", csg.obop[iob], ios);
                  ds_put0i(unit, "obndat", csg.obndat[iob], ios);
                  ds_put0b(unit, "noat", csg.noat[iob], ios);
                  ds_put1r(unit, "at", MAXD, csg.at[iob], ios);
                  ds_put0b(unit, "noor", csg.noor[iob], ios);
                  ds_put0r(unit, "azim", csg.azim[iob], ios);
                  ds_put0r(unit, "altd", csg.altd[iob], ios);
                  ds_put0r(unit, "twistd", csg.twistd[iob], ios);
                  //ds_put1r(unit, "ar", MAXD4*MAXD4, csg.ar[iob], ios);
                  //ds_put1r(unit, "arinv", MAXD4*MAXD4, csg.arinv[iob],
                  //  ios);
                  ds_put1r(unit, "ar", MAXD4*MAXD4, &csg.ar[iob][0][0], ios);	// BM
                  ds_put1r(unit, "arinv", MAXD4*MAXD4, &csg.arinv[iob][0][0],
                    ios);	// BM
               }
            }
            ds_put0i(unit, "ndat", csg.ndat, ios);
            if (csg.ndat > 0) {
               ds_put1r(unit, "dat", csg.ndat, csg.dat, ios);
            }
            ds_put_end(unit, "csg", ios);
         }
         ds_close(unit, ios);
      }
   }
   return 0;
}

int csg_load(file, ios)
   char file[];
   int *ios;
{
   int unit,id,icl,iob,offset;
   float version;
   char nameb[SZNAM],typeb[SZNAM];

   if (!*ios) {
      ds_open(file, &unit, ios);
      if (*ios) {
         printf("error:  opening file '%s' for CSG ", file);
         printf("csg_load)\n");
      } else {
         if (!ds_get(unit, nameb, typeb, ios)) {
            printf("error:  loading data type from CSG file '%s' ",
              file);
            printf("(csg_load)\n");
         } else if (!csame(typeb,"csg")) {
            printf("error:  expecting csg dataset\n");
            *ios = ERR;
         } else {
            ds_get0r(unit, "version", &version, ios);
            ds_get0s(unit, "type", csg.type, ios);
            ds_get0s(unit, "title", csg.title, ios);
            for (id=0; id < MAXD; id++) {
               ds_get0b(unit, "noxr", &(csg.noxr[id]), ios);
               if (!csg.noxr[id]) {
                  ds_get0r(unit, "xmin", &(csg.xmin[id]), ios);
                  ds_get0r(unit, "xmax", &(csg.xmax[id]), ios);
               }
               csg.noxmin[id] = csg.noxr[id];
               csg.noxmax[id] = csg.noxr[id];
            }
            ds_get0i(unit, "ncl", &(csg.ncl), ios);
            if (*ios) {
            } else if (csg.ncl > MAXCL) {
               printf("error:  too many object clusters\n");
               *ios = ERR;
            } else if (csg.ncl > 0) {
               for (icl=0; icl < csg.ncl; icl++) {
                  ds_get0s(unit, "cltag", csg.cltag[icl], ios);
                  ds_get0i(unit, "clob0", &(csg.clob0[icl]), ios);
                  ds_get0i(unit, "clnob", &(csg.clnob[icl]), ios);
               }
            }
            offset = 0;
            ds_get0i(unit, "nob", &(csg.nob), ios);
            if (*ios) {
            } else if (csg.nob > MAXOB) {
               printf("error:  too many primitive csg objects\n");
               *ios = ERR;
            } else if (csg.nob > 0) {
               for (iob=0; iob < csg.nob; iob++) {
                  ds_get0s(unit, "obnam", csg.obnam[iob], ios);
                  ds_get0i(unit, "obop", &(csg.obop[iob]), ios);
                  ds_get0i(unit, "obndat", &(csg.obndat[iob]), ios);
                  ds_get0b(unit, "noat", &(csg.noat[iob]), ios);
                  ds_get1r(unit, "at", MAXD, csg.at[iob], ios);
                  ds_get0b(unit, "noor", &(csg.noor[iob]), ios);
                  ds_get0r(unit, "azim", &(csg.azim[iob]), ios);
                  ds_get0r(unit, "altd", &(csg.altd[iob]), ios);
                  ds_get0r(unit, "twistd", &(csg.twistd[iob]), ios);
                  //ds_get1r(unit, "ar", MAXD4*MAXD4, csg.ar[iob], ios);
                  //ds_get1r(unit, "arinv", MAXD4*MAXD4, csg.arinv[iob],
                  //  ios);
                  ds_get1r(unit, "ar", MAXD4*MAXD4, &csg.ar[iob][0][0], ios);		// BM
                  ds_get1r(unit, "arinv", MAXD4*MAXD4, &csg.arinv[iob][0][0],
                    ios);	// BM
                  if (!*ios) {
                     csg.obdat[iob] = &(csg.dat[offset]);
                     offset = offset + csg.obndat[iob];
                  }
               }
            }
            ds_get0i(unit, "ndat", &(csg.ndat), ios);
            if (*ios) {
            } else if (csg.ndat > MAXDAT) {
               printf("error:  too much primitive object data\n");
               *ios = ERR;
            } else if (csg.ndat > 0) {
               ds_get1r(unit, "dat", csg.ndat, csg.dat, ios);
            }
            ds_get_end(unit, typeb, ios);
            nocsg = *ios == ERR;
         }
         ds_close(unit, ios);
      }
   }
   return 0;
}

int csg_begin(tag, ios)
   char tag[];
   int *ios;
{
   if (*ios) {
   } else if (csg_tagged) {
      printf("error:  already inside tag-group\n");
      *ios = ERR;
   } else {
      csg_tagged = TRUE;
      if (csg.ncl >= MAXCL) {
         printf("error:  too many tag clusters\n");
         *ios = ERR;
      } else {
         csg.ncl = csg.ncl + 1;
         ccopy(tag, 0, SZNAM, csg.cltag[csg.ncl-1]);
         csg.clob0[csg.ncl-1] = csg.nob;
         csg.clnob[csg.ncl-1] = 0;
      }
   }
   return 0;
}

int csg_end(ios)
   int *ios;
{
   if (*ios) {
      csg_tagged = FALSE;
   } else if (!csg_tagged) {
      printf("error:  not in tag-group\n");
      *ios = ERR;
   } else {
      csg_tagged = FALSE;
   }
   return 0;
}

int csg_at(x, y, z, ios)
   float x,y,z;
   int *ios;
{
   if (!*ios) {
      csg_x = x;
      csg_y = y;
      csg_z = z;
      nocsg_at = FALSE;
   }
   return 0;
}

int csg_orient(azim, altd, twistd, ios)
   float azim,altd,twistd;
   int *ios;
{
   if (!*ios) {
      csg_azim = azim;
      csg_altd = altd;
      csg_twistd = twistd;
      nocsg_or = FALSE;
   }
   return 0;
}

int csg_op(op, shape, npar, par, ios)
   char op[],shape[];
   int npar,*ios;
   float par[];
{
   int i;

   if (*ios) {
   } else if (!csg_tagged) {
      printf("error:  noty in tagged group\n");
      *ios = ERR;
   } else if ((!csame("add",op)) && (!csame("sub",op))) {
      printf("error:  bad csg-op\n");
      *ios = ERR;
   } else if (!csg_check_shape(shape, npar, par)) {
      printf("error:  bad csg-shape specification\n");
      *ios = ERR;
   } else if (csg.nob >= MAXOB) {
      printf("error:  too many primitive objects\n");
      *ios = ERR;
   } else if (csg.ndat+npar > MAXDAT) {
      printf("error:  too much csg-object data\n");
      *ios = ERR;
   } else {
/*   increment number of primitive objects */
      csg.nob = csg.nob + 1;
      csg.clnob[csg.ncl-1] = csg.clnob[csg.ncl-1] + 1;
/*   copy operator as coded integer */
      if (csame("add",op)) {
         csg.obop[csg.nob-1] = ADD;
      } else if (csame("sub",op)) {
         csg.obop[csg.nob-1] = SUB;
      } else {
         csg.obop[csg.nob-1] = NOP;
      }
/*   copy shape */
      ccopy(shape, 0, SZNAM, csg.obnam[csg.nob-1]);
/*   copy at-orient values */
      csg.noat[csg.nob-1] = nocsg_at;
      csg.at[csg.nob-1][0] = csg_x;
      csg.at[csg.nob-1][1] = csg_y;
      csg.at[csg.nob-1][2] = csg_z;
      csg.noor[csg.nob-1] = nocsg_or;
      csg.azim[csg.nob-1] = csg_azim;
      csg.altd[csg.nob-1] = csg_altd;
      csg.twistd[csg.nob-1] = csg_twistd;
/*   calc matrix to transform from world to object coordinates */
      csg_ar(csg_x, csg_y, csg_z, csg_azim, csg_altd,
        csg_twistd, csg.ar[csg.nob-1]);
/*   calc matrix to transform from object to world coordinates */
      csg_arinv(csg_x, csg_y, csg_z, csg_azim, csg_altd,
        csg_twistd, csg.arinv[csg.nob-1]);
/*   reset at-orient values */
      nocsg_at = TRUE;
      csg_x = 0.0;
      csg_y = 0.0;
      csg_z = 0.0;
      nocsg_or = TRUE;
      csg_azim = 0.0;
      csg_altd = 0.0;
      csg_twistd = 0.0;
/*   number of data values */
      csg.obndat[csg.nob-1] = npar;
/*   address of object data array */
      csg.obdat[csg.nob-1] = &(csg.dat[csg.ndat]);
      for (i=0; i < npar; i++) {
         csg.obdat[csg.nob-1][i] = par[i];
      }
      csg.ndat = csg.ndat + npar;
   }
   return 0;
}

int csg_npar(shape, npar, ios)
   char shape[];
   int *npar,*ios;
{
   if (*ios) {
   } else if (csame("sphere",shape)) {
      *npar = 1;
   } else if (csame("cylinder",shape)) {
      *npar = 2;
   } else if (csame("block",shape)) {
      *npar = 3;
   } else if (csame("rect",shape)) {
      *npar = 2;
   } else if (csame("cone",shape)) {
      *npar = 2;
   } else {
      printf("error:  unknown shape=%s\n", shape);
      *ios = ERR;
   }
   return 0;
}

int csg_ar(at_x, at_y, at_z, azim, altd, twistd, ar)
   float azim,altd,twistd,ar[MAXD4][MAXD4],
     at_x,at_y,at_z;
{
   int id,jd;

/* zero out matrix */
   for (id=0; id < MAXD4; id++) {
      for (jd=0; jd < MAXD4; jd++) {
         ar[id][jd] = 0.0;
      }
      ar[id][id] = 1.0;
   }
/* translation */
   csg_tran(at_x, at_y, at_z, ar);
/* azimuthal rotation about z-axis */
   csg_rot(2, azim, ar);
/* altitudinal rotation about y-axis */
   csg_rot(1, altd, ar);
/* twist about resultant z-axis */
   csg_rot(2, twistd, ar);
   return 0;
}

int csg_arinv(at_x, at_y, at_z, azim, altd, twistd, ar)
   float azim,altd,twistd,ar[MAXD4][MAXD4],
     at_x,at_y,at_z;
{
   int id,jd;

/* zero out matrix */
   for (id=0; id < MAXD4; id++) {
      for (jd=0; jd < MAXD4; jd++) {
         ar[id][jd] = 0.0;
      }
      ar[id][id] = 1.0;
   }
/* twist about resultant z-axis */
   csg_rot(2, -twistd, ar);
/* altitudinal rotation about y-axis */
   csg_rot(1, -altd, ar);
/* azimuthal rotation about z-axis */
   csg_rot(2, -azim, ar);
/* translation */
   csg_tran(-at_x, -at_y, -at_z, ar);
   return 0;
}

int csg_tran(at_x, at_y, at_z, a)
   float at_x,at_y,at_z,a[MAXD4][MAXD4];
{
   float t[MAXD4][MAXD4],b[MAXD4][MAXD4];
   int i,j;

/* construct T-matrix */
   for (i=0; i < MAXD4; i++) {
      for (j=0; j < MAXD4; j++)  {
         t[i][j] = 0.0;
      }
      t[i][i] = 1.0;
   }
   t[MAXD][0] = -at_x;
   t[MAXD][1] = -at_y;
   t[MAXD][2] = -at_z;
/* apply to A-matrix */
   csg_mult(t, a, b);
/* copy back into a-matrix */
   csg_copy(b, a);
   return 0;
}

int csg_mult(a, b, c)
   float a[MAXD4][MAXD4],b[MAXD4][MAXD4],c[MAXD4][MAXD4];
{
   int i,j,k;
   float d;

   for (i=0; i < MAXD4; i++) {
      for (j=0; j < MAXD4; j++) {
         d = 0.0;
         for (k=0; k < MAXD4; k++) {
            d = d + a[k][i]*b[j][k];
         }
         c[j][i] = d;
      }
   }
   return 0;
}

int csg_copy(a, b)
   float a[MAXD4][MAXD4],b[MAXD4][MAXD4];
{
   int i,j;

   for (i=0; i < MAXD4; i++) {
      for (j=0; j < MAXD4; j++) {
         b[i][j] = a[i][j];
      }
   }
   return 0;
}

int csg_rot(id, angled, a)
   int id;
   float angled,a[MAXD4][MAXD4];
{
   float r[MAXD4][MAXD4],b[MAXD4][MAXD4],cosa,sina;
   int i,j,jd,kd;

/* construct r-matrix */
   for (i=0; i < MAXD4; i++) {
      for (j=0; j < MAXD4; j++) {
         r[i][j] = 0.0;
      }
   }
   jd = rem(id+1,3);
   kd = rem(id+2,3);
   cosa = cos(angled*DEGRAD);
   sina = sin(angled*DEGRAD);
   r[MAXD][MAXD] = 1.0;
   r[id][id] = 1.0;
   r[jd][jd] = cosa;
   r[kd][kd] = cosa;
   r[kd][jd] = sina;
   r[jd][kd] = -sina;
/* apply to A-matrix,result - B-matrix */
   csg_mult(r, a, b);
/* copy back into A */
   csg_copy(b, a);
   return 0;
}

int csg_vmult(a, x, y)
   float a[MAXD4][MAXD4],x[MAXD4],y[MAXD4];
{
   float d;
   int i,j;

   for (i=0; i < MAXD4; i++) {
      d = 0.0;
      for (j=0; j < MAXD4; j++) {
         d = d + a[j][i]*x[j];
      }
      y[i] = d;
   }
   return 0;
}

int csg_vmult3(a, x, y)
   float a[MAXD4][MAXD4],x[MAXD4],y[MAXD4];
{
   float d;
   int i,j;

   for (i=0; i < MAXD; i++) {
      d = 0.0;
      for (j=0; j < MAXD; j++) {
         d = d + a[j][i]*x[j];
      }
      y[i] = d;
   }
   return 0;
}

int csg_vis(ios)
   int *ios;
{

   if (!*ios) {
      vis_begin();
      if (!nocsg) {
         vis_csg(&csg, ios);
      }
      vis_end();
   }
   return 0;
}

int csg_set0s(name, value, ios)
   char name[];
   int *ios;
   char value[];
{

   if (*ios) {
   } else if (csame("title",name)) {
      ccopy(value, 0, SZLIN, csg.title);
   } else {
      printf("error:  unknown string attribute '%s' (csg_set0s)\n",
        name);
      *ios = ERR;
   }
   return 0;
}

int csg_set0i(name, value, ios)
   char name[];
   int *ios;
   int value;
{

   if (*ios) {
   } else {
      printf("error:  unknown integer attribute '%s' (csg_set0i)\n",
        name);
   }
   return 0;
}

int csg_set0r(name, value, ios)
   char name[];
   int *ios;
   float value;
{

   if (*ios) {
   } else if (csame("xmin",name)) {
      csg.noxmin[0] = FALSE;
      csg.xmin[0] = value;
      csg.noxr[0] = csg.noxmax[0];
   } else if (csame("xmax",name)) {
      csg.noxmax[0] = FALSE;
      csg.xmax[0] = value;
      csg.noxr[0] = csg.noxmin[0];
   } else if (csame("ymin",name)) {
      csg.noxmin[1] = FALSE;
      csg.xmin[1] = value;
      csg.noxr[1] = csg.noxmax[1];
   } else if (csame("ymax",name)) {
      csg.noxmax[1] = FALSE;
      csg.xmax[1] = value;
      csg.noxr[1] = csg.noxmin[1];
   } else if (csame("zmin",name)) {
      csg.noxmin[2] = FALSE;
      csg.xmin[2] = value;
      csg.noxr[2] = csg.noxmax[2];
   } else if (csame("zmax",name)) {
      csg.noxmax[2] = FALSE;
      csg.xmax[2] = value;
      csg.noxr[2] = csg.noxmin[2];
   } else {
      printf("error:  unknown real attribute '%s' (csg_set0r)\n",
        name);
   }
   return 0;
}

int csg_set0l(name, value, ios)
   char name[];
   int *ios;
   bool value;
{

   if (*ios) {
   } else {
      printf("error:  unknown logical attribute '%s' (csg_set0l)\n",
        name);
   }
   return 0;
}

int csg_seta(name, ios)
   char name[];
   int *ios;
{

   if (*ios) {
   } else {
      printf("error:  unknown assertion '%s' (csg_seta)\n",
        name);
   }
   return 0;
}

int csg_zerom(a)
   float a[MAXD4][MAXD4];
{
   int id,jd;

   for (id=0; id < MAXD4; id++) {
      for (jd=0; jd < MAXD4; jd++) {
         a[id][jd] = 0.0;
      }
      a[id][id] = 1.0;
   }
   return 0;
}

