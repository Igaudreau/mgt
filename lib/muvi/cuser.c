/* cuser.c:  C-Fortran USER mappings

     This module contains routines which ease the calling of
     Fortran 77 user routines from C.

     In the case of SUN/SGI these routines call Fortran 77 directly.
     In the FX case, we pass through an additional Fortran layer
     (GFX.F), which declares private CHARACTER strings and copies
     into them before calling actual user routines.

     In the VAX case, something similar to FX is required (but
     different in its details).

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"
#include "model.h"

#ifdef FX
#define usetm_ usetm__
#define uchks_ uchks__
#define ugets_ ugets__
#ifdef CURV
#define geomet_ geomet__
#define geox_ geox__
#endif
#endif

void user_setup_model(char geom[], int nd, char co[], char un[], 
  char phys[], int nq, char vl[], char vu[], int vmlt[][3], int *ios);
void user_check_scalar(char name[], int *ios);
void user_get_scalar(float q[], char name[], float *value);
bool geo_met(char geom[], int nd);
void geo_x(int lumd, float *z, float *x);	// Changed to pointers BM

/* EXPLICIT DECLARATIONS ADDED BY BM */
extern void cfopy(char a[], int lenb, char b[]);




void user_setup_model(char geom[], int nd, char co[], char un[], 
  char phys[], int nq, char vl[], char vu[], int vmlt[][3], int *ios)
{
   /*
   char geom_f[SZNAM],co_f[SZCO],un_f[SZCO],phys_f[SZNAM],
     vl_f[SZCQ],vu_f[SZCQ];
   int geom_len=SZNAM,co_len=SZCO,un_len=SZCO,phys_len=SZNAM,
     vl_len=SZCQ,vu_len=SZCQ;
   cfopy(geom, SZNAM, geom_f);
   cfopy(co, SZCO, co_f);
   cfopy(un, SZCO, un_f);
   cfopy(phys, SZNAM, phys_f);
   cfopy(vl, SZCQ, vl_f);
   cfopy(vu, SZCQ, vu_f);
   usetm_(geom_f, &nd, co_f, un_f, phys_f, &nq, vl_f, vu_f, vmlt, ios,
     ADRV(geom_len), ADRV(co_len), ADRV(un_len), ADRV(phys_len),
     ADRV(vl_len), ADRV(vu_len));
   */
}

void user_check_scalar(char name[], int *ios)
{
   /*
   char name_f[SZNAM];
   int name_len=SZNAM;
   cfopy(name, SZNAM, name_f);
   uchks_(name_f, ios, ADRV(name_len));
   */
}

void user_get_scalar(float q[], char name[], float *value)
{
   /*
   char name_f[SZNAM];
   int name_len=SZNAM;
   cfopy(name, SZNAM, name_f);
   ugets_(q, name_f, value, ADRV(name_len));
   */
}

#ifdef CURV
bool geo_met(char geom[], int nd)
{
   char geom_f[SZNAM];
   //int geom_len=SZNAM,ios;

   /*   printf("geo_met:  geom=%s,nd=%d\n", geom, nd); */

   cfopy(geom, SZNAM, geom_f);

   /* I don't think this routine exists anywhere so now we just return true; */
   /* geomet_(geom_f, &nd, &ios, ADRV(geom_len)); */
   /* printf("geo_met:  GEOMET returns ios=%d\n", ios); */
   //return(ios==0);
   return 1;
}

void geo_x(int lumd, float *z, float *x)
{
   //int limcon=1,limd=MAXD,ncon=1;
   /*geox_(&limcon, &limd, &ncon, z, x); */
   /*printf("geo_x:  z=[%f,%f] => x=[%f,%f]\n", z[0], z[1], x[0], x[1]); */
}

#endif
