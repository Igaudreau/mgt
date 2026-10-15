#define MAXD4 (MAXD+1)
#ifdef OLD
#define MAXOB 2000
#define MAXCL 2000
#define MAXDAT 5000
#else
#ifdef STL
#define MAXOB 200
#define MAXCL 200
#define MAXDAT 200000
#else
#define MAXOB 200
#define MAXCL 200
#define MAXDAT 500
#endif
#endif
#define NOP 0
#define ADD 1
#define SUB 2
#define DEGRAD 1.745329E-2

typedef struct csg_type {
   char name[SZNAM];
   char type[SZNAM];
   char title[SZLIN];
   bool noxr[MAXD];
   float xmin[MAXD];
   float xmax[MAXD];
   bool noxmin[MAXD];
   bool noxmax[MAXD];
   int ncl;
   char cltag[MAXCL][SZNAM];
   int clob0[MAXCL];
   int clnob[MAXCL];
   int nob;
   bool noat[MAXOB];
   float at[MAXOB][MAXD];
   bool noor[MAXOB];
   float azim[MAXOB];
   float altd[MAXOB];
   float twistd[MAXOB];
   float ar[MAXOB][MAXD4][MAXD4];
   float arinv[MAXOB][MAXD4][MAXD4];
   int obop[MAXOB];
   char obnam[MAXOB][SZNAM];
   int obndat[MAXOB];
   float *obdat[MAXOB];
   int ndat;
   float dat[MAXDAT];
} csg_struct;
