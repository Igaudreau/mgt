#undef MAIN
#define USER
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "user.h"
#include "funcs.h"
#include "para.h"
/*Fortran Routines */
extern "C"
{
void aprox13_init_(float aion[], float zion[], float bion[]);
void timreact13_(int *nd, float pa[], float ydot[]);
void timrat13_(int *nd, float pa[], float RR[]);
void aprox21_init_(float aion[], float zion[], float bion[]);  
void timreact21_(int *nd, float pa[], float ydot[]);
void timrat21_(int *nd, float pa[], float RR[]);      
void torch_init_(int *niso, int *nrat, float aion[], float zion[],
  float bion[], char netname[]);
  void torch_react_(int *nd, int *doweak, int *doneut, float pa[],
  float ydot[]);
void torch_rat_(int *nd, int *doweak, int *doneut,float pa[],
  float RR[]);      
}
/* GSL headers*/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
extern bool get_r(float &value);
extern bool get_i(int &value);
extern bool get_n(char key[]);
extern bool get_cmd(char cmd[]);
extern void assign_part(struct part_list *start);
extern void rm_part(struct part *pptr);
extern void report();
extern void write_var(const char item[], float value);
extern void write_var(const char item[], int size, float value[]); 
extern void write_var(const char item[], int value);
extern void write_var(const char item[], int size, int value[]);
extern void write_var(const char item[], bool value);
extern void write_var(const char item[], char value[]);
extern void read_var(float *value);
extern void read_var(int *value);
extern void read_var(bool *value);
extern void read_var(char *value);
extern void del_allpart(struct cell *cptr);
extern void merge_array(float buf[], int nvalues[]);
extern struct cell *find_fcell(int &ilev, float xc[]);
extern float global_sum(float value);
extern float global_max(float value);
extern void broadcast(bool *value);
extern void broadcast(float buf[], int nvalues);
extern void mk_partlist(struct part_list *&start, 
  struct part_list *&end);
extern struct part_list *mov_part(bool all, struct part_list *pmov);
//extern int write_silo(int flev, char silofile[]);
bool usrqis(float xc[], float qa[]);
float usrss(float pa[]);
float usrssit(float pa[]);
void usrqp(float pa[], float qa[]);
void usrpq(float qa[], float pa[]);
void usrfxq(cell *cptr);
void usrfxp(float pa[]);
bool usrps(cell *cptr, float &s);
void usrsaf(int norder, float dt, cell *cptr, float so[]);
void usrref(cell* cptr, float qref[]);
float usrqsol( cell *cptr);
float usrte(float pa[]);
void usrmu(float pa[], float mu[]);
void usrnu(float pa[], float nu[]);
void usrbd(int id, int bound, join *jptr, float pl[], float pr[]);
bool usrcmd(char key[]);
void usrprt();
void usrwr();
bool usrrd(char vname[]);
float usrsum();
bool usropenxq(char monmod[], int mode);
bool usrrewxq(char monmod[]);
void usrclosexq(const char monmod[], int mode);		// TODO ADDED const BM
void usrmon(char monmod[]);
void usrmalconv();
bool usrmaldatconv();
void usrsmooth(int nav);
void usrphaseplan();
void usrsoncross(float xmax);
void usrmaxab();
void usrtimconv();
bool usrsetint(char monpref[], int mode, float size);
void usrinter(float x, float pa[]);
void usrfindinter(int &n1, int &n2, float xi);
void usrdef();
void usrkeps();
bool usroff();
float usrvel(int id, struct cell *cptr);
float usrbfield(int id, struct cell *cptr);
float usrgfield(int id, struct cell *cptr);
void usrga(int id, struct cell *cptr, struct cell *left,
  struct cell *right);
void usrprompt(char prompt[]);
void usrquit();
void usrpath(char path[]);
float usrdmdt(float pmass, float pa[]);
struct part *usrmkpart(float xp[]);
void usrprate(struct part *pptr, float dt, float pa[], float dvp[],
  float &rho, float &dm);
void usrpartpop(int ilev, float tcur);
bool usrpartps(char pitem[], struct part *pptr, float &s); 
/*-----sootfoil Addition------*/
void usrSootFoil(cell *cptr);
/* Supernova Routines */
float usrBE(float pa[]);
void usrX2Y(float pa[], float X[]);
void usrY2X(float pa[], float X[]);
float usrYE(float pa[]);
float usrAbar(float pa[]);
void usrYAZ(float pa[], float &Ye, float &Abar, float &Zbar,
  float &Z2bar);
float usrein(float temp, float rho, float Abar, float Ye);
float usrcv(float pa[]);
float usrpg(float temp, float rho, float Abar, float Ye);
void usrion(float temp, float rho, float Abar, float he[]);
void usrrad(float temp, float rho, float he[]);
float usrein2te(float tempi, float rho, float ein, float Abar,
  float Ye);
float usrpg2te(float tempi, float rho, float pg, float Abar,
  float Ye);
void usrnadyo(float temp, float rho, float Ye, float &pres,
  float &ein, float &sound);
void usrinitnet();
void usrshock(float pa[], float rho_r, float T_r, float v_r,
  float &rho_l, float &T_l, float &v_l, float &p_l);
void usrsetionam();
int usrgetionind(char name[]);
//file hermite.c
float psi0(float z);
float dpsi0(float z);
float ddpsi0(float z);
float psi1(float z);
float dpsi1(float z);
float ddpsi1(float z);
float psi2(float z);
float dpsi2(float z);
float ddpsi2(float z);
/* bicubic hermite polynomial statement function */
float herm5(int i, int j, float w0t, float w1t, float w2t,
  float w0mt, float w1mt, float w2mt, float w0d, float w1d,
  float w2d, float w0md, float w1md, float w2md);
//end hermite.c
bool usrinithelm();
void usrhelm(float temp, float rho, float Ye, float he[]);
void usrrate(float pa[], float RR[]);
void usrreact(float pa[], float ydot[]);
void usrinitscreen(float z1, float a1, float z2, float a2, int jscreen);
void usrscreen(float temp, float rho, float Zbar, float Abar,
  float z2bar, float z1, float a1, float z2, float a2, int jscreen,
  float& scor);
void usrrieman(float pvl[], float pvr[], float pres[]);
void  rara(float dl, float ul, float pl, float cel, float pvl[],
  float dr, float ur, float pr, float cer, float pvr[], float &d,
  float &u, float &p);
void shosho(float dl, float ul, float pl, float cl, float pvl[],
  float dr, float ur, float pr, float cr, float pvr[], float &d,
  float &u, float &p);
void shora(int is, float ds, float us, float ps, float cs, float pvs[],
  float dr, float ur, float pr, float cer, float pvr[], float &d,
  float &u, float &p);
void usrrki(int is,float du, float &d, float &u, float &p, float &c,
  float pvr[]);
float wave(float p, float di, float pi, float ci, float pa[]);
void usrsteady(float pa[], float so[], float dt);
bool usratmos(float rho_0, float rho_min, float T_0, float M_0,
  float R_0, float M_shell, char mode[], float pa[]);
bool usrcore(float rho_0, float rho_min, float T_0, char mode[],
  float pa[]);
float usrshockpos(bool verbose, int dir, float xn, float qa[],
  float pa[]);
float usrsonpos(float xs, float detv);
void usrpfcell(cell *cptr);
void usrpref(mcell *mcptr, float pref[]);
bool usrinitintp(char filename[]);
void intp(float quan, float qmin, float qmax, float nq, float dq,
  int &i, float &fraci);
// mod.c functions
bool modsol(struct cell *pcptr, bool verbose, bool bisect);
bool modgsl(int nsub, float pa[],float dt);
int modgsldiff(float t,const float y[],float ydot[],
  void* param);
int modjac(float t,const float y[],float* dfdy,float dfdt[],
  void* param);
void modsrc(float pa[],float ydot[]);
void modfix(float pa[]);
extern void sys_cmd(char cmd[]);
#include <stdarg.h>	// BM
int safe_fscanf(FILE *stream, const char *format, ...);	// BM
// TODO Update all fscanf commands to use safe_fscanf


/* GSL structure  */
struct gsl_param {
  float pa[MAXQ];
  int nvar;
};

/* gsl_matrix 
typedef struct 
{
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
  gsl_block * block;
  int owner;
} gsl_matrix;
*/
/* GSL wrappers */
bool usrgsl(cell *cptr, float pa[], float dt);
int usrgslreact(float t, const float y[], float dydt[],
  void *param);
/* GSL Right hand side and Jacobian  */
int usrrhs(float t, const float y[], float ydot[], void *param);
int usrjac(float t, const float y[], float *ydotdy, float dfdt[],
  void *param);
/* Monitor files */
FILE *xatfile[MAXSPECIES], *xaxfile[MAXSPECIES], *tempfilet, *machfilet,
  *rhofilet, *enfilet, *phifilet, *kapfilet, *bindfilet, *pgfilet,
  *ufilet, *abarfilet, *tempfilex, *machfilex, *rhofilex, *enfilex,
   *enifilex, *phifilex, *kapfilex, *bindfilex,*pgfilex, *ufilex,
  *abarfilex, *phimachfilex;
/*Shock speed, sonic pt, Malcolm convert filenames etc*/
float oshockpos, oshockt;
FILE  *shockposfile, *shockvelfile, *sonposfile, *sonwidthfile,
  *shockin, *shockout, *convin, *dkout, *dlcout, *trho, *engen;
/* Atmosphere initial conditions */
int intpts; 
float intrho[10000], inttemp[10000], intdx, intxmax;
FILE  *rhofile, *tempfile;
/* TIM net name */
char netname[SZNAM];
/* Resimulation files */
FILE *resim[MAXQ];
/* Interpolation structure and npoints*/
  struct intdata{
    float x;
    float pa[MAXQ];
  };   
struct intdata **intarray;
int npoints;
float xmin, xmax, ldx, lxmin, lxmax;

//file usrqis.c
/*
+ usrqis - classify and initialise cell, given coordinates.
     Description:
        This initialises qa in the  cell with coordinates xc according
        to cas.
*/
bool usrqis(float xc[], float qa[])
{
  int i, id;
  float pa[MAXQ], temp, rho, ye, X[MAXSPECIES];
  (void)ye;
  float sine_value;
  /* Defaults (safety) */
  pa[iqd] = 1.0;
  for (id = 1; id <= nd; id++)
    pa[iqu0+id] = 0.0;
  pa[iqe] = 1.0;
  for (i = iqe + 1; i < nq; i++)
    pa[i] = 0.0;
  /* Abundances */
  for (i = 0; i < niso; i++)
    X[i] = 0.0;  
  /* Zero all mass fractions for safety */
  for (i = 0; i < niso; i++)
    X[i] = 0.0;  
  if (strcmp(cas, "KHOKHLOV") == 0){
    /* A.M. Khokhlov, MNRAS, 239, 785, 1989. */
    X[ic12] = 0.5;
    X[io16] = 0.5;
    if (xc[0] <= 0.0){
    //    if (xc[0] <= 3.95e2){     // Cylindrical
      rho = 7.310070e+08;
      temp = 3.938624e+09;
      pa[iqu0+1] = 1.21e9 - 4.965753e+08;
    }
    else{
      rho = 3.0e8;
      temp = 2.0e8;
      pa[iqu0+1] = 0.0;      
    }      
  }
  else if (strcmp(cas, "KHO_13_HE") == 0){
    /* A.M. Khokhlov, MNRAS, 239, 785, 1989.  Helium detonation*/
    X[ihe4] = 1.0;
    if (xc[0] <= 0.0){
      rho = 5.747439e+05;
      temp = 2.145964e+09;
      //      pa[iqu0+1] = 1.5e9 - 3.154055e+08;
      pa[iqu0+1] = 0.0;
      if (nal > 13)
	pa[iqal0+nal] = 1.0;
    }
    else{
      rho = rho_0;
      temp = T_0;
      pa[iqu0+1] = 0.0;
      if (nal > 13)
	pa[iqal0+nal] = 0.0;      
    }      
  }  
  else if (strcmp(cas, "TIMMIX") == 0){
    temp = T_0;
    rho =  rho_0;
    if (tim7 || tim13)    
      X[ihe4] = 0.9;
    else{
      X[ihe4] = 0.891;
      X[in14] = 0.009;
    }
    X[ic12] = 0.05;          
    X[io16] = 0.05;      
  }
  else if (strcmp(cas, "MIX1D") == 0){
    if (xc[0] < 0.1e6){
      if (niso == 27)
	X[iti44] = 1.0;
      else
	X[ini56] = 1.0;	
      temp = T_s;
      rho = rho_s;
      //      pa[iqu0+1] = 1.0281311e9 - 2.113653e+08;
    }
    else{
      if (tim7 || tim13)    
	X[ihe4] = 0.9;
      else{
	X[ihe4] = 0.891;
	X[in14] = 0.009;
      }
      X[ic12] = 0.05;          
      X[io16] = 0.05;  
      temp = T_0;
      rho = rho_0;      
    }
  }
  else if (strcmp(cas, "ATMOS_MIX_INIT_X") == 0){
    if ((niso == 7) || (niso == 13))
      X[ihe4] = 0.9;
    else{
      X[ihe4] = 0.891;
      X[in14] = 0.009;
    }
    X[ic12] = 0.05;          
    X[io16] = 0.05;
    if (nd == 1)
      usrinter(xc[0], pa);
    else
      usrinter(xc[1], pa);
    temp = T_0;    
    if (xc[0] > xmax)
      rho = rho_0;
    else
      rho = pa[iqd];
    float tv;
    tv = xc[0]*xc[0] + xc[1]*xc[1];
    if (tv < depth*depth){
      temp = T_s;
      rho = rho_s;
    }    
  }    
  else if (strcmp(cas, "ATMOS_MIX_XY") == 0){
    if ((niso == 7) || (niso == 13))
      X[ihe4] = 0.9;
    else{
      X[ihe4] = 0.891;
      X[in14] = 0.009;
    }
    X[ic12] = 0.05;          
    X[io16] = 0.05;
    if ((xc[0] < 0.0) && (xc[1] < depth)){
      temp = T_s;
      rho = rho_s;
    }
    else{
      temp = T_0;
      rho = rho_0*exp(-xc[1]/depth);      
    }
  }  
  else if (strcmp(cas, "MIX_13_XY") == 0){
    if (xc[0] < 0.1e6){
      temp = T_s;
      rho = rho_s;
      //      pa[iqu0+1] = 1.0281311e9 - 2.113653e+08;
      //      X[ini56] = 1.0; //Ni_56
      X[ini56] = 1.0;
    }
    else{
      temp = 3.0e8;
      //      rho = rho_0*(1.0 + 0.05*sin(gkpi*xc[1]*2.0/(xr[1]- xl[1])));
      rho = rho_0;
      X[ihe4] = 0.90;  //He_4
      X[ic12] = 0.05;  //C_12
      X[io16] = 0.05;  //O_16  
    }
  }
  else if (strcmp(cas, "MIX_46") == 0){
    if (xc[0] < 0.1e6){
      temp = T_s;
      rho = rho_s;
      //      pa[iqu0+1] = 1.0281311e9 - 2.113653e+08;
      //      X[ini56] = 1.0; //Ni_56
      X[ini56] = 1.0;
    }
    else{
      temp = 3.0e8;
      //      rho = rho_0*(1.0 + 0.05*sin(gkpi*xc[1]*2.0/(xr[1]- xl[1])));
      rho = rho_0;
      X[ihe4] = 0.891;  //He_4
      X[in14] = 0.009;  //N_14      
      X[ic12] = 0.05;  //C_12
      X[io16] = 0.05;  //O_16  
    }
  }  
  else if (strcmp(cas, "MIX_46_INIT") == 0){
    usrinter(xc[0], pa);
    if (xc[0] > xmax){
      rho = rho_0;
      temp = T_0;
    }
    else{
      rho = pa[iqd];
      temp = pa[iqte];
    }
    for (i = 0; i < niso; i++)
      X[i] = pa[iqal0+1+ i];
  } 
  else if (strcmp(cas, "TIMMESHE") == 0){
    temp = T_0;
    rho = rho_0;
    X[ihe4] = 1.0;
  }
  else if (strcmp(cas, "STEADYMIX") == 0){
    temp = T_0;
    rho = rho_0;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    if ((niso ==13) || (niso == 7))
      X[ihe4] = 0.9;      
    else{
      X[ihe4] = 0.891;
      X[in14] = 0.009;
    }
  }  
  else if (strcmp(cas, "STEADYNON46") == 0){
    temp = T_0;
    rho = rho_0;
    X[ihe4] = 0.9;
    X[ic12] = 0.05;
    X[io16] = 0.05;
  }  
  else if (strcmp(cas, "STEADYC46") == 0){
    temp = T_0;
    rho = rho_0;
    X[ihe4] = 0.9;
    X[ic12] = 0.1;
  }
  else if (strcmp(cas, "STEADYCO46") == 0){
    temp = T_0;
    rho = rho_0;
    X[io16] = 0.5;
    X[ic12] = 0.5;
  }    
  else if (strcmp(cas, "MIX") == 0){
    if (xc[0] < 0.0){
      temp = T_s;
      rho = rho_s;
      X[ini56] = 1.0;
      pa[iqu0+1] = 1.0406e+09 - 6.320852e+08;
      pa[iqal0+47] = 1.0;
    }
    else{
      X[ihe4] = 0.891;
      X[in14] = 0.009;
      X[ic12] = 0.05;
      X[io16] = 0.05;
      temp = T_0;
      rho = rho_0;      
    }
  }
  else if (strcmp(cas, "MIX_SHELL") == 0){
    float fraci;
    if (xc[0] < 0.0){
      temp = T_s;
      rho = rho_s;
      X[ini56] = 1.0;
    }
    else{
      X[ihe4] = 0.891;
      X[in14] = 0.009;
      X[ic12] = 0.05;
      X[io16] = 0.05;
      intp(xc[1], 0, intxmax, intpts, intdx, i, fraci);
      rho = intrho[i] + fraci*(intrho[i+1] - intrho[i]);
      temp = inttemp[i] + fraci*(inttemp[i+1] - inttemp[i]);            
    }
  }     
  else if (strcmp(cas, "NOMOTOA") == 0){
    temp = 2.45471e7;
    rho = 1.202264e8;
    X[ic12] = 1.0;
    if (xc[0] < 10.0){
      rho = 3.465556e+08;
      temp = 4.278478e+09;
      pa[iqu0+1] = 1.2e9 - 4.163017e+08;
    }
  }
  else if (strcmp(cas, "SHARPE1") == 0){
    X[ic12] = 0.5;
    X[io16] = 0.5;
    if (xc[0] < 0.0){
      rho = 2.887717e+08;
      temp = 4.106040e+09;
      pa[iqu0+1] = 1.161e9 - 4.020477e+08;
    }
    else{
      rho = 1.0e8;
      temp = 1.0e8;
      pa[iqu0+1] = 0;
    }
  }
  else if (strcmp(cas, "GAMEZO") == 0){
    int ix, iy;
    /* Gamezo, APJ, 512, 827, 1999. */
    X[ic12] = 0.5;
    X[io16] = 0.5;
    // carbon scale
    //double A = 25.0;
    //double lam = 0.05 * M_PI;
    // Si scale
    //sine_value = amplitude * sin(period * M_PI * xc[1]);
    // First do 1D evolution
    //if (xc[0] <= sine_value){
    //if (xc[0] <= 0.0){
    if (nd < 2) {
      if (xc[0] < 0.0){
        // Post shock from Timmes program
        rho = 2.255E+07;
        temp = 3.974E+09;
        pa[iqu0+1] = 1.176E+09 - 2.607E+08;
      }
      else{
        rho = 5.0e6;
        temp = 2.0e8;
        pa[iqu0+1] = 0.0;
      }
    }
    else {
      /* Set the orientation of the domain */
      if (strcmp(coord[0], "X") == 0) {
        ix = 0;
        iy = 1;
      }
      else {
        ix = 1;		// if x is the second coordinate then reset the index.
        iy = 0;
      }
      //if (xc[0] <= 0.0 && xc[1] <= 0.0){
      sine_value = amplitude * sin(period * M_PI * xc[iy]);
      if (xc[ix] <= sine_value){
        rho = 2.255E+07;
        temp = 3.974E+09;
        pa[iqu0+1+ix] = 1.176E+09 - 2.607E+08;
      }
      else{
        rho = 5.0e6;
        temp = 2.0e8;
        pa[iqu0+1+ix] = 0.0;
      }
    }

  }
  else if (strcmp(cas, "HELIUMMIX") == 0){
    /* Gamezo, APJ, 512, 827, 1999. */
    X[ihe4] = 0.6;
    X[ic12] = 0.2;
    X[io16] = 0.2;
    double A = 50.0;
    double lam = 0.1 * M_PI;
    double sine_value = A * sin(lam * xc[1]);
    /*if (xc[0] <= 0.0){*/
    if (xc[0] <= sine_value){
    //    if (xc[0] <= 3.95e2){     // Cylindrical
      rho = 5.323E+06;
      temp = 3.530E+09;
      pa[iqu0+1] = 1.396E+09 - 2.623E+08;
    }
    else{
      rho = 1.0e6;
      temp = 2.0e8;
      pa[iqu0+1] = 0.0;      
    }      
  }
/*==========================HeCON==================================*/
  else if (strcmp(cas, "HECOZND13") == 0){
    X[ihe4] = 0.9;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    rho = rho_0;
    temp = T_0;
    pa[iqu0+1] = 0.0;
  }
  else if (strcmp(cas, "HECONZND") == 0){
    X[ihe4] = 0.891;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    X[in14] = 0.009;
    rho = rho_0;
    temp = T_0;
    pa[iqu0+1] = 0.0;
  }
  else if (strcmp(cas, "HECO_13ISO_K1E-9") == 0){
    X[ihe4] = 0.9;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    sine_value = amplitude * sin(period * M_PI * xc[1]);
    if (xc[0] <= sine_value){
      rho = 9.729e5;
      temp = 2.044e9;
      pa[iqu0+1] = 1.028e9 - 2.114e8;
    }
    else{
      rho = 2e5;
      temp = 3e8;
      pa[iqu0+1] = 0.0;
    }
  }
  else if (strcmp(cas, "HECO_13ISO_K1E-8") == 0){
    X[ihe4] = 0.9;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    sine_value = amplitude * sin(period * M_PI * xc[1]);
    if (xc[0] <= sine_value){
      rho = 8.102e5;
      temp = 1.468e9;
      pa[iqu0+1] = 7.455e8 - 1.840e8;
    }
    else{
      rho = 2e5;
      temp = 3e8;
      pa[iqu0+1] = 0.0;
    }
  }
  else if (strcmp(cas, "HECO_13ISO_K105E-8") == 0){
    X[ihe4] = 0.9;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    sine_value = amplitude * sin(period * M_PI * xc[1]);
    if (xc[0] <= sine_value){
      rho = 6.010e5;
      temp = 8.473e8;
      pa[iqu0+1] = 4.852e8 - 1.615e8;
    }
    else{
      rho = 2e5;
      temp = 3e8;
      pa[iqu0+1] = 0.0;
    }
  }
  else if (strcmp(cas, "HECO_27ISO_K1E-8") == 0){
    X[ihe4] = 0.891;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    X[in14] = 0.009; 
    sine_value = amplitude * sin(period * M_PI * xc[1]);
    if (xc[0] <= sine_value){
      rho = 8.811e5;
      temp = 1.717e9;
      pa[iqu0+1] = 8.585e8 - 1.949e8;
    }
    else{
      rho = 2e5;
      temp = 3e8;
      pa[iqu0+1] = 0.0;
    }
  }
  else if (strcmp(cas, "HECO_27ISO_K1E-7") == 0){
    X[ihe4] = 0.891;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    X[in14] = 0.009;
    sine_value = amplitude * sin(period * M_PI * xc[1]);
    if (xc[0] <= sine_value){
      rho = 8.034e5;
      temp = 1.443e9;
      pa[iqu0+1] = 7.343e8 - 1.828e8;
    }
    else{
      rho = 2e5;
      temp = 3e8;
      pa[iqu0+1] = 0.0;
    }
  }
  else if (strcmp(cas, "HECO_34ISO_K1E-8") == 0){
    X[ihe4] = 0.891;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    X[in14] = 0.009; 
    sine_value = amplitude * sin(period * M_PI * xc[1]);
    if (xc[0] <= sine_value){
      rho = 8.811e5;
      temp = 1.717e9;
      pa[iqu0+1] = 8.585e8 - 1.949e8;
    }
    else{
      rho = 2e5;
      temp = 3e8;
      pa[iqu0+1] = 0.0;
    }
  }
  else if (strcmp(cas, "HECO_34ISO_K1E-7") == 0){
    X[ihe4] = 0.891;
    X[ic12] = 0.05;
    X[io16] = 0.05;
    X[in14] = 0.009;
    sine_value = amplitude * sin(period * M_PI * xc[1]);
    if (xc[0] <= sine_value){
      rho = 8.040e5;
      temp = 1.445e9;
      pa[iqu0+1] = 7.351e8 - 1.829e8;
    }
    else{
      rho = 2e5;
      temp = 3e8;
      pa[iqu0+1] = 0.0;
    }
  }
/*==========================HeCON End==============================*/
  else{
    /* Case not found return false */
    return(false);
  }  
  /* Put abundances in pa using mass fraction X */
  usrX2Y(pa, X);
  pa[iqd] = rho;
  pa[iqte] = temp;    
  pa[iqe] = usrpg(temp, rho, usrAbar(pa), usrYE(pa));
  usrqp(pa, qa);
  return(true);
}
//end usrqis.c

//file usrss.c
/*
  + usrss - Frozen sound speed
*/
float usrss(float pa[])
{
  float rho, temp, Abar, Ye, dpdd, dpdt, dsdd, dsdt, he[12];
  Abar = usrAbar(pa);
  Ye = usrYE(pa);
  rho = pa[iqd];
  temp = usrpg2te(pa[iqte], rho, pa[iqe], Abar, Ye);    
  dpdd = 0.0;
  dpdt = 0.0;      
  dsdd = 0.0;
  dsdt = 0.0;
  usrhelm(temp, rho, Ye, he);
  dpdd += he[3];
  dpdt += he[4];      
  dsdd += he[7];
  dsdt += he[8];
  usrion(temp, rho, Abar, he);
  dpdd += he[3];
  dpdt += he[4];      
  dsdd += he[7];
  dsdt += he[8];    
  usrrad(temp, rho, he);
  dpdd += he[3];
  dpdt += he[4];      
  dsdd += he[7];
  dsdt += he[8];
  /* Sound speed from equation 5 in sound.pdf */
  return(sqrt(dpdd - dpdt*dsdd/dsdt));
}
//end usrss.c

//file usrssit.c
/*
+ usrssit - Isothermal sound speed
*/
float usrssit(float pa[])
{
  float rho, temp, Abar, Ye, dpdd, he[12];
  Abar = usrAbar(pa);
  Ye = usrYE(pa);
  rho = pa[iqd];
  temp = pa[iqte];
  dpdd = 0.0;
  usrhelm(temp, rho, Ye, he);
  dpdd += he[3];
  usrion(temp, rho, Abar, he);
  dpdd += he[3];
  usrrad(temp, rho, he);
  dpdd += he[3];
  /* Sound speed */
  return(sqrt(dpdd));
}
//end usrssit.c

//file usrqp.c
/*
+ usrqp - conserved  variables from primitive
*/
void usrqp(float pa[], float qa[])
{
  float rho, divrho, ke, u, tv;
  (void)divrho, (void)tv;
  int i, id;
  rho = pa[iqd];
  qa[iqd] = rho;
  ke = 0.0;
  for (id = 1; id <= nd; id++){
    u = pa[iqu0+id];
    qa[iqu0+id] = u*rho;
    ke += u*u;
  }
  /* Supernova EOS */
  float Ye;
  Ye = usrYE(pa);
  //    qa[iqe] = rho*usrein(pa[iqte], rho, usrAbar(pa), usrYE(pa))
  qa[iqe] = rho*usrein(pa[iqte], rho, usrAbar(pa), Ye)    
    + 0.5*rho*ke;
  for (i = iqal0 + 1; i <= iqal0 + scad; i++)
    qa[i] = rho*pa[i];
  for (i = iqal0 + scad + 1; i < nq; i++)
    qa[i] = pa[i];
  if (inten)
    qa[iqie] = rho*pa[iqie];
  return;
}
//end usrqp.c

//file usrpq.c
/*
+ usrpq - primitive variables from conserved
*/
void usrpq(float qa[], float pa[])
{
  float rho, divrho, ke, u, Abar, Ye, ein, temp;
  int i, id;
  rho = max(qa[iqd], dfloor);
  divrho = 1.0/rho;
  pa[iqd] = rho;
  ke = 0.0;
  for (id = 1; id <= nd; id++){
    u = qa[iqu0+id]*divrho;
    pa[iqu0+id] = u;
    ke += u*u;
  }
  ke *= 0.5*rho;
  /* Do advecting scalars for abundances */
  for (i = iqe + 1; i <= iqal0 + scad; i++)
    pa[i] = qa[i]*divrho;
  /* Do non-advecting scalars (temperature etc) */  
  for (i = iqal0 + scad + 1; i < nq; i++)
    pa[i] = qa[i];
  if (inten)
    pa[iqie] = qa[iqie]*divrho;
  /* Supernova EOS */
  ein = (qa[iqe] - ke)*divrho;
  Abar = usrAbar(pa);
  Ye = usrYE(pa);
  /* Note does not recalculate pa[iqte] (temperature) */
  pa[iqe] = usrpg(qa[iqte], rho, Abar, Ye);
  if (fxp)
    usrfxp(pa);
  else
    pa[iqe] = max(pfloor,pa[iqe]);  
  return;
}
//end usrpq.c

//file usrfxq.c
/*
+ usrfxq - Fix up solution.
     Description :
        This fixes the  solution in a cell to  ensure energy, pressure
        and density are positive. It  is given the cell coordinates in
        cptr->xc[], so  that it  can also be  used for  other purposes
        e.g. forcing  the solution  to take a  given value  at certain
        points.
*/
void usrfxq(struct cell * cptr)
{
  int i, id, iymax;
  float ke, rho, divrho, temp, mu, ein, Abar, Ye, qa[MAXQ], pa[MAXQ],
    xc[MAXQ], tv, sum;
  rho = max(cptr->qa[iqd], dfloor);    
  divrho = 1.0/rho;
  cptr->qa[iqd] = rho;
  /* Thermal energy */
  ke = 0.0;
  for (id = 1; id <= nd; id++){
    mu = cptr->qa[iqu0+id];
    ke += mu*mu;
  }
  ke *= 0.5*divrho;
  ein = (cptr->qa[iqe] - ke)*divrho;
  /* Scalars times atomic weight in 0-rho */
  for (i = 0; i < niso; i++)
    cptr->qa[iqal0+1+i] = min(cptr->qa[iqd]/Ai[i],
      max(cptr->qa[iqal0+1+i], 0.0));
  if (cons){    
    /* Find max mass */
    pa[iqd] = cptr->qa[iqd];
    pa[iqte] = cptr->qa[iqte];
    for (i = 0; i < niso; i++)
      pa[iqal0+1+i] = cptr->qa[iqal0+1+i]*divrho;
    tv = -1.0;
    for (i = 0; i < niso; i++){
      if (pa[iqal0+1+i]*Ai[i] > tv){
	tv = pa[iqal0+1+i]*Ai[i];
	iymax = i;
      }
    }
    /* Max y[i] from conservation */
    sum = 0.0;
    for (i = 0; i < niso; i++){
      if (i != iymax)
	sum += cptr->qa[iqal0+1+i]*Ai[i];
    }
    cptr->qa[iqal0+1+iymax] = (rho - sum)/Ai[iymax];
  }
  /* Put iymax abundance in pa */
  pa[iqal0+1+iymax] = cptr->qa[iqal0+1+iymax]*divrho;
  if (!isothermal){
    Abar = usrAbar(pa);
    Ye = usrYE(pa);    
    temp =  usrein2te(cptr->qa[iqte], rho, ein, Abar, Ye);
    cptr->qa[iqte] = temp;
  }
  /*
  cptr->qa[iqal0+scad] = min(cptr->qa[iqd],
    max(cptr->qa[iqal0+scad], 0.0));
  */
  if (keps){
    cptr->qa[iqk] = max(cptr->qa[iqk], rho*kmin);
    cptr->qa[iqeps] = max(cptr->qa[iqeps], rho*epsmin);
  }
  return;
}
//end usrfxq.c

//file usrfxp.c
/*
+ usrfxp - Fix up primitive solution.
     Description :
        This  fixes the  solution in  pa to  ensure that  pressure and
        density are positive.
*/
void usrfxp(float pa[])
{
  int i;
  /* Fix density */
  pa[iqd] = max(pa[iqd], dfloor);
  /* Fix pressure */
  pa[iqe] = max(pa[iqe], pfloor);
  if (keps){
    pa[iqk] = max(pa[iqk], kmin);
    pa[iqeps] = max(pa[iqeps], epsmin);
  }
  /* Scalars */
  for (i = 1; i <= niso; i++)
    pa[iqal0+i] = min(1.0/Ai[i-1], max(pa[iqal0+i], 0.0));
  //  pa[iqal0+scad] = min(1.0, max(pa[iqal0+scad], 0.0));   
  return;
}
//end usrfxp.c

//file usrps.c
/* 
+ usrps - construct plot scalar
     Description :
        This  returns  the value  of  a plot  scalar  in  s given  the
        conserved solution  in qa. The  scalar is defined by  ITEM and
        returns TRUE if the plot item is found, FALSE otherwise.
*/
bool usrps(cell *cptr, float &s)
{
  int i, j, id, ilev;
  bool found;
  float qref[MAXQ], pa[MAXQ], tv[MAXQ], temp;
  struct cell *left, *right;
  found = false;
  if (strncmp(item, "QA", 2) == 0){
    sscanf(&item[2], "%d", &i);
    if ((i < 1) || (i > nq))
      found = false;
    else{
      s = cptr->qa[i-1];
      found = true;
    }
  }
  else if (strncmp(item, "QR", 2) ==0){
    sscanf(&item[2], "%d", &i);
    if ((i < 1) || (i > nq))
      found = false;
    else{
      s = cptr->dqdt[i-1]/cptr->vol;
      found = true;
    }
  }
  else if (strcmp(item, "RHO") == 0){
    s = cptr->qa[iqd];
    found = true;
  }
  else if (strcmp(item, "U1") == 0){
    s = cptr->qa[iqu0+1]/cptr->qa[iqd];
    found = true;
  }
  else if (strcmp(item, "U2") == 0){
    if (nd >= 2){
      s = cptr->qa[iqu0+2]/cptr->qa[iqd];
      found = true;
    }
    else
      found = false;
  }
  else if (strcmp(item, "U3") == 0){
    if (nd == 3){
      s = cptr->qa[iqu0+3]/cptr->qa[iqd];
      found = true;
    }
    else
      found = false;
  }
  else if (strcmp(item, "UT") == 0){
    s = 0.0;
    for (id = 1; id <= nd; id++){
      temp = cptr->qa[iqu0+id]/cptr->qa[iqd];
      s += temp*temp;
    }
    s = sqrt(s);
    found = true;
  }
  else if (strcmp(item, "POT") == 0){
    if (grav){
      s = gconst*cptr->qa[iqpot];
      found = true;
    }
    else
      found = false;
  }
  else if (strcmp(item, "DIVV") == 0){
    float tv;
    s = 0.0;
    for (id = 1; id <= nd; id++){
      i = 0;
      tv = 0.0;
      left = cptr->ljoin[id-1]->left;
      if (left != NULL){
	tv += (cptr->qa[iqu0+id] - left->qa[iqu0+id])
	  /(cptr->xc[id-1] - left->xc[id-1]);
	i++;
      }
      right = cptr->rjoin[id-1]->right;
      if (right != NULL){
	tv += (right->qa[iqu0+id] - cptr->qa[iqu0+id])
	  /(right->xc[id-1] - cptr->xc[id-1]);
	i++;
      }
      s += tv/(float)i;
    }
    found = true;
  }  
  else if (strncmp(item, "GA", 2) ==0){
    sscanf(&item[2], "%d", &id);
    if ((!grav) || (id < 1) || (id > nd))
      found = false;
    else{
      i = 0;
      s = 0.0;
      left = cptr->ljoin[id-1]->left;
      if (left != NULL){
	s += cptr->qa[iqpot] - left->qa[iqpot];
	i++;
      }
      right = cptr->rjoin[id-1]->right;
      if (right != NULL){
	s +=  right->qa[iqpot] - cptr->qa[iqpot];
	i++;
      }
      if (i == 0)
	found = false;
      else
	s /= - (float)(i)*
	  (cptr->rjoin[id-1]->xj[id-1] - cptr->ljoin[id-1]->xj[id-1])
	  /gconst;
      found = true;
    }
  }
  else if (strcmp(item, "QSOL") == 0){
    s = cptr->qsol;
    found = true;
  }
  else if (strcmp(item, "ICELL") == 0){
    s = cptr->icell;
    found = true;
  }
  else if (strcmp(item, "PROC") == 0){
    s = (float)mype;
    found = true;
  }
  else if (strcmp(item, "WORK") == 0){
    i = (int)((cptr->xc[nd-1] - xl[nd-1])/lev[nlevs-1]->dxbar[nd-1]);
    s = (float)wka[i];
    found = true;
  }
  else{
    usrpq(cptr->qa, pa);
    if (strncmp(item, "PA", 2) == 0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nq))
	found = false;
      else{
	s = pa[i-1];
	found = true;
      }
    }
    // abundance
    else if (strncmp(item, "AL", 2) == 0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nal))
	found = false;
      else{
	s = pa[iqal0+i];
	found = true;
      }
    }
    // mass fraction
    else if (strncmp(item, "XAI", 3) == 0){
      sscanf(&item[3], "%d", &i);
      if ((i < 1) || (i > nal))
	found = false;
      else{
	s = pa[iqal0+i]*Ai[i-1];
	found = true;
      }
    }
    else if (strncmp(item, "XAN", 3) == 0){
      char name[SZNAM], tail[SZNAM];
      int len, ind;
      for (i = 0; i < 10; i++)
	tail[i] = item[i+3];
      len = strlen(tail);
      for (i = 0; i < len; i++)
	tail[i] = tolower(tail[i]);            
      ind = usrgetionind(tail);      
      s = pa[iqal0+ind+1]*Ai[ind];
      found = true;
    }
    else if (strcmp(item, "PG") == 0){
      s = pa[iqe];
      found = true;
    }
    else if (strcmp(item, "INTEN") == 0){
      if (inten){
	s = pa[iqie];
	found = true;
      }
      else{
	float Abar, Ye;
	Abar = usrAbar(pa);
	Ye = usrYE(pa);	
	s = usrein(pa[iqte], pa[iqd], Abar, Ye);
	found = true;
      }
    }
    else if (strcmp(item, "C") == 0){
      s = usrss(pa);
      found = true;
    }
    else if (strcmp(item, "WS") == 0){
      s = mabs(pa[iqu0+1]) + usrss(pa);
      found = true;
    }    
    else if (strcmp(item, "GG") == 0){
      s = usrss(pa);
      s *= s;
      s = pa[iqd]*s/pa[iqe];
      found = true;
    }
    else if (strcmp(item, "MACH") == 0){
      s = 0.0;
      for (id = 1; id <= nd; id++)
	s += pa[iqu0+id]*pa[iqu0+id];
      s = sqrt(s)/usrss(pa);
      found = true;
    }
    else if (strcmp(item, "TE") == 0){
      s = usrte(pa);
      found = true;
    }
    /*-----sootfoil Addition------*/
    else if (strcmp(item, "SF") == 0) {
      if (sootfoil) {
        s = cptr->qa[iqsf];
        found = true;
      } else {
        found = false;
      }
    }
    // Refinement level
    else if (strcmp(item, "REF") == 0) {
      struct cell *fcptr;
      fcptr = find_fcell(ilev, cptr->xc);
      s = (float)ilev;
      found = true;
    }
    /*----------------------------*/
    else if (strcmp(item, "K") == 0){
      if (!keps)
	found = false;
      else{
	s = pa[iqk];
	found = true;
      }
    }
    else if (strcmp(item, "EPS") == 0){
      if (!keps)
	found = false;
      else{
	s = pa[iqeps];
	found = true;
      }
    }
    else if (strcmp(item, "TL") == 0){
      if (!keps)
	found = false;
      else{
	s = ceps*pa[iqk]*sqrt(pa[iqk])/max(pa[iqeps],epsmin);
	found = true;
      }
    }
    else if (strncmp(item, "MU", 2) ==0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nq))
	found = false;
      else{
	usrmu(pa, tv);
	s = tv[i];
	found = true;
      }
    }
    else if (strncmp(item, "SA", 2) ==0){
      char old_ode[SZNAM];
      strcpy(old_ode, ode);
      strcpy(ode, "NONE");      
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nq))
	found = false;
      else{
	for (j = 0; j < nq; j++)
	  cptr->pa[j] = pa[j];
	usrsaf(1, 1.0e-8, cptr, tv);
	s = tv[i-1];
	found = true;
      }
      strcpy(ode, old_ode);          
    }
    else if (strcmp(item, "PROGSO") ==0){
      char old_ode[SZNAM];
      strcpy(old_ode, ode);
      strcpy(ode, "NONE");      
      for (i = 0; i < nq; i++)
	cptr->pa[i] = pa[i];
      usrsaf(1, 1.0e-8, cptr, tv);
      s = tv[iqal0+nal]/pa[iqd];
      found = true;
      strcpy(ode, old_ode);
    }    
    else if (strcmp(item, "SUMXA") == 0){
      s = 0.0;
      for (i = 1; i <= niso; i++)
	s += pa[iqal0+i]*Ai[i-1];
      found = true;
    }    
    else if (strcmp(item, "SUMXSO") == 0){
      float ydot[MAXSPECIES+1];
      usrreact(pa, ydot);
      s = 0.0;
      for (i = 0; i < niso; i++)
	s += ydot[i]*Ai[i];
      found = true;
    }
    else if (strcmp(item, "IYMAX") == 0){
      int iymax;
      float tv, ydot[MAXSPECIES+1];
      if (!cons){
	if (mype == 0)
	  cout << "Error : Cons not on \n";
      }
      else{      
	/* Find max rate  */
	usrreact(pa, ydot);
	tv = -1.0;
	for (i = 0; i < niso; i++){
	  if (ydot[i]*Ai[i] > tv){
	    tv = ydot[i]*Ai[i];
	    iymax = i;
	  }
	}      
	s = iymax;
	found = true;
      }
    }
    else if (strcmp(item, "XSOW") == 0){
      float ydot[MAXSPECIES+1];
      bool old_isothermal;
      old_isothermal = isothermal;
      isothermal = false;
      usrreact(pa, ydot);
      isothermal = old_isothermal;
      s = 0.0;
      for (i = 0; i < niso; i++)
	s += (mn*(Ai[i] - Zi[i]) + mp*Zi[i])*ydot[i];
      s *= enuc_conv2*usrcv(pa);
      found = true;
    }
    else if (strcmp(item, "XSOEMAXI") == 0){
      float ydot[MAXSPECIES+1], tv;
      usrreact(pa, ydot);
      s = -1.0;
      tv = -1.0;
      for (i = 0; i < niso; i++){
	if (fabs(ydot[i]*Qi[i]) > tv){
	  tv = fabs(ydot[i]*Qi[i]);
	  s = (float)i;
	}
      }      
      found = true;
    }
    else if (strcmp(item, "XSOEMAX") == 0){
      float ydot[MAXSPECIES+1];
      usrreact(pa, ydot);
      s = -1.0;
      for (i = 0; i < niso; i++)
	s = max(s, ydot[i]*Qi[i]);
      found = true;
    }        
    else if (strncmp(item, "XSO", 3) == 0){
      sscanf(&item[3], "%d", &i);
      if ((i < 1) || (i > niso + 1))
	found = false;
      else{
	float ydot[MAXSPECIES+1];
	bool old_isothermal;
	old_isothermal = isothermal;
	isothermal = false;
	usrreact(pa, ydot);
	isothermal = old_isothermal;
	if (i == niso + 1)
	  /* Energy source per unit mass */
	  s = ydot[i-1]*usrcv(pa);
	else
	  /* Isotope source per unit mass */
	  s = ydot[i-1]*Ai[i-1];
	found = true;
      }
    }
    else if (strcmp(item, "SRATE") ==0){
      float ydot[MAXSPECIES+1], srate;
      bool old_isothermal;
      old_isothermal = isothermal;
      isothermal = false;
      usrreact(pa, ydot);
      isothermal = old_isothermal;
      srate = -1.0;
      /* Isotope rate */      
      for (i = 0; i < niso; i++)
	srate = max(srate, fabs(Ai[i]*ydot[i]));
      /* Energy rate */
      srate = max(srate, fabs(ydot[niso]/pa[iqte]));
      s = srate;
      found = true;
    }
    else if (strcmp(item, "DTRAT") ==0){
      s = cptr->dtrat;
      found = true;
    }      
    else if (strcmp(item, "NSUB") ==0){
      s = cptr->nsub;
      found = true;
    }
    else if (strcmp(item, "BINDEN") == 0){
      int i;
      s = 0.0;
      for (i = 0; i < niso; i++)
	s -= pa[iqal0+1+i]*Qi[i];	
      found = true;
    }
    else if (strcmp(item, "NNEUT") == 0){
      int i;
      s = 0.0;
      for (i = 0; i < niso; i++)
	s += pa[iqal0+1+i]*(Ai[i] - Zi[i]);
      found = true;
    }
    else if (strcmp(item, "NPROT") == 0){
      int i;
      s = 0.0;
      for (i = 0; i < niso; i++)
	s += pa[iqal0+1+i]*Zi[i];
      found = true;
    }        
    else if (strcmp(item, "ABAR") == 0){
      int i;
      s = usrAbar(pa);
      found = true;
    }
    else if (strcmp(item, "ZBAR") == 0){
      int i;
      float Ye, Abar, Zbar, Z2bar;
      usrYAZ(pa, Ye, Abar, Zbar, Z2bar);
      s = Zbar;
      found = true;
    }
    else if (strcmp(item, "DTRAT") ==0){
      s = cptr->dtrat;
      found = true;
    }      
    else if (strcmp(item, "NSUB") ==0){
      s = cptr->nsub;
      found = true;
    }
    else if (strcmp(item, "PHI") == 0){
      if (strncmp(cas, "STEADY", 6) == 0){
	float ss, u, pa[MAXQ], ydot[MAXSPECIES+1];      
	cptr = carr[1];
	usrpq(cptr->qa, pa);
	modsrc(pa, ydot);
	ss = usrss(pa);
	u = pa[iqu0+1];
	s = 1e9*ydot[nal]*(ss*ss - u*u)/u;
	found = true;	
      }
    }
  }
  cptr->mat = ffluid;
  return(found);
}
//end usrps.c

//file usrsaf.c
/*
+ usrsaf - source terms.
     Description :
        Calculates source  terms given the cell  pointer, cptr. norder
        is the order of the step (First 2nd). Calls to gsl routines to
        integrate over interval  dt. Note that cptr->qa is  at the old
        time for  both 1st and  2nd order, pa is  at old time  for 1st
        order, half-time  for 2nd order.   bool second sets  2nd order
        density for gsl integrator.
*/
void usrsaf(int norder, float dt, cell *cptr, float so[])
{
  int i, j;
  float rc, vcof, tv, divv, pp,  dudx[MAXD][MAXD], pa[MAXQ], qao[MAXQ],
    qan[MAXQ], mu[MAXQ], ydot[MAXSPECIES+1], rho, temp;
  bool rzp, zrp, rtp;
  /* Zero so,  copy primitives and old conservatives */
  for (i = 0; i < nq; i++){
    so[i] = 0.0;
    pa[i] = cptr->pa[i];
    qao[i] = cptr->qa[i];    
  }
  /* Geometry */
  rzp = ((strcmp(geom, "RZP") == 0) && (strcmp(coord[0], "R") == 0));
  zrp = ((strcmp(geom, "RZP") == 0) && (strcmp(coord[1], "R") == 0));
  rtp = (strcmp(geom, "RTP") == 0);
  /* Geometric source term (if any) */
  if (rzp){
    /* Cylindrical polars xc[0] = R */
    rc = cptr->xc[0];
    so[iqu0+1] += pa[iqe]/rc;
  }
  else if (zrp){
    /* Cylindrical polars xc[1] = R */
    rc = cptr->xc[1];
    so[iqu0+2] += pa[iqe]/rc;
  }
  else if (rtp){
    /* Spherical polars */
    rc = cptr->vol/(cptr->rjoin[0]->area - cptr->ljoin[0]->area);
    so[iqu0+1] += pa[iqe]/rc;
  }
  if (viscous){
    /* Diffusion coefficients etc */
    usrmu(pa, mu);
    vcof = mu[iqu0+1];
    divv = 0.0;
    /*Cartesian divergence */
    for (i = 0; i < nd; i++){
      divv += cptr->dudx[i][i];
    }
    /* Geometric source term (if any) */
    if (rzp){
      /* Cylindrical polars xc[0] = R */
      divv += pa[iqu0+1]/rc;
      so[iqu0+1] += 2.0*vcof*(divv/3.0 - pa[iqu0+1]/rc)/rc;
    }
    else if (zrp){
      /* Cylindrical polars xc[1] = R */
      divv += pa[iqu0+2]/rc;
      so[iqu0+2] += 2.0*vcof*(divv/3.0 - pa[iqu0+2]/rc)/rc;
    }
    else if (rtp){
      /* Spherical polars */
      divv += 2.0*pa[iqu0+1]/rc;
      so[iqu0+1] += 4.0*vcof*(divv/3.0 - pa[iqu0+1]/rc)/rc;
    }
    if (keps){
      /* Velocity gradients */
      for (i = 0; i < nd ; i++)
	for (j = 0; j < nd; j++)
	  dudx[i][j] = cptr->dudx[i][j];
      /* Production term */
      if (rzp){
	/* Turbulent pressure r-momentum source */
	so[iqu0+1] += 2.0*pa[iqd]*pa[iqk]/(3.0*rc);
	/* Geometric production term */
	tv = pa[iqu0+1]/rc;
	pp = 2.0*tv*tv;
      }
      else if (zrp){
	/* Turbulent pressure r-momentum source */
	so[iqu0+2] += 2.0*pa[iqd]*pa[iqk]/(3.0*rc);
	/* Geometric production term */
	tv = pa[iqu0+2]/rc;
	pp = 2.0*tv*tv;
      }
      else if (rtp){
	/* Turbulent pressure r-momentum source */
	so[iqu0+1] += 4.0*pa[iqd]*pa[iqk]/(3.0*rc);
	/* Geometric production term */
	tv = pa[iqu0+1]/rc;
	pp = 4.0*tv*tv;
      }
      else
	pp = 0.0;
      for (i = 0; i < nd; i++){
	for (j = 0; j < nd; j++)
	  pp += dudx[i][j]*(dudx[i][j] + dudx[j][i]);
      }
      pp = vcof*pp - 2.0*divv*(pa[iqd]*pa[iqk] + divv*vcof)/3.0;
      /* k source */
      so[iqk] += pp - pa[iqd]*pa[iqeps];
      /* Energy source */
      so[iqe] += (pa[iqd]*pa[iqeps] - pp);
      /* Epsilon source */
      so[iqeps] += pa[iqeps]*(ce1*pp - pa[iqd]*ce2*pa[iqeps])/pa[iqk]; 
    }
  }
  /*
     Physical  sources  (per unit  volume).   Note  primitives are  in
     pa. If the inten is used, you also need to add in its sources
  */
  if (react){ 
    if (strncmp(cas, "STEADY", 6) == 0){
      /* Get pa from solution at old timestep */
      usrpq(qao, pa);
      /* Put current solution in cptr->diff */
      for (i = 0; i < nq; i++)
	cptr->diff[i] = pa[i];      
      usrsteady(pa, so, dt);
      return;
    }   
    if (nal > niso){
      /* Reaction progress variable */      
      if (pa[iqd] > 1.1*rho_0){
	tv = pa[iqal0+nal];
	so[iqal0+nal] = rat*pa[iqd]*(tv + 0.01)
	  *(1.0 - tv);	
      }
    }
    if (curv > tiny){
      float qat[MAXQ];
      usrqp(pa, qat);		// get qa based on 2nd order pa
      tv = - curv*qat[iqu0+1];
      so[iqd] += tv;
      so[iqu0+1] += tv*pa[iqu0+1];
      so[iqe] += (pa[iqe] + qat[iqe])*tv/pa[iqd];
      for (i = iqal0+1; i <= iqal0+nal; i++)
	so[i] += tv*pa[i];
    }
    if (strcmp(ode, "NONE") == 0){
      /* Simple solver */
      usrreact(pa, ydot);
      for (i = 0; i < niso; i++)
	so[iqal0+1+i] = pa[iqd]*ydot[i];
      if (!isothermal){      
	so[iqe] = pa[iqd]*usrcv(pa)*ydot[niso];
	so[iqte] = ydot[niso];
      }
    }
    else{
      /* Get pa from solution at old timestep */
      usrpq(qao, pa);
      /* gsl step (update abundances and temperature in pa) */
      usrgsl(cptr, pa, dt);
      /* Get conserved variables after reaction */	
      usrqp(pa, qan);
      /* Nuclear source terms */
      for (i = iqal0+1; i < iqal0+1+niso; i++)
	so[i] += (qan[i] - qao[i])/dt;
      if (!isothermal){
	/* Energy  and temperature */
	so[iqe] += (qan[iqe] - qao[iqe])/dt;
	so[iqte] += (qan[iqte] - qao[iqte])/dt;
      }
    }
    if (inten)
      so[iqie] += (qan[iqie] - qao[iqie])/dt;
  }
  if (strncmp(cas, "ATMOS", 5) == 0){
    tv = par[8]*cptr->pa[iqd];
    so[iqu0+nd] = tv;
    so[iqe] += tv*cptr->pa[iqu0+nd];
  }
  return;
}
//end usrsaf.c

//file usrref.c
/*
+ usrref - reference solution.
     Description:
        Provides a reference solution  for cell cptr. 
*/
void usrref(cell *cptr, float qref[])
{
  int i, id;
  float pa[MAXQ], rho, ss;
  for (i = 0; i < nq; i++)
    qref[i] = mabs(cptr->qa[i]);
  ss = sqrt(cptr->qa[iqe]*cptr->qa[iqd]);
  for (id = 1; id <= nd; id++)
    qref[iqu0+id] = mabs(cptr->qa[iqu0+id]) + ss;
  rho = cptr->qa[iqd];
  if (keps){
    qref[iqk] = 0.01*qref[iqe]/rho;
    qref[iqeps] = qref[iqk];
  }
  for (i = iqal0 + 1; i <= iqal0 + scad; i++)
    qref[i] = rho;
  for (i = iqal0 + scad + 1; i <= iqal0 + nal; i++)
    qref[i] = 1.0;
  if (grav)
    qref[iqpot] = big;
  return;
}
//end usrref.c

//file usrqsol.c
/*
+ usrqsol - Solution quality.
     Description:
        This returns  qsol = badcell if  the solution is  poor, qsol =
        good if  not. This is  determined by comparing  cptr->qa (fine
        cell   value),  with   cptr->pa  (mapped   down   coarse  cell
        value). usrref provides a reference solution.
*/
float usrqsol(int ilev, struct cell *cptr)
{
  int i;
  float qsol, qref[MAXQ];
  usrref(cptr, qref);
  qsol = good;
  if (keps){
    for (i = 0; i < iqk; i++)
      if (mabs(cptr->qa[i] - cptr->pa[i]) > rtol*qref[i])
	qsol = badcell;
    for (i = iqal0+1; i < nq; i++)
      if (mabs(cptr->qa[i] - cptr->pa[i]) > rtol*qref[i])
	qsol = badcell;
  }
  else{
    for (i = 0; i <= iqal0+nal; i++){
      if (mabs(cptr->qa[i] - cptr->pa[i]) > rtol*qref[i])
	qsol = badcell;
    }
  }
  /*
  if ((cptr->xc[1] > par[0]) &&  (ilev > int(par[1])))
    qsol = good;
  if ((cptr->xc[1] > par[2]) &&  (ilev > int(par[3])))
    qsol = good;
  if ((cptr->xc[1] > par[4]) &&  (ilev > int(par[5])))
    qsol = good;
  if ((cptr->xc[1] > par[6]) &&  (ilev > int(par[7])))
    qsol = good;
  */
  return(qsol);
}
//end usrqsol.c

//file usrte.c
/*
+ usrte - Temperature.
*/
float usrte(float pa[])
{
  /* Supernova EOS */
  return(pa[iqte]);
}
//end usrte.c

//file usrmu.c
/*
+ usrmu - Dynamic diffusion coefficients.
     Description :
        Returns dynamic diffusion coefficients.
*/
void usrmu(float pa[], float mu[])
{
  float rho, nu[MAXQ];
  int i;
  /* Zero for safety */
  for (i = 0; i < nq; i++)
    mu[i] = 0.0;
  rho = pa[iqd];
  /* Get kinematic diffsuion coefficients */
  usrnu(pa, nu);
  /* Everything up to advecting scalars multiplied by density */
  for (i = 0; i <= iqal0 + scad; i++)
    mu[i] = rho*nu[i];
  /* Non-advecting scalars as is
  for (i = iqal0 + scad + 1; i <= iqal0 + nal; i++)
    mu[i] = nu[i];*/
  return;
}
//end usrmu.c

//file usrnu.c
/*
+ usrmu - Kinematic diffusion coefficients.
     Description :
        Returns     kinematic     diffusion     coefficients.      SEE
        http://cococubed.asu.edu/code_pages/deflagration.shtml     for
        FLASH coeffs?
*/
void usrnu(float pa[], float nu[])
{
  float tv, tv1, tl, nut, rho, g_eff;
  int i;
  /* Zero for safety */
  for (i = 0; i < nq; i++)
    nu[i] = 0.0;
  if (viscous){
    /* Laminar diffusion coefficients */
    for (i = iqu0 + 1; i <= iqu0 + nd; i++)
      nu[i] = vis;
    nu[iqe] = tcond;
    tv = sdiff;
    for (i = iqal0 + 1; i <= iqal0 + nal; i++)
      nu[i] = tv;
  }
  if (keps){
    /* k-epsilon */
    rho = pa[iqd];
    tv = max(pa[iqk], kmin);
    tv1 = sqrt(tv);
    tl = tv*tv1/max(pa[iqeps],epsmin);
    tl = min(ceps*tl, mtl);
    nut = tv1*tl;
    /* Turbulent viscosity */
    tv = cs*nut;
    for (i = iqu0 + 1; i <= iqu0 + nd; i++)
      nu[i] += tv;
    /* Turbulent thermal conduction */
    nu[iqe] += ct*nut;
    //nu[iqe] += ct*nut;
    /* k diffusion */
    nu[iqk] = tv;
    /* epsilon diffusion */
    nu[iqeps] = ce*nut;
    /* Turbulent scalar diffusion (advecting scalars only) */
    tv = cd*nut;
    for (i = iqal0 + 1; i <= iqal0 + scad; i++)
      nu[i] += tv;
  }
  return;
}
//end usrnu.c

//file usrbd.c
/*
+ usrbd - Boundary conditions.
     Description:
        Sets  the primitive arrays,  pl(left), pr(right)  according to
        the boundary condition specified by the boundary flag, bound.
*/
void usrbd(int id, int bound, join *jptr, float pl[], float pr[])
{
  int i;
  float qa[MAXQ], pa[MAXQ], xc[MAXQ];
  if (bound == lebound){
    /* Left boundary -- default is pl = pr */
    for (i = 0; i < nq; i++)
      pl[i] = pr[i];
    if (bcl[id-1] == bcref)
      /* Symmetry -- reverse normal velocity */
      pl[iqu0+id] = - pr[iqu0+id];
    else if (bcl[id-1] == bcwall)
      /* Wall -- reverse all velocities */
      for (i = 1; i <= nd; i++)
	pl[iqu0+i] = - pr[iqu0+i];
    else if (bcl[id-1] == bcfix){
      /* Fix -- fix solution */
      for (i = 0; i < nd; i++)
	xc[i] = jptr->xj[i];
      //      xc[id-1] = 2.0*xc[id-1] - jptr->right->xc[id-1];
      xc[id-1] = jptr->xj[id-1];
      usrqis(xc, qa);
      usrpq(qa, pa);
      for (i = 0; i < nq; i++)
	pl[i] = pa[i];
    }
    if ((grav) &&(bcl[id-1] != bcref))
    /* Gravitational potential is zero unless symmetry */
      pl[iqpot] = -pr[iqpot];
  }
  else{
    /* Right boundary -- default is pr = pl */
    for (i = 0; i < nq; i++)
      pr[i] = pl[i];
    if (bcr[id-1] == bcref)
      /* Symmetry -- reverse normal velocity */
      pr[iqu0+id] = - pl[iqu0+id];
    else if (bcr[id-1] == bcwall)
      /* Wall -- reverse all velocities */
      for (i = 1; i <= nd; i++)
	pr[iqu0+i] = - pl[iqu0+i];
    else if (bcr[id-1] == bcfix){
      /* Fix -- fix solution */
      for (i = 0; i < nd; i++)
	xc[i] = jptr->xj[i];
      xc[id-1] = 2.0*xc[id-1] - jptr->left->xc[id-1];
      usrqis(xc, qa);
      usrpq(qa, pa);
      for (i = 0; i < nq; i++)
	pr[i] = pa[i];
    }
    if ((grav) &&(bcr[id-1] != bcref))
    /* Gravitational potential is zero unless symmetry */
      pr[iqpot] = -pl[iqpot];
  }
  return;
}
//end usrbd.c

// file usrcmd.c
/*
+ usrcmd - User commands.
    Description :
     This checks to see if it knows the command string in key, and and
     if necessary  performs the associated actions,  returning true to
     stop other command menus from seeing it.
*/
bool usrcmd(char key[])
{
  int i, len;
  float tv;
  char cmd[SZNAM], tname[SZNAM], filename[SZNAM], tail[3];
  bool found, fail;
  struct cell *cptr;
  found = true;
  if (strcmp(key, "PFLOOR") == 0){
    if (!get_r(pfloor)){
      cout << "Error : reading PFLOOR" << '\n';
      report();
    }
  }
  else if (strcmp(key, "DFLOOR") == 0){
    if (!get_r(dfloor)){
      cout << "Error : reading DFLOOR" << '\n';
      report();
    }
  }
  else if (strcmp(key, "AVISP") == 0){
    if (!get_r(avisp)){
      cout << "Error : reading AVISP" << '\n';
      report();
    }
  }
  else if (strcmp(key, "AVISE") == 0){
    if (!get_r(avise)){
      cout << "Error : reading AVISE" << '\n';
      report();
    }
  }
  else if (strcmp(key, "GAMMA") == 0){
    if (!get_r(g)){
      cout << "Error : reading G" << '\n';
      report();
    }
  }
  else if (strncmp(key, "PAR", 3) == 0){
    int i;
    char temp[SZNAM];    
    if (strlen(key) == 4)
      sscanf(&key[3], "%d", &i);
    else{
      temp[0] = key[3];
      temp[1] = key[4];
      temp[2] = '\0';
      sscanf(temp, "%d", &i);
    }
    if ((i < 0) || (i > MAXPAR-1)){
      cout << "Error: reading parameter number " << '\n';
      report();
    }
    else{
      if (get_r(par[i]) == 0){
         cout << "Error : reading par[" << i << "]" << '\n';
         report();
      }
    }
  }
  else if (strncmp(key, "PPAR", 4) == 0){
    int i;
    char temp[SZNAM];
    if (strlen(key) == 5)
      sscanf(&key[4], "%d", &i);
    else{
      temp[0] = key[4];
      temp[1] = key[5];
      temp[2] = '\0';
      sscanf(temp, "%d", &i);
    }
    if (mype == 0){
      if ((i < 0) || (i > MAXPAR-1)){
	cout << "Error: reading parameter number " << '\n';
	report();
      }
      else
	cout << "PAR" << i << " = " << par[i] << '\n';
    }
  }  
  else if (strcmp(key, "REACT") == 0){
    /* Read network name */
    if (!get_n(netname)){
      if (mype == 0)
	cout << "Error : reading network name(.net) \n";
      report();
    }
    else if (!usrinithelm()){
      if(mype == 0)
	cout << "Error: Not using super eos \n";
    }
    else{
      tim7 = false;		
      tim13 = false;
      tim21 = false;
      timtorch = false;
      fail = false;
      len = strlen(netname);
      tail[0] = netname[len-3]; tail[1] = netname[len-2];
      tail[2] = netname[len-1];
      if (strncmp(netname, "tim07", 5) == 0){
	tim7 = true;	
	niso = 7;
	nrat = 16;
	rrscreen = true;      
	revscreen = true;
      }
      else if (strncmp(netname, "tim13", 5) == 0){
	tim13 = true;
	niso = 13;
	nrat = 67;
	if (!init){
	  /* Initialise 13 network */
	  aprox13_init_(aion, zion, bion);
	  init = true;      
	}		
        aprox13_init_(aion, zion, bion);
      }
      else if (strncmp(netname, "tim21", 5) == 0){
	tim21 = true;
	niso = 21;
	nrat = 112;
	if (!init){
	  /* Initialise 21 network */
	  aprox21_init_(aion, zion, bion);
	  init = true;      
	}	
      }
      else if (strncmp(tail, "net", 3) == 0){
	timtorch = true;	
	/* Get network from .net file */	
	if (!init){
	  /* Initialise MESA_TIM network */
	  torch_init_(&niso, &nrat, aion, zion, bion, netname);
	  init = true;      
	}
      }
      else{
	fail = true;
	if (mype ==0)
	  cout << "Error unknown network \n";
      }
      if (!fail){
	cons = true;      
	usrinitnet();      
	react = true;
	if (nom){
	  nd = 1;
	  nal = niso+1;
	  usroff();
	}
      }
    }
  }
  else if (strcmp(key, "NOREACT") == 0){
    init = false;
    react = false;
  }
  else if (strcmp(key, "TESTEOS") == 0){
    float temp, rho, pres, ptot, ein, eintot, sound, Abar, Ye,
      pa[MAXQ], he[12];
    if (mype == 0)
      cout << "Input density \n";
    if (!get_r(rho)){
      if (mype == 0)
	cout << "Error : reading density \n";
      report();
    }
    else{
      if (mype == 0)
	cout << "Input temperature \n";
      if(!get_r(temp)){
	if (mype == 0)
	  cout << "Error : reading temperature \n";
	report();
      }
      else{
	usrpq(carr[1]->qa, pa);
	Abar = usrAbar(pa);
	Ye = usrYE(pa);
	pa[iqd] = rho;
	pa[iqte] = temp;
	usrnadyo(temp, rho, Ye, pres, ein, sound);
	cout << std::scientific;
	cout << std::setprecision(6);
	cout << "Nadyo" << '\n';
	cout << "pressure =  " << pres << " internal energy = " << ein
	     << " sound = " << sound << '\n';
	cout << "Super" << '\n';
	ptot = 0.0;
	eintot  = 0.0;	  
	usrhelm(temp, rho, Ye, he);
	pres = he[0];
	ein = he[1];
	ptot += pres;
	eintot += ein;
	cout << "Helm: pressure =  " << pres << " internal energy = "
	     << ein << '\n';
	usrion(temp, rho, Abar, he);
	pres = he[0];
	ein = he[1];
	ptot += pres;
	eintot += ein;
	cout << "Ion: pressure =  " << pres << " internal energy = "
	     << ein << '\n';	  
	usrrad(temp, rho, he);
	pres = he[0];
	ein = he[1];
	ptot += pres;
	eintot += ein;
	cout << "Rad: pressure =  " << pres << " internal energy = "
	     << ein << '\n';	  	  
	cout << "Total: pressure =  " << ptot << " internal energy = "
	     << eintot << " sound = " << usrss(pa) << '\n';
      }
    }
  }
  else if (strcmp(key, "CONS") == 0)
    cons = true;
  else if (strcmp(key, "NOCONS") == 0)
    cons = false;
  else if (strcmp(key, "CNS") == 0){
    if (!get_r(cns)){
      cout << "Error : reading CNS" << '\n';
      report();
    }
  }
  else if (strcmp(key, "CNS1") == 0){
    if (get_r(cns1) == 0){
      cout << "Error : reading CNS1" << '\n';
      report();
    }
  }
  else if (strcmp(key, "CNS2") == 0){
    if (get_r(cns2) == 0){
      cout << "Error : reading CNS2" << '\n';
      report();
    }
    else if (cns2 > cns1)
      cout << "Warning : Must have CNS2 <= cns1" << '\n';      
  }    
  else if (strcmp(key, "SOLEV") == 0){
    if (get_i(solev) == 0){
      cout << "Error : reading SOLEV" << '\n';
      report();
    }
  }
  else if (strcmp(key, "DTP") == 0){
    if (!get_r(dtp)){
      cout << "Error : reading DTP" << '\n';
      report();
    }
  }
  else if (strcmp(key, "TPOLD") == 0){
    if (!get_r(tp_old)){
      cout << "Error : reading TPOLD" << '\n';
      report();
    }
  }  
  else if (strcmp(key, "GOD") == 0){
    god = true;
    hll = false;
    kt = false;
  }
  else if (strcmp(key, "HLL") == 0){
    god = false;
    hll = true;
    kt = false;
  }
  else if (strcmp(key, "KT") == 0){
    god = false;
    hll = false;
    kt = true;
  }
  else if (strcmp(key, "KEPS") == 0){
    keps = true;
    usrkeps();
  }
  else if (strcmp(key, "KMIN") == 0){
    if (get_r(tv) == 0){
      cout << "Error : reading KMIN" << '\n';
      report();
    }
    else{
      kmin = tv;
      if (!keps)
      cout << "Warning : not turbulent" << '\n';
    }
  }
  else if (strcmp(key, "EPSMIN") == 0){
    if (get_r(tv) == 0){
      cout << "Error : reading EPSMIN" << '\n';
      report();
    }
    else{
      epsmin = tv;
      if (!keps)
      cout << "Warning : not turbulent" << '\n';
    }
  }
  else if (strcmp(key, "MTL") == 0){
    if (get_r(tv) == 0){
      cout << "Error : reading MTL" << '\n';
      report();
    }
    else{
      mtl = tv;
      if (!keps)
      cout << "Warning : not turbulent" << '\n';
    }
  }
  else if (strcmp(key, "NONE") == 0)
    /* no ode solver */
    strcpy(ode, "NONE");  
  else if (strcmp(key, "RK2") == 0)
    /* rk2 ode solver */
    strcpy(ode, "RK2");
  else if (strcmp(key, "RK2IMP") == 0)
    /* rk2imp ode solver */
    strcpy(ode, "RK2IMP");
  else if (strcmp(key, "MSBDF") == 0)
    /* msbdf ode solver */
    strcpy(ode, "MSBDF");
  else if (strcmp(key, "NSUB") == 0){
    if (!get_i(nsub)){    
      if (mype == 0)
	cout << "Error : reading NSUB (number of substeps) \n";
      report();
    }
  }  
  else if (strcmp(key, "SCREEN") == 0)
    rrscreen = true;
  else if (strcmp(key, "NOSCREEN") == 0)
    rrscreen = false;
  else if (strcmp(key, "REVSCREEN") == 0)
    revscreen = true;
  else if (strcmp(key, "NOREVSCREEN") == 0)
    revscreen = false;    
  else if (strcmp(key, "DONEUT") == 0)
    doneut = 1;
  else if (strcmp(key, "NODONEUT") == 0)
    doneut = 0;
  else if (strcmp(key, "DOWEAK") == 0)
    doweak = 1;
  else if (strcmp(key, "NODOWEAK") == 0)
    doweak = 0;    
  else if (strcmp(key, "HYDROSTATIC") == 0)
    /* Hydrostatic (constant density, temperature) */
    hydrostatic = true;
  else if (strcmp(key, "DYNAMIC") == 0)
    /* Not hydrostatic */
    hydrostatic = false;
  else if (strcmp(key, "ISOTHERMAL") == 0)
    /* Isothermal */
    isothermal = true;
  else if (strcmp(key, "NOISOTHERMAL") == 0)
    /* Not isothermal */    
    isothermal = false;    
  else if (strcmp(key, "PRINT") == 0){
    /* Print pressure, internal energy, etc */
    int i, ilev;
    float temp, rho, opg, otemp, orho, Abar, Ye, xc[MAXD], pa[MAXQ],
      ydot[MAXSPECIES+1], X[MAXSPECIES], pg, ein, pos, ss, binden;
    bool old_isothermal;      
    struct cell *cptr;   
    if (mype == 0)
      cout << "Input density, temperature and position\n";
    if (!get_r(rho)){
      if (mype == 0)
	cout << "Error : reading density \n";
      report();
    }
    else if (!get_r(temp)){
      if (mype == 0)
	cout << "Error : reading temperature \n";
      report();
    }
    else if (!get_r(pos)){
      if (mype == 0)
	cout << "Error : reading position \n";
      report();
    }      
    else{
      for (i = 0; i < nd; i++)
	xc[i] = 0.0;
      xc[0] = pos;
      cptr = find_fcell(ilev, xc);
      if (cptr == NULL){
	if (mype == 0)
	  cout << "Error : no such cell \n";
	report();
      }
      else{
	usrpq(cptr->qa, pa);
	pg = usrpg(temp, rho, usrAbar(pa), usrYE(pa));
	ein = usrein(temp, rho, usrAbar(pa), usrYE(pa)); 
	cout << std::scientific;
	cout << std::setprecision(4);
	cout << "pressure " << pg << " internal energy " << ein <<
	  '\n';
	orho = pa[iqd];
	opg = pa[iqe];
	otemp = pa[iqte];	  
	pa[iqd] = rho;
	pa[iqe] = pg;
	pa[iqte] = temp;
	ss = usrss(pa);
	cout << "frozen sound speed " << ss <<
	  " frozen Mach " << mabs(pa[iqu0+1])/ss << '\n';
	ss = usrssit(pa);
	cout << "isothermal sound speed " << ss <<
	  " isothermal Mach " << mabs(pa[iqu0+1])/ss << '\n';	  
	old_isothermal = isothermal;
	isothermal = false;
	usrreact(pa, ydot);
	isothermal = old_isothermal;
	binden = 0.0;
	for (i = 0; i < niso; i++)
	  binden += Qi[i]*pa[iqal0+1+i];
	cout << "energy source (unit mass) " <<
	  ydot[niso]*usrcv(pa) << " reaction timescale "
	     << ein/(usrcv(pa)*fabs(ydot[niso])) << " binden " <<
	  binden << '\n';
      }
    }
  }
  else if (strcmp(key, "PFCELL") == 0){
    int id, i, ilev;
    bool fail = false;
    float xc[MAXD], pa[MAXQ];
    struct cell *cptr;
    for (id = 1; id <= nd; id++){
      if (get_r(xc[id-1]) == 0){
	cout << "Error : reading " << id << " coordinate \n";
	report();
	fail = true;
	break;
      }
    }
    if (!fail){
      cptr = find_fcell(ilev, xc);
      if (cptr == NULL){
	if(mype == 0)
	  cout << "Error : no such cell \n";
	report();
      }
      else{
	usrpfcell(cptr);
      }
    }
  }
  else if (strcmp(key, "SHOCK") == 0){
    int i, ilev;
    float rho_r, T_r, v_r, rho_l, T_l, v_l, p_l, xc[MAXD], pa[MAXQ],
      pos;
      
    if (mype == 0)
      cout << "Input upstream density, temperature, velocity" <<
	" and position\n";
    if (!get_r(rho_r)){
      if (mype == 0)
	cout << "Error : reading rho_r" << '\n';
      report();
    }
    else if (!get_r(T_r)){
      if (mype == 0)
	cout << "Error : reading T_r" << '\n';
      report();
    }
    else if (!get_r(v_r)){
      if (mype == 0)
	cout << "Error : reading v_r" << '\n';
      report();
    }
    else if (!get_r(pos)){
      if (mype == 0)
	cout << "Error : reading position" << '\n';
      report();
    }      
    else{
      for (i = 0; i < nd; i++)
	xc[i] = 0.0;
      xc[0] = pos;
      cptr = find_fcell(ilev, xc);
      if (cptr == NULL){
	if (mype == 0)
	  cout << "Error : no such cell \n";
	report();
      }
      else{
	usrpq(cptr->qa, pa);
	usrshock(pa, rho_r, T_r, v_r, rho_l, T_l, v_l, p_l);
	cout << std::scientific;
	cout << std::setprecision(6);
	cout << "rho_l = " << rho_l << " T_l = " << T_l << " v_l = "
	     << v_l << " p_l = " << p_l << '\n';
      }
    }
  }
  else if (strcmp(key, "SETSTEADY") == 0){
    int i, ilev, icell;
    float rho_r, T_r, v_r, rho_l, T_l, v_l, p_l, xc[MAXD], pa[MAXQ],
      pos;
    if (mype == 0)
      cout << "Input upstream density, temperature, velocity" <<
	" and position\n";
    if (!get_r(rho_r)){
      if (mype == 0)
	cout << "Error : reading rho_r" << '\n';
      report();
    }
    else if (!get_r(T_r)){
      if (mype == 0)
	cout << "Error : reading T_r" << '\n';
      report();
    }
    else if (!get_r(v_r)){
      if (mype == 0)
	cout << "Error : reading v_r" << '\n';
      report();
    }
    else if (!get_r(pos)){
      if (mype == 0)
	cout << "Error : reading position" << '\n';
      report();
    }      
    else{
      for (i = 0; i < nd; i++)
	xc[i] = 0.0;
      xc[0] = pos;
      cptr = find_fcell(ilev, xc);
      if (cptr == NULL){
	if (mype == 0)
	  cout << "Error : no such cell \n";
	report();
      }
      else{
	usrpq(cptr->qa, pa);
	usrshock(pa, rho_r, T_r, v_r, rho_l, T_l, v_l, p_l);
	/* Impose post shock solution everywhere */
	pa[iqd] = rho_l;
	pa[iqu0+1] = v_l;
	pa[iqte] = T_l;
	pa[iqe] = p_l;
	/* Set  */
	D = v_r;
	/* Distance from shock (Note nal must be no of species +1) */
	pa[iqal0+nal] = 0.0;
	for (icell = 1; icell <= celltop; icell++){
	  cptr = carr[icell];
	  if (!cptr->halo){
	    for (i = 0; i < nq; i++)
	      cptr->pa[i] = pa[i];
	    usrqp(cptr->pa, cptr->qa);
	  }
	}
	hydrostatic = true;
      }
    }
  }
  else if (strcmp(key, "DOSTEADY") == 0){
    int mon;
    char verbose[SZNAM], bisect[SZNAM];
    if (nom){
      if (mype == 0)
	cout << "Error : no model" << '\n';
      report();
    }
    else{
      if (mype == 0)
	cout << "Input lower, upper, initial D, verbose and bisect \n" ;
      if (!get_r(D_l)){
	if (mype == 0)
	  cout << "Error : reading D_l" << '\n';
	report();
      }
      else if (!get_r(D_r)){
	if (mype == 0)
	  cout << "Error : reading D_r" << '\n';
	report();
      }
      else if (!get_r(D)){
	if (mype == 0)
	  cout << "Error : reading D" << '\n';
	report();
      }
      else if (!get_n(verbose)){      
	if (mype == 0)
	  cout << "Error : reading verbose" << '\n';
	report();
      }
      else if (!get_n(bisect)){      
	if (mype == 0)
	  cout << "Error : reading bisect" << '\n';
	report();
      }            
      else{
	cptr = carr[1];
	if(!modsol(cptr, strcmp(verbose, "verbose") == 0,
		   strcmp(bisect, "bisect") == 0)){
	  if (mype == 0)
	    cout << "Error : modsol failed" << '\n';
	  report();
	}     
      }
    }
  }
  else if (strcmp(key, "ATMOS") == 0){
    int i, ilev;
    float rho_0, rho_min, T_0, M_0, R_0, M_shell, pa[MAXQ], xc[MAXD];
    struct cell *cptr;
    cout << "base density, min density, temperature, stellar mass, ";
    cout << "radius, shell mass, adiabatic/isothermal, \n";
    cout << "output filename \n";
    if (!get_r(rho_0)){
      if (mype == 0)
	cout << "Error : reading density" << '\n';
      report();
    }
    else if (!get_r(rho_min)){
      if (mype == 0)
	cout << "Error : reading min density" << '\n';
      report();
    }    
    else if (!get_r(T_0)){
      if (mype == 0)
	cout << "Error : reading temperature" << '\n';
      report();
    }
    else if (!get_r(M_0)){
      if (mype == 0)
	cout << "Error : reading stellar mass" << '\n';
      report();
    }
    else if (!get_r(R_0)){
      if (mype == 0)
	cout << "Error : reading stellar radius" << '\n';
      report();
    }
    else if (!get_r(M_shell)){
      if (mype == 0)
	cout << "Error : reading shell mass" << '\n';
      report();
    }
    else if (!get_cmd(cmd)){
      if (mype == 0)
	cout << "Error : reading ADIABATIC/ISOTHERMAL" << '\n';
      report();
    }
    else if (!get_n(monpref)){
      if (mype == 0)
	cout << "Error : reading filename" << '\n';
      report();
    }    
    else {
      strcpy(monmod, "POS");            
      usropenxq(monmod, 3);     
      for (i = 0; i < nd; i++)
	xc[i] = 0.0;
      cptr = find_fcell(ilev, xc);
      usrpq(cptr->qa, pa);
      if (!usratmos(rho_0, rho_min, T_0, M_0, R_0, M_shell, cmd, pa)){
	if (mype == 0)
	  cout << "Error : usratmos failed" << '\n';
	report();
      }
      else
	usrclosexq(monmod, 3);
    }
  }
  else if (strcmp(key, "CORE") == 0){
    int i, ilev;
    float rho_0, rho_min, T_0, pa[MAXQ], xc[MAXD];
    struct cell *cptr;
    cout << "central, min density, central temperature, ";
    cout << "adiabatic/isothermal output filename \n";
    if (!get_r(rho_0)){
      if (mype == 0)
	cout << "Error : reading density" << '\n';
      report();
    }
    else if (!get_r(rho_min)){
      if (mype == 0)
	cout << "Error : reading min density" << '\n';
      report();
    }    
    else if (!get_r(T_0)){
      if (mype == 0)
	cout << "Error : reading temperature" << '\n';
      report();
    }
    else if (!get_cmd(cmd)){
      if (mype == 0)
	cout << "Error : reading ADIABATIC/ISOTHERMAL" << '\n';
      report();
    }
    else if (!get_n(monpref)){
      if (mype == 0)
	cout << "Error : reading filename" << '\n';
      report();
    }    
    else {
      strcpy(monmod, "POS");            
      usropenxq(monmod, 3);     
      for (i = 0; i < nd; i++)
	xc[i] = 0.0;
      cptr = find_fcell(ilev, xc);
      usrpq(cptr->qa, pa);
      if (!usrcore(rho_0, rho_min, T_0, cmd, pa)){
	if (mype == 0)
	  cout << "Error : usrcore failed" << '\n';
	report();
      }
      else
	usrclosexq(monmod, 3);
    }
  }  
  else if (strcmp(key, "ATMINIT") == 0){
    cout << "Input filename \n";
    if (!get_n(filename)){
      if (mype == 0)
	cout << "Error : reading atmosphere filename \n";
      report();
    }
    else
      usrinitintp(filename);
  }
  else if (strcmp(key, "TESTATMINIT") == 0){
    float x;
    cout << "Input x \n";
    if (!get_r(x)){
      if (mype == 0)
	cout << "Error : reading atmosphere x \n";
      report();
    }
    else{
      float fraci, value;
      int i;
      intp(x, 0, intxmax, intpts, intdx, i, fraci);
      value = intrho[i] + fraci*(intrho[i+1] - intrho[i]);
      cout << "value " << value << '\n';
    }
  }
  else if (strcmp(key, "SHOCKPOS")== 0){
    int dir;
    float tv, nx, qa[MAXQ], pa[MAXQ];
    fail = false;
    nx =0.0;
    if (!get_i(dir)){
      if (mype == 0)
	cout << "Error : reading id \n";
      report();
    }      
    else if ((dir < 1) || (dir > nd)){
      if (mype == 0)
	cout << "Error : dir is not in 1..nd \n";
      report();
      fail = true;
    }
    if (!fail){
      if (!get_r(nx)){
	if (mype == 0)
	  cout << "Error : reading other direction x \n";
	report();
      }
      else
	tv = usrshockpos(true, dir, nx, qa, pa);
    }
  }
  else if (strcmp(key, "SONPOS") == 0){
    int dir, nx;
    float xs, pos, width, detv, qa[MAXQ], pa[MAXQ];
    struct cell *cptr;
    char name[SZNAM];
    xs = usrshockpos(false, 1, 0, qa, pa);
    detv = pa[iqu0+1]/(1.0 - qa[iqd]/pa[iqd]);
    pos = usrsonpos(xs, detv);
    if (mype == 0)
      cout << "Sonic pos " << pos << " length " << xs - pos <<
	" D " << detv <<'\n'; 
  }
  else if (strcmp(key, "SMOOTH") == 0){
    int nav;
    char name[SZNAM], filein[SZNAM], fileout[SZNAM];
    if (limpe >1){
      if (mype == 0)
	cout << "Error : only works for single processor \n";
      report();
    }
    cout << "Input nav and .xq file \n";    
    if (!get_i(nav)){    
      if (mype == 0)
	cout << "Error : reading NAV (number of averaging pts \n";
      report();
    }
    else if (!get_n(name)){
      if (mype == 0)
	cout << "Error : reading input filename" << '\n';
      report();
    }
    else{
      strcpy(filein, name);
      strcat(filein, ".xq");
      strcpy(fileout, name);
      strcat(fileout, "_smooth.xq");            
      if (!(shockin=fopen(filein, "r"))){
	if (mype ==0)
	  cout << "Error : opening input file \n";
	report();
      }
      else if (!(shockout=fopen(fileout, "w"))){
	if (mype ==0)
	  cout << "Error : opening output file \n";
	report();
      }      
      else{
	fprintf(shockout, "%d %d \n", 1, 10);
	fprintf(shockout, "D\n");
	fprintf(shockout, "t D \n");
	usrsmooth(nav);
      }
    }
  }
  else if (strcmp(key, "PHASEPLAN") == 0){
    int nav;
    char name1[SZNAM], name2[SZNAM], filein1[SZNAM], filein2[SZNAM],
      fileout[SZNAM];
    if (limpe >1){
      if (mype == 0)
	cout << "Error : only works for single processor \n";
      report();
    }
    cout << "Input abscissa and ordinate .xq files \n";    
    if (!get_n(name1)){
      if (mype == 0)
	cout << "Error : reading abscissa filename" << '\n';
      report();
    }
    else if (!get_n(name2)){
      if (mype == 0)
	cout << "Error : reading ordinate filename" << '\n';
      report();
    }
    else{
      strcpy(filein1, name1);
      strcat(filein1, ".xq");
      strcpy(filein2, name2);
      strcat(filein2, ".xq");      
      strcpy(fileout, name1);
      strcat(fileout, "_phase.xq");            
      if (!(shockin=fopen(filein1, "r"))){
	if (mype ==0)
	  cout << "Error : opening first input file \n";
	report();
      }
      else if (!(shockout=fopen(filein2, "r"))){
	if (mype ==0)
	  cout << "Error : opening second file \n";
	report();
      }
      else if (!(convin=fopen(fileout, "w"))){
	if (mype ==0)
	  cout << "Error : opening output file file \n";
	report();
      }            
      else{
	fprintf(convin, "%d %d \n", 1, 10);
	fprintf(convin, "phaseplan\n");
	fprintf(convin, "var1 var2 \n");
	usrphaseplan();
      }
    }
  }  
  else if (strcmp(key, "SONCROSS") == 0){
    float xmax;
    char name[SZNAM], filein[SZNAM];
    if (limpe >1){
      if (mype == 0)
	cout << "Error : only works for single processor \n";
      report();
    }
    cout << "Input mach_x.xq file  and xmax \n";    
    if (!get_n(name)){
      if (mype == 0)
	cout << "Error : reading input filename" << '\n';
      report();
    }
    else if (!get_r(xmax)){
      if (mype == 0)
	cout << "Error : reading max x" << '\n';
      report();
    }    
    else{
      strcpy(filein, name);
      strcat(filein, "_mach_x.xq");
      if (!(shockin=fopen(filein, "r"))){
	if (mype ==0)
	  cout << "Error : opening input mach file \n";
	report();
      }
      else{
	strcpy(filein, name);
	strcat(filein, "_u_x.xq");
	shockout=fopen(filein, "r");
	usrsoncross(xmax);
      }
    }
  }
  else if (strcmp(key, "MAXXA") == 0){
    char name[SZNAM], filein[SZNAM];
    if (limpe >1){
      if (mype == 0)
	cout << "Error : only works for single processor \n";
      report();
    }
    cout << "Input iso .xq file \n";    
    if (!get_n(name)){
      if (mype == 0)
	cout << "Error : reading input filename" << '\n';
      report();
    }
    else{
      strcpy(filein, name);
      strcat(filein, "_x.xq");
      if (!(shockin=fopen(filein, "r"))){
	if (mype ==0)
	  cout << "Error : opening input file \n";
	report();
      }
      else{
	usrmaxab();
      }
    }
  }  
  else if (strcmp(key, "MALDATCONV") == 0){
    if (get_n(monpref) == 0){
      if (mype == 0)
	cout << "Error : reading infile prefix" << '\n';
      report();
    }
    else{
      mon = true;
      strcpy(monmod, "POS");      
      usropenxq(monmod, 1);
      if (!usrmaldatconv()){
	if (mype == 0)
	  cout << "Error : reading opening react/prim file \n";
      }
    }    
  }  
  else if (strcmp(key, "MALCONV") == 0){
    char name[SZNAM], filein[SZNAM], fileout[SZNAM];
    if (limpe >1){
      if (mype == 0)
	cout << "Error : only works for single processor \n";
      report();
    }
    cout << "Input .dat file name \n";
    if (!get_n(name)){    
      if (mype == 0)
	cout << "Error : reading .dat file \n";
      report();
    }
    else{
      strcpy(filein, name);
      strcat(filein, ".dat");
      if (!(convin=fopen(filein, "r"))){
	if (mype ==0)
	  cout << "Error : opening input file \n";
	report();
      }
      else if (!(dkout=fopen("dkout.xq", "w"))){
	if (mype ==0)
	  cout << "Error : opening dk output file \n";
	report();
      }
      else if (!(dlcout=fopen("dlcout.xq", "w"))){
	if (mype ==0)
	  cout << "Error : opening dlc output file \n";
	report();
      }            
      else{
	fprintf(dkout, "%d %d \n", 1, 10);
	fprintf(dkout, "D \n");
	fprintf(dkout, "K D \n");
	fprintf(dlcout, "%d %d \n", 1, 10);
	fprintf(dlcout, "l_c \n");
	fprintf(dlcout, "K l_c \n");	
	usrmalconv();
      }
    }
  }
  else if (strcmp(key, "TIMCONV") == 0){
    char name[SZNAM], filein[SZNAM], fileout[SZNAM];
    if (limpe >1){
      if (mype == 0)
	cout << "Error : only works for single processor \n";
      report();
    }
    cout << "Input .dat file name \n";
    if (!get_n(name)){    
      if (mype == 0)
	cout << "Error : reading .dat file \n";
      report();
    }
    else{
      strcpy(filein, name);
      strcat(filein, ".dat");
      if (!(convin=fopen(filein, "r"))){
	if (mype ==0)
	  cout << "Error : opening input file \n";
	report();
      }
      else if (!(dkout=fopen("se.xq", "w"))){
	if (mype ==0)
	  cout << "Error : opening se output file \n";
	report();
      }
      else{
	fprintf(dkout, "%d %d \n", 1, 10);
	fprintf(dkout, "se \n");
	fprintf(dkout, "t se \n");
	usrtimconv();
      }
    }
  }
  else if (strcmp(key, "RIEMAN") == 0){
    int i, ilev;
    float rho_r, u_r, p_r, rho_l, u_l, p_l, pos, xc[MAXD], pa[MAXQ],
      pvl[MAXQ], pvr[MAXQ];
    if (mype == 0)
      cout << "Input left state: rho, u, p \n";
    if (!get_r(rho_l)){
      if (mype == 0)
	cout << "Error : reading rho_l \n";
      report();
    }
    else if (!get_r(u_l)){
      if (mype == 0)
	cout << "Error : reading u_l \n";
      report();
    }
    else if (!get_r(p_l)){
      if (mype == 0)
	cout << "Error : reading p_l \n";
      report();
    }
    if (mype == 0)
      cout << "Input right state: rho, u, p \n";
    if (!get_r(rho_r)){
      if (mype == 0)
	cout << "Error : reading rho_r \n";
      report();
    }
    else if (!get_r(u_r)){
      if (mype == 0)
	cout << "Error : reading u_r \n";
      report();
    }
    else if (!get_r(p_r)){
      if (mype == 0)
	cout << "Error : reading p_r \n";
      report();
    }
    else{
      if (mype == 0)
	cout << "input position \n";
      if (!get_r(pos)){
	if (mype == 0)
	  cout << "Error : reading position \n";	  
      }
      for (i = 0; i < nd; i++)
	xc[i] = 0.0;
      xc[0] = pos;
      cptr = find_fcell(ilev, xc);
      if (cptr == NULL){
	if (mype == 0)
	  cout << "Error : no such cell \n";
	report();
      }
      else{
	usrpq(cptr->qa, pa);
	for (i = iqal0+1; i < nq; i++){
	  pvr[i] = pa[i];
	  pvl[i] = pa[i];
	}
	pvr[iqd] = rho_r;
	pvr[iqu0+1] = u_r;
	pvr[iqe] = p_r;
	pvl[iqd] = rho_l;
	pvl[iqu0+1] = u_l;
	pvl[iqe] = p_l;
	usrrieman(pvl, pvr, pa);
	cout << std::scientific;
	cout << std::setprecision(6);
	cout << "Resolved state \n";
	cout << "rho = " << pa[iqd] << " u = " << pa[iqu0+1] <<
	  " p = " << pa[iqe] << '\n';
      }
    }
  }
  else if (strcmp(key, "RAT") == 0){
    if (get_r(rat) == 0){
      cout << "Error : reading RAT" << '\n';
      report();
    }
  }
  else if (strcmp(key, "CURV") == 0){
    if (get_r(curv) == 0){
      cout << "Error : reading CURV" << '\n';
      report();
    }
  }
  else if (strcmp(key, "DEPTH") == 0){
    if (get_r(depth) == 0){
      cout << "Error : reading DEPTH" << '\n';
      report();
    }
  }  
  else if (strcmp(key, "TEQUI") == 0){
    if (get_r(T_e) == 0){
      if (mype == 0)
	cout << "Error : reading TEQUI" << '\n';
      report();
    }
  }
  else if (strcmp(key, "RHO_0") == 0){
    if (get_r(rho_0) == 0){
      if (mype == 0)
	cout << "Error : reading rho_0" << '\n';
      report();
    }
  }    
  else if (strcmp(key, "T_0") == 0){
    if (get_r(T_0) == 0){
      if (mype == 0)
	cout << "Error : reading T_0" << '\n';
      report();
    }
  }  
  else if (strcmp(key, "RHO_S") == 0){
    if (get_r(rho_s) == 0){
      if (mype == 0)
	cout << "Error : reading rho_s" << '\n';
      report();
    }
  }
  else if (strcmp(key, "T_S") == 0){
    if (get_r(T_s) == 0){
      if (mype == 0)
	cout << "Error : reading T_s" << '\n';
      report();
    }
  }  
  else if (strcmp(key, "SETMASSFRAC") == 0){
    int i, icell;
    float XN, yn, xl, xr;
    struct cell *cptr;
    if (mype == 0)
      cout << "Input isotope number (1..7), mass fraction xl, xr \n";    
    if (get_i(i) == 0){
      if (mype == 0)      
	cout << "Error : reading i \n";
      report();
    }
    else if ((i < 1) || (i > niso)){
      if (mype == 0)      
	cout << "Error : i out of range" << '\n';
      report();
    }        
    else if (get_r(XN) == 0){
      if (mype == 0)      
	cout << "Error : reading XN \n";
      report();
    }
    else if (get_r(xl) == 0){
      if (mype == 0)      
	cout << "Error : reading left limit \n";
      report();
    }
    else if (get_r(xr) == 0){
      if (mype == 0)      
	cout << "Error : reading right limit \n";
      report();
    }
    else{
      yn = XN/Ai[i-1];
      for (icell = 1; icell <= celltop; icell++){
	cptr = carr[icell];
	if ((!cptr->halo) && (cptr->xc[0] > xl) && (cptr->xc[0] < xr)){
	  usrpq(cptr->qa, cptr->pa);
	  cptr->pa[iqal0+i] = yn;
	  usrqp(cptr->pa, cptr->qa);
	}
      }
    }
  }
  else if (strcmp(key, "SETRPROG") == 0){
    int icell;
    float xl, xr, al;
    struct cell *cptr;
    if (mype == 0)
      cout << "Input al, xl, xr \n";
    if (get_r(al) == 0){
      if (mype == 0)      
	cout << "Error : reading AL \n";
      report();
    }    
    if (get_r(xl) == 0){
      if (mype == 0)      
	cout << "Error : reading left limit \n";
      report();
    }
    else if (get_r(xr) == 0){
      if (mype == 0)      
	cout << "Error : reading right limit \n";
      report();
    }
    else{
      for (icell = 1; icell <= celltop; icell++){
	cptr = carr[icell];
	if ((!cptr->halo) && (cptr->xc[0] > xl) && (cptr->xc[0] < xr)){
	  usrpq(cptr->qa, cptr->pa);
	  cptr->pa[iqal0+nal] = al;
	  usrqp(cptr->pa, cptr->qa);
	}
      }
    }
  }
  else if (strcmp(key, "SETPRIMV") == 0){
    int i=-1, icell;
    float tv, xl, xr;
    struct cell *cptr;
    if (mype == 0)
      cout << "Input \'rho\' or \'ux\' or \'temp\' value, xl, xr \n";
    if (!get_cmd(cmd)){    
      if (mype == 0)      
	cout << "Error : reading variable \n";
      report();
    }
    else if (strcmp(cmd, "RHO") == 0)
      i = iqd;
    else if (strcmp(cmd, "UX") == 0)
      i = iqu0 + 1;
    else if (strcmp(cmd, "TEMP") == 0)
      i = iqte;
    else{
      if (mype == 0)      
	cout << "Error : unknown variable \n";
      report();
    }
    if (i >=0){
      if (get_r(tv) == 0){
	if (mype == 0)      
	  cout << "Error : reading value \n";
	report();
      }
      else if (get_r(xl) == 0){
	if (mype == 0)      
	  cout << "Error : reading left limit \n";
	report();
      }
      else if (get_r(xr) == 0){
	if (mype == 0)      
	  cout << "Error : reading right limit \n";
	report();
      }    
      else{
	for (icell = 1; icell <= celltop; icell++){
	  cptr = carr[icell];
	  if ((!cptr->halo) && (cptr->xc[0] > xl)
	      && (cptr->xc[0] < xr)){
	    usrpq(cptr->qa, cptr->pa);
	    cptr->pa[i] = tv;
	    usrqp(cptr->pa, cptr->qa);
	  }
	}
      }
    }
  }
  else if (strncmp(key, "CLEANS", 6) == 0){
    struct cell *cptr;
    int dir,  icell, side;
    float xc[MAXD];
    fail = false;
    sscanf(&key[6], "%s", tname);
    if (strncmp(tname, "L", 1) == 0)
      side = -1;
    else if (strncmp(tname, "R", 1) == 0)
      side = 1;
    else{
      if (mype == 0)
	cout << "Error : Unknown side \n";
      report();
      fail = true;
    }
    if (!fail){
      if (mype == 0)
	cout << "Input cleans coord id and position \n";
      if (!get_i(dir)){ 
	if (mype == 0)
	  cout << "Error : reading coord id \n";
	fail = true;
	report();
      }
    }
    if (!fail)
      fail =(dir < 1) || (dir > nd);
    if (fail){  
      if (mype == 0)
	cout << "Error : dir is not in 1..nd \n";	  
      report();
    }
    if (!fail){
      if (get_r(tv) == 0){
	if (mype == 0)
	  cout << "Error : reading clean position" << '\n';
	report();
      }
      else{
	for (i = 0; i < nd; i++)
	  xc[i] = 0.0;
	xc[dir-1] = xl[dir-1] + 0.1*lev[dir-1]->dxbar[dir-1]; 
	for (icell = 1; icell <= celltop; icell++){
	  cptr = carr[icell];
	  if (!cptr->halo){
	    if (side*(cptr->xc[dir-1] - tv) > 0)
	      usrqis(cptr->xc, cptr->qa);
	  }
	}
      }
    }
  }  
  else if (strcmp(key, "VSHIFT") == 0){
    if (mype == 0)
      cout << "Input velocity shift \n";
    if (get_r(tv) == 0){
      if (mype == 0)
	cout << "Error : reading v \n";
      report();
    }
    else{
      int icell;
      float pt[MAXQ];
      for (icell = 1; icell <= celltop; icell++){
	cptr = carr[icell];	
	if (!cptr->halo){	
	  usrpq(cptr->qa, pt);
	  pt[iqu0+1] += tv;
	  usrqp(pt, cptr->qa);
	}
      }
    }
  }  
  else if (strcmp(key, "MON") == 0){
    if (get_n(monpref) == 0){
      if (mype == 0)
	cout << "Error : reading monfile prefix" << '\n';
      report();
    }
    else{
      mon = true;
      strcpy(monmod, "POS");      
      usropenxq(monmod, 1);
    }
  }
  else if (strcmp(key, "NOMON") == 0){
    if (mon){
      mon = false;
      strcpy(monmod, "POS");      
      usrclosexq(monmod, 1);
    }
    else
      cout << "no monitor\n";
  }
  else if (strcmp(key, "MONV") == 0){
    if (get_n(monpref) == 0){
      if (mype == 0)
	cout << "Error : reading monfile prefix" << '\n';
      report();
    }
    else{
      mon = true;
      strcpy(monmod, "VEL");      
      usropenxq(monmod, 1);
    }
  }
  else if (strcmp(key, "NOMONV") == 0){
    if (mon){
      mon = false;
      strcpy(monmod, "VEL");
      usrclosexq(monmod, 1);
    }
    else
      cout << "no monitor\n";
  }
  else if (strcmp(key, "MONP") == 0){
    mon = true;
    strcpy(monmod, "PART");      
    usropenxq(monmod, 1);
  }
  else if (strcmp(key, "NOMONP") == 0){
    if (mon){
      mon = false;
      strcpy(monmod, "PART");
      usrclosexq(monmod, 1);
    }
    else
      cout << "no monitor\n";
  }
  else if (strcmp(key, "MONE") == 0){
    mon = true;
    strcpy(monmod, "ENGEN");      
    usropenxq(monmod, 1);
  }
  else if (strcmp(key, "NOMONE") == 0){
    if (mon){
      mon = false;
      strcpy(monmod, "ENGEN");
      usrclosexq(monmod, 1);
    }
    else
      cout << "no monitor\n";
  }
  else if (strcmp(key, "SETINT") == 0){
    int mode;
    float size;
    if (mype == 0)
      cout << "Input filename, mode and size \n";
    if (get_n(monpref) == 0){
      if (mype == 0)
	cout << "Error : reading interpolation file  prefix" << '\n';
      report();
    }
    else if (get_i(mode) == 0){
      if (mype == 0)
	cout << "Error : reading mode " << '\n';
      report();
    }          
    else if (get_r(size) == 0){
      if (mype == 0)
	cout << "Error : reading size " << '\n';
      report();
    }    
    else if (!usrsetint(monpref, mode, size)){
      if (mype == 0)
	cout << "Error : reading opening interpolation file" << '\n';
      report();
    }
  }
  else if (strcmp(key, "GSLATOL") == 0){
    if (get_r(gslatol) == 0){
      cout << "Error : reading GSLATOL" << '\n';
      report();
    }
  }
  else if (strcmp(key, "GSLRTOL") == 0){
    if (get_r(gslrtol) == 0){
      cout << "Error : reading GSLRTOL" << '\n';
      report();
    }
  }
  else if (strcmp(key, "GSLHMIN") == 0){
    if (get_r(gslhmin) == 0){
      cout << "Error : reading GSLHMIN" << '\n';
      report();
    }
  }
  else if (strcmp(key, "GSLHMAX") == 0){
    if (get_r(gslhmax) == 0){
      cout << "Error : reading GSLHMAX" << '\n';
      report();
    }
  }
  else if (strcmp(key, "FXP") == 0)
    fxp = true;
  else if (strcmp(key, "NOFXP") == 0)
    fxp = false;
  else if (strcmp(key, "VISCOUS") == 0)
    viscous = true;
  else if (strcmp(key, "NOVISCOUS") == 0)
    viscous = false;
  else if (strcmp(key, "INTEN") == 0)
    inten = true;
  else if (strcmp(key, "NOINTEN") == 0)
    inten = false;
  else if (strcmp(key, "TCOND") == 0){
    if (get_r(tcond) == 0){
      cout << "Error : reading TCOND" << '\n';
      report();
    }
  }
  else if (strcmp(key, "SDIFF") == 0){
    if (get_r(sdiff) == 0){
      cout << "Error : reading SDIFF" << '\n';
      report();
    }
  } 
  else if (strcmp(key, "VIS") == 0){
    if (get_r(vis) == 0){
      cout << "Error : reading VIS" << '\n';
      report();
    }
  }
  else if (strcmp(key, "GCONST") == 0){
    if (get_r(gconst) == 0){
      cout << "Error : reading GCONST" << '\n';
      report();
    }
  }
 /*-----IG and BM Additions------*/
  else if (strcmp(key, "SOOTFOIL") == 0) {
    sootfoil = true;
  }
  else if (strcmp(key, "NOSOOTFOIL") == 0)
    sootfoil = false;
  else if (strcmp(key, "TOTALS") == 0) {
    // calculate total mass and energy in domain
    struct cell *cptr; 
    float pa[MAXQ];
    float mass = 0.0;
    float en = 0.0;
    float be = 0.0;
    cptr = lev[0]->fcell; // coarsest
    while (cptr != NULL) {
      usrpq(cptr->qa, pa);
      mass += cptr->qa[iqd] * cptr->vol;
      en += cptr->qa[iqe] * cptr->vol;
      be += cptr->qa[iqd] * usrBE(pa)* cptr->vol;
      cptr = cptr->nxt;
    }
    cout << "Total domain mass = " << mass << endl;
    cout << "Total domain int energy = " << en << endl;
    cout << "Total domain binding energy = " << be << endl;
    cout << "Total domain energy = " << en + be << endl;
  }
  else if (strcmp(key, "AMPLITUDE") == 0){
    if (get_r(amplitude) == 0){
      cout << "Error : reading CNS1" << '\n';
      report();
    }
  }
  else if (strcmp(key, "PERIOD") == 0){
    if (get_r(period) == 0){
      cout << "Error : reading PERIOD" << '\n';
      report();
    }
  }
 /*----------------------------*/
  else if (strcmp(key, "PARTICLE") == 0)
    particle = true;
  else if (strcmp(key, "NOPARTICLE") == 0){
    int icell;
    for (icell = 1; icell <= celltop; icell++)
      del_allpart(carr[icell]);
    particle = false;
  }
  else if (strcmp(key, "MKPART") == 0){
    int id, ilev;
    float xc[MAXD], pa[MAXQ];
    bool done = true;
    struct part *pptr;
    struct cell *cptr;
    if (mype == 0)
      cout << "Input position, \n";
    for (id = 1; id <= nd; id++){
      if (!get_r(xc[id-1])){
	if (mype == 0)
	  cout << "Error : reading pos for dim" << id << '\n';
	done = false;
	break;
	report();
      }
    }
    if (done){
      pptr = usrmkpart(xc);
      if (pptr != NULL){
	cptr = pptr->cptr;
	usrpq(cptr->qa, pa);
	pptr->props[0] = -1.0;
	pptr->props[1] = pa[iqd];	
	pptr->props[2] = pa[iqu0+1];
	pptr->props[3] = pa[iqte];	
	for (id = 1; id <= nd;id++)	
	  pptr->vp[id-1] = 0.0;
	pptr->vp[0] = pa[iqu0+1];
	cout << "particle created on " << mype <<'\n';
      }
    }
  }
  else if (strcmp(key, "PRINTPART") == 0){
    int npart;
    struct part *pptr;
    struct part_list *plist, *tplist, *start, *end, *ostart;
    start = NULL;
    ostart = NULL;    
    /* Get particles from other processors */
    ostart = mov_part(true, NULL);
    if (mype == 0){
      /* Add particles on master */
      mk_partlist(start, end);
      if (start == NULL)
	start = ostart;
      else
	end->nxt = ostart;
      if (start == NULL)
	cout << "note:  no particles" << '\n';
      else{
      /* Print particle data */
	npart = 1;
	plist = start;
	  cout << std::scientific;
	  cout << std::setprecision(10);	
	while (plist != NULL){
	  pptr = plist->pptr;
	  cout << "no " << npart << " pos " << pptr->xp[0] << " t_j " <<
	    t - pptr->props[0] << '\n';
	  cout << "rho " << pptr->props[1] <<  " vel " <<
	    pptr->props[2] << " te " << pptr->props[3] << '\n';
	  cout << "pptr pos " << pptr->xp[0] <<  " vel " <<
	    pptr->vp[0] << '\n';
	  npart++;
	  plist = plist->nxt;
	}
      }
    }
  }
  else
    /* Command not known here, pass it back to caller */
    found = false;

  return(found);
}
//end usrcmd.c

//file usrprt.c
/*
+ usrprt - Print user information.
    Description :
     Prints user information. It is called by the USERS command.
*/
void usrprt()
{
  cout.setf(ios::scientific);
  cout << setprecision(4);
  cout << "no of scalars " << nal << " no of advecting scalars "
       << scad << '\n';
  cout << "dfloor " << dfloor;
  cout << " pfloor " << pfloor;
  cout << " avisp " << avisp;
  cout << " avise " << avise << '\n';
  cout << "gslatol " << gslatol << '\n';
  cout << "gslrtol " << gslrtol << '\n';
  cout << "gslhmin " << gslhmin << '\n';
  cout << "gslhmax " << gslhmax << '\n';
  cout << "nsub " << nsub << '\n';

  cout << "nsub " << nsub << " solev " << solev << " cns " << cns
       << " cns1 " << cns1 << " cns2 " << cns2 << '\n';  
  cout << nlevs << " levels " << "Finest mesh spacing " <<
    lev[nlevs-1]->dxbar[0] << '\n';
  cout << "g " << g <<" rat " << rat << " curv " << curv <<
    " depth " << depth <<" D " << D << '\n';
  cout << "T_0 " << T_0 << " rho_0 " << rho_0;
  cout << " T_s " << T_s << " rho_s " << rho_s <<  '\n';
  cout << "tp_old " << tp_old << " dtp " << dtp << '\n';    
  if (god)
    cout << "Exact solver ";
  if (hll)
    cout << "HLL solver ";
  if (kt)
    cout << "Kurganov-Tadmor solver ";
  if (viscous)
    cout << " viscous " << " vis " << vis << " tcond " << tcond << 
      " sdiff " << sdiff <<'\n';
  else 
    cout << " inviscid " << '\n';
  if (keps)
    cout << "k-epsilon: kmin " << kmin << "  epsmin " << epsmin <<
      " max turb length scale mtl " << mtl << '\n';
  if (inten)
    cout << "calculating inten " << '\n';
  if (particle)
    cout << "particles " << " deleted " << delpart << " particles" << 
      '\n';
  else
    cout << "noparticles " << '\n';
  if (grav){
    cout << "gravity nrel " << nrel << " normrat " << 
      normrat << " gconst " << gconst;
    if (boundpot)
      cout << " boundpot levb " <<  levb << '\n';
    else
      cout << " noboundpot " << '\n';
  }
  else
    cout << "nogravity " << '\n';
  if (mon)
    cout << "monitoring " << '\n';
  else 
    cout << "nomonitoring " << '\n';
  if (react){
    if (rrscreen)
      cout << "screen ";
    else
      cout << "noscreen ";
    if (revscreen)
      cout << "revscreen ";
    else
      cout << "norevscreen ";
    if (doneut == 1)
      cout << "neut loss ";
    else
      cout << " ";
    if (doweak == 1)
      cout << "weak ";
    else
      cout << "no weak ";        
    cout << "iso" << niso << " rates " << nrat << " " <<ode;
    if (cons)
      cout << " cons ";
    else
      cout << " nocons ";
  }
  else
    cout << "No reaction" << '\n';
  if (hydrostatic)
    cout << "hydrostatic";
  else
    cout << "dynamic";
  if (isothermal)
    cout << " isothermal" << '\n';
  else
    cout << " not isothermal" << '\n';               
  return;
}
//end usrprt.c

//file usrwr.c
/*
+ usrwr - Write user information.
    Description :
     Writes  user information to  a model  file. It  is called  by the
     WRITE command.
*/
void usrwr()
{
  write_var("iqd", iqd);
  write_var("iqu0", iqu0);
  write_var("iqe", iqe);
  write_var("keps", keps);
  if (keps){
    write_var("iqk", iqk);
    write_var("iqeps", iqeps);
    write_var("ce1", ce1);
    write_var("ce2", ce2);
    write_var("cs", cs);
    write_var("ct", ct);
    write_var("cd", cd);
    write_var("ce", ce);
    write_var("ceps", ceps);
    write_var("kmin", kmin);
    write_var("epsmin", epsmin);
    write_var("mtl", mtl);
  }
  write_var("iqal0", iqal0);
  write_var("iqie", iqie);
  write_var("delpart", delpart); 
  write_var("iqsf", iqsf);
  write_var("iqpot", iqpot);
  write_var("iqte", iqte);
  write_var("niso", niso);
  write_var("g", g);
  write_var("pfloor", pfloor);
  write_var("dfloor", dfloor);
  write_var("avisp", avisp);
  write_var("avise", avise);
  write_var("mon", mon);
  write_var("fxp", fxp);
  write_var("inten", inten);
  write_var("god", god);
  write_var("hll", hll);
  write_var("kt", kt);
  write_var("react", react);
  write_var("doneut", doneut);
  write_var("doweak", doweak);    
  write_var("cons", cons);    
  write_var("ode", ode);
  write_var("nsub", nsub);  
  write_var("gslatol", gslatol);
  write_var("gslrtol", gslrtol);
  write_var("gslhmin", gslhmin);
  write_var("gslhmax", gslhmax);
  write_var("hydrostatic", hydrostatic);
  write_var("isothermal", isothermal);    
  write_var("rrscreen", rrscreen);
  write_var("revscreen", revscreen);    
  write_var("viscous", viscous);
  write_var("vis", vis);
  write_var("tcond", tcond);
  write_var("sdiff", sdiff);
  write_var("gconst", gconst);
  write_var("rat", rat);
  write_var("curv", curv);
  write_var("depth", depth);    
  write_var("T_e", T_e);
  write_var("T_0", T_0);
  write_var("rho_0", rho_0);
  write_var("T_s", T_s);
  write_var("rho_s", rho_s);
  write_var("par", MAXPAR, par);
  write_var("solev", solev);
  write_var("cns1", cns1);
  write_var("cns2", cns2);    
  write_var("netname", netname);
  write_var("cns", cns);
  write_var("sootfoil", sootfoil);
  write_var("amplitude", amplitude);
  write_var("period", period);

  return;
}
//end usrwr.c

//file usrrd.c
/*
+ usrrd - Read user information.
    Description :

      Reads user  information to a  model file.   It is called  by the
      READ command.
*/
bool usrrd(char vname[])
{
  bool found;
  found = true;
  if (strcmp(vname, "iqd") == 0)
    read_var(&iqd);
  else if (strcmp(vname, "iqu0") == 0)
    read_var(&iqu0);
  else if (strcmp(vname, "iqe") == 0)
    read_var(&iqe);
  else if (strcmp(vname, "keps") == 0)
    read_var(&keps);
  else if (strcmp(vname, "iqk") == 0)
    read_var(&iqk);
  else if (strcmp(vname, "iqeps") == 0)
    read_var(&iqeps);
  else if (strcmp(vname, "iqal0") == 0)
    read_var(&iqal0); 
  else if (strcmp(vname, "iqie") == 0)
    read_var(&iqie);
  else if (strcmp(vname, "delpart") == 0)
    read_var(&delpart);
  else if (strcmp(vname, "iqsf") == 0)
    read_var(&iqsf);
  else if (strcmp(vname, "iqpot") == 0)
    read_var(&iqpot);
  else if (strcmp(vname, "iqte") == 0)
    read_var(&iqte);
  else if (strcmp(vname, "niso") == 0)
    read_var(&niso);
  else if (strcmp(vname, "ce1") == 0)
    read_var(&ce1);
  else if (strcmp(vname, "ce2") == 0)
    read_var(&ce2);
  else if (strcmp(vname, "cs") == 0)
    read_var(&cs);
  else if (strcmp(vname, "ct") == 0)
    read_var(&ct);
  else if (strcmp(vname, "cd") == 0)
    read_var(&cd);
  else if (strcmp(vname, "ce") == 0)
    read_var(&ce);
  else if (strcmp(vname, "ceps") == 0)
    read_var(&ceps);
  else if (strcmp(vname, "kmin") == 0)
    read_var(&kmin);
  else if (strcmp(vname, "epsmin") == 0)
    read_var(&epsmin);
  else if (strcmp(vname, "mtl") == 0)
    read_var(&mtl);
  else if (strcmp(vname, "g") == 0)
      read_var(&g);
  else if (strcmp(vname, "pfloor") == 0)
    read_var(&pfloor);
  else if (strcmp(vname, "dfloor") == 0)
    read_var(&dfloor);
  else if (strcmp(vname, "avisp") == 0)
    read_var(&avisp);
  else if (strcmp(vname, "avise") == 0)
    read_var(&avise);
  else if (strcmp(vname, "mon") == 0)
    read_var(&mon);
  else if (strcmp(vname, "fxp") == 0)
    read_var(&fxp);
  else if (strcmp(vname, "inten") == 0)
    read_var(&inten);
  else if (strcmp(vname, "hll") == 0)
    read_var(&hll);  
  else if (strcmp(vname, "god") == 0)
    read_var(&god);
  else if (strcmp(vname, "kt") == 0)
    read_var(&kt);
  else if (strcmp(vname, "react") == 0)
    read_var(&react);
  else if (strcmp(vname, "doneut") == 0)
    read_var(&doneut);
  else if (strcmp(vname, "doweak") == 0)
    read_var(&doweak);    
  else if (strcmp(vname, "cons") == 0)
    read_var(&cons);  
  else if (strcmp(vname, "ode") == 0)
    read_var(ode);
  else if (strcmp(vname, "nsub") == 0)
    read_var(&nsub);  
  else if (strcmp(vname, "gslatol") == 0)
    read_var(&gslatol);
  else if (strcmp(vname, "gslrtol") == 0)
    read_var(&gslrtol);
  else if (strcmp(vname, "gslhmin") == 0)
    read_var(&gslhmin);
  else if (strcmp(vname, "gslhmax") == 0)
    read_var(&gslhmax);
  else if (strcmp(vname, "hydrostatic") == 0)
    read_var(&hydrostatic);
  else if (strcmp(vname, "isothermal") == 0)
    read_var(&isothermal);    
  else if (strcmp(vname, "rrscreen") == 0)
    read_var(&rrscreen);
  else if (strcmp(vname, "revscreen") == 0)
    read_var(&revscreen);  
  else if (strcmp(vname, "viscous") == 0)
    read_var(&viscous);
  else if (strcmp(vname, "vis") == 0)
    read_var(&vis);
  else if (strcmp(vname, "tcond") == 0)
    read_var(&tcond);
  else if (strcmp(vname, "sdiff") == 0)
    read_var(&sdiff);
  else if (strcmp(vname, "gconst") == 0)
    read_var(&gconst);
  else if (strcmp(vname, "rat") == 0)
    read_var(&rat);
  else if (strcmp(vname, "curv") == 0)
    read_var(&curv);
  else if (strcmp(vname, "depth") == 0)
    read_var(&depth);    
  else if (strcmp(vname, "T_e") == 0)
    read_var(&T_e);
  else if (strcmp(vname, "T_0") == 0)
    read_var(&T_0);
  else if (strcmp(vname, "rho_0") == 0)
    read_var(&rho_0);
  else if (strcmp(vname, "T_s") == 0)
    read_var(&T_s);
  else if (strcmp(vname, "rho_s") == 0)
    read_var(&rho_s);
  else if (strcmp(vname, "par") == 0)
    read_var(par);
  else if (strcmp(vname, "solev") == 0)
    read_var(&solev);
  else if (strcmp(vname, "cns1") == 0)
    read_var(&cns1);
  else if (strcmp(vname, "cns2") == 0)
    read_var(&cns2);
  else if (strcmp(vname, "netname") == 0)
    read_var(netname);        
  else if (strcmp(vname, "cns") == 0)
    read_var(&cns);
  else if (strcmp(vname, "sootfoil") == 0)
    read_var(&sootfoil); 
  else if (strcmp(vname, "amplitude") == 0)
    read_var(&amplitude);
  else if (strcmp(vname, "period") == 0)
    read_var(&period);
  else
    /* Unknown here  */
    found = false;
  /* Check for turbulence */
  keps = iqk > 0;
  if (strcmp(vname, "cns") == 0){
    if (strcmp(netname, "netname") == 0)
      init = true;
    else{
      /* Initialise networks */
      int len;
      char tail[SZNAM];
      bool fail = true;
      tim7 = false;		
      tim13 = false;
      tim21 = false;
      timtorch = false;
      fail = false;
      len = strlen(netname);
      tail[0] = netname[len-3]; tail[1] = netname[len-2];
      tail[2] = netname[len-1]; tail[3] = netname[len];              
      if (strncmp(netname, "tim07", 5) == 0){
	tim7 = true;		
	niso = 7;
	nrat = 16;
	rrscreen = true;      
	revscreen = true;
      }
      else if (strncmp(netname, "tim13", 5) == 0){
	tim13 = true;
	niso = 13;
	nrat = 67;
	if (!init)
	  /* Initialise 13 network */
	  aprox13_init_(aion, zion, bion);
      }
      else if (strncmp(netname, "tim21", 5) == 0){
	tim21 = true;
	niso = 21;
	nrat = 112;
	if (!init)
	  /* Initialise 21 network */
	  aprox21_init_(aion, zion, bion);
      }
      else if (strncmp(tail, "net", 3) == 0){      
	timtorch = true;	
	/* Get network from netname */
	if (!init)
	  /* Initialise MESA_TIM network */
	  torch_init_(&niso, &nrat, aion, zion, bion, netname);
      }
      else{
	fail = true;
	if (mype == 0)
	  cout << "Error unknown network; use REACT to fix \n";
      }
      if (!fail){
	usrinitnet();      
	init = true;
      }
    }
  }      
  /* Initialise supernova EOS */
  usrinithelm();
  /* nomon for safety */
  mon = false;
  return(found);
}
//end usrrd.c

//file usrsum.c
/*
+ usrsum - sum variables over the grid.
     Description :
        Calculates the sum of variables over the grid. The sum is over
        levels  plmin  to  plmax.  If theer  is  a  reaction  progress
        variable, it returnd the detonation speed.
*/
float usrsum()
{
  int i, ilev;
  float cvol, tot[MAXQ], pa[MAXQ], ydot[MAXSPECIES+1], so[MAXQ], sumso,
    xc[MAXD];
  char old_ode[SZNAM];  
  struct cell *cptr, *tcptr;
  tcptr = new cell;
  /* Return if no model */
  if (nom)
    return(-1.0);
  /* Use NONE */
  strcpy(old_ode, ode);    
  strcpy(ode, "NONE");  
  /* Zero sum, */
  for (i = 0; i < nq; i++)
    tot[i] = 0.0;
  sumso = 0.0;
  /* Sum of all variables */
  for (ilev = plmin; ilev <= plmax; ilev++){
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      if (!cptr->refined){
	cvol = cptr->vol;
	for (i = 0; i < nq; i++){
	  tcptr->qa[i] = cptr->qa[i];
	  tot[i] += cvol*cptr->qa[i];
	}
	usrpq(tcptr->qa, tcptr->pa);
	if (nal > niso){
	  usrsaf(1, 1.0e-8, tcptr, so);
	  sumso += cvol*(curv*tcptr->qa[iqu0+1]*tcptr->pa[iqal0+nal]
	    + so[iqal0+nal]);
	}
      }
      cptr = cptr->nxt;
    }
  }
  delete tcptr;
  /* Restore old ode solver */
  strcpy(ode, old_ode);  
  /* Get global sum */
  for (i = 0; i < nq; i++)
    tot[i] = global_sum(tot[i]);
  sumso = global_sum(sumso);
  for (i = 1; i <= nd; i++)
    xc[i] = 0.0; 
  xc[0] = xr[0] - 0.01*(xr[0] - xl[0]);
  cptr = find_fcell(ilev, xc);
  if (cptr != NULL){
    usrqis(xc, cptr->qa);
    if (nal > niso)
      sumso = sumso/cptr->qa[iqd];
    else
      sumso = -1.0;
  }
  else
    sumso = -1.0;
  sumso = global_max(sumso);
  /* Output sum of all variables */
  if (nd == 2)
    sumso /=(xr[1] - xl[1]);  
  if (mype == 0){
    cout.setf(ios::scientific);
    cout << "det speed (ignore -f -1 " << setprecision(8) << sumso <<
      '\n';  
    if (shockvelfile != NULL)
      fprintf(shockvelfile, "%e %e \n", t, sumso);
  }
  return(sumso);
}
//end usrsum.c  

//file usropenxq.c
/*
+ usropenxq - Opens xq files
     Description :
        This  is called  by  the command  MON and  sets  up files  for
        outputting abundances  etc.  Can be used  for constructing .xq
        files in usrmaldatconv.c mode =  1 generates all files, mode =
        3 generates only files for atmosphere.
*/
bool usropenxq(char monmod[], int mode)
{
  char filename[SZNAM], tail[SZNAM], tempfn[SZNAM];
  int i;
  if (mype != 0)
    return(false);
  if (strcmp(monmod, "POS") == 0){
    if (mode == 1){
      for (i = 0; i < niso; i++){
	strcpy(filename, monpref);
	//      sprintf(tail, "_xa_%d_t.xq", i+1);
	sprintf(tail, "_xa_%s_t.xq", ionam[i]);      
	strcat(filename,tail);      
	if (!(xatfile[i]=fopen(filename, "w"))){
	  break;
	  return(false);
	}
	fprintf(xatfile[i], "%d %d \n", 1, 10);
	fprintf(xatfile[i], "xat\n");
	fprintf(xatfile[i], "t xa\n");      
	strcpy(filename, monpref);
	//      sprintf(tail, "_xa_%d_x.xq", i+1);
	sprintf(tail, "_xa_%s_x.xq", ionam[i]);      
	strcat(filename,tail);            
	if (!(xaxfile[i]=fopen(filename, "w"))){
	  break;
	  return(false);
	}
	fprintf(xaxfile[i], "%d %d \n", 1, 10);
	fprintf(xaxfile[i], "xax \n");
	fprintf(xaxfile[i], "x xa\n");
      }
      strcpy(filename, monpref);
      sprintf(tail, "_temp_t.xq");    
      strcat(filename,tail);
      if (!(tempfilet=fopen(filename, "w")))
	return(false);    
      fprintf(tempfilet, "%d %d \n", 1, 10);
      fprintf(tempfilet, "temp \n");
      fprintf(tempfilet, "t temp\n");
      strcpy(filename, monpref);
      sprintf(tail, "_mach_t.xq");        
      strcat(filename,tail);
      if (!(machfilet=fopen(filename, "w")))
	return(false);
      fprintf(machfilet, "%d %d \n", 1, 10);
      fprintf(machfilet, "mach \n");
      fprintf(machfilet, "t mach\n");
      strcpy(filename, monpref);
      sprintf(tail, "_mach_x.xq");        
      strcat(filename,tail);
      if (!(machfilex=fopen(filename, "w")))
	return(false);
      fprintf(machfilex, "%d %d \n", 1, 10);
      fprintf(machfilex, "mach \n");
      fprintf(machfilex, "x mach\n");
      strcpy(filename, monpref);
      sprintf(tail, "_rho_t.xq");
      strcat(filename,tail);
      if (!(rhofilet=fopen(filename, "w")))
	return(false);    
      fprintf(rhofilet, "%d %d \n", 1, 10);
      fprintf(rhofilet, "rho \n");
      fprintf(rhofilet, "t rho\n");
      strcpy(filename, monpref);
      sprintf(tail, "_en_t.xq");
      strcat(filename,tail);    
      if (!(enfilet=fopen(filename, "w")))
	return(false);    
      fprintf(enfilet, "%d %d \n", 1, 10);
      fprintf(enfilet, "en \n");
      fprintf(enfilet, "t en\n");
      strcpy(filename, monpref);
      sprintf(tail, "_en_x.xq");
      strcat(filename,tail);    
      if (!(enfilex=fopen(filename, "w")))
	return(false);    
      fprintf(enfilex, "%d %d \n", 1, 10);
      fprintf(enfilex, "en \n");
      fprintf(enfilex, "x en\n");
      strcpy(filename, monpref);
      sprintf(tail, "_phi_t.xq");
      strcat(filename,tail);    
      if (!(phifilet=fopen(filename, "w")))
	return(false);    
      fprintf(phifilet, "%d %d \n", 1, 10);
      fprintf(phifilet, "phi \n");
      fprintf(phifilet, "t phi\n");
      strcpy(filename, monpref);
      sprintf(tail, "_phi_x.xq");
      strcat(filename,tail);    
      if (!(phifilex=fopen(filename, "w")))
	return(false);    
      fprintf(phifilex, "%d %d \n", 1, 10);
      fprintf(phifilex, "phi \n");
      fprintf(phifilex, "x phi\n");
      strcpy(filename, monpref);
      sprintf(tail, "_kap_t.xq");
      strcat(filename,tail);    
      if (!(kapfilet=fopen(filename, "w")))
	return(false);    
      fprintf(kapfilet, "%d %d \n", 1, 10);
      fprintf(kapfilet, "kap \n");
      fprintf(kapfilet, "t kap\n");
      strcpy(filename, monpref);
      sprintf(tail, "_kap_x.xq");
      strcat(filename,tail);    
      if (!(kapfilex=fopen(filename, "w")))
	return(false);    
      fprintf(kapfilex, "%d %d \n", 1, 10);
      fprintf(kapfilex, "kap \n");
      fprintf(kapfilex, "x kap\n");
      strcpy(filename, monpref);
      sprintf(tail, "_bind_t.xq");
      strcat(filename,tail);    
      if (!(bindfilet=fopen(filename, "w")))
	return(false);    
      fprintf(bindfilet, "%d %d \n", 1, 10);
      fprintf(bindfilet, "bind \n");
      fprintf(bindfilet, "t bind\n");
      strcpy(filename, monpref);
      sprintf(tail, "_bind_x.xq");
      strcat(filename,tail);    
      if (!(bindfilex=fopen(filename, "w")))
	return(false);    
      fprintf(bindfilex, "%d %d \n", 1, 10);
      fprintf(bindfilex, "bind \n");
      fprintf(bindfilex, "x bind\n");            
      strcpy(filename, monpref);
      sprintf(tail, "_pg_t.xq");
      strcat(filename,tail);    
      if (!(pgfilet=fopen(filename, "w")))
	return(false);    
      fprintf(pgfilet, "%d %d \n", 1, 10);
      fprintf(pgfilet, "pg \n");
      fprintf(pgfilet, "t pg\n");
      strcpy(filename, monpref);
      sprintf(tail, "_pg_x.xq");
      strcat(filename,tail);    
      if (!(pgfilex=fopen(filename, "w")))
	return(false);    
      fprintf(pgfilex, "%d %d \n", 1, 10);
      fprintf(pgfilex, "pg \n");
      fprintf(pgfilex, "x pg\n");    
      strcpy(filename, monpref);
      sprintf(tail, "_u_t.xq");
      strcat(filename,tail);    
      if (!(ufilet=fopen(filename, "w")))
	return(false);    
      fprintf(ufilet, "%d %d \n", 1, 10);
      fprintf(ufilet, "u \n");
      fprintf(ufilet, "t u\n");
      strcpy(filename, monpref);
      sprintf(tail, "_u_x.xq");
      strcat(filename,tail);    
      if (!(ufilex=fopen(filename, "w")))
	return(false);    
      fprintf(ufilex, "%d %d \n", 1, 10);
      fprintf(ufilex, "u \n");
      fprintf(ufilex, "x u\n");
      strcpy(filename, monpref);
      sprintf(tail, "_abar_t.xq");
      strcat(filename,tail);    
      if (!(abarfilet=fopen(filename, "w")))
	return(false);    
      fprintf(abarfilet, "%d %d \n", 1, 10);
      fprintf(abarfilet, "abar \n");
      fprintf(abarfilet, "t abar\n");
      strcpy(filename, monpref);
      sprintf(tail, "_abar_x.xq");
      strcat(filename,tail);    
      if (!(abarfilex=fopen(filename, "w")))
	return(false);
      fprintf(abarfilex, "%d %d \n", 1, 10);
      fprintf(abarfilex, "abar \n");
      fprintf(abarfilex, "x abar\n");      
      strcpy(filename, monpref);
      sprintf(tail, "_phimach.xq");
      strcat(filename,tail);        
      if (!(phimachfilex=fopen(filename, "w")))
	return(false);    
      fprintf(phimachfilex, "%d %d \n", 1, 10);
      fprintf(phimachfilex, "phimach\n");
      fprintf(phimachfilex, " mach phi\n");      
    }
    strcpy(filename, monpref);
    sprintf(tail, "_temp_x.xq");    
    strcat(filename,tail);
    if (!(tempfilex=fopen(filename, "w")))
      return(false);    
    fprintf(tempfilex, "%d %d \n", 1, 10);
    fprintf(tempfilex, "temp \n");
    fprintf(tempfilex, "x temp\n");         
    strcpy(filename, monpref);
    sprintf(tail, "_rho_x.xq");
    strcat(filename,tail);
    if (!(rhofilex=fopen(filename, "w")))
      return(false);    
    fprintf(rhofilex, "%d %d \n", 1, 10);
    fprintf(rhofilex, "rho \n");
    fprintf(rhofilex, "x rho\n");
    strcpy(filename, monpref);
    sprintf(tail, "_eni_x.xq");
    strcat(filename,tail);
    if (!(enifilex=fopen(filename, "w")))
      return(false);    
    fprintf(enifilex, "%d %d \n", 1, 10);
    fprintf(enifilex, "eni \n");
    fprintf(enifilex, "x eni\n");    
    return (true);            
  }
  else if (strcmp(monmod, "VEL") == 0){
    sprintf(tail, "sumvel.xq");
    strcpy(filename, monpref);    
    strcat(filename,tail);
    if (!(shockvelfile=fopen(filename, "w")))
      return(false);
    fprintf(shockvelfile, "%d %d \n", 1, 10);
    fprintf(shockvelfile, "vel \n");
    fprintf(shockvelfile, "t vel \n");
  }
  else if (strcmp(monmod, "PART") == 0){
    if (!(trho=fopen("trho.raj", "w")))
      return(false);
    fprintf(trho, "%d %d \n", 1, 10);
    fprintf(trho, "Dens \n");;        
    fprintf(trho,"Temp Dens \n");
  }
  else if (strcmp(monmod, "ENGEN") == 0){
    if (!(engen=fopen("engen.xq", "w")))
      return(false);
    fprintf(engen, "%d %d \n", 1, 10);
    fprintf(engen, "Engen \n");;        
    fprintf(engen,"t engen \n");
  }  
  else
    return(false);
  return(true);
}    
//end usropenxq.c

//file usrrewxq.c
/*
+ usrrewxq - Reqind  xq files
     Description :
        This rewinds xq files that  have been opened by usropenxq.c It
        only does  mode = 1  files. Used  in modsol.c to  rewind files
        after each detonation iteration.
*/
bool usrrewxq(char monmod[])
{
  int i;
  if (mype != 0)
    return(false);
  for (i = 0; i < niso; i++){
    rewind(xatfile[i]);
    fprintf(xatfile[i], "%d %d \n", 1, 10);
    fprintf(xatfile[i], "xat\n");
    fprintf(xatfile[i], "t xa\n");
    rewind(xaxfile[i]);	
    fprintf(xaxfile[i], "%d %d \n", 1, 10);
    fprintf(xaxfile[i], "xax \n");
    fprintf(xaxfile[i], "x xa\n");
  }
  rewind(tempfilet);
  fprintf(tempfilet, "%d %d \n", 1, 10);
  fprintf(tempfilet, "temp \n");
  fprintf(tempfilet, "t temp\n");
  rewind(tempfilex);
  fprintf(tempfilex, "%d %d \n", 1, 10);
  fprintf(tempfilex, "temp \n");
  fprintf(tempfilex, "x temp\n");         
  rewind(machfilet);
  fprintf(machfilet, "%d %d \n", 1, 10);
  fprintf(machfilet, "mach \n");
  fprintf(machfilet, "t mach\n");
  rewind(machfilex);
  fprintf(machfilex, "%d %d \n", 1, 10);
  fprintf(machfilex, "mach \n");
  fprintf(machfilex, "x mach\n");
  rewind(rhofilet);
  fprintf(rhofilet, "%d %d \n", 1, 10);
  fprintf(rhofilet, "rho \n");
  fprintf(rhofilet, "t rho\n");
  rewind(rhofilex);
  fprintf(rhofilex, "%d %d \n", 1, 10);
  fprintf(rhofilex, "rho \n");
  fprintf(rhofilex, "x rho\n");
  rewind(enfilet);
  fprintf(enfilet, "%d %d \n", 1, 10);
  fprintf(enfilet, "en \n");
  fprintf(enfilet, "t en\n");
  rewind(enfilex);
  fprintf(enfilex, "%d %d \n", 1, 10);
  fprintf(enfilex, "en \n");
  fprintf(enfilex, "x en\n");
  rewind(phifilet);
  fprintf(phifilet, "%d %d \n", 1, 10);
  fprintf(phifilet, "phi \n");
  fprintf(phifilet, "t phi\n");
  rewind(phifilex);
  fprintf(phifilex, "%d %d \n", 1, 10);
  fprintf(phifilex, "phi \n");
  fprintf(phifilex, "x phi\n");
  rewind(kapfilet);
  fprintf(kapfilet, "%d %d \n", 1, 10);
  fprintf(kapfilet, "kap \n");
  fprintf(kapfilet, "t kap\n");
  rewind(kapfilex);
  fprintf(kapfilex, "%d %d \n", 1, 10);
  fprintf(kapfilex, "kap \n");
  fprintf(kapfilex, "x kap\n");
  rewind(bindfilet);
  fprintf(bindfilet, "%d %d \n", 1, 10);
  fprintf(bindfilet, "bind \n");
  fprintf(bindfilet, "t bind\n");
  rewind(bindfilex);
  fprintf(bindfilex, "%d %d \n", 1, 10);
  fprintf(bindfilex, "bind \n");
  fprintf(bindfilex, "x bind\n");            
  rewind(pgfilet);
  fprintf(pgfilet, "%d %d \n", 1, 10);
  fprintf(pgfilet, "pg \n");
  fprintf(pgfilet, "t pg\n");
  rewind(pgfilex);
  fprintf(pgfilex, "%d %d \n", 1, 10);
  fprintf(pgfilex, "pg \n");
  fprintf(pgfilex, "x pg\n");    
  rewind(ufilet);
  fprintf(ufilet, "%d %d \n", 1, 10);
  fprintf(ufilet, "u \n");
  fprintf(ufilet, "t u\n");
  rewind(ufilex);
  fprintf(ufilex, "%d %d \n", 1, 10);
  fprintf(ufilex, "u \n");
  fprintf(ufilex, "x u\n");
  rewind(abarfilet);
  fprintf(abarfilet, "%d %d \n", 1, 10);
  fprintf(abarfilet, "abar \n");
  fprintf(abarfilet, "t abar\n");
  rewind(abarfilex);
  fprintf(abarfilex, "%d %d \n", 1, 10);
  fprintf(abarfilex, "abar \n");
  fprintf(abarfilex, "x abar\n");      
  rewind(phimachfilex);
  fprintf(phimachfilex, "%d %d \n", 1, 10);
  fprintf(phimachfilex, "phimach\n");
  fprintf(phimachfilex, " mach phi\n");      
  fprintf(enifilex, "%d %d \n", 1, 10);
  fprintf(enifilex, "eni \n");
  fprintf(enifilex, "x eni\n");
  return(true);
}
//end usrrewxq.c  
  
//file usrmyclose.c
/*
+ usrmyclose - Closes files
     Description :
        This is  closes file  and sets  its pointer  to NULL.  Used in
        usrclosexq.
*/
void usrmyclose(FILE *file)
{
  if (file != NULL){
    fclose(file);
    file = NULL;
  }
  return;
}
//file usrmyclose.c

//file usrclosexq.c
/*
+ usrclosexq - Closes xq files
     Description :
        This is called  by the command NOMON and closes  the files set
        up by  usropenxq. mode = 1  closes all files, mode  = 3 closes
        only files for atmosphere.
*/
void usrclosexq(const char monmod[], int mode)
{
  int i;
  if (mype != 0)
    return;  
  if (strcmp(monmod, "POS") == 0){
    if (mode == 1){
      for (i = 0; i < niso; i++){
	usrmyclose(xatfile[i]);
	usrmyclose(xaxfile[i]);  
      }
      usrmyclose(rhofilet);
      usrmyclose(rhofilex);
      usrmyclose(enifilex);      
      usrmyclose(tempfilet);
      usrmyclose(tempfilex);    
      usrmyclose(machfilet);
      usrmyclose(machfilex);    
      usrmyclose(enfilet);
      usrmyclose(enfilex);    
      usrmyclose(phifilet);
      usrmyclose(phifilex);
      usrmyclose(kapfilet);
      usrmyclose(kapfilex);
      usrmyclose(bindfilet);
      usrmyclose(bindfilex);            
      usrmyclose(pgfilet);
      usrmyclose(pgfilex);    
      usrmyclose(ufilet);
      usrmyclose(ufilex);
      usrmyclose(abarfilet);
      usrmyclose(abarfilex);
      usrmyclose(phimachfilex);      
    }
    else{
      usrmyclose(tempfilex);
      usrmyclose(rhofilex);
    }
  }
  else if (strcmp(monmod, "VEL") == 0)
    usrmyclose(shockvelfile);
  else if (strcmp(monmod, "PART") == 0)
    usrmyclose(trho);
  else if (strcmp(monmod, "ENGEN") == 0)
    usrmyclose(engen);
  return;
}
//end usrclosexq.c

//file usrmon.c
/*
+ usrmon - Monitors solution.
     Description :
        This is called  every timestep if mon=true. It  can be used to
        output quantities to stdout or a file at selected intervals.
*/
void usrmon(char monmod[])
{
  int ilev, i, icell;
  float xc[MAXD], pa[MAXQ], tj, rho, temp, mach, ss, phi, u, kterm,
    binden, abar, ydot[MAXSPECIES+1];
  bool old_isothermal;
  struct cell *cptr;
  if (mon && (t - tp_old > dtp)){
    tp_old = t;
    if (strcmp(monmod, "POS") == 0){
      if (mype == 0){
	cptr = carr[1];
	usrpq(cptr->qa, pa);
	for (i = 0; i < niso; i++){
	  fprintf(xatfile[i], "%e %e \n", t, pa[iqal0+1+i]*Ai[i]);
	  fprintf(xaxfile[i], "%e %e \n", pa[iqal0+nal],
		  pa[iqal0+1+i]*Ai[i]);
	}
	fprintf(tempfilet, "%e %e \n", t, pa[iqte]);
	fprintf(tempfilex, "%e %e \n", pa[iqal0+nal], pa[iqte]);
	ss = usrss(pa);
	u = pa[iqu0+1];		
	mach = u/ss;
	modsrc(pa, ydot);
	phi = 1.0e9*ydot[nal]*(ss*ss - u*u)/u;
	binden = 0.0;
	for (i = 0; i < niso; i++)
	  binden -= pa[iqal0+1+i]*Qi[i];
        kterm = curv*pa[iqd]*u*(D - u);	
	fprintf(machfilet, "%e %e \n", t, mach);
	fprintf(machfilex, "%e %e \n", pa[iqal0+nal], mach);	
	fprintf(rhofilet, "%e %e \n", t, pa[iqd]);
	fprintf(rhofilex, "%e %e \n", pa[iqal0+nal], pa[iqd]);
	old_isothermal = isothermal;
	isothermal = false;
	usrreact(pa, ydot);
	isothermal = old_isothermal;
	fprintf(enfilet, "%e %e \n", t, ydot[niso]*usrcv(pa));
	fprintf(enfilex, "%e %e \n", pa[iqal0+nal],
	  ydot[niso]*usrcv(pa));
	fprintf(phifilet, "%e %e \n", t, phi);
	fprintf(phifilex, "%e %e \n", pa[iqal0+nal], phi);
	fprintf(kapfilet, "%e %e \n", t, fabs(kterm/phi));
	fprintf(kapfilex, "%e %e \n", pa[iqal0+nal],
	  fabs(kterm/(phi - kterm)));
	fprintf(bindfilet, "%e %e \n", t, binden);
	fprintf(bindfilex, "%e %e \n", pa[iqal0+nal], binden);	
	fprintf(pgfilet, "%e %e \n", t, pa[iqe]);
	fprintf(pgfilex, "%e %e \n", pa[iqal0+nal], pa[iqe]);	
	fprintf(ufilet, "%e %e \n", t, u);
	fprintf(ufilex, "%e %e \n", pa[iqal0+nal], pa[iqu0+1]);
	abar = usrAbar(pa);
	fprintf(abarfilet, "%e %e \n", t, abar);
	fprintf(abarfilex, "%e %e \n", pa[iqal0+nal], abar);
      }
    }
    else if (strcmp(monmod, "VEL") == 0)
      usrsum();
    else if (strcmp(monmod, "PART") == 0){
      for (icell = celltop; icell >= 1; icell--){
	cptr = carr[icell];
	if (!cptr->halo){
	  if (cptr->pptr != NULL)
	    break;
	}
      }
      tj = -1.0;
      rho = -1.0;
      temp = -1.0;
      if (cptr->pptr != NULL){
	tj = t - cptr->pptr->props[0];
	rho = cptr->pptr->props[1];
	temp = cptr->pptr->props[3];
      }
      tj = global_max(tj);
      rho = global_max(rho);
      temp = global_max(temp);
      if (mype == 0)
	fprintf(trho, "%e %e %e \n", tj, temp, rho);
    }
    if (strcmp(monmod, "ENGEN") == 0){
      float ydot[MAXSPECIES+1];
      bool old_isothermal;
      if (mype == 0){
	for (i = 0; i < nd; i++)
	  xc[i] = 0.0;
	cptr = find_fcell(ilev, xc);
	usrpq(cptr->qa, pa);
	old_isothermal = isothermal;
	isothermal = false;
	usrreact(pa, ydot);
	isothermal = old_isothermal;	
	fprintf(engen, "%e %e \n", t, ydot[niso]*usrcv(pa));
      }
    }
  }
  return;
  }
//end usrmon.c

//file usrmaldatconv.c
/*
+ usrmaldatconv - Converts MALCOLM .dat file
     Description :
        This converts  a file  from Malcolm's SNTI  to .xq  files with
        prefix monpref input in the MALDATCONV command.
*/
bool usrmaldatconv()
{
  int i, ndata, imax;
  bool doit = true, old_isothermal;
  float tv, x, t, p, rho, T, ss, u, mach, phi, binden, kterm, pa[MAXQ],
    abar, ydot[MAXSPECIES+1], maxen;
  char str[2000], filename[SZNAM];
  //  FILE *reactfile, *primfile;  
  struct datap{
    float x, t;
    float pa[MAXQ];
    struct datap *nxt;
  };
  struct datap *start, *dptr, *tdptr;    
  strcpy(filename, monpref);
  strcat(filename,"_react.data");
  //  if (!(reactfile=fopen(filename, "r")))
  ifstream reactfile(filename);
  if (!reactfile)
    return(false);
  strcpy(filename, monpref);
  strcat(filename,"_prim.data");
  ifstream primfile(filename);
  if (!primfile)  
    return(false);
  /* Read reactfile*/
  for (i = 1; i <= niso +1; i++)
    reactfile >> str;    
  start = new datap;
  dptr = start;
  ndata = 0;
  while (reactfile >> tv){
    dptr->x = tv;      
    for (i = 1; i <= niso; i++){     
      reactfile >> tv;
      dptr->pa[iqal0+i] = tv/Ai[i-1];
    }
    dptr->nxt = new datap;
    tdptr = dptr;
    dptr = dptr->nxt;      
  }
  /* Terminate list */
  delete tdptr->nxt;
  tdptr->nxt = NULL;  
  /* Read primfile*/
  for (i = 1; i <= 17; i++)
    primfile >> str;      
  dptr = start;  
  while (dptr != NULL){
    primfile >> x;
    primfile >> dptr->pa[iqe];
    primfile >> dptr->pa[iqd];
    primfile >> dptr->pa[iqte];
    primfile >> tv;
    primfile >> tv;    
    primfile >> dptr->pa[iqu0+1];
    for (i = 1; i <= 10; i++)
      primfile >> tv;	
    dptr->t = tv;
    dptr = dptr->nxt;
  }
  dptr = start;
  while (dptr != NULL){
    t = dptr->t;
    x = dptr->x;
    for (i = 0; i < nq; i++)
      pa[i] = dptr->pa[i];
    for (i = 0; i < niso; i++){
      fprintf(xatfile[i], "%e %e \n", t, pa[iqal0+1+i]*Ai[i]);
      fprintf(xaxfile[i], "%e %e \n", x,
	      pa[iqal0+1+i]*Ai[i]);
    }
    fprintf(tempfilet, "%e %e \n", t, pa[iqte]);
    fprintf(tempfilex, "%e %e \n", x, pa[iqte]);
    ss = usrss(pa);
    u = pa[iqu0+1];		
    mach = u/ss;
    modsrc(pa, ydot);
    phi = 1.0e9*ydot[nal]*(ss*ss - u*u)/u;
    binden = 0.0;
    for (i = 0; i < niso; i++)
      binden -= pa[iqal0+1+i]*Qi[i];
    kterm = curv*pa[iqd]*u*(D - u);	
    fprintf(machfilet, "%e %e \n", t, mach);
    fprintf(machfilex, "%e %e \n", x, mach);	
    fprintf(rhofilet, "%e %e \n", t, pa[iqd]);
    fprintf(rhofilex, "%e %e \n", x, pa[iqd]);
    old_isothermal = isothermal;
    isothermal = false;
    usrreact(pa, ydot);
    maxen = - big;
    for (i = 0; i < niso; i++){
      tv = fabs(ydot[i]*Qi[i]);
      if (tv  > maxen){
        maxen =  tv;
	imax = i;
      }
    }
    fprintf(enifilex, "%e %e \n", x, float(imax));        
    isothermal = old_isothermal;
    fprintf(enfilet, "%e %e \n", t, ydot[niso]*usrcv(pa));
    fprintf(enfilex, "%e %e \n", x,
	    ydot[niso]*usrcv(pa));
    fprintf(phifilet, "%e %e \n", t, phi);
    fprintf(phifilex, "%e %e \n", x, phi);
    fprintf(kapfilet, "%e %e \n", t, fabs(kterm/phi));
    fprintf(kapfilex, "%e %e \n", x,
	    fabs(kterm/(phi - kterm)));
    fprintf(bindfilet, "%e %e \n", t, binden);
    fprintf(bindfilex, "%e %e \n", x, binden);	
    fprintf(pgfilet, "%e %e \n", t, pa[iqe]);
    fprintf(pgfilex, "%e %e \n", x, pa[iqe]);	
    fprintf(ufilet, "%e %e \n", t, u);
    fprintf(ufilex, "%e %e \n", x, pa[iqu0+1]);
    abar = usrAbar(pa);
    fprintf(abarfilet, "%e %e \n", t, abar);
    fprintf(abarfilex, "%e %e \n", x, abar);
    fprintf(phimachfilex, "%e %e \n", mach, phi);    
    dptr = dptr->nxt;
  }
  /* Delete data point list */
  dptr = start;
  while (dptr != NULL){
    tdptr = dptr;
    dptr = dptr->nxt;
    delete tdptr;
  }     
  usrclosexq("POS", 1);
  return(true);
}
//end usrmaldatconv.c

//file usrsmooth.c
/*
+ usrsmooth - Smooths shock velocity file.
     Description :
        This reads  the file shockin,  smooths it over nav  points and
        outputs to shockout.
*/
void usrsmooth(int nav)
{
  int i, j, count;
  float tin, vin;
  float avtin = 0.0, avvin = 0.0;
  char str1[10], str2[10], str3[10];
  safe_fscanf(shockin, "%d %d", &i, &j);
  safe_fscanf(shockin, "%s", str1);
  safe_fscanf(shockin, "%s %s", str2, str3);
  count = 0;
  /* Count lines */
  while (!(safe_fscanf(shockin,"%lf %lf", &tin, &vin) == EOF)){
    count++;
  }
  rewind(shockin);
  safe_fscanf(shockin, "%d %d", &i, &j);
  safe_fscanf(shockin, "%s", str1);
  safe_fscanf(shockin, "%s %s", str2, str3);
  i = 0;  
  for (j = 0; j < count; j++){
    safe_fscanf(shockin,"%lf %lf", &tin, &vin);
    i++;
    avtin += tin;
    avvin += vin;
    if (j > count - 2*nav)
      nav = 1;
    if (i >= nav){
      fprintf(shockout, "%e %e \n", avtin/float(i), avvin/float(i));
      i = 0;
      avtin = 0.0;
      avvin = 0.0;
    }
  }
  usrmyclose(shockout);
  return;
}
//end usrsmooth.c

//file usrphaseplan.c
/*
+ usrphaseplan - Generate phase plane file
     Description :
        This reads  the files shockin,  shockout .xq file  that should
        contain steady solution  .xq files. It then  outputs the phase
        plane .xg file in convin.
*/
void usrphaseplan()
{
  int i, j, count;
  float x, tv1, tv2;
  char str1[10], str2[10], str3[10];
  safe_fscanf(shockin, "%d %d", &i, &j);
  safe_fscanf(shockin, "%s", str1);
  safe_fscanf(shockin, "%s %s", str2, str3);
  safe_fscanf(shockout, "%d %d", &i, &j);
  safe_fscanf(shockout, "%s", str1);
  safe_fscanf(shockout, "%s %s", str2, str3);  
  count = 0;
  /* Count lines */
  while (!(safe_fscanf(shockin,"%lf %lf", &tv1, &tv2) == EOF)){
    count++;
  }
  rewind(shockin);
  safe_fscanf(shockin, "%d %d", &i, &j);
  safe_fscanf(shockin, "%s", str1);
  safe_fscanf(shockin, "%s %s", str2, str3);
  for (j = 0; j < count; j++){
    safe_fscanf(shockin,"%lf %lf", &x, &tv1);
    safe_fscanf(shockout,"%lf %lf", &x, &tv2);    
    fprintf(convin, "%e %e \n", tv1, tv2);
  }
  usrmyclose(shockout);
  usrmyclose(shockin);
  usrmyclose(convin);
  return;
}
//end usrphaseplan.c

//file usrsoncross.c
/*
+ usrsoncross - Computes upstream sound crossing time
     Description :
        This reads  the files shockin (mach_x.xq),  shockout (u_x.xq),
        uses  them to  compute the  upstream sound  crossing time  and
        prints it.
*/
void usrsoncross(float xmax)
{
  int i, j, count;
  float xin, machin, uin, tson, xo;
  char str1[10], str2[10], str3[10];
  safe_fscanf(shockin, "%d %d", &i, &j);
  safe_fscanf(shockin, "%s", str1);
  safe_fscanf(shockin, "%s %s", str2, str3);
  safe_fscanf(shockout, "%d %d", &i, &j);
  safe_fscanf(shockout, "%s", str1);
  safe_fscanf(shockout, "%s %s", str2, str3);  
  count = 0;
  /* Count lines */
  while (!(safe_fscanf(shockin,"%lf %lf", &xin, &machin) == EOF)){
    count++;
  }
  rewind(shockin);
  safe_fscanf(shockin, "%d %d", &i, &j);
  safe_fscanf(shockin, "%s", str1);
  safe_fscanf(shockin, "%s %s", str2, str3);
  tson = 0;
  for (i = 0; i < count; i++){
    safe_fscanf(shockin,"%lf %lf", &xin, &machin);
    safe_fscanf(shockout,"%lf %lf", &xin, &uin);
    if (xin > xmax)
      break;
    if (i > 0)  
      tson += (xin - xo)/(uin*(1.0- machin));
    xo = xin;
  }
  usrmyclose(shockin);    
  usrmyclose(shockout);
  cout << std::scientific;
  cout << std::setprecision(6);  
  cout << "Uptream sound crossing time = " << tson << '\n';
  return;
}
//end usrsoncross.c

//file usrmaxab.c
/*
+ usrmaxab - Computes max mass fraction.
     Description :
        This reads the  files shockin (iso_x) and uses it  to find max
        mass fraction and its location and then prints them,
*/
void usrmaxab()
{
  int i, j, count;
  float X, x, xin, Xin;
  char str1[10], str2[10], str3[10];
  safe_fscanf(shockin, "%d %d", &i, &j);
  safe_fscanf(shockin, "%s", str1);
  safe_fscanf(shockin, "%s %s", str2, str3);
  X = 0.0;
  x = 0.0;
  while (!(safe_fscanf(shockin,"%lf %lf", &xin, &Xin) == EOF)){
    if (Xin > X){
      X = Xin;
      x = xin;
    }    
  }
  usrmyclose(shockin);    
  cout << std::scientific;
  cout << std::setprecision(6);  
  cout << "Max X " << X << " at " << x <<'\n';
  return;
}
//end usrmaxab.c

//file usrmalconv.c
/*
+ usrmalconv - Converts Malcolm's D - kappa files.
     Description :
         This reads the file convin, converts it xq files and outputs
        outputs to dkout and dlcout.
*/
void usrmalconv()
{
  float D, kap, xloc;
  char str1[10], str2[10], str3[10];
  safe_fscanf(convin, "%s %s %s", str1, str2, str3);
  while (!(safe_fscanf(convin,"%lf %lf %lf", & kap, &D, &xloc) == EOF)){
    fprintf(dkout, "%e %e \n", kap, D);
    fprintf(dlcout, "%e %e \n", kap, xloc);
  }
  usrmyclose(dkout);
  usrmyclose(dlcout);
  return;
}
//end usrmalconv.c

//file usrtimconv.c
/*
+ usrtimconv - Converts Timmmes .dat files to xq file.
     Description :
         This reads the file convin, converts it xq files and outputs
        outputs to dkout and dlcout.
*/
void usrtimconv()
{
  bool done = false;
  int it, i;
  float tim, enuc, tv;
  char str[10];
  while (!(safe_fscanf(convin, "%d %s %lf %s %s %lf", &it, str, &tim,    
    str, str, &tv) == EOF)){//n, time, dt,  
    safe_fscanf(convin, "%s %lf", str, &tv);//1st abund
    safe_fscanf(convin, "%s %lf", str, &tv);//2nd abund
    safe_fscanf(convin, "%s %lf", str, &tv);//3rd abund
    safe_fscanf(convin, "%s %lf", str, &tv);//temp
    safe_fscanf(convin, "%s %lf", str, &tv);//den
    safe_fscanf(convin, "%s %lf", str, &enuc);//enuc
    safe_fscanf(convin, "%s %lf", str, &tv);//ye
    safe_fscanf(convin, "%s %lf", str, &tv);//sum 
    fprintf(dkout, "%e %e \n", tim, enuc);
  }
  usrmyclose(dkout);
  return;
}
//end usrtimconv.c

//file usrsetint.c
/*
+ usrsetint - Set up interpolation arrays for initial conditions.
     Description:
        Reads in  the interpolation arrays in  primitive variables for
        the  initial  state from  .xq  files  generated by  export  or
        usrmaldatconv.c.    name   is   the  prefix   in   export   or
        usrmaldatconv.    It  uses   the   same  file   names  as   in
        usropenxq.c. mode = 1 - from dns, mode = 2 from steady, mode =
        3  from  atmosphere.  It  uses  a list  of  structures  as  in
        usrmaldatconv  and readxq.  It  must have  setup the  network,
        either from a model or REACT.
*/
bool usrsetint(char monpref[], int mode, float size)
{
  int i, n, nt;
  float x, tv, pa[MAXQ];
  char str1[10], str2[10], str3[10], filename[SZNAM], tail[SZNAM];
  struct intdata *dptr, *tdptr;
  /* Check mode */
  /* Check if already done */
  if (intarray != NULL){
    /* Delete interpolation array */
    if (npoints > 0){
      for (n = 0; n < npoints - 1; n++)
	if (intarray[n] != NULL)
	  delete intarray[n];
    }
    delete [] intarray;
    intarray = NULL;
    npoints = 0;
  }
  /* Open files for reading and remove headers */
  if (mode != 3){
    for (i = 0; i < niso; i++){
      strcpy(filename, monpref);
      //    sprintf(tail, "_xa_%d_x.xq", i+1);
      sprintf(tail, "_xa_%s_x.xq", ionam[i]);          
      strcat(filename,tail);            
      if (!(xaxfile[i]=fopen(filename, "r"))){
	break;
	return(false);
      }
      safe_fscanf(xaxfile[i], "%d %d", &nt, &nt);
      safe_fscanf(xaxfile[i], "%s", str1);
      safe_fscanf(xaxfile[i], "%s %s", str2, str3);    
    }
    strcpy(filename, monpref);
    sprintf(tail, "_temp_x.xq");    
    strcat(filename,tail);
    if (!(tempfilex=fopen(filename, "r")))
      return(false);
    safe_fscanf(tempfilex, "%d %d", &nt, &nt);
    safe_fscanf(tempfilex, "%s", str1);
    safe_fscanf(tempfilex, "%s %s", str2, str3);
    strcpy(filename, monpref);
    sprintf(tail, "_rho_x.xq");
    strcat(filename,tail);
    if (!(rhofilex=fopen(filename, "r")))
      return(false);    
    safe_fscanf(rhofilex, "%d %d", &nt, &nt);
    safe_fscanf(rhofilex, "%s", str1);
    safe_fscanf(rhofilex, "%s %s", str2, str3);
    strcpy(filename, monpref);
    sprintf(tail, "_pg_x.xq");
    strcat(filename,tail);    
    if (!(pgfilex=fopen(filename, "r")))
      return(false);    
    safe_fscanf(pgfilex, "%d %d", &nt, &nt);
    safe_fscanf(pgfilex, "%s", str1);
    safe_fscanf(pgfilex, "%s %s", str2, str3);
    strcpy(filename, monpref);
    sprintf(tail, "_u_x.xq");
    strcat(filename,tail);    
    if (!(ufilex=fopen(filename, "r")))
      return(false);    
    safe_fscanf(ufilex, "%d %d", &nt, &nt);
    safe_fscanf(ufilex, "%s", str1);
    safe_fscanf(ufilex, "%s %s", str2, str3);
  }
  else{
    strcpy(filename, monpref);
    sprintf(tail, "_temp_x.xq");    
    strcat(filename,tail);
    if (!(tempfilex=fopen(filename, "r")))
      return(false);
    safe_fscanf(tempfilex, "%d %d", &nt, &nt);
    safe_fscanf(tempfilex, "%s", str1);
    safe_fscanf(tempfilex, "%s %s", str2, str3);
    strcpy(filename, monpref);
    sprintf(tail, "_rho_x.xq");
    strcat(filename,tail);
    if (!(rhofilex=fopen(filename, "r")))
      return(false);    
    safe_fscanf(rhofilex, "%d %d", &nt, &nt);
    safe_fscanf(rhofilex, "%s", str1);
    safe_fscanf(rhofilex, "%s %s", str2, str3);
  }
  /* Get true number of lpoints, xmax, xmin, ldx */
  npoints = 0;
  xmax = -big;
  xmin = big;
  while (!(safe_fscanf(rhofilex,"%lf %lf", &x, &tv) == EOF)){
    npoints++;
    xmax = max(xmax, x);
    xmin = min(xmin, x);
    /* Break if x too large */
    if ((size > 0) && (x > size))
      break;
  }
  rewind(rhofilex);
  safe_fscanf(rhofilex, "%d %d", &nt, &nt);
  safe_fscanf(rhofilex, "%s", str1);
  safe_fscanf(rhofilex, "%s %s", str2, str3);  
  /* Create and fill array of structures */
  intarray = new intdata *[npoints];  
  for (n = 0; n < npoints; n++)
    intarray[n] = new intdata;
  for (n = npoints - 1; n >= 0; n--){
    if (mode == 2)
      /* From Steady */
      nt = n;	 
    else
      /* From DNS or ATMOS*/
      nt = npoints - 1 - n;
    if (mode != 3){
      for (i = 0; i < niso; i++){
	safe_fscanf(xaxfile[i], "%lf %lf", &x, &tv);
	pa[iqal0+i+1] = tv;
      }
      safe_fscanf(ufilex, "%lf %lf", &x, &tv);
      pa[iqu0+1] = tv;
      safe_fscanf(pgfilex, "%lf %lf", &x, &tv);
      pa[iqe] = tv;
    }      
    safe_fscanf(rhofilex, "%lf %lf", &x, &tv);
    pa[iqd] = tv;
    safe_fscanf(tempfilex, "%lf %lf", &x, &tv);
    pa[iqte] = tv;
    dptr = new intdata;
    intarray[nt] = dptr;
    if (mode != 2)
      dptr->x = x;
    else
      dptr->x = xmax - x;
    for (i = 0; i < nq; i++)
      dptr->pa[i] = pa[i];
  }
  return (true);   
}
//end usrsetint.c

//file usrinter.c
/*
+ usrinter - Interpolate primitive variables.
     Description:
        Uses  the   interpolation  array  of  structures   created  by
        usrsetint.c calculate the primitive variables at xc.
*/
void usrinter(float x, float pa[])
{
  int i, nql, n1, n2;
  float fracx, tv, pl[MAXQ];
  struct intdata *dptr;
  /* Check for 2D (interpolantt is 1D) */
  nql = nq - nd + 1;
  if (x >= xmax){
    for (i = 0; i < nql; i++)
      pl[i] = intarray[npoints - 1]->pa[i];
  }
  else if (x <= xmin){
    for (i = 0; i < nql; i++)
      pl[i] = intarray[0]->pa[i];
  }
  else{
    n1 = 0;
    n2 = npoints-1;
    usrfindinter(n1, n2, x);
    dptr = intarray[n1];
    fracx = (x - intarray[n1]->x)/(intarray[n2]->x - intarray[n1]->x);
    for (i = 0; i < nql; i++){
      tv = dptr->pa[i];
      pl[i] = tv + fracx*(intarray[n2]->pa[i] - tv);
    }
  }
  if (nd == 1){
    for (i = 0; i < nq; i++){
      pa[i] = pl[i];
    }
  }
  else{
    pa[iqd] = pl[iqd];
    pa[iqu0+1] = pl[iqu0+1];
    pa[iqu0+2] =  0.0;
    for (i = iqe; i < nq; i++){
      pa[i] = pl[i-1];
    }
  }
  return;
}
//end usrinter.c

//file usrfindinter.c
/*
+ usrfindinter - Find interpolation position
     Description:
        Uses recursive bisection to find interpolation index.
*/
void usrfindinter(int &n1, int &n2, float xi)
{
  int ic;
  ic = (n1 + n2)/2;
  if (xi < intarray[ic]->x)
    n2 = ic;
  else
    n1 = ic;
  if (n2 - n1 < 2)
    return;
  else
    usrfindinter(n1, n2, xi);
}
//end usrfindinter.c

//file usrdef.c
/*
+ usrdef - set defaults.
     Description :
        Sets the user defaults at start up.
*/
void usrdef()
{
  g = 5.0/3.0;
  mon = false;
  fxp = true;
  inten = false;
  god = true;
  hll = false;
  kt = false;
  tim7 = false;
  tim13 = false;
  tim21 = false;
  timtorch = false;  
  init = false;
  doneut = 0;
  doweak = 0;
  rat = 0.0;
  curv = 0.0;
  depth = 1.0;
  cns1 = 1e-3;
  cns2 = 1e-3;  
  solev = -1;
  T_0 = 1.0e8;
  rho_0 = 1.0e8;  
  strcpy(ode, "MSBDF");
  nsub = 100;
  gslatol = 1E-8;
  gslrtol = 1E-8;
  gslhmin = -1.0;
  gslhmax = -1.0;
  cns = 0.1;
  hydrostatic = false;
  isothermal = false;
  dfloor = 0.001;
  pfloor = 0.001;
  avisp = 0.2;
  avise = 0.2;
  viscous = false;
  grav = false;
  gconst = 1.0;
  delpart = 0;
  sonposfile = NULL;
  sonwidthfile = NULL;
  dtp = -1.0;
  tp_old = 0;
  sootfoil = false;
  amplitude = 25.0;
  period = 0.05;
}
//end usrdef.c

//file usrkeps.c
/*
+ usrkeps - set k-epsilon constants.
     Description :
        Sets the k-epsilon constants.
*/
void usrkeps()
{
  ce1 = 1.4;
  ce2 = 1.94;
  cs = 0.09;
  ct = 0.225;
  cd = 0.09;
  ce = 0.0692;
  ceps = 1.0;
  kmin = 0.0;
  epsmin = 0.0;
  mtl = 1.0;
  viscous = true;
  return;
}
//end usrkeps.c

//file usroff.c
/*
+ usroff - set ofsets.
     Description :
        Sets the offsets for the variables at start up.
*/
bool usroff()
{
  iqd = 0;
  iqu0 = 0;
  iqe = 1 + nd;
  /* k- epsilon */
  if (keps){
    iqk = iqe + 1;
    iqeps = iqk + 1;
    iqal0 = iqeps;
  }
  else{
    iqk = -1;
    iqal0 = iqe;
  }
  nq = iqal0 + nal + 1;
  if (inten){
    iqie = nq;
    nq = iqie + 1;
  }
  iqte = nq;
  nq = iqte + 1;
  if (grav){
    /* Gravity */
    iqpot = nq;
    nq++;
  }
  else
    iqpot = 0;
  if (sootfoil) {
    iqsf = nq;
    nq = iqsf + 1;
  }
  if (nq > MAXQ){
    cout << "Error : nq = " << nq << " > MAXQ = " << MAXQ <<
      " Change MAXQ in struct.h" << '\n';
    return(false);
  }
  /* Check for reaction */
  if (react && (niso > nal)){
    cout << "Error : niso = " << niso << " > nal = " << nal << '\n';
    return(false);    
  }  
  return(true);
}
//end file usroff.c

//file usrvel.c
/*
+ usrvel - Returns velocity.
     Description :
        This returns the id direction  velocity in the cell pointed to
        by cptr. Used for plotting velocity vectors and streamlines.
*/
float usrvel(int id, struct cell *cptr)
{
  return(cptr->qa[iqu0+id]/cptr->qa[iqd]);
}
//end usrvel.c

//file usrbfield.c
/*
+ usrbfield - Returns velocity.
     Description :
        This returns the id direction  velocity in the cell pointed to
        by   cptr.   Used    for   plotting   velocity   vectors   and
        streamlines. This is just a copy of usrvel.
*/
float usrbfield(int id, struct cell *cptr)
{
  return(cptr->qa[iqu0+id]/cptr->qa[iqd]);
}
//end usrbfield.c

//file usrgfield.c
/*
+ usrgfield - Returns gravitational acceleration.
     Description :

        This  returns the id  direction gravitational  acceleration in
        the cell  pointed to by cptr.  Used  for plotting acceleration
        vectors and streamlines.
*/
float usrgfield(int id, struct cell *cptr)
{
  return(cptr->pa[iqu0+id]);
}
//end usrgfield.c

//file usrga.c
/*
+ usrga - Returns gravitational acceleration.
     Description :
        This calculates the id direction gravitational acceleration in
        the cell  pointed to by cptr,  with cells left,  right on left
        and right.  The result is  put in cptr->pa[iqu0+id].
*/
void usrga(int id, struct cell *cptr, struct cell *left,
  struct cell *right)
{
  int count;
  float ga;
  /* Initialise */
  count = 0;
  ga = 0.0;
  if (left != NULL){
    ga += cptr->qa[iqpot] - left->qa[iqpot];
    count++;
  }
  if (right != NULL){
    ga +=  right->qa[nq-1] - cptr->qa[nq-1];
    count++;
  }
  if (count == 0)
    cptr->pa[iqu0+id] = 0.0;
  else{
    ga *= - gconst/((float)(count)*
      (cptr->rjoin[id-1]->xj[id-1] - cptr->ljoin[id-1]->xj[id-1]));
    cptr->pa[iqu0+id] = ga;
  }
}
//end usrga.c

/*-----sootfoil Addition------*/
/* Creates sootfoil using the maximum pressure */
void usrSootFoil(cell *cptr) {
  float pres;
  pres = cptr->pa[iqe]; // uses 2nd order step pressure
  if (pres > cptr->qa[iqsf]) {
    cptr->qa[iqsf] = pres;
    cptr->pa[iqsf] = pres;
  }

  return;
}
/*----------------------------*/


//file usrprompt.c
/*
+ usrprompt - set prompt.
     Description:
        This sets the prompt at start up.
*/
void usrprompt(char prompt[])
{
  strcpy(prompt, "MG_TORCH>");
}
//end usrprompt.c

//file usrquit.c
/*
+ usrquit - clean up user dynamic arrays.
     Description:
        This  deletes any dynamic arrays defined by the user.
*/
void usrquit()
{
  return;
}
//end usrquit.c

//file usrpath.c
/*
+ usrpath - set path.
     Description:
        This sets the path of the model file directory at start up.
*/
void usrpath(char path[])
{
  strcpy(path, "mg_g_data");
}
//end usrpath.c

//file usrdmdt.c
/*
+ usrdmdt - Rate of change of particle mass.
     Description :
        This returns  the rate  of change of  particle mass  with mass
        pmass with the fluid state  in pa[]. It is used for lagrangian
        particles.
*/
float usrdmdt(float pmass, float pa[])
{
  float tv;
  /* Zero for now */
  tv = 0.0;
  return(tv);
}
//end usrdmdt.c

//file usrmkpart.c
/*
+ usrmkpart - Make a lagrangian particle.
     Description :
        This  creates a  particle  at xp,  assigns  it to  a cell  and
        returns a  pointer to the  particle. The particle  position is
        set, but  none of  the other properties.  The particle  is not
        made if xp is not on the processor or in the domain.
*/
struct part *usrmkpart(float xp[])
{
  int id, ilev, iedge;
  float dx;
  struct cell *cptr;
  struct part *pptr;
  /* Check whether close to nd cell edge */
  dx = lev[nlevs-1]->dxbar[nd-1];
  iedge = (int)((xp[nd-1] - xl[nd-1])/dx);
  if (mabs(xp[nd-1] - iedge*dx - xl[nd-1]) < 1.0e-3*dx)
    /* Near cell cell */
    xp[nd-1] = iedge*dx + xl[nd-1] + 1.0e-3*dx;
  /* Check right boundary */  
  if (xr[nd-1] - xp[nd-1] < 1.0e-3*dx)
    /* Near right boundary */
    xp[nd-1] = xr[nd-1] - 1.0e-3*dx;
  /* Find cell */
  cptr = find_fcell(ilev, xp);
  if (cptr == NULL) 
    /* Not on processor - do nothing */
    return(NULL);
  else if ((limpe > 1) && (cptr->halo)) 
    /* Halo cell - do nothing */
    return(NULL);
  /* Move particle off nd cell boundary if necessary */
  if (mabs(xp[nd-1] - cptr->ljoin[nd-1]->xj[nd-1]) < 
      1.0e-3*lev[ilev]->dxbar[nd-1]){
    xp[nd-1] = cptr->ljoin[nd-1]->xj[nd-1] + 
      1.0e-3*lev[ilev]->dxbar[nd-1];
  }
  if (mabs(xp[nd-1] - cptr->rjoin[nd-1]->xj[nd-1]) < 
      1.0e-3*lev[ilev]->dxbar[nd-1]){
    xp[nd-1] = cptr->rjoin[nd-1]->xj[nd-1] - 
      1.0e-3*lev[ilev]->dxbar[nd-1];
  }
  /* Create particle */
  pptr = new part;
  /* Initial NULL pointers (safety) */
  pptr->prev = NULL;
  pptr->nxt = NULL;
  /* Assign to cell */
  pptr->cptr = cptr;
  if (cptr->pptr == NULL)
    cptr->pptr = pptr;
  else{
    cptr->pptr->prev = pptr;
    pptr->nxt = cptr->pptr;
    cptr->pptr = pptr;
  }
  /* Set particle time, position and velocty */
  pptr->t = t;
  for (id = 1; id <= nd; id++){
    pptr->xp[id-1] = xp[id-1];
    pptr->vp[id-1] = cptr->qa[iqu0+id]/cptr->qa[iqd];    
  }
  return(pptr);
}
//end usrmkpart.c

//file usrprate.c
/*
+ usrprate - Rate of change of lagrangian particle variables.
     Description :
        This  should  return  the  time derivatives  of  the  particle
        variables. It is given
          pprt    -- particle  pointer
          pa[]    -- fluid state at particle position (pptr->xp)
        It should return
          dvp[]   -- particle acceleration
          den     -- rate of change of particle energy
          dm      -- rate of change of particle mass
*/
void usrprate(struct part *pptr, float dt, float pa[], float dvp[],
 float &den, float &dm)
{
  float ql, qr;
  int id;
  struct cell *cptr;
  struct mcell *mcptr, *left, *right;
  /* Zero for now */
  for (id = 0; id < nd; id++)
    dvp[id] = 0.0;
  /* Energy and mass source */
  den = 0.0;
  dm = 0.0;
  /* Particle has fluid velocity (1D) */
  dvp[0] = (pa[iqu0+1] - pptr->vp[0])/dt;
  /*
    cout << "usrate vp " << pptr->vp[0] << " dt " << dt << " u1 " <<
      pa[iqu0+1] << " dvp " << dvp[0] << '\n';
  */
  /* Set shock time */
  if (pptr->props[0] < 0.0)
    if (pa[iqd] > 1.1*rho_0)
      pptr->props[0] = t;
  /* Set density, velocity and temperature */
  pptr->props[1] = pa[iqd];
  pptr->props[2] = pa[iqu0+1];
  pptr->props[3] = pa[iqte];
  if (grav){
    /* Calculate left/right gravitational potentials */
    cptr = pptr->cptr;
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
    /* Apply gravitational acceleration */
      dvp[id-1] -= 0.5*gconst*(qr - ql)/
	(cptr->rjoin[id-1]->xj[id-1] - cptr->ljoin[id-1]->xj[id-1]);
    }
  } 
  return;
}
//end usrprate.c

//file usrpartpop.c
/*
+ usrpartpop - Update lagrangian particle population.
     Description :
        This  should  be  used   to  update  the  lagrangian  particle
        population. tcur is the  current model time.
*/
void usrpartpop(int ilev, float tcur)
{
  /* Nothing for now */
  return;
}
//end usrpartpop.c

//file usrpartps.c
/* 
+ usrpartps - construct plot scalar for particle.
     Description :
        s returns  the value  of a particle  property in  the particle
        pointed to  by pptr. The  choice of property is  determined by
        pitem true if the plot pitem is found, false otherwise.
*/
bool usrpartps(char pitem[], struct part *pptr, float &s)
{
  bool found;
  if (strcmp(item, "PPARTMASS") == 0){
    s = pptr->mass;
    found = true;
  }
  else if (strcmp(item, "PPARTRHO") == 0){
    s = pptr->props[0];
    found = true;
  }
  else if (strcmp(item, "PPARTEN") == 0){
    s = pptr->en;
    found = true;
  }
  else if (strcmp(item, "PPARTU1") == 0){
    s = pptr->vp[0];
    found = true;
  }
  else if (strcmp(item, "PPARTU2") == 0){
    if (nd < 2)
      found = false;
    else{
      s = pptr->vp[1];
      found = true;
    }
  }
  else if (strcmp(item, "PPARTU3") == 0){
    if (nd < 3)
      found = false;
    else{
      s = pptr->vp[2];
      found = true;
    }
  }
  else if (strcmp(item, "PPARTUT") == 0){
    int id;
    float tv1;
    s = 0.0;
    for (id = 1; id <= nd; id++){
      tv1 = pptr->vp[id-1];
      s += tv1*tv1;
    }
    s = sqrt(s);
    found = true;
  }
  else if (strcmp(item, "PPARTT") == 0){
    s = pptr->t;
    found = true;
  }
  return(found);
}
//end usrpartps.c

//file usrBE.c
/* 
+ usrBE - binding energy routine
     Description :
        Returns  the  value of  the  total  binding energy  given  the
        primitive variables.
*/
float usrBE(float pa[])
{
  float sum, sumX;
  /* Zero energy and mass fraction sum */
  sum = 0.0;
  sumX = 0.0;
  for (int i = 0; i < niso; i++){
    sum += pa[iqal0+i+1]*Qi[i];
    sumX += pa[iqal0+i+1]*Ai[i];
  }
  return sum;
}
//end usrBE.c

//file usrY2X.c
/* 
+ usrY2X - mass fractions
     Description :
        Returns array with mass  fractions, given nuclei abundance for
        each species.
*/
void usrY2X(float pa[], float X[])
{
  float sum = 0.0;
  for (int i = 0; i < niso; i++){
    X[i] = pa[iqal0+i+1]*Ai[i];
    sum += X[i];
  }
  return;
}
//end usrY2X.c

//file usrX2Y.c
/* 
+ usrX2Y - nuclei abundances
     Description :
        Returns  nuclei  abundance  for  each species  (Y)  from  mass
        fraction (X).
*/
void usrX2Y(float pa[], float X[])
{
  for (int i = 0; i < niso; i++)
    pa[iqal0+i+1] = X[i]/Ai[i];
  return;
}
//end usrX2Y.c

//file usrYE.c
/* 
+ usrYE - Electrons per nucleon (Ye)
     Description :
        Returns  Ye given  temperature and  nuclei abundance  for each
        species (Y). Requires Yi for  each species (pa) Can be modified
        to return Abar and Zbar.
*/
float usrYE(float pa[])
{
  float Abar, Zbar, sumA, sumZ, Ye;
  sumA = 0.0;
  sumZ = 0.0;
  for (int i = 0; i < niso; i++) {
    sumA += pa[iqal0+i+1];
    sumZ += Zi[i]*pa[iqal0+i+1];
  }
  Abar = 1.0/sumA;
  Zbar = Abar*sumZ;
  Ye = Zbar/Abar;
  return(Ye);
}
//end usrYE.c

//file usrAbar.c
/* 
+ usrAbar - Mean atomic mass
     Description :
        Returns mean atomic mass, Abar,  used in ion().
*/
float usrAbar(float pa[])
{
  float Abar, sumA;
  /* Zero */
  sumA = 0.0;
  for (int i = 0; i < niso; i++)
    sumA += pa[iqal0+i+1];
  Abar = 1.0/sumA;
  return(Abar);
}
//end usrAbar.c

//file usrYAZ.c
/* 
+ usrYAZ - Alternative electron mass fraction.
     Description :
        Returns  Ye, Abar, Zbar, Z2bar  given nuclei abundance  for each
        species (pa[]).
*/
void usrYAZ(float pa[], float &Ye, float &Abar, float &Zbar,
  float &Z2bar)
{
  float sumA, sumZ, sumZ2;
  sumA = 0.0;
  sumZ = 0.0;
  sumZ2 = 0.0;
  for (int i = 0; i < niso; i++){
    sumA += pa[iqal0+i+1];
    sumZ += Zi[i]*pa[iqal0+i+1];
    sumZ2 += Zi[i]*Zi[i]*pa[iqal0+i+1];
  }
  Abar = 1.0/sumA;
  Zbar = Abar*sumZ;
  Z2bar = Abar*sumZ2;
  Ye = Zbar/Abar;
  return;
}
//end usrYAZ.c

//file usrein.c
/*
+ usrein - specific internal energy
     Description :
     Returns  the  specific  internal energy  given  the  temperature,
     density,   Abar    and   Ye.   This   is    the   total   energy:
     electron-positron + radiation + ion.
*/
float usrein(float temp, float rho, float Abar, float Ye)
{
  float ein, he[12];
  ein = 0.0;
  usrhelm(temp, rho, Ye, he);
  ein += he[1];
  usrion(temp, rho, Abar, he);
  ein += he[1];
  usrrad(temp, rho, he);
  ein += he[1];  
  return(ein);
}
//end usrein.c

//file usrcv.c
/*
+ usrcv - specific heat at constant volume (de/dT)
     Description :
     Returns the specific heat  energy given the temperature, density,
     Abar and Ye.  This is  the total: electron-positron + radiation +
     ion.  */
float usrcv(float pa[])
{
  float rho, temp, he[12], cv;
  rho = pa[iqd];
  temp = pa[iqte];
  usrhelm(temp, rho, usrYE(pa), he);
  cv = he[6];
  usrion(temp, rho, usrAbar(pa), he);
  cv += he[6];
  usrrad(temp, rho, he);
  cv += he[6];  
  return(cv);
}
//end usrcv.c

//file usrpg.c
/*
+ usrpg - Total pressure
     Description :
     Returns the  total pressure given the  temperature, density, Abar
     and  Ye.   This  is   the  total  pressure:  electron-positron  +
     radiation + ion.
*/
float usrpg(float temp, float rho, float Abar, float Ye)
{
  float pg, he[12];
  pg = 0.0;
  usrhelm(temp, rho, Ye, he);
  pg += he[0];
  usrion(temp, rho, Abar, he);
  pg += he[0];
  usrrad(temp, rho, he);
  pg += he[0];  
  return(pg);
}
//end usrpg.c

//file usrion.c
/*
+ usrion - Ion thermodynamics
     Description :
        Calculates ion thermodynamics
input :
temp = temperature(in K)
rho = density (in g cm~3)
Abar = average atomic mass of nucleons
output :
he[0] = pressure
he[1] = specific internal energy
he[2] = specific entropy
he[3] = partial derivative of pressure with density
he[4] = partial derivative of pressure with temperature
he[5] = partial derivative of energy with density
he[6] = partial derivative of energy with temperature
he[7] = partial derivative of entropy with density
he[8] = partial derivative of entropy with temperature
he[9] = adiabatic sound speed
all in c.g.s
*/
void usrion(float temp, float rho, float Abar, float he[])
{
  float tv1, tv2, pres, etot, dpdt, dpdd, dedt, dedd, stot, dsdt, dsdd;
  tv1 = k_b*N_A/Abar;
  tv2 = tv1*rho;
  pres = tv2*temp;
  dpdt = tv2;
  dpdd = tv1*temp;
  etot = 1.5*tv1*temp;
  dedt = 1.5*tv1;   
  dedd = 0.0;
  stot = 0.0; /* Warning do not use */
  dsdt = 1.5*tv1/temp; /* From C_v */
  dsdd = - tv1/rho;    /* From Maxwell relation */
  he[0] = pres;
  he[1] = etot;
  he[2] = stot;
  he[3] = dpdd;
  he[4] = dpdt;
  he[5] = dedd;
  he[6] = dedt;
  he[7] = dsdd;
  he[8] = dsdt;
  return;
}
//end usrion.c

//file usrrad.c
/*
+ usrrad - Radiation thermodynamics
     Description :
        Calculates radiation thermodynamics
input :
temp = temperature(in K)
rho = density (in g cm~3)
output :
he[0] = pressure
he[1] = specific internal energy
he[2] = specific entropy
he[3] = partial derivative of pressure with density
he[4] = partial derivative of pressure with temperature
he[5] = partial derivative of energy with density
he[6] = partial derivative of energy with temperature
he[7] = partial derivative of entropy with density
he[8] = partial derivative of entropy with temperature
he[9] = adiabatic sound speed
all in c.g.s
*/
void usrrad(float temp, float rho, float he[])
{
  float tv1, tv2, pres, etot, dpdt, dpdd, dedt, dedd, stot, dsdt, dsdd;
  tv1 = temp*temp;
  tv1 = tv1*tv1;
  tv2 = stefan*tv1;
  pres = tv2/3.0;
  dpdt = 4.0*pres/temp;
  dpdd = 0.0;
  etot = tv2/rho;
  dedt = 4.0*etot/temp;   
  dedd = -etot/rho;
  stot = 4.0*etot/(3.0*temp);
  dsdt = 3.0*stot/temp;
  dsdd = - stot/rho;
  he[0] = pres;
  he[1] = etot;
  he[2] = stot;
  he[3] = dpdd;
  he[4] = dpdt;
  he[5] = dedd;
  he[6] = dedt;
  he[7] = dsdd;
  he[8] = dsdt;    
  return;
}
//end usrrad.c

//file usrein2te.c
/*
+ usrein2te  - Temperature from specific internal  energy
     Description :
     Returns the temperature given the total specific internal energy,
     density,   Abar   and   Ye.     This   is   the   total   energy:
     electron-positron + radiation + ion. It uses a secant method with
     tempi as the initial guess.
*/
float usrein2te(float tempi, float rho, float ein, float Abar, float Ye)
{
  int nit;
  float tol = 1e-6, temp, fo, f, dt, dto;
  nit = 0;
  fo = usrein(tempi, rho, Abar, Ye) - ein;
  dto = 1.0e-4*tempi;
  temp = tempi + dto;
  do {
    f = usrein(temp, rho, Abar, Ye) - ein;
    dt = - f*dto/(f - fo);
    temp += dt;
    nit++;
    if (nit > 100){
      cout << "Convergence failure in usrpg2te" << '\n';
      break;
    }    
    fo = f;
    dto = dt;
  }
  while (fabs(dt) > tol*temp);
  return(temp);
}
//end usrein2te.c

//file usrpg2te.c
/*
+ usrpg2te  - Temperature from pressure
     Description :
     Returns the temperature given  the total total pressure, density,
     Abar  and Ye.   This  is the  total  energy: electron-positron  +
     radiation  + ion.  It  uses a  secant method  with  tempi as  the
     initial guess.
*/
float usrpg2te(float tempi, float rho, float pg, float Abar, float Ye)
{
  int nit;
  float tol = 1e-8, temp, fo, f, dt, dto;
  nit = 0;
  fo = usrpg(tempi, rho, Abar, Ye) - pg;  
  dto = 1.0e-4*tempi;
  temp = tempi + dto;
  do {
    f = (usrpg(temp, rho, Abar, Ye) - pg);
    dt = - f*dto/(f - fo);
    if (fabs(dt) > 0.1*temp)
      dt *= 0.1*fabs(temp/dt);    
    temp += dt;
    nit++;
    if (nit > 50){
      //      cout << "Convergence failure in usrpg2te" << '\n';
      break;
    }    
    fo = f;
    dto = dt;
  }
  while (fabs(dt) > tol*temp);
  return(temp);
}
//end usrpg2te.c

//file usrnadyo.c
/*
+ usrnadyo - Nadyozhin approx (apjs106)
     Description :
        Implements Nayoshin approximation
*/
void usrnadyo(float temp, float rho, float Ye, float &pres,
  float &ein, float &sound)
{
  const float Q = 4.82744e-8, B = 2.68192e13, a1 = 45.0/16.0,
    a2 = 2865.0/512.0, a3 = 15.0/4.0, a4 = 165.0/32, a5 = 5.0/4.0,
    a6 = 9.0/16.0, a7 = 15.0/16.0, a8 = 955.0/512.0, a9 = 21.0/4.0,
    a10 = 405.0/32;
  float al, tv1, tv2, tv3;
  float t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15,
    t17, t18, t19, t20, t21, t22, t23, t24, t25, t27, t29, t31, t33,
    t38, t39, t40, t45, t49, t63;
  al = 5.929896509e9/temp;
  t4 = pow(al,1.5);
  t5 = 1/al;
  t7 = al*al;
  t8 = 1/t7;
  t13 = rho*rho;
  t19 = exp(-2.0*al);
  /* Pressure */
  pres = k_b*temp*rho*N_A*(1.0+Q*rho*t4*(1.0-a1*t5+a2*t8)+
    B/t13/t7/al*t19*(1.0+a3*t5+a4*t8));
  /*
  std::cout << "pres = " << pres << '\n';
  std::cout << "k_b = " << k_b << '\n';
  std::cout << "N_A = " << N_A << '\n';
  std::cout << "Q = " << Q << '\n';
  std::cout << "t4 = " << t4 << '\n';
  std::cout << "t5 = " << t5 << '\n';
  std::cout << "t8 = " << t8 << '\n';
  std::cout << "a1 = " << a1 << '\n';
  std::cout << "a2 = " << a2 << '\n';
  std::cout << "a3 = " << a3 << '\n';
  std::cout << "a4 = " << a4 << '\n';
  std::cout << "B = " << B << '\n';
  std::cout << "t7 = " << t7 << '\n';
  std::cout << "t13 = " << t13 << '\n';
  std::cout << "t19 = " << t19 << '\n';
  std::cout << "al = " << al << '\n';
  */
  /* Internal Energy */
  t2 = al;
  t3 = 1/t2;
  t11 = pow(t2,1.5);
  t18 = rho*rho;
  t21 = t2*t2;
  t24 = exp(-2.0*t2);
  ein = 1.5*k_b*temp*N_A*(1.0+a5*(0.1E1-(1.0-a6*t3)*t3)*t3+
    Q*rho*t11*(1.0-(a7+a8*t3)*t3)+
    2.0*B/3.0/t18/t21*t24*(1.0+(a9+a10*t3)*t3));
  float  dpdd, dpdt, dedd, dedt, dedd_p, dedp_d;
  /* (dp/dd)_T */
  t1 = k_b*temp;
  t3 = pow(al,1.5);
  t4 = 1/al;
  t6 = al*al;
  t7 = 1/t6;
  t9 = 1.0-a1*t4+a2*t7;
  t12 = rho*rho;
  t18 = exp(-2.0*al);
  t23 = 1/t6/al*t18*(1.0+a3*t4+a4*t7);
  dpdd = t1*N_A*(1.0+Q*rho*t3*t9+B/t12*t23)+t1*rho*N_A*(Q*t3*t9
    -2.0*B/t12/rho*t23);
  /* (dp/dT)_d */
  t2 = Q*rho;
  t3 = al;
  t4 = pow(t3,1.5);
  t5 = 1/t3;
  t7 = t3*t3;
  t8 = 1/t7;
  t10 = 1.0-a1*t5+a2*t8;
  t13 = rho*rho;
  t15 = B/t13;
  t17 = 1/t7/t3;
  t19 = exp(-2.0*t3);
  t20 = t17*t19;
  t23 = 1.0+a3*t5+a4*t8;
  t31 = pow(t3,0.5);
  t33 = -al/temp;
  t45 = t7*t7;
  t49 = t19*t23*t33;
  dpdt = k_b*rho*N_A*(1.0+t2*t4*t10+t15*t20*t23)+
    k_b*temp*rho*N_A*(0.15E1*t2*t31*t10*t33+t2*t4*(a1*t8*t33-
    2.0*a2*t17*t33)-3.0*t15/t45*t49-2.0*t15*t17*t49+t15*t20*(-a3*t8*t33
    -2.0*a4*t17*t33));
  /* (de/dd)_T */
  t2 = al;
  t3 = pow(t2,1.5);
  t5 = 1/t2;
  t11 = rho*rho;
  t15 = t2*t2;
  t18 = exp(-2.0*t2);
  dedd = 1.5*k_b*temp*N_A*(Q*t3*(1.0-(a7+a8*t5)*t5)-4.0/3.0*B/t11/
    rho/t15*t18*(1.0+(a9+a10*t5)*t5));
  /* (de/dT)_d */
  t2 = al;
  t3 = 1/t2;
  t5 = 1.0-a6*t3;
  t8 = a5*(1.0-t5*t3);
  t10 = Q*rho;
  t11 = pow(t2,1.5);
  t13 = a7+a8*t3;
  t15 = 1.0-t13*t3;
  t18 = rho*rho;
  t20 = B/t18;
  t21 = t2*t2;
  t22 = 1/t21;
  t24 = exp(-2.0*t2);
  t25 = t22*t24;
  t27 = a9+a10*t3;
  t29 = 1.0+t27*t3;
  t38 = 1/t21/t2;
  t40 = -al/temp;
  t49 = sqrt(t2);
  t63 = t24*t29*t40;
  dedt = 1.5*k_b*N_A*(1.0+t8*t3+t10*t11*t15+2.0/3.0*t20*t25*t29)+
    1.5*k_b*temp*N_A*(a5*(-a6*t38*t40+t5*t22*t40)*t3-t8*t22*t40+
    1.5*t10*t49*t15*t40+t10*t11*(a8*t38*t40+t13*t22*t40)-
    4.0/3.0*t20*t38*t63-4.0/3.0*t20*t22*t63+2.0/3.0*t20*t25*
    (-a10*t38*t40-t27*t22*t40));
  /* (de/dd)_p */
  dedd_p = dedd - dedt*dpdd/dpdt;
  /* (de/dp)_d */
  dedp_d = dedt/dpdt;
  /* sound speed */
  sound = sqrt((pres/(rho*rho) - dedd_p)/dedp_d);
  return;
}
//end usrnadyo.c

//file usrshock.c
/*
+ usrshock - Post shock state
     Description :
        Returns the downstream state rho_l,  T_l, v_l, p_l for a shock
        with speed v_r  and upstream state rho_r, T_r,  v_r, p_r. This
        is in the shock frame. Gets the abundances from pa.
*/
void usrshock(float pa[], float rho_r, float T_r, float v_r,
  float &rho_l, float &T_l, float &v_l, float &p_l)
{
  int nit;
  float tol = 1e-8, Abar, Ye, p_r, Q, Q2, P, E, fo, f, drho, drhoo;
  /* Get Abar and Ye */
  Abar = usrAbar(pa);
  Ye = usrYE(pa);
  /* Upstream pressure */
  p_r = usrpg(T_r, rho_r, Abar, Ye);
  /* Mass flux */
  Q = rho_r*v_r;
  Q2 = Q*Q;
  /* Momentum flux */
  P = v_r*Q + p_r;
  /* Energy (Bernoulli) */
  E = usrein(T_r, rho_r, Abar, Ye) + 0.5*v_r*v_r + p_r/rho_r;
  nit = 0;
  rho_l = 2.0*rho_r;
  p_l = P - Q2/rho_l;
  v_l = Q/rho_l;
  T_l = usrpg2te(T_r, rho_l, p_l, Abar, Ye); 
  fo = usrein(T_l, rho_l, Abar, Ye) + 0.5*v_l*v_l + p_l/rho_l - E;
  drhoo = 1.0e-4*rho_l;
  rho_l += drhoo;
  do {
    p_l = P - Q2/rho_l;
    v_l = Q/rho_l;    
    T_l = usrpg2te(T_l, rho_l, p_l, Abar, Ye);
    f = usrein(T_l, rho_l, Abar, Ye) + 0.5*v_l*v_l + p_l/rho_l - E;    
    drho = - f*drhoo/(f - fo);
    if (fabs(drho) > 0.05*rho_l)
      drho *= 0.05*fabs(rho_l/drho);
    rho_l += drho;
    nit++;
    if (nit > 100){
      cout << "Convergence failure in usrshock" << '\n';
      break;
    }
    fo = f;
    drhoo = drho;
  }
  while(fabs(drho) > tol*rho_l);
  p_l = P - Q2/rho_l;
  v_l = Q/rho_l;    
  T_l = usrpg2te(T_l, rho_l, p_l, Abar, Ye);
}  
//end usrshock.c

//file usrsetionam.c
/*
  + usrsetionam - Set ion names.
  Description :     
  Uses aion, zion to set the names in ionam.
*/	
void usrsetionam()
{
  int i, ia, iz;
  for (i = 0; i < niso; i++){
    if (zion[i] == 0){
      strcpy(ionam[i], "neut");
      ineut = i;
    }
    else if (zion[i] == 1){
      if (aion[i] == 1){
	strcpy(ionam[i], "h1");
	ih1 = i;
      }
      else if (aion[i] == 2){
	strcpy(ionam[i], "h2");
	ih2 = i;
      }
      else if (aion[i] == 3){
	strcpy(ionam[i], "h3");
	ih3 = i;
      }
    }
    else if (zion[i] == 2){
      if (aion[i] == 3){
	strcpy(ionam[i], "he3");
	ihe3 = i;
      }
      else if (aion[i] == 4){
	strcpy(ionam[i], "he4");
	ihe4 = i;
      }
    }      
    else if (zion[i] == 3){
      if (aion[i] == 6){
	strcpy(ionam[i], "li6");
	ili6 = i;
      }
      else if (aion[i] == 7){
	strcpy(ionam[i], "li7");
	ili7 = i;
      }
      else if (aion[i] == 8){
	strcpy(ionam[i], "li8");
	ili8 = i;
      }
    }
    else if (zion[i] == 4){
      if (aion[i] == 7){
	strcpy(ionam[i], "be7");
	ibe7 = i;
      }      
      else if (aion[i] == 8){
	strcpy(ionam[i], "be8");
	ibe8 = i;
      }      
      else if (aion[i] == 9){
	strcpy(ionam[i], "be9");
	ibe9 = i;
      }
      else if (aion[i] == 10){
	strcpy(ionam[i], "be10");
	ibe10 = i;
      }
    }
    else if (zion[i] == 5){
      if (aion[i] == 8){
	strcpy(ionam[i], "b8");
	ib8 = i;
      }
      else if (aion[i] == 10){
	strcpy(ionam[i], "b10");
	ib10 = i;
      }      
      else if (aion[i] == 11){
	strcpy(ionam[i], "b11");
	ib11 = i;
      }
    }
    else if (zion[i] == 6){
      if (aion[i] == 11){
	strcpy(ionam[i], "c11");
	ic11 = i;
      }      
      if (aion[i] == 12){
	strcpy(ionam[i], "c12");
	ic12 = i;
      }
      else if (aion[i] == 13){
	strcpy(ionam[i], "c13");
	ic13 = i;
      }
      else if (aion[i] == 14){
	strcpy(ionam[i], "c14");
	ic14 = i;
      }      
    }
    else if (zion[i] == 7){
      if (aion[i] == 13){
	strcpy(ionam[i], "n13");
	in13 = i;
      }      
      else if (aion[i] == 14){
	strcpy(ionam[i], "n14");
	in14 = i;
      }
      else if (aion[i] == 15){
	strcpy(ionam[i], "n15");
	in15 = i;
      }
    }
    else if (zion[i] == 8){
      if (aion[i] == 14){
	strcpy(ionam[i], "o14");
	io14 = i;
      }            
      else if (aion[i] == 15){
	strcpy(ionam[i], "o15");
	io15 = i;
      }      
      else if (aion[i] == 16){
	strcpy(ionam[i], "o16");
	io16 = i;
      }
      else if (aion[i] == 17){
	strcpy(ionam[i], "o17");
	io17 = i;
      }
      else if (aion[i] == 18){
	strcpy(ionam[i], "o18");
	io18 = i;
      }
      else if (aion[i] == 19){
	strcpy(ionam[i], "o19");
	io19 = i;
      }      
    }
    else if (zion[i] == 9){
      if (aion[i] == 17){
	strcpy(ionam[i], "f17");
	if17 = i;
      }
      else if (aion[i] == 18){
	strcpy(ionam[i], "f18");
	if18 = i;
      }
      else if (aion[i] == 19){
	strcpy(ionam[i], "f19");
	if19 = i;
      }                   
    }
    else if (zion[i] == 10){
      if (aion[i] == 18){
	strcpy(ionam[i], "ne18");
	ine18 = i;
      }
      else if (aion[i] == 19){
	strcpy(ionam[i], "ne19");
	ine19 = i;
      }
      else if (aion[i] == 20){
	strcpy(ionam[i], "ne20");
	ine20 = i;
      }
      else if (aion[i] == 21){
	strcpy(ionam[i], "ne21");
	ine21 = i;
      }      
      else if (aion[i] == 22){
	strcpy(ionam[i], "ne22");
	ine22 = i;
      }
      else if (aion[i] == 23){
	strcpy(ionam[i], "ne23");
	ine23 = i;
      }      
    }
    else if (zion[i] == 11){
      if (aion[i] == 18){
	strcpy(ionam[i], "na18");
	ina18= i;
      }
      else if (aion[i] == 19){
	strcpy(ionam[i], "na19");
	ina19= i;
      }
      else if (aion[i] == 20){
	strcpy(ionam[i], "na20");
	ina20= i;
      }      
      else if (aion[i] == 21){
	strcpy(ionam[i], "na21");
	ina21= i;
      }
      else if (aion[i] == 22){
	strcpy(ionam[i], "na22");
	ina22= i;
      }
      else if (aion[i] == 23){
	strcpy(ionam[i], "na23");
	ina23= i;
      }
      else if (aion[i] == 24){
	strcpy(ionam[i], "na24");
	ina24= i;
      }
    }
    else if (zion[i] == 12){
      if (aion[i] == 22){
	strcpy(ionam[i], "mg22");
	img22= i;
      }
      else if (aion[i] == 23){
	strcpy(ionam[i], "mg23");
	img23= i;
      }
      else if (aion[i] == 24){
	strcpy(ionam[i], "mg24");
	img24= i;
      }
      else if (aion[i] == 25){
	strcpy(ionam[i], "mg25");
	img25= i;
      }
      else if (aion[i] == 26){
	strcpy(ionam[i], "mg26");
	img26= i;
      }
      else if (aion[i] == 27){
	strcpy(ionam[i], "mg27");
	img27= i;
      }      
    }
    else if (zion[i] == 13){
      if (aion[i] == 25){
	strcpy(ionam[i], "al25");
	ial25= i;
      }
      else if (aion[i] == 26){
	strcpy(ionam[i], "al26");
	ial26= i;
      }
      else if (aion[i] == 27){
	strcpy(ionam[i], "al27");
	ial27= i;
      }
      else if (aion[i] == 28){
	strcpy(ionam[i], "al28");
	ial27= i;
      }      
    }
    else if (zion[i] == 14){
      if (aion[i] == 27){
	strcpy(ionam[i], "si27");
	isi27= i;
      }
      else if (aion[i] == 28){
	strcpy(ionam[i], "si28");
	isi28= i;
      }
      else if (aion[i] == 29){
	strcpy(ionam[i], "si29");
	isi29= i;
      }
      else if (aion[i] == 30){
	strcpy(ionam[i], "si30");
	isi30= i;
      }
      else if (aion[i] == 31){
	strcpy(ionam[i], "si31");
	isi31= i;
      }
      else if (aion[i] == 32){
	strcpy(ionam[i], "si32");
	isi32= i;
      }      
    }
    else if (zion[i] == 15){
      if (aion[i] == 29){
	strcpy(ionam[i], "p29");
	ip29= i;
      }
      else if (aion[i] == 30){
	strcpy(ionam[i], "p30");
	ip30= i;
      }
      else if (aion[i] == 31){
	strcpy(ionam[i], "p31");
	ip31= i;
      }
      else if (aion[i] == 32){
	strcpy(ionam[i], "p32");
	ip32= i;
      }
      else if (aion[i] == 33){
	strcpy(ionam[i], "p33");
	ip33= i;
      }
      else if (aion[i] == 34){
	strcpy(ionam[i], "p34");
	ip34= i;
      }      
    }
    else if (zion[i] == 16){
      if (aion[i] == 30){
	strcpy(ionam[i], "s30");
	is30= i;
      }
      else if (aion[i] == 31){
	strcpy(ionam[i], "s31");
	is31= i;
      }
      else if (aion[i] == 32){
	strcpy(ionam[i], "s32");
	is32= i;
      }
      else if (aion[i] == 33){
	strcpy(ionam[i], "s33");
	is33= i;
      }
      else if (aion[i] == 34){
	strcpy(ionam[i], "s34");
	is34= i;
      }
      else if (aion[i] == 35){
	strcpy(ionam[i], "s35");
	is35= i;
      }
      else if (aion[i] == 36){
	strcpy(ionam[i], "s36");
	is36= i;
      }
      else if (aion[i] == 37){
	strcpy(ionam[i], "s37");
	is37= i;
      }      
    }
    else if (zion[i] == 17){
      if (aion[i] == 33){
	strcpy(ionam[i], "cl33");
	icl33= i;
      }
      else if (aion[i] == 34){
	strcpy(ionam[i], "cl34");
	icl34= i;
      }      
      if (aion[i] == 35){
	strcpy(ionam[i], "cl35");
	icl35= i;
      }      
      if (aion[i] == 36){
	strcpy(ionam[i], "cl36");
	icl36= i;
      }      
      if (aion[i] == 37){
	strcpy(ionam[i], "cl37");
	icl37= i;
      }
      if (aion[i] == 38){
	strcpy(ionam[i], "cl38");
	icl38= i;
      }      
    }
    else if (zion[i] == 18){
      if (aion[i] == 35){
	strcpy(ionam[i], "ar35");
	iar35= i;
      }
      else if (aion[i] == 36){
	strcpy(ionam[i], "ar36");
	iar36= i;
      }
      else if (aion[i] == 37){
	strcpy(ionam[i], "ar37");
	iar37= i;
      }            
      else if (aion[i] == 38){
	strcpy(ionam[i], "ar38");
	iar38= i;
      }
      else if (aion[i] == 39){
	strcpy(ionam[i], "ar39");
	iar39= i;
      }
      else if (aion[i] == 40){
	strcpy(ionam[i], "ar40");
	iar40= i;
      }
      else if (aion[i] == 41){
	strcpy(ionam[i], "ar41");
	iar41= i;
      }            
    }
    else if (zion[i] == 19){
      if (aion[i] == 37){
	strcpy(ionam[i], "k37");
	ik37= i;
      }
      else if (aion[i] == 38){
	strcpy(ionam[i], "k38");
	ik38= i;
      }      
      else if (aion[i] == 39){
	strcpy(ionam[i], "k39");
	ik39= i;
      }
      else if (aion[i] == 40){
	strcpy(ionam[i], "k40");
	ik40= i;
      }
      else if (aion[i] == 41){
	strcpy(ionam[i], "k41");
	ik41= i;
      }
      else if (aion[i] == 42){
	strcpy(ionam[i], "k42");
	ik42= i;
      }      
    }
    else if (zion[i] == 20){
      if (aion[i] == 39){
	strcpy(ionam[i], "ca39");
	ica39= i;
      }
      else if (aion[i] == 40){
	strcpy(ionam[i], "ca40");
	ica40= i;
      }
      else if (aion[i] == 41){
	strcpy(ionam[i], "ca41");
	ica41= i;
      }
      else if (aion[i] == 42){
	strcpy(ionam[i], "ca42");
	ica42= i;
      }
      else if (aion[i] == 43){
	strcpy(ionam[i], "ca43");
	ica43= i;
      }
      else if (aion[i] == 44){
	strcpy(ionam[i], "ca44");
	ica44= i;
      }
      else if (aion[i] == 45){
	strcpy(ionam[i], "ca45");
	ica45= i;
      }
      else if (aion[i] == 46){
	strcpy(ionam[i], "ca46");
	ica46= i;
      }
      else if (aion[i] == 47){
	strcpy(ionam[i], "ca47");
	ica47= i;
      }
      else if (aion[i] == 48){
	strcpy(ionam[i], "ca48");
	ica48= i;
      }
      else if (aion[i] == 49){
	strcpy(ionam[i], "ca49");
	ica49= i;
      }
    }
    else if (zion[i] == 21){
      if (aion[i] == 41){
	strcpy(ionam[i], "sc41");
	isc41= i;
      }
      else if (aion[i] == 42){
	strcpy(ionam[i], "sc42");
	isc42= i;
      }
      else if (aion[i] == 43){
	strcpy(ionam[i], "sc43");
	isc43= i;
      }      
      else if (aion[i] == 44){
	strcpy(ionam[i], "sc44");
	isc45= i;
      }
      else if (aion[i] == 45){
	strcpy(ionam[i], "sc45");
	isc45= i;
      }
      else if (aion[i] == 46){
	strcpy(ionam[i], "sc46");
	isc46= i;
      }
      else if (aion[i] == 47){
	strcpy(ionam[i], "sc47");
	isc47= i;
      }
      else if (aion[i] == 48){
	strcpy(ionam[i], "sc48");
	isc48= i;
      }
      else if (aion[i] == 49){
	strcpy(ionam[i], "sc49");
	isc49= i;
      }
      else if (aion[i] == 50){
	strcpy(ionam[i], "sc46");
	isc46= i;
      }            
    }
    else if (zion[i] == 22){
      if (aion[i] == 43){
	strcpy(ionam[i], "ti43");
	iti43= i;
      }      
      else if (aion[i] == 44){
	strcpy(ionam[i], "ti44");
	iti44= i;
      }
      else if (aion[i] == 45){
	strcpy(ionam[i], "ti45");
	iti45= i;
      }
      else if (aion[i] == 46){
	strcpy(ionam[i], "ti46");
	iti46= i;
      }
      else if (aion[i] == 47){
	strcpy(ionam[i], "ti47");
	iti47= i;
      }
      else if (aion[i] == 48){
	strcpy(ionam[i], "ti48");
	iti48= i;
      }
      else if (aion[i] == 49){
	strcpy(ionam[i], "ti49");
	iti49= i;
      }
      else if (aion[i] == 50){
	strcpy(ionam[i], "ti50");
	iti50= i;
      }
      else if (aion[i] == 51){
	strcpy(ionam[i], "ti51");
	iti51= i;
      }
      else if (aion[i] == 52){
	strcpy(ionam[i], "ti52");
	iti52= i;
      }
      else if (aion[i] == 53){
	strcpy(ionam[i], "ti52");
	iti52= i;
      }
      else if (aion[i] == 54){
	strcpy(ionam[i], "ti54");
	iti54= i;
      }
      else if (aion[i] == 55){
	strcpy(ionam[i], "ti55");
	iti55= i;
      }
      else if (aion[i] == 56){
	strcpy(ionam[i], "ti56");
	iti56= i;
      }
      else if (aion[i] == 57){
	strcpy(ionam[i], "ti57");
	iti57= i;
      }
      else if (aion[i] == 58){
	strcpy(ionam[i], "ti58");
	iti58= i;
      }      
    }
    else if (zion[i] == 23){
      if (aion[i] == 45){
	strcpy(ionam[i], "v45");
	iv45= i;
      }
      else if (aion[i] == 46){
	strcpy(ionam[i], "v46");
	iv46= i;
      }      
      else if (aion[i] == 47){
	strcpy(ionam[i], "v47");
	iv47= i;
      }
      else if (aion[i] == 48){
	strcpy(ionam[i], "v48");
	iv48= i;
      }
      else if (aion[i] == 49){
	strcpy(ionam[i], "v49");
	iv49= i;
      }
      else if (aion[i] == 50){
	strcpy(ionam[i], "v50");
	iv50= i;
      }
      else if (aion[i] == 51){
	strcpy(ionam[i], "v51");
	iv51= i;
      }
      else if (aion[i] == 52){
	strcpy(ionam[i], "v52");
	iv52= i;
      }
      else if (aion[i] == 53){
	strcpy(ionam[i], "v53");
	iv53= i;
      }
      else if (aion[i] == 54){
	strcpy(ionam[i], "v54");
	iv54= i;
      }
      else if (aion[i] == 55){
	strcpy(ionam[i], "v55");
	iv55= i;
      }
      else if (aion[i] == 56){
	strcpy(ionam[i], "v56");
	iv56= i;
      }
      else if (aion[i] == 57){
	strcpy(ionam[i], "v57");
	iv57= i;
      }
      else if (aion[i] == 58){
	strcpy(ionam[i], "v58");
	iv58= i;
      }      
    }
    else if (zion[i] == 24){
      if (aion[i] == 47){
	strcpy(ionam[i], "cr47");
	icr47= i;
      }
      else if (aion[i] == 48){
	strcpy(ionam[i], "cr48");
	icr48= i;
      }
      else if (aion[i] == 49){
	strcpy(ionam[i], "cr49");
	icr49= i;
      }
      else if (aion[i] == 50){
	strcpy(ionam[i], "cr50");
	icr50= i;
      }
      else if (aion[i] == 51){
	strcpy(ionam[i], "cr51");
	icr51= i;
      }
      else if (aion[i] == 52){
	strcpy(ionam[i], "cr52");
	icr52= i;
      }
      else if (aion[i] == 53){
	strcpy(ionam[i], "cr53");
	icr53= i;
      }
      else if (aion[i] == 54){
	strcpy(ionam[i], "cr54");
	icr54= i;
      }
      else if (aion[i] == 55){
	strcpy(ionam[i], "cr55");
	icr55= i;
      }      
      else if (aion[i] == 56){
	strcpy(ionam[i], "cr56");
	icr56= i;
      }
      else if (aion[i] == 57){
	strcpy(ionam[i], "cr57");
	icr57= i;
      }
      else if (aion[i] == 58){
	strcpy(ionam[i], "cr58");
	icr58= i;
      }
      else if (aion[i] == 59){
	strcpy(ionam[i], "cr59");
	icr59= i;	
      }
      else if (aion[i] == 60){
	strcpy(ionam[i], "cr60");
	icr60= i;
      }      
    }
    else if (zion[i] == 25){
      if (aion[i] == 51){
	strcpy(ionam[i], "mn51");
	imn51= i;
      }
      else if (aion[i] == 52){
	strcpy(ionam[i], "mn52");
	imn52= i;
      }
      else if (aion[i] == 53){
	strcpy(ionam[i], "mn53");
	imn53= i;
      }
      else if (aion[i] == 54){
	strcpy(ionam[i], "mn54");
	imn54= i;
      }
      else if (aion[i] == 55){
	strcpy(ionam[i], "mn55");
	imn55= i;
      }
      else if (aion[i] == 56){
	strcpy(ionam[i], "mn56");
	imn56= i;
      }
      else if (aion[i] == 57){
	strcpy(ionam[i], "mn57");
	imn57= i;
      }
      else if (aion[i] == 58){
	strcpy(ionam[i], "mn58");
	imn58= i;
      }
      else if (aion[i] == 59){
	strcpy(ionam[i], "mn59");
	imn59= i;
      }
      else if (aion[i] == 60){
	strcpy(ionam[i], "mn60");
	imn60= i;
      }
      else if (aion[i] == 61){
	strcpy(ionam[i], "mn61");
	imn61= i;
      }            
    }
    else if (zion[i] == 26){
      if (aion[i] == 52){
	strcpy(ionam[i], "fe52");
	ife52= i;
      }
      else if (aion[i] == 53){
	strcpy(ionam[i], "fe53");
	ife53= i;
      }
      else if (aion[i] == 54){
	strcpy(ionam[i], "fe54");
	ife54= i;
      }
      else if (aion[i] == 55){
	strcpy(ionam[i], "fe54");
	ife55= i;
      }      
      else if (aion[i] == 56){
	strcpy(ionam[i], "fe56");
	ife56= i;
      }
      else if (aion[i] == 57){
	strcpy(ionam[i], "fe57");
	ife57= i;
      }
      else if (aion[i] == 58){
	strcpy(ionam[i], "fe58");
	ife58= i;
      }
      else if (aion[i] == 59){
	strcpy(ionam[i], "fe59");
	ife59= i;
      }
      else if (aion[i] == 60){
	strcpy(ionam[i], "fe60");
	ife60= i;
      }
      else if (aion[i] == 61){
	strcpy(ionam[i], "fe61");
	ife61= i;
      }      
    }
    else if (zion[i] == 27){
      if (aion[i] == 55){
	strcpy(ionam[i], "co55");
	ico55= i;
      }
      else if (aion[i] == 56){
	strcpy(ionam[i], "co56");
	ico56= i;
      }
      else if (aion[i] == 57){
	strcpy(ionam[i], "co57");
	ico57= i;
      }
      else if (aion[i] == 58){
	strcpy(ionam[i], "co58");
	ico58= i;
      }
      else if (aion[i] == 59){
	strcpy(ionam[i], "co59");
	ico59= i;
      }
      else if (aion[i] == 60){
	strcpy(ionam[i], "co60");
	ico60= i;
      }
      else if (aion[i] == 61){
	strcpy(ionam[i], "co61");
	ico61= i;
      }
      else if (aion[i] == 62){
	strcpy(ionam[i], "co62");
	ico62= i;
      }
      else if (aion[i] == 63){
	strcpy(ionam[i], "co63");
	ico63= i;
      }
      else if (aion[i] == 64){
	strcpy(ionam[i], "co64");
	ico64= i;
      }
      else if (aion[i] == 65){
	strcpy(ionam[i], "co65");
	ico65= i;
      }
    }
    else if (zion[i] == 28){
      if (aion[i] == 56){
	strcpy(ionam[i], "ni56");
	ini56= i;
      }
      else if (aion[i] == 57){
	strcpy(ionam[i], "ni57");
	ini57= i;
      }
      else if (aion[i] == 58){
	strcpy(ionam[i], "ni58");
	ini58= i;
      }
      else if (aion[i] == 59){
	strcpy(ionam[i], "ni59");
	ini59= i;
      }
      else if (aion[i] == 60){
	strcpy(ionam[i], "ni60");
	ini60= i;
      }
      else if (aion[i] == 61){
	strcpy(ionam[i], "ni61");
	ini61= i;
      }
      else if (aion[i] == 62){
	strcpy(ionam[i], "ni62");
	ini62= i;
      }
      else if (aion[i] == 63){
	strcpy(ionam[i], "ni63");
	ini63= i;
      }
      else if (aion[i] == 64){
	strcpy(ionam[i], "ni64");
	ini64= i;
      }
      else if (aion[i] == 65){
	strcpy(ionam[i], "ni65");
	ini65= i;
      }            
    }
    else if (zion[i] == 29){
      if (aion[i] == 57){
	strcpy(ionam[i], "cu57");
	icu57= i;
      }
      else if (aion[i] == 58){
	strcpy(ionam[i], "cu58");
	icu58= i;
      }      
      else if (aion[i] == 59){
	strcpy(ionam[i], "cu59");
	icu59= i;
      }
      else if (aion[i] == 60){
	strcpy(ionam[i], "cu60");
	icu60= i;
      }
      else if (aion[i] == 61){
	strcpy(ionam[i], "cu61");
	icu61= i;
      }
      else if (aion[i] == 62){
	strcpy(ionam[i], "cu62");
	icu62= i;
      }      
      else if (aion[i] == 63){
	strcpy(ionam[i], "cu63");
	icu63= i;
      }
      else if (aion[i] == 64){
	strcpy(ionam[i], "cu64");
	icu64= i;
      }
      else if (aion[i] == 65){
	strcpy(ionam[i], "cu65");
	icu65= i;
      }
      else if (aion[i] == 66){
	strcpy(ionam[i], "cu66");
	icu66= i;
      }      
    }
    else if (zion[i] == 30){
      if (aion[i] == 60){
	strcpy(ionam[i], "zn60");
	izn60= i;
      }
      else if (aion[i] == 61){
	strcpy(ionam[i], "zn61");
	izn61= i;
      }
      else if (aion[i] == 62){
	strcpy(ionam[i], "zn62");
	izn62= i;
      }
      else if (aion[i] == 63){
	strcpy(ionam[i], "zn63");
	izn63= i;
      }
      else if (aion[i] == 64){
	strcpy(ionam[i], "zn64");
	izn64= i;
      }      
      else if (aion[i] == 64){
	strcpy(ionam[i], "zn64");
	izn64= i;
      }
      else if (aion[i] == 65){
	strcpy(ionam[i], "zn65");
	izn65= i;
      }
      else if (aion[i] == 66){
	strcpy(ionam[i], "zn66");
	izn66= i;
      }
      else if (aion[i] == 67){
	strcpy(ionam[i], "zn67");
	izn67= i;
      }
      else if (aion[i] == 68){
	strcpy(ionam[i], "zn68");
	izn68= i;
      }
      else if (aion[i] == 69){
	strcpy(ionam[i], "zn69");
	izn69= i;
      }            
    }      
  }
  return;
}
// end usrsetionam.c

//file usrgetionind.c
/*
  + usrgetionind - Get ion offset.
  Description :
  Uses ionam to get the ion offset for the ion "name".
*/	
int usrgetionind(char name[])
{
  int i, ind;
  for (i = 0; i < niso; i++){
    if (strcmp(ionam[i], name) == 0){
      ind = i;
      break;
    }
  }
  return(ind);
}
// end usrgetionind.c

//file usrinitnet.c
/*
+ usrinitnet - Initialise isotope nucleons and proton numbers
     Description :
        Sets the isotope nucleons and proton numbers for each species.
        Does 7, 13, 21 and torch networks.
Species networks
7-species network (actually 9)
He4, C12, O16, Ne20, Mg24, Si28, Ni56, (Ca40, Ti44)
13-species network
He4, C12, O16, Ne20, Mg24, Si28, S32, Ar36, Ca40, Ti44, Cr48, Fe52,
Ni56
21 Species network
H1, He3, He4, C12, O16, Ne20, Mg24, Si28, S32, Ar36, Ca40, Ti44, Cr48,
Cr56, Fe52, Fe54, Fe56, Ni56, neut, prot

Ai = number of nucleons (protons+neutrons) per isotope
Zi = number of protons per isotope
Qi = Binding energy (> 0)
binding energy in mev/isotope

*/
void usrinitnet()
{
  float mev2erg  = 1.60217648740e-6;
  float mn = 1.67492721184e-24;	// mass of neutron
  float mp = 1.67262163783e-24;  // mass of protron
  float clight = 2.99792458e10;
  float mev2gr = mev2erg/(clight*clight);
  /* Species names and indices */
  if (tim7){
    /* 7 isotope doesn't use Timmes F90 routines */
    ihe4 = 0;
    ic12 = 1;
    io16 = 2;
    ine20 = 3;
    img24 = 4;
    isi28 = 5;
    ini56 = 6;
    /* Reaction Indices */
    ir3a = 0;     /* 3He_4 -> C_12 */
    irg3a = 1;    /* C_12 -> 3He_4  */
    ircag = 2;    /* C_12(a,g)O_16 */ 
    iroga = 3;    /* O16(g,a)C_12 */ 
    iroag = 4;    /* O_16(a,g)Ne_20 */
    irnega = 5;   /* Ne_20(g,a)O_16 */
    irneag = 6;   /* Ne_20(a,g)Mg_24 */
    irmgga = 7;   /* Mg_24(g,a)Ne_20 */
    irmgag = 8;   /* Mg_24(a,g)Si_28 */
    irsiga = 9;   /* Si_28(g,a)Mg_24 */
    ircaag = 10;  /* Ca_40(a,g)Ti_44 */
    irtiga = 11;  /* Ti_44(g,a)Ca_40 */
    ir1212 = 12;  /* C_12 + C_12-> Ne_20 + He_4 */
    ir1216 = 13;  /* C_12 + O_16-> Ne_20 + He_4 */
    ir1616 = 14;  /* O_16 + O_16 -> Si_28 + He_4 */
    /* Composite Si + 7He_4 <-> Ni_56 reactions */
    irsi2ni = 15; /* Forward */
    irni2si = 16; /* Reverse */
    /* Atomic mass, Charge and Binding energy */
    /* He_4 */
    aion[ihe4] = 4;
    zion[ihe4] = 2;
    bion[ihe4] = 28.29603;
    /* C */
    aion[ic12] = 12;
    zion[ic12] = 6;
    bion[ic12] = 92.16294;
    /* O */
    aion[io16] = 16;
    zion[io16] = 8;
    bion[io16] = 127.62093;
    /* Ne */ 
    aion[ine20] = 20;
    zion[ine20] = 10;
    bion[ine20] = 160.64788;
    /* Mg */
    aion[img24] = 24;
    zion[img24] = 12;
    bion[img24] = 198.25790;
    /* Si */
    aion[isi28] = 28;
    zion[isi28] = 14;
    bion[isi28] = 236.53790;
    /* Ni */
    aion[ini56] = 56;
    zion[ini56] = 28;
    bion[ini56] = 484.00300;
    if (rrscreen) {
      /* Initialise screening */
      usrinitscreen(zion[ihe4], aion[ihe4], zion[ihe4], aion[ihe4],1);
      usrinitscreen(zion[ihe4], aion[ihe4], 4.0, 8.0, 2);
      usrinitscreen(zion[ic12], aion[ic12], zion[ihe4], aion[ihe4], 3);
      usrinitscreen(zion[ic12], aion[ic12], zion[ic12], aion[ic12], 4);
      usrinitscreen(zion[ic12], aion[ic12], zion[io16], aion[io16], 5);
      usrinitscreen(zion[io16], aion[io16], zion[io16], aion[io16], 6);
      usrinitscreen(zion[io16], aion[io16], zion[ihe4], aion[ihe4], 7);
      usrinitscreen(zion[ine20], aion[ine20], zion[ihe4], aion[ihe4], 8);
      usrinitscreen(zion[img24], aion[img24], zion[ihe4], aion[ihe4], 9);
      usrinitscreen(zion[ica40], aion[ica40], zion[ihe4], aion[ihe4], 10);
    }
  }
  /* Set ion names */
  usrsetionam();
  /* Atomic mass, Charge and Binding energy alternative  arrays */
  for (int i = 0; i < niso; i++){
    Ai[i] = aion[i];
    Zi[i] = zion[i];  
    Qi[i] = bion[i];
  }
  /* Use binding energy including weak changes */
  if (doweak == 1){
    for (int i = 0; i < niso; i++)
      Qi[i]  = enuc_conv2*((Ai[i] - Zi[i])*mn + Zi[i]*mp
			   - mev2gr*Qi[i]);
  }
  else{
    for (int i = 0; i < niso; i++)
      Qi[i] *= - mev2erg*N_A;    
  }
  return;
}
//end usrinitnet.c

//file hermite.c
/*
     Description :
        quintic hermite basis statement functions
*/
  /* psi0 and its derivatives */
float psi0(float z)
{
  return(z*z*z*(z*(-6.0*z + 15.0) - 10.0) + 1.0);
}
float dpsi0(float z)
{
  return(z*z*(z*(-30.0*z + 60.0) - 30.0));
}
float ddpsi0(float z)
{
  return(z*(z*(-120.0*z + 180.0) - 60.0));
}
  /* psi1 and its derivatives */
float psi1(float z)
{
  return(z*(z*z*(z*(-3.0*z + 8.0) - 6.0) + 1.0));
}
float dpsi1(float z)
{
  return(z*z*(z*(-15.0*z + 32.0) - 18.0) + 1.0);
}
float ddpsi1(float z)
{
  return(z*(z*(-60.0*z + 96.0) - 36.0));
}
  /* psi2 and its derivatives */
float psi2(float z)
{
  return(0.5*z*z*(z*(z*(-z + 3.0) - 3.0) + 1.0));
}
float dpsi2(float z)
{
  return(0.5*z*(z*(z*(-5.0*z + 12.0) - 9.0) + 2.0));
}
float ddpsi2(float z)
{
  return(0.5*(z*(z*(-20.0*z + 36.0) - 18.0) + 2.0));
}
/* biquintic  hermite polynomial statement function */
float herm5(int i, int j, float w0t, float w1t, float w2t,
  float w0mt, float w1mt, float w2mt, float w0d, float w1d,
  float w2d, float w0md, float w1md, float w2md)
{
  float tv;
  tv = fl[1]*w0d*w0t + fl[2]*w0md*w0t
    + fl[3]*w0d*w0mt+ fl[4]*w0md*w0mt
    + fl[5]*w0d*w1t + fl[6]*w0md*w1t 
    + fl[7]*w0d*w1mt+ fl[8]*w0md*w1mt 
    + fl[9]*w0d*w2t + fl[10] *w0md*w2t 
    + fl[11]*w0d*w2mt+ fl[12]*w0md*w2mt 
    + fl[13]*w1d*w0t + fl[14]*w1md*w0t 
    + fl[15]*w1d*w0mt+ fl[16]*w1md*w0mt 
    + fl[17]*w2d*w0t + fl[18]*w2md*w0t 
    + fl[19]*w2d*w0mt+ fl[20]*w2md*w0mt 
    + fl[21]*w1d*w1t + fl[22]*w1md*w1t 
    + fl[23]*w1d*w1mt+ fl[24]*w1md*w1mt 
    + fl[25]*w2d*w1t + fl[26]*w2md*w1t 
    + fl[27]*w2d*w1mt+ fl[28]*w2md*w1mt 
    + fl[29]*w1d*w2t + fl[30]*w1md*w2t 
    + fl[31]*w1d*w2mt+ fl[32]*w1md*w2mt
    + fl[33]*w2d*w2t + fl[34]*w2md*w2t
    + fl[35]*w2d*w2mt + fl[36]*w2md*w2mt;
  return(tv);
}
//end hermite.c

//file usrinithelm.c
/*
+ usrinithelm - Initialise EOS (helm)
     Description :
        Initializes EOS  table
*/
bool usrinithelm()
{
  FILE *input;
  int i, j;
  float tf, tfd, tft, tfdd, tftt, tfdt, tfddt, tfdtt, tfddtt;
  if (!(input = fopen("helm_table.dat","r"))){
    cout << "Error: no helm_table.dat \n";
    return(false);
  }  
  for (j = 1; j <= JMAX; j++){
    for (i = 1; i <= IMAX; i++){
      safe_fscanf(input,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &tf, &tfd,
        &tft, &tfdd, &tftt, &tfdt, &tfddt, &tfdtt, &tfddtt);
      f[j][i] = tf;
      fd[j][i] = tfd;
      ft[j][i] = tft;
      fdd[j][i] = tfdd;
      ftt[j][i] = tftt;
      fdt[j][i] = tfdt;
      fddt[j][i] = tfddt;
      fdtt[j][i] = tfdtt;
      fddtt[j][i] = tfddtt;
    }
  }
 usrmyclose(input);

 /* Density and temperature arrays */
 for (i = 1; i <= IMAX; i++)
   d[i] = powf(10.0, dlo + (i-1)*dstp);
 for (j = 1; j <= JMAX; j++)
   te[j] = powf(10.0, tlo + (j-1)*tstp);
 return(true);
}
//end usrinithelm.c

//file usrhelm.c
/*
+ usrhelm - Electron-positron equation of state
     Description :
        This   routine   performs   a   thermodynamically   consistent
        interpolation in a tabular electron-positron equation of state
        using  biquintic  hermite  basis  functions.   Based  on  F.X.
        Timmes, D.F. Swesty,  APjS, 126, 501, 2000. Data  is stored in
        helm_table.dat

input :
temp = temperature(in K)
rho = density (in g cm~3)
Ye = electrons per baryon = Zbar/Abar

helm_table.dat contains the  table of helmholtz free  energy and eight
of  its partial  derivatives :  f,  df_d, df_t,  df_dd, df_tt,  df_dt,
df_ddt, df_dtt and df_ddtt

output :
he[0] = pressure
he[1] = specific internal energy
he[2] = specific entropy
he[3] = partial derivative of pressure with density
he[4] = partial derivative of pressure with temperature
he[5] = partial derivative of energy with density
he[6] = partial derivative of energy with temperature
he[7] = partial derivative of entropy with density
he[8] = partial derivative of entropy with temperature
he[9] = adiabatic sound speed
all in c.g.s

  Range:  1e4 < T <1e11
          1e-3 < d < 1e11
*/
void usrhelm(float temp, float rho, float Ye, float he[])
{
  int i, j, iat, jat;
  float free, df_d, df_t, df_dd, df_tt, df_dt, dt, dt2, dt2i, dd,
  dd2, ddi, dd2i, xt, xd, mxt, mxd, si0t, si1t, si2t, si0mt, si1mt,
  si2mt, si0d, si1d, si2d, si0md, si1md, si2md, dsi0t, dsi1t, dsi2t,
  dsi0mt, dsi1mt, dsi2mt, dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md,
  ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, ddsi0d, ddsi1d,
  ddsi2d, ddsi0md, ddsi1md, ddsi2md, z, w0t, w1t, w2t, w0mt, w1mt, w2mt,
  w0d, w1d, w2d, w0md, w1md, w2md, din, ptot, dpdt, dpdd, stot, dsdt,
    dsdd, etot, dedt, dedd, dedd_p, dedp_d;

  /* f and its derivatives
  float f[JMAX][IMAX], fd[JMAX][IMAX], ft[JMAX][IMAX], fdd[JMAX][IMAX],
  ftt[JMAX][IMAX], fdt[JMAX][IMAX], fddt[JMAX][IMAX],
  fdtt[JMAX][IMAX],fddtt[JMAX][IMAX]; */

  /* Check if temperature is in range */
  if ((log10(temp) <= tlo) || (log10(temp) >= thi)){
    printf("Temperature out of range = %11.4e \n", log10(temp));
    //temp = 1e11;		// temporarily set a max temp
    *(int*)0 = 0;
    exit(EXIT_FAILURE);
  }
  
  /* Find the table locations, entering with Ye * density */

  jat = int((log10(temp) - tlo)/tstp) + 1;
  din = Ye*rho;
  iat = int((log10(din) - dlo)/dstp) + 1;
  /* Various differences */
  dt = te[jat+1] - te[jat];
  dt2 = dt*dt;
  dd = d[iat+1] - d[iat];
  dd2 = dd*dd;
  xt = max((temp - te[jat])/dt, 0.0);
  xd = max((din - d[iat])/dd,0.0);
  mxt = 1.0 - xt;
  mxd = 1.0 - xd;
  /* Access the table locations only once */
  fl[1] = f[jat][iat];
  fl[2] = f[jat][iat+1];
  fl[3] = f[jat+1][iat];
  fl[4] = f[jat+1][iat+1];
  fl[5] = ft[jat][iat];
  fl[6] = ft[jat][iat+1];
  fl[7] = ft[jat+1][iat];
  fl[8] = ft[jat+1][iat+1];
  fl[9] = ftt[jat][iat];
  fl[10] = ftt[jat][iat+1];
  fl[11] = ftt[jat+1][iat];
  fl[12] = ftt[jat+1][iat+1];
  fl[13] = fd[jat][iat];
  fl[14] = fd[jat][iat+1];
  fl[15] = fd[jat+1][iat];
  fl[16] = fd[jat+1][iat+1];
  fl[17] = fdd[jat][iat];
  fl[18] = fdd[jat][iat+1];
  fl[19] = fdd[jat+1][iat];
  fl[20] = fdd[jat+1][iat+1];
  fl[21] = fdt[jat][iat];
  fl[22] = fdt[jat][iat+1];
  fl[23] = fdt[jat+1][iat];
  fl[24] = fdt[jat+1][iat+1];
  fl[25] = fddt[jat][iat];
  fl[26] = fddt[jat][iat+1];
  fl[27] = fddt[jat+1][iat];
  fl[28] = fddt[jat+1][iat+1];
  fl[29] = fdtt[jat][iat];
  fl[30] = fdtt[jat][iat+1];
  fl[31] = fdtt[jat+1][iat];
  fl[32] = fdtt[jat+1][iat+1];
  fl[33] = fddtt[jat][iat];
  fl[34] = fddtt[jat][iat+1];
  fl[35] = fddtt[jat+1][iat];
  fl[36] = fddtt[jat+1][iat+1];
  /* Evaluate the basis functions */
  si0t = psi0(xt);
  si1t = psi1(xt)*dt;
  si2t = psi2(xt)*dt2;
  si0mt = psi0(mxt);
  si1mt = -psi1(mxt)*dt;
  si2mt = psi2(mxt)*dt2;
  si0d = psi0(xd);
  si1d = psi1(xd)*dd;
  si2d = psi2(xd)*dd2;
  si0md = psi0(mxd);
  si1md = -psi1(mxd)*dd;
  si2md = psi2(mxd)*dd2;
  /* their first derivatives */
  dsi0t = dpsi0(xt)/dt;
  dsi1t = dpsi1(xt);
  dsi2t = dpsi2(xt)*dt;
  dsi0mt = -dpsi0(mxt)/dt;
  dsi1mt = dpsi1(mxt);
  dsi2mt = -dpsi2(mxt)*dt;
  dsi0d = dpsi0(xd)/dd;
  dsi1d = dpsi1(xd);
  dsi2d = dpsi2(xd)*dd;
  dsi0md = -dpsi0(mxd)/dd;
  dsi1md = dpsi1(mxd);
  dsi2md = -dpsi2(mxd)*dd;
  /*their second derivatives */
  ddsi0t = ddpsi0(xt)/dt2;
  ddsi1t = ddpsi1(xt)/dt;
  ddsi2t = ddpsi2(xt);
  ddsi0mt = ddpsi0(mxt)/dt2;
  ddsi1mt = -ddpsi1(mxt)/dt;
  ddsi2mt = ddpsi2(mxt);
  ddsi0d = ddpsi0(xd)/dd2;
  ddsi1d = ddpsi1(xd)/dd;
  ddsi2d = ddpsi2(xd);
  ddsi0md = ddpsi0(mxd)/dd2;
  ddsi1md = -ddpsi1(mxd)/dd;
  ddsi2md = ddpsi2(mxd);

  /* free energy */
  free = herm5(iat, jat, si0t, si1t, si2t, si0mt, si1mt, si2mt, si0d,
	       si1d, si2d, si0md, si1md, si2md);

/* derivative of the free energy wrt density */
  df_d = herm5(iat, jat, si0t, si1t, si2t, si0mt, si1mt, si2mt, dsi0d,
	       dsi1d, dsi2d, dsi0md, dsi1md, dsi2md);

  /* derivative of the free energy wrt temperature */
  df_t = herm5(iat, jat, dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt,
	       si0d, si1d, si2d, si0md, si1md, si2md);

  /* second derivative of the free energy wrt to density */
  df_dd = herm5(iat, jat, si0t, si1t, si2t, si0mt, si1mt, si2mt,
		ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md);

  /* second derivative of the free energy wrt temperature */
  df_tt = herm5(iat, jat, ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt,
		ddsi2mt, si0d, si1d, si2d, si0md, si1md, si2md);

  /* second derivative of the free energy wrt temperature and density */
  df_dt = herm5(iat, jat, dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt,
		dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md);

  /* set  the  return  arguments  ;  the  electron-positron  pressure,
     specific  entropy,  and  internal  energy  along  with  their
     partial derivatives */

  ptot = din*din*df_d;
  dpdt = din*din*df_dt;
  dpdd = Ye*(din*din*df_dd + 2.0*din*df_d);
  stot = -df_t*Ye;
  dsdt = -df_tt*Ye;
  dsdd = -df_dt*Ye*Ye;

  etot = Ye*free + temp*stot;
  dedt = temp*dsdt;
  dedd = Ye*Ye*df_d + temp*dsdd;
  he[0] = ptot;
  he[1] = etot;
  he[2] = stot;
  he[3] = dpdd;
  he[4] = dpdt;
  he[5] = dedd;
  he[6] = dedt;
  he[7] = dsdd;
  he[8] = dsdt;
  return;
}
//end usrhelm.c

//file usrrate.c
/*
+ usrrate - Nuclear reaction rate coefficients
     Description :
        Reaction rate coefficients. They  are calculated here for tim7
        and from  Timmes' F90  routines for  other networks.  They are
        returned in RR[].
*/
void usrrate(float pa[], float RR[])
{
  int i;
  float scor = 1.0;
  float scor1 = 1.0;
  float scor2 = 1.0;
  float rho, temp;
  float t9, t92, t93, t94, t95, t96, t912, t932, t952, t972, t913,
    t923, t943, t953, t973, t9113, t914, t934, t954, t974, t915, t935,
    t945, t965, t917, t927, t947, t918, t938, t958, t9i, t9i2, t9i3,
    t9i12, t9i32, t9i52, t9i72, t9i13, t9i23, t9i43, t9i53, t9i14,
    t9i34, t9i54, t9i15, t9i35, t9i45, t9i65, t9i17, t9i27, t9i47,
    t9i18, t9i38, t9i58, t916, t976, t9i76, t9l, t92l, t93l, t9a, t9a13,
    t9a56, t9a23, rc148;  
  float r2abe, rbeac, tmp, tmp2, term, rev, term1, term2, term3,
    rsi2ni, rni2si;
  float q1, q2, rc28, aa, bb, cc, dd, ee, ff, gg, hh, xx, uu, yy, vv,
    zz, f1, f2, yeff_ca40, yeff_ti44, t9i92, denom, t992, hhi, rc121,
    rc102;  
  /* Zero for safety */
  for (i = 0; i < niso; i++)
    RR[i] = 0.0;
  if (tim7){
    /* 7 isotope doesn't use Timmes F90 routines */    
    rho = pa[iqd];
    temp = pa[iqte];
    t9 = temp*1.0e-9;
    /* Direct from Maxwell sntI usr.c with limits removed */
    t9l = min(t9, 10.0);
    t92l = t9l*t9l;
    t93l = t92l*t9l;
    /* Unlimited */
    t92 = t9*t9;
    t93 = t9*t92;
    t94 = t9*t93;
    t95 = t9*t94;
    t96 = t9*t95;
    t912 = sqrt(t9);
    t932 = t9*t912;
    t952 = t9*t932;
    t972 = t9*t952;
    t913 = pow(t9, 1.0/3.0);
    t923 = t913*t913;
    t943 = t9*t913;
    t953 = t9*t923;
    t973 = t953*t923;
    t9113 = t973*t943;
    t914 = pow(t9, 0.25);
    t934 = t914*t914*t914;
    t954 = t9*t914;
    t974 = t9*t934;
    t915 = pow(t9, 0.2);
    t935 = t915*t915*t915;
    t945 = t915*t935;
    t965 = t9*t915;
    t916 = pow(t9, 1.0/6.0);
    t976 = t9*t916;
    t9i76 = 1.0/t976;
    t917 = pow(t9, 1.0/7.0);
    t927 = t917*t917;
    t947 = t927*t927;
    t918 = sqrt(t914);
    t938 = t918*t918*t918;
    t958 = t938*t918*t918;
    t9i  = 1.0/t9;
    t9i2 = t9i*t9i;
    t9i3 = t9i2*t9i;
    t9i12 = 1.0/t912;
    t9i32 = t9i*t9i12;
    t9i52 = t9i*t9i32;
    t9i72 = t9i*t9i52;
    t9i13 = 1.0/t913;
    t9i23 = t9i13*t9i13;
    t9i43 = t9i*t9i13;
    t9i53 = t9i*t9i23;
    t9i14 = 1.0/t914;
    t9i34 = t9i14*t9i14*t9i14;
    t9i54 = t9i*t9i14;
    t9i15 = 1.0/t915;
    t9i35 = t9i15*t9i15*t9i15;
    t9i45 = t9i15*t9i35;
    t9i65 = t9i*t9i15;
    t9i17 = 1.0/t917;
    t9i27 = t9i17*t9i17;
    t9i47 = t9i27*t9i27;
    t9i18 = 1.0/t918;
    t9i38 = t9i18*t9i18*t9i18;
    t9i58 = t9i38*t9i18*t9i18;
    /* Rates from http://cococubed.asu.edu/*/
    /* 3He_4 -> C_12 reaction and reverse */
    q1 = 1.0 / 0.009604;
    q2 = 1.0 / 0.055225;
    rc28 = 0.1;
    /* a(a,g)Be_8 */
    aa = 7.40e5*t9i32*exp(-1.0663*t9i);
    bb = 4.164e9*t9i23*exp(-13.49*t9i13 - t92*q1);
    cc = 1.0 + 0.031*t913 + 8.009*t923 + 1.732*t9 + 49.883*t943
      + 27.426*t953;
    r2abe = aa + bb*cc;
    /* Be_8(a,g)C_12 */
    dd = 130.0*t9i32*exp(-3.3364*t9i);
    ee = 2.510e7*t9i23*exp(-23.57*t9i13 - t92*q2);
    ff = 1.0 + 0.018*t913 + 5.249*t923 + 0.650*t9 + 19.176*t943
      + 6.034*t953;
    rbeac = dd + ee*ff;
    /* Alpha factor */
    xx = rc28*1.35e-7*t9i32*exp(-24.811*t9i);
    if (t9 > 0.08)
      /* High temperature */
      term = 2.90e-16*r2abe*rbeac + xx;
    else{
      /* Low temperature */
      uu = 0.8*exp(-pow(0.025*t9i, 3.263));
      yy = 0.2 + uu;
      vv = 4.0*exp(-pow(t9/0.025, 9.227));
      zz = 1.0 + vv;
      aa = 1.0/zz;
      f1 = 0.01 + yy*aa;
      term = 2.90e-16*r2abe*rbeac*f1 + xx;
    }
    rev = 2.00e20*t93*exp(-84.424*t9i);
    RR[ir3a] = term*rho*rho;
    RR[irg3a] = rev*term;
    /* C_12(a,g)O_16 and reverse */
    q1 = 1.0 / 12.222016;
    aa = 1.0 + 0.0489*t9i23;
    bb = t92*aa*aa;
    cc = exp(-32.120*t9i13 - t92*q1);
    dd = 1.0 + 0.2654*t9i23;
    ee = t92*dd*dd;
    ff = exp(-32.120*t9i13);
    gg = 1.25e3*t9i32*exp(-27.499*t9i);
    hh = 1.43e-2*t95*exp(-15.541*t9i);
    zz = 1.0 / bb;
    f1 = cc*zz;
    zz = 1.0 / ee;
    f2 = ff*zz;
    term = 1.04e8*f1 + 1.76e8*f2 + gg + hh;
    term = 1.7*term;
    rev = 5.13e10*t932*exp(-83.111*t9i);  
    /* C_12(a,g)O_16 and reverse */
    RR[ircag] = term*rho;
    RR[iroga] = rev*term;
    /* O_16(a,g)Ne_20 and reverse */
    q1 = 1.0/2.515396;
    term1 = 9.37e9*t9i23*exp(-39.757*t9i13 - t92*q1);
    aa = 62.1*t9i32*exp(-10.297*t9i);
    bb = 538.0*t9i32*exp(-12.226*t9i);
    cc = 13.0*t92*exp(-20.093*t9i);
    term2 = aa + bb + cc;
    term = term1 + term2;
    rev = 5.65e10*t932*exp(-54.937*t9i);  
    RR[iroag] = term*rho;
    RR[irnega] = rev*term;
    /* Ne_20(a,g)Mg_24 + g and reverse*/
    rc102 = 0.1;
    q1 = 1.0 / 4.923961;
    aa = 4.11e11*t9i23*exp(-46.766*t9i13 - t92*q1);
    bb = 1.0 + 0.009*t913 + 0.882*t923 + 0.055*t9 + 0.749*t943
      + 0.119*t953;
    term1 = aa*bb;
    aa = 5.27e3*t9i32*exp(-15.869*t9i);
    bb = 6.51e3*t912*exp(-16.223*t9i);
    term2 = aa + bb;
    aa = 42.1*t9i32*exp(-9.115*t9i);
    bb = 32.0*t9i23*exp(-9.383*t9i);
    term3 = rc102*(aa + bb);
    aa = 5.0*exp(-18.960*t9i);
    bb = 1.0 + aa;
    zz = 1.0 / bb;
    term = (term1 + term2 + term3)*zz;
    rev = 6.01e10*t932*exp(-108.059*t9i);
    RR[irneag] = term*rho;
    RR[irmgga] = rev*term;
    /* Mg_24(a,g)Si_28 + g and reverse */
    rc121 = 0.1;
    aa = 4.78e1*t9i32*exp(-13.506*t9i);
    bb = 2.38e3*t9i32*exp(-15.218*t9i);
    cc = 2.47e2*t932*exp(-15.147*t9i);
    dd = rc121*1.72e-09*t9i32*exp(-5.028*t9i);
    ee = rc121*1.25e-03*t9i32*exp(-7.929*t9i);
    ff = rc121*2.43e1*t9i*exp(-11.523*t9i);
    gg = 5.0*exp(-15.882*t9i);
    hh = 1.0 + gg;
    hhi = 1.0 / hh;
    term = (aa + bb + cc + dd + ee + ff)*hhi;
    rev = 6.27e10*t932*exp(-115.862*t9i);
    RR[irmgag] = term*rho;
    RR[irsiga] = rev*term;
    /* Ca_40(a,g)Ti_44 and the reverse - use t9l, t92l, t93l */
    aa = 1.0 + 1.650e-2*t9l + 5.973e-3*t92l - 3.889e-4*t93l;
    term = 4.66e24*t9i23*exp(-76.435*t9i13*aa);
    rev = 6.843e10*t932*exp(-59.510*t9i);
    RR[ircaag] = term*rho;
    RR[irtiga] = rev*term;
    /* C_12 + C_12-> Ne_20 + He_4 */
    aa = 1.0 + 0.0396*t9;
    zz = 1.0/aa;
    t9a = t9*zz;
    t9a13 = pow(t9a, 1.0/3.0);
    t9a56 = pow(t9a, 5.0/6.0);
    term = 4.27e26*t9a56*t9i32*exp(-84.165/t9a13 - 2.12e-3*t93);
    RR[ir1212] = term*rho;
    if (t9 >= 0.5) {
      aa = 1.0 + 0.055*t9;
      zz = 1.0/aa;
      t9a = t9*zz;
      t9a13 = pow(t9a, 1.0/3.0);
      t9a23 = t9a13*t9a13;
      t9a56 = pow(t9a, 5.0/6.0);
      aa = exp(-0.18*t9a*t9a);
      bb = 1.06e-3*exp(2.562*t9a23);
      cc = aa + bb;
      zz = 1.0/cc;
      term = 1.72e31*t9a56*t9i32*exp(-106.594/t9a13)*zz;
    }
    else
      term = 0.0;
    /* C_12 + O_16 -> Mg_24 + He_4 see cf88 references 47-4 */
    RR[ir1216] = term*rho;
    /* O_16 + O_16 -> Si_28 + He_4 */
    term = 7.10e36*t9i23
      *exp(-135.93*t9i13 - 0.629*t923 - 0.445*t943 + 0.0103*t9*t9);
    RR[ir1616] = term*rho;
    if (rrscreen){
      /* Calculate Zbar, Abar, Z2bar */
      float Ye, Abar, Zbar, Z2bar;
      usrYAZ(pa, Ye, Abar, Zbar, Z2bar);    
      /* screening functions */
      /* Triple alpha */
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[ihe4], Ai[ihe4],
		Zi[ihe4], Ai[ihe4], 1, scor1);
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[ihe4], Ai[ihe4], 4.0,
		8.0, 2, scor2);
      scor = scor1*scor2;
      RR[ir3a] *= scor;
      /* C_12 -> O_16 */
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[ic12], Ai[ic12],
		Zi[ihe4], Ai[ihe4], 3, scor);
      RR[ircag] *= scor;
      /* C_12 + C_12 */     
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[ic12], Ai[ic12],
		Zi[ic12], Ai[ic12], 4, scor);
      RR[ir1212] *= scor;
      /* C_12 + O_16 */
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[ic12], Ai[ic12],
		Zi[io16], Ai[io16], 5, scor);
      RR[ir1216] *= scor;
      /* O_16 + O_16 */
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[io16], Ai[io16],
		Zi[io16], Ai[io16], 6, scor);
      RR[ir1616] *= scor;
      /* O_16 to Ne_20 */ 
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[io16], Ai[io16],
		Zi[ihe4], Ai[ihe4], 7, scor);
      RR[iroag] *= scor;
      /* Ne_20 to Mg_24 */
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[ine20], Ai[ine20],
		Zi[ihe4], Ai[ihe4], 8, scor);
      RR[irneag] *= scor;
      /* Ca_40 to Ti_44 (note jscreen=10 is different for iso7 */      
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[ica40], Ai[ica40],
		Zi[ihe4], Ai[ihe4], 10, scor);
      RR[ircaag] *= scor;
      /* Mg_24 to Si_28 */
      usrscreen(temp, rho, Zbar, Abar, Z2bar, Zi[img24], Ai[img24],
		Zi[ihe4], Ai[ihe4], 9, scor);
      RR[irmgag] *= scor;
      /* Si_28 + 7He_4 -> Ni_56 and reverse*/
      RR[irni2si] = 0.0;
      RR[irsi2ni] = 0.0;
      if (t9 > 2.5) {
	if ((pa[iqal0+1+ic12] + pa[iqal0+1+io16]) <= 0.004) {
	  t992 = t972*t9;
	  t9i92 = 1.0 / t992;
	  yeff_ca40 = t9i92*exp(239.42*t9i - 74.741);
	  yeff_ti44 = t992*exp(-274.12*t9i + 74.914);
	  denom = rho*rho*rho*pa[iqal0+1+ihe4]*pa[iqal0+1+ihe4]*
	    pa[iqal0+1+ihe4];
	  /* Si_28 + 7He_4 -> Ni_56 */
	  RR[irsi2ni] = yeff_ca40*denom*RR[ircaag]*pa[iqal0+1+isi28];
	  /* Ni_56 -> Si_28 + 7He_4 */
	  if (denom != 0.0) {
	    zz = 1.0 / denom;
	    RR[irni2si] = yeff_ti44*RR[irtiga]*zz;
	    if (RR[irni2si] > 1.0e10)
	      RR[irni2si] = 1.0e10;
	  }
	}
      }
    }
  }
  else if (tim13)
    /* 13 iso from Timmes F90 */
    timrat13_(&nd, pa, RR);
  else if (tim21)
    /* 21 iso from Timmes F90 */
    timrat21_(&nd, pa, RR);
  else if (timtorch)
    /* Get reaction rate coefficients from torch */
    torch_rat_(&nd, &doweak, &doneut, pa, RR);
  return;
}
//end usrrate.c

//file usrreact.c
/*
+ usrreact - Nuclear reaction source terms
     Description :
        Returns  the  rate  of  change of  the  molar  abundances  and
        temperature, dydt[i], given the  abundances and temperature in
        pa[]. Uses from routines from Timmes aprox21 code
*/
void usrreact(float pa[], float ydot[])
{
  int i, iymax;
  float y[MAXSPECIES+1], s[MAXREACT], RR[MAXREACT], tv, enuc;
  reactcall++;  
  /* Set molar abundances, zero ydot */
  for (i = 0; i < niso; i++){
    y[i] = pa[iqal0+i+1];
    ydot[i] = 0.0;
  }
  if (tim7){
    /* Get 7 iso reaction rates */
    usrrate(pa, RR);
    /* Set 7 iso source terms */
    s[ir3a] = 0.5*y[ihe4]*y[ihe4]*y[ihe4]*RR[ir3a];
    s[irg3a] = 3.0*y[ic12]*RR[irg3a];
    s[ircag] = y[ic12]*y[ihe4]*RR[ircag];
    s[iroga] = y[io16]*RR[iroga];
    s[iroag] = y[io16]*y[ihe4]*RR[iroag];
    s[irnega] = y[ine20]*RR[irnega];
    s[irneag] = y[ine20]*y[ihe4]*RR[irneag];
    s[irmgga] = y[img24]*RR[irmgga];
    s[irmgag] = y[img24]*y[ihe4]*RR[irmgag];
    s[irsiga] = y[isi28]*RR[irsiga];
    s[ir1212] = 0.5*y[ic12]*y[ic12]*RR[ir1212];
    s[ir1216] = 0.5*y[ic12]*y[io16]*RR[ir1216];
    s[ir1616] = 0.5*y[io16]*y[io16]*RR[ir1616];
    s[irsi2ni] = 7.0*y[ihe4]*RR[irsi2ni];
    s[irni2si] = 7.0*y[ini56]*RR[irni2si];
    /* He_4 source */
    ydot[ihe4] = - s[ir3a] + s[irg3a] - s[ircag] + s[iroga] - s[iroag]
      + s[irnega] - s[irneag]  + s[irmgga] - s[irmgag] + s[irsiga]
      + s[ir1212] + s[ir1216] + 0.56*s[ir1616];
    ydot[ihe4] += - s[irsi2ni] + s[irni2si];
    /* C_12 source */
    ydot[ic12] = (s[ir3a] - s[irg3a])/3.0 - s[ircag] + s[iroga]
      - 2.0*(s[ir1212] + s[ir1216]);
    /* O_16 source */
    ydot[io16] = s[ircag] - s[iroga] - s[iroag]  + s[irnega]
      - 2.0*(s[ir1216] + s[ir1616]);
    /* Ne_20 source */
    ydot[ine20] = s[iroag] - s[irnega] - s[irneag] + s[irmgga]
      + s[ir1212];
    /* Mg_24 source */
    ydot[img24] = s[irneag] - s[irmgga] - s[irmgag] + s[irsiga]
      + s[ir1216];
    /* Si_28 source */
    ydot[isi28] = s[irmgag] - s[irsiga] + s[ir1216] + 0.56*s[ir1616];  
    /* Si_28 source */
    ydot[isi28] += -(s[irsi2ni] - s[irni2si])/7.0;
    /* Ni_56 source */
    ydot[ini56] = (s[irsi2ni] - s[irni2si])/7.0;
  }
  else if (tim13)
    /* 13 iso source terms from Timmes 13 iso F90 */
    timreact13_(&nd, pa, ydot);
  else if (tim21)
    /* 21 iso from Timmes 21 iso F90 */
    timreact21_(&nd, pa, ydot);
  else if (timtorch)
    /* Get reaction rates from Timmes torch */
    torch_react_(&nd, &doweak, &doneut, pa, ydot);
  if (cons){
    /* Find max mass fraction */
    tv = -1.0;
    for (i = 0; i < niso; i++){
      if (y[i]*Ai[i] > tv){
	tv = y[i]*Ai[i];
	iymax = i;
      }
    }      
    /* Max y[i] from conservation */
    tv = 0.0;
    for (i = 0; i < niso; i++){
      if (i != iymax)
	tv += ydot[i]*Ai[i];
    }
    ydot[iymax] = - tv/Ai[iymax];
  }
  /* ydot Temperature K */
  if (isothermal)
    ydot[niso] = 0.0;
  else{
    if (!timtorch){
      /* Reaction energy generation rate */
      enuc = 0.0;
      for (i = 0; i < niso; i++)
	enuc += ydot[i]*Qi[i];
      ydot[niso] = - enuc;
    }
    /* Temperature K */      
    ydot[niso] /= usrcv(pa);
  }
  return;
}
//end usrreact.c

//file usrinitscreen.c
/*
+ usrinitscreen - Set up screening factor calculation
     Description :
        Initialisation  for  screening  factors used  in  calculating
        iso 7 reaction rates
*/
void usrinitscreen(float z1, float a1, float z2, float a2, int jscreen)
{
  float x13, x53, x512;
  x13 = 0.333333333;
  x53 = 1.666666667;
  x512 = 0.416666667;

  zs13[jscreen] = pow((z1 + z2), x13);
  zs13inv[jscreen] = 1.0 / zs13[jscreen];
  zhat[jscreen] = pow((z1 + z2), x53) - pow(z1, x53) - pow(z2, x53);
  zhat2[jscreen] = pow((z1 + z2), x512) - pow(z1, x512) - pow(z2, x512);
  lzav[jscreen] = x53*log(z1*z2 / (z1 + z2));
  aznut[jscreen] = pow((z1*z1*z2*z2*a1*a2 / (a1 + a2)), x13);
  return;
}
//end usrinitscreen.c

//file usrscreen.c
/*
+ usrscreen - Calculate screening factors
     Description :
        Returns screening factors for iso7 reaction rates.
*/
void usrscreen(float temp, float rho, float Zbar, float Abar,
  float z2bar, float z1, float a1, float z2, float a2, int jscreen,
  float& scor)
{
  /*
common
    jscreen = 1	triple alpha (He)
    = 2		its inverse
    = 3		c12 to o16
    = 4		c12 + c12
    = 5		c12 + o16
    = 6		o16 + o16
    = 7		o16 to ne20
iso 7
    = 8		o16 to mg24
    = 9		mg24 to si28
    = 10	ca40 to ti44
iso13
    = 8         ne20 to mg24
    = 9         mg24 to si28, mg24 to al27
    = 10        al27 to si28
    = 11        si28 to s32, si28 to p31
    = 12        p31 to s32
    = 13        s32 to ar36, s32 to cl35
    = 14        cl35 to ar36
    = 15        ar36 to ca40, ar36 to k39
    = 16        k39 to ca40
    = 17        ca40 to ti44, ca40 to sc43
    = 18        sc43 to ti44
    = 19        ti44 to cr48. ti44 to v47
    = 20        v47 to cr48
    = 21        cr48 to fe52, cr48 to mn51
    = 22        mn51 to fe52
    = 23        fe52 to ni56, fe52 to co55
    = 24        co55 to ni56
  */
  float aa, bb, cc, pp, qq, rr, ss, tt, uu, vv, a3, tempi, rhoi, qlam0z,
    h12w, h12, h12x, x1, denom, taufac, gamp, gamef, tau12, alph12,
    xlgfac, gamp14, xni, ytot;
  float x13, x14, x532, fact, co2, gamefx, gamefs, dgamma, blend_frac;
  x13 = 0.333333333;
  x14 = 0.25;
  /* x53    = 1.666666667; */
  x532 = 0.15625;
  /* x512   = 0.416666667; */
  fact = 1.25992104989487;
  co2 = x13*4.248719e3;
  gamefx = 0.1;
  gamefs = 0.4;
  dgamma = 1.0/(gamefs - gamefx);
  blend_frac = 0.05;
  ytot = 1.0/Abar;
  rr = rho*ytot;
  tempi = 1.0/temp;
  rhoi = 1.0/rho;
  /* equation 4 of dewitt et al becomes ztilde = z2bar + Zbar */
  pp = sqrt(rr*tempi*(z2bar + Zbar));
  /* equation 6 of dewitt et al */
  qlam0z = 1.88e8*tempi*pp;
  taufac = co2*pow(tempi, x13);
  qq = rr*Zbar;
  xni = pow(qq, x13);
  aa = 2.27493e5*tempi*xni;
  /*calculate individual screening factors */
  bb = z1*z2;
  gamp = aa;
  qq = fact*bb*zs13inv[jscreen];
  gamef = qq*gamp;
  /*equation 6 of itoh et al */
  tau12 = taufac*aznut[jscreen];
  qq = 1.0/tau12;
  alph12 = gamef*qq;
  /* limit alph12 to 1.6 to prevent unphysical behavior.*/
  /* this should really be replaced by a pycnonuclear reaction rate
     formula */
  if (alph12 > 1.6) {
    alph12 = 1.6;
    gamef = 1.6*tau12;
    qq = zs13[jscreen]/(fact*bb);
    gamp = gamef*qq;
  }
  /* weak screening regime */
  h12w = bb*qlam0z;
  h12 = h12w;
  /* intermediate and strong sceening regime */
  if (gamef > gamefx) {
    gamp14 = pow(gamp, x14);
    /* equation 25 of alastuey */
    /* coefficients from equation 16 of jancovici */
    cc = 0.896434*gamp*zhat[jscreen];
    cc -= 3.44740*gamp14*zhat2[jscreen];
    cc -= 0.5551*(log(gamp) + lzav[jscreen]);
    cc -= 2.996;
    a3 = alph12*alph12*alph12;
    /* equation 28 of alastuey */
    qq = 0.014 + 0.0128*alph12;
    rr = x532 - alph12*qq;
    ss = tau12*rr;
    /* equation 31 of alastuey */
    tt = -0.0098 + 0.0048*alph12;
    uu = 0.0055 + alph12*tt;
    vv = gamef*alph12*uu;
    /* equation 32 of alastuey */
    h12 = cc - a3*(ss + vv);
    /* equation 34 of alastuey */
    rr = 1.0 - 0.0562*a3;
    if (rr >= 0.77) {
      xlgfac = rr;
    }
    else{
      xlgfac = 0.77;
    }
    /* quantum correction of jancovici 1977 */
    h12 = log(xlgfac) + h12;
    /*intermediate regime */
    if (gamef <= gamefs) {
      rr = dgamma*(gamefs - gamef);
      ss = dgamma*(gamef - gamefx);
      /* store current values for blending with weak regime */
      h12x = h12;
      /* real linear ; gives same result as below */
      h12 = h12w*rr + h12x*ss;
    }
  }
  /*  machine limit the output */
  /*further limit to avoid the pyconuclear regime */
  if (h12 > 30.0)
    h12 = 30.0;
  if (h12 < 0.0)
    h12 = 0.0;
  scor = exp(h12);
  return;
}
//end usrscreen.c
 
//file usrrieman.c
/*
+ usrrieman - Riemann solver
     Description :
        Given  the left  state in  pvl and  right state  in pvr  of an
        arbitrary  discontinuity this  returns the  resolved state  in
        pres.
*/
void usrrieman(float pvl[], float pvr[], float pres[])
{
  int i;
  float dl, pl, ul, tel, dr, pr, ur, ter, cel, cer, cl, cr, p, u, d;
  /* Extract densities, pressures, velocities and temperatures */
  dl = pvl[iqd];
  pl = pvl[iqe];
  ul = pvl[iqu0+1];
  tel = pvl[iqte];
  dr = pvr[iqd];
  pr = pvr[iqe];
  ur = pvr[iqu0+1];
  ter = pvr[iqte];
  /* Linear guess */
  cel = usrss(pvl);
  cer = usrss(pvr);
  cr = dr*cer;
  cl = dl*cel;
  p = (cr*pl + cl*pr - cr*cl*(ur - ul))/(cr + cl);
  u = (cr*ur + cl*ul - (pr - pl))/(cr + cl);
  if (u > 0.0)
    d = cl*dl/(cl + (u - ul)*dl);
  else
    d = cr*dr/(cr - (u - ur)*dr);
  if (p > pl){
    /* Left wave is a shock */
    if (p > pr)
      /* Both waves are shocks (case 1) */
      shosho(dl, ul, pl, cl, pvl, dr, ur, pr, cr, pvr, d, u, p);
    else
      /* Left wave is shock, right wave is rarefaction (case 2) */
      shora(-1, dl, ul, pl, cl, pvl ,dr , ur, pr, cer, pvr, d, u,
	    p);
  }
  else{
    if (p > pr)
      /* Left wave is rarefaction, right wave is shock (case 3)*/
      shora(1, dr, ur, pr, cr, pvr, dl, ul, pl, cel, pvl, d, u, p);
    else
      /* Both waves are rarefactions (case 4) */
      rara(dl, ul, pl, cel, pvl, dr, ur, pr, cer, pvr, d, u, p);
  }
  /* Resolved state */
  pres[iqd] = d;
  pres[iqu0+1] = u;
  pres[iqe] = p;
  if (u > 0.0){
    for (i = 1; i <= niso; i++)
      pres[iqal0+i] = pvl[iqal0+i];
  }
  else{
    for (i = 1; i <= niso; i++)
      pres[iqal0+i] = pvr[iqal0+i];
  }
  return;
}
//end rieman.c

//file usrrara.c
/*
+ usrrara - Two rarefactions.
     Description :
        Returns the resolved state when both waves are rarefactions.
        dl  - density in left state
        ul  - velocity in left state
	pl  - pressure in left state
	cel - eulerian sound speed in left state
	pvl - primitive variables in left state
	dr  - density in right state
	ur  - velocity in right state
	pr  - pressure in right state
	cer - eulerian sound speed in right state
	pvr - primitive variables in left state
	d   - density in resolved state]
	u   - velocity in resolved state
	p   - pressure in resolved state
*/
void  rara(float dl, float ul, float pl, float cel, float pvl[],
  float dr, float ur, float pr, float cer, float pvr[], float &d,
  float &u, float &p)
{
  int iterl,iterr;
  float step = 0.01;
  float  dtl, utl, ptl, ctl, dil, uil, pil, dul, wl, dtr, utr, ptr, ctr,
    dir, uir, pir, dur, wr, cl, cr, pi, wli, wri;
  wl = ul - cel;
  if (wl > 0.0){
    /* Disturbance is on right of interface */
    d = dl;
    u = ul;
    p = pl;
    return;
  }
  wr = ur + cer;
  if (wr < 0.0){
    /* Disturbance is on left of interface */
    d = dr;
    u = ur;
    p = pr;
    return;
  }
  /* Initialise r-k integrator for left rarefaction */
  dtl = dl;
  utl = ul;
  ptl = pl;
  ctl = cel;
  iterl = 0;
  /* Initialise r-k integrator for right rarefaction */  
  dtr = dr;
  utr = ur;
  ptr = pr;
  ctr = cer;
  iterr = 0;
  /* r-k integrator for left and right rarefactions*/
  while (true){ 
    if (ptl > ptr){
      /* Integrate left rarefaction */
      dil = dtl;
      uil = utl;
      pil = ptl;
      wli = wl;
      dul = step*min(ptl/(dtl*ctl),ctl);
      usrrki(-1, dul, dtl, utl, ptl, ctl, pvl);
      iterl = iterl + 1;
      if (iterl >= 1000){
	printf("left state %11.4e %11.4e %11.4e %11.4e \n", dl, ul, pl,
	  cel);
	printf("right state %11.4e %11.4e %11.4e %11.4e \n", dr, ur, pr,
	  cer);      
	exit(EXIT_FAILURE);
      }
      /* Check for left sonic point */
      wl = utl - ctl;
      if (wl > 0.0){
	/*  Interpolate and exit */
	dul = wli/(wli - wl);
	d = dil + dul*(dtl - dil);
	u = uil + dul*(utl - uil);
	p = pil + dul*(ptl - pil);
	return;
      }
      else{
	/*  Integrate right rarefaction */
	dir = dtr;
	uir = utr;
	pir = ptr;
	wri = wr;
	dur = - step*min(ptr/(dtr*ctr),ctr);
	usrrki(1, dur, dtr, utr, ptr, ctr,pvr);
	iterr = iterr + 1;
	if (iterr >= 1000){
	  printf("left state %11.4e %11.4e %11.4e %11.4e \n", dl, ul, pl,
	    cel);
	  printf("right state %11.4e %11.4e %11.4e %11.4e \n", dr, ur,
	    pr, cer);      
	  exit(EXIT_FAILURE);
	}
	/* Check for right sonic point */
	wr = utr + ctr;
	if (wr < 0.0){
	  /* Interpolate and exit */
	  dur = wri/(wri - wr);
	  d = dir + dur*(dtr - dir);
	  u = uir + dur*(utr - uir);
	  p = pir + dur*(ptr - pir);
	  return;
	}
      }
      /* Check for cavitation */
      if ((ptl <= pfloor) && (ptr <= pfloor)){
	/* Cavitation - take current state and exit */
	if (utl > 0.0){
	  /* Take left state */
	  d = dtl;
	  u = utl;
	  p = ptl;
	}
	else if (utr < 0.0){
	  /* Take right state */
	  d = dtr;
	  u = utr;
	  p = ptr;
	}
	else{
	  /* Total cavitation */
	  d = dfloor;
	  u = 0.0;
	  p = pfloor;
	}
	return;
      }
      /* Check for intersection */
      if (utr <= utl){
	if (iterl < 1){
	  /* Shock - rarefaction*/
	  cl = dl*cel;
	  shora(-1, dl, ul, pl, cl, pvl, dr, ur, pr, cer, pvr, d, u, p);
	  return;
	}
	if (iterr < 1){
	  /* Rarefaction - shock */
	  cr = dr*cer;
	  shora(1, dr, ur, pr, cr, pvr, dl, ul, pl, cel, pvl, d, u, p);
	  break;
	}
	/* Interpolate and exit */
	pi = (ptr - pir)/(utr - uir);
	dul = (pir - pil + (uil - uir)*pi)/(ptl - pil -
	  (utl - uil)*pi);
	u = uil + dul*(utl - uil);
	p = pil + dul*(ptl - pil);
	if (u > 0.0){
	  /* Contact is on right of interface */
	  d = dil + dul*(dtl - dil);
	  break;
	}
	else{
	  /* Contact is on left of interface */
	  dur = (uil - uir + dul*(utl - uil))/
	    (utr - uir);
	  d = dir + dur*(dtr - dir);
	  return;
	}
      }
    }
  }
  return;
}
// end usrrara.c

//file usrshosho.c
/*
+ usrshosho - Two shocks.
     Description :
        Returns the resolved state when both waves are shocks.
        dl  - density in left state
        ul  - velocity in left state
	pl  - pressure in left state
	cl  - lagrangian sound speed in left state
	pvl - primitive variables in left state
	dr  - density in right state
	ur  - velocity in right state
	pr  - pressure in right state
	cr  - lagrangian sound speed in right state
	pvr - primitive variables in left state
	d   - density in resolved state]
	u   - velocity in resolved state
	p   - pressure in resolved state
*/
void shosho(float dl, float ul, float pl, float cl, float pvl[],
  float dr, float ur, float pr, float cr, float pvr[], float &d,
  float &u, float &p)
{
  int mxit = 100, iter;
  float linear = 0.2, tol = 1.0e-3;
  float pi, wl, wr, frac;
  pi = p;
  wl = - wave(p, dl, pl, cl, pvl);
  wr = wave(p, dr, pr, cr, pvr);
  iter = 0;
  do {
    pi = p;    
    p = (wr*pl - wl*pr + wl*wr*(ur - ul))/(wr - wl);
    wl = - wave(p, dl, pl, cl, pvl);
    wr = wave(p, dr, pr, cr, pvr);    
    iter++;
    if (iter > mxit){
      printf("shosho having trouble \n");
      printf("iter %d  pi %11.4e p %11.4e wl %11.4e %11.ee wr \n", iter,
        pi, p, wl, wr);
      break;
    }
  }
  while (fabs(p - pi) > tol*p);
  u = (wr*ur - wl*ul - (pr - pl))/(wr - wl);
  if (u > 0.0){
    /* Contact is on right of interface */
    if (ul + wl/dl > 0.0){
      /* Disturbance is on right of interface */
      d = dl;
      u = ul;
      p = pl;
    }
    else
      /* Left wave is on left of interface; */
      d = wl*dl/(wl - (u - ul)*dl);
  }
  else{
    /* Contact is on left of interface */
    if (ur + wr/dr < 0.0){
      /* Disturbance is on left of interface; */
      d = dr;
      u = ur;
      p = pr;
    }
    else
      /* Right wave is on right of interface; */
      d = wr*dr/(wr - (u - ur)*dr);
  }
  return;
}
//end usrshosho.c

//file usrshora.c
/*
+ usrshora - Shock and rarefaction.
     Description :
        Returns the resolved  state when one waves is a  shock and the
        other is a rarefaction.
	is  - left shock is = - 1, right shock, is = 1
	ds  - density ahead of shock
	us  - velocity ahead of shock
	ps  - pressure ahead of shock
	cs  - lagrangian sound speed ahead of shock
	pvs - primitive variables in shock state
	dr  - density ahead of rarefaction
	ur  - velocity ahead of rarefaction
	pr  - pressure ahead of rarefaction
	cer - eulerian sound speed ahead of rarefaction
	pvr - primitive variables in left state
	d   - density in resolved state]
	u   - velocity in resolved state
	p   - pressure in resolved state
*/
void shora(int is, float ds, float us, float ps, float cs, float pvs[],
  float dr, float ur, float pr,float cer, float pvr[], float &d,
  float &u, float &p)
{
  int iter;
  float step = 0.01;
  float cr, dt, ut, pt, ct, di, ui, pi, du, ws, wr, res, wri, resi;
  /* Check if shock - shock */
  if (is*(ur - us) > 0.0){
    ws = wave(pr, ds, ps, cs, pvs);
    res = is*(pr - ps)/(ur - us) - ws;
    if (res < 0.0){
      /* Shock - shock */
      cr = cer*dr;
      if (is == -1)
	shosho(ds, us, ps, cs, pvs, dr, ur, pr, cr, pvr, d, u, p);
      else
	shosho(dr, ur, pr, cr, pvr, ds, us, ps, cs, pvs, d, u, p);
    }
    return;
  }
  wr = ur - is*cer;
  if (is*wr >= 0.0){
    /* Disturbance is on rarefaction side of interface */
    d = dr;
    u = ur;
    p = pr;
    return;
  }
  /* Initialise r-k integrator for rarefaction */
  dt = dr;
  ut = ur;
  pt = pr;
  ct = cer;
  iter = 0;
  while (true){
    /*  r-k integrator for rarefaction */
    di = dt;
    ui = ut;
    pi = pt;
    wri = wr;
    resi = res;
    du = is*step*min(pt/(dt*ct),ct);
    usrrki(-is, du, dt, ut, pt, ct, pvr);
    /* Check for sonic point */
    wr = ut - is*ct;
    if (is*wr > 0.0){
      /* Sonic point - interpolate and exit */
      du = wri/(wri - wr);
      d = di + du*(dt - di);
      u = ui + du*(ut - ui);
      p = pi + du*(pt - pi);
      return;
    }
    ws = wave(pt, ds, ps, cs, pvs);
    if (is*(us + is*ws/ds) < 0.0){
      /* Take shock state and exit */
      d = ds;
      u = us;
      p = ps;
      return;
    }
    if (is*(ut - us) > 0.0){
      iter = iter + 1;
      if (iter >= 1000){
	if (is < 0)
	  printf("left shock failure %11.4e %11.4e %11.4e %11.4e \n",
		 ds, us, ps, cs);
	else
	  printf("right shock failure %11.4e %11.4e %11.4e %11.4e \n",
		 ds, us, ps, cs);
	return;
      }
      res = is*(pt - ps)/(ut - us) - ws;
      if (res < 0.0){
	/* Interpolate and exit */
	if (iter < 2)
	  /* State close to switch in velocity */
          du = (is*(pi - ps) + ws*(us - ui))/(ws*(ut -ui) -
	    is*(pt - pi));
	else
	  /* State not close to switch in velocity */
	  du = resi/(resi - res);
	u = ui + du*(ut - ui);
	p = pi + du*(pt - pi);
	if (is*u > 0.0)
	  /* Take rarefaction density */
	  d = di + du*(dt - di);
	else{
	  /* Take shock density */
	  ws = wave(p, ds, ps, cs, pvs);
	  d = ws*ds/(ws - is*(u - us)*ds);
	}
	return;
      }
    }
    /* Check pressure */
    if (pt < ps){
      /* Pressure too low - take shock pressure and exit */
      du = (ps - pi)/(pt - pi);
      u = ui + du*(ut - ui);
      p = ps;
      if (is*u < 0.0)
	/* Take shock state */
	d =  ds;
      else
	/* Take rarefaction state */
	d = di + du*(dt - di);
      return;
    }
  }
}
//end usrshora.c

//file usrrki.c
/*
+ usrrki - rk integrator.
     Description :
        This routine performs one  second order runge-kutta step along
        the rarefaction curve.
	is  - left rarefaction is = -1, right rarefaction is = 1
	du  - step size in velocity
	d   - density
	u   - velocity
	p   - pressure
	c   - eulerian sound speed
	pvr - array of primitive variables for rarefaction state

*/
void usrrki(int is, float du, float &d, float &u, float &p, float &c,
  float pvr[])
{
  int i;
  float pt[MAXQ], kd1, kd2, kp1, kp2, tv;
  /* Copy of pvr */
  for (i = 0; i < nq; i++)
    pt[i] = pvr[i];
  kp1 = is*du*d*c;
  kd1 = kp1/(c*c);
  tv = d + kd1;  
  pt[iqd] = tv;
  pt[iqe] = p + kp1;
  c = usrss(pt);
  kp2 = is*du*tv*c;
  kd2 = kp2/(c*c);
  d += 0.5*(kd1 + kd2);
  u += du;
  p += 0.5*(kp1 + kp2);
  pt[iqd] = d;
  pt[iqe] = p;  
  c = usrss(pt);
  return;
}
//end usrrki.c
 
//file usrwave.c
/*
+ usrwave - slope of p - u curve.
     Description :
        Returns the "wave speed" required  by the Riemann solver. This
        is the slope of the p -  u curve for shock states and is equal
        to mass flux sthrough the wave.
	p  - pressure in resolved state
	pi - pressure in adjacent gas
	di - density in resolved state
	ci - lagrangian sound speed in resolved state
	pa - array of "background" primitive variables
*/
float wave(float p, float di, float pi, float ci, float pa[])
{
  int nit, mxit = 20;
  float eps = 0.1;
  float Abar, Ye, temp, tv1, tv2, ei, vi, v1, v2, fo, f, dv, dvo;
  /* Linear case:  return lagrangian sound speed */
  if (fabs(p - pi) < eps*pi)
    return(ci);
  Abar = usrAbar(pa);
  Ye = usrYE(pa);  
  temp = usrpg2te(pa[iqte], di, pi, Abar, Ye);
  ei = usrein(temp, di, Abar, Ye);
  vi = 1.0/di;
  tv1 = 2.0/(p + pi);
  tv2 =  vi + tv1*ei;
  v1 = vi;
  v2 = 0.75*vi;
  temp = usrpg2te(temp, di, p, Abar, Ye);
  fo = vi - tv2 + tv1*usrein(temp, 1.0/vi, Abar, Ye);
  nit = 0;
  dvo = v2 - v1;
  do {
    temp = usrpg2te(temp, 1.0/v2, p, Abar, Ye);    
    f = v2 - tv2 + tv1*usrein(temp, 1.0/v2, Abar, Ye);
    dv = - f*dvo/(f - fo);
    v2 += dv;
    nit++;
    if (nit > mxit){
      cout << "Convergence failure in wave" << '\n';
      break;
    }      
    fo = f;
    dvo = dv;
  }
  while (fabs(dv) > 1.0e-04*v2);
  return(sqrt((p - pi)/(vi - v2)));
}
// usrwave.c

//file usrsteady.c
/*
+ usrsteady - Sources for steady solution
     Description :
        Returns  the  sources for  a  steady  detonation using  a  gsl
        integrator.
*/
void usrsteady(float pa[], float so[], float dt)
{
  int i, nsubl;
  float pat[MAXQ], qao[MAXQ], qan[MAXQ];
  bool old_isothermal;
  old_isothermal = isothermal;
  isothermal = false;
  /* Get old conservatives */ 
  usrqp(pa, qao);
  /* Temporary copy of pa */
  for (i = 0; i < nq; i++)
    pat[i] = pa[i];
  /* gsl step (update pat) */
  nsubl = nsub;
  modgsl(nsubl, pat, dt);
  /* Get conserved variables after step */	
  usrqp(pat, qan);
  /* Source terms */
  for (i = 0; i < nq; i++)
    so[i] += (qan[i] - qao[i])/dt;
  isothermal = old_isothermal;
  return;
}
//end usrsteady.c

//file usratmos.c
/*
+ usratmos - Calculate hydrostatic atmosphere solution
     Description :
        Calculates  the  density  in  a  hydrostatic  atmosphere  with
        surface density rho_0, temperature T_0 on a star with mass M_0
        and  radius R_0.   The results  is written  as an  xq file  in
        filename.xq.   M_shell is  the  desired mass  of the  shell.It
        returns false if it fails.
*/
bool usratmos(float rho_0, float rho_min, float T_0, float M_0,
  float R_0, float M_shell, char mode[], float pa[])
{
  FILE *outM, *outg;
  char atmos[SZNAM];
  int i;
  bool done, adiabatic;
  float z, dz, rho, grav, M, M_n, dM_1, dM_2, Abar, Ye, drho_1, rho_n,
    drho_2, pg, pg_n, dpg_1, dpg_2, ss2, area, R, R_n, T, T_n, pt[MAXQ];
  if (strcmp(mode, "ADIABATIC") == 0)
    adiabatic = true;
  else
    adiabatic = false;  
  strcpy(atmos, monpref);
  strcat(atmos, "_M_x.xq");
  if (!(outM=fopen(atmos, "w")))
    return(false);
  strcpy(atmos, monpref);
    strcat(atmos, "_grav_x.xq");
  if (!(outg=fopen(atmos, "w")))
    return(false);  
  fprintf(outM, "%d %d \n", 1, 10);
  fprintf(outM, "M \n");
  fprintf(outM, "z M \n");
  fprintf(outg, "%d %d \n", 1, 10);
  fprintf(outg, "g \n");
  fprintf(outg, "z g \n");  
  for (i = 0; i < nq; i++)
    pt[i] = pa[i];
  Abar = usrAbar(pt);
  Ye = usrYE(pt);  
  rho = rho_0;
  T = T_0;
  pg = usrpg(T_0, rho, Abar, Ye);
  M = 0.0;
  R = R_0;  
  z = 0.0;
  grav = phygconst*(M_0 + M)/(R*R);  
  dz = R_0*1e-4;
  done = false;
  float zp = 0;
  while (!done){
    fprintf(rhofilex, "%e %e \n", z, rho);
    fprintf(tempfilex, "%e %e \n", z, T);    
    fprintf(outM, "%e %e \n", z, M);
    fprintf(outg, "%e %e \n", z, grav);        
    pt[iqd] = rho;
    pt[iqe] = pg;
    pt[iqte] = T;
    grav = phygconst*(M_0 + M)/(R*R);
    area = 4*gkpi*R*R;
    if (adiabatic)
      ss2 = usrss(pt);
    else
      ss2 = usrssit(pt);
    ss2 *= ss2;
    dpg_1 = - grav*rho*dz;
    drho_1 = dpg_1/ss2;    
    pg_n = pg + dpg_1;
    rho_n  = rho + drho_1;
    if (adiabatic)
      T_n = usrpg2te(T, rho_n, pg_n, Abar, Ye);
    else
      T_n = T_0;
    pt[iqd] = rho_n;
    pt[iqe] = pg_n;
    pt[iqte] = T_n;
    dM_1 = rho*area*dz;
    R_n = R + dz;    
    if (adiabatic)
      ss2 = usrss(pt);
    else
      ss2 = usrssit(pt);    
    ss2 *= ss2;
    grav = phygconst*(M_0 + M + dM_1)/(R_n*R_n);
    area = 4*gkpi*R_n*R_n;    
    dpg_2 = - grav*rho_n*dz;
    drho_2 = dpg_2/ss2;
    dM_2 = rho_n*area*dz;
    pg += 0.5*(dpg_1 + dpg_2);    
    rho += 0.5*(drho_1 + drho_2);
    M += 0.5*(dM_1 + dM_2);
    if (adiabatic)
      T = usrpg2te(T_n, rho, pg, Abar, Ye);
    else
      T = T_0;
    z += dz;
    R += dz;
    done = (M > M_shell) || (rho < rho_min);
  }
  usrmyclose(outM);
  usrmyclose(outg);  
  cout << std::scientific;
  cout << std::setprecision(6);
  cout << "shell height " << z << "  Mass " << M << '\n';
  cout << "density " << rho << " temp " << T << " sound speed "
       << sqrt(ss2) << '\n';
  return(true);
}
//end usratmos.c

//file usrcore.c
/*
+ usratmos - Calculate hydrostatic core  solution
     Description :
        Calculates  the  density   as  a  function  of   radius  in  a
        hydrostatic WD  core with  central density  rho_0, temperature
        T_0.The results  is written as  an xq file in  filename.xq. It
        returns false if it fails.
*/
bool usrcore(float rho_0, float rho_min, float T_0, char mode[],
  float pa[])
{
  FILE *outM;
  char core[SZNAM];
  int i;
  bool done, adiabatic;
  float pt[MAXQ], Abar, Ye, rho, T, pg, z, dz, grav, area, ss2, dpg_1,
    drho_1, pg_n, rho_n, T_n, M, dM_1, dpg_2, drho_2, dM_2;
  if (strcmp(mode, "ADIABATIC") == 0)
    adiabatic = true;
  else
    adiabatic = false;  
  strcpy(core, monpref);
  strcat(core, "_M_x.xq");
  if (!(outM=fopen(core, "w")))
    return(false);
  fprintf(outM, "%d %d \n", 1, 10);
  fprintf(outM, "M \n");
  fprintf(outM, "z M \n");
  for (i = 0; i < nq; i++)
    pt[i] = pa[i];
  Abar = usrAbar(pt);
  Ye = usrYE(pt);  
  rho = rho_0;
  T = T_0;
  pg = usrpg(T_0, rho, Abar, Ye);
  M = 0.0;
  z = 0.0;
  dz = 1.0e4;
  grav = 0.0;
  done = false;
  while (!done){
    fprintf(rhofilex, "%e %e \n", z, rho);
    fprintf(tempfilex, "%e %e \n", z, T);    
    fprintf(outM, "%e %e \n", z, M);
    pt[iqd] = rho;
    pt[iqe] = pg;
    pt[iqte] = T;
    if (z < 2.0*dz)
      grav = phygconst*4*gkpi*z*rho_0/3.0;
    else
      grav = phygconst*M/(z*z);      
    area = 4*gkpi*z*z;
    if (adiabatic)
      ss2 = usrss(pt);
    else
      ss2 = usrssit(pt);
    ss2 *= ss2;
    dpg_1 = - grav*rho*dz;
    drho_1 = dpg_1/ss2;
    pg_n = pg + dpg_1;
    rho_n  = rho + drho_1;
    if (adiabatic)
      T_n = usrpg2te(T, rho_n, pg_n, Abar, Ye);
    else
      T_n = T_0;
    pt[iqd] = rho_n;
    pt[iqe] = pg_n;
    pt[iqte] = T_n;
    dM_1 = rho*area*dz;
    if (adiabatic)
      ss2 = usrss(pt);
    else
      ss2 = usrssit(pt);    
    ss2 *= ss2;
    z +=dz;
    grav = phygconst*(M + dM_1)/(z*z);
    area = 4*gkpi*z*z;    
    dpg_2 = - grav*rho_n*dz;
    drho_2 = dpg_2/ss2;
    dM_2 = rho_n*area*dz;
    pg += 0.5*(dpg_1 + dpg_2);    
    rho += 0.5*(drho_1 + drho_2);
    M += 0.5*(dM_1 + dM_2);
    if (adiabatic)
      T = usrpg2te(T_n, rho, pg, Abar, Ye);
    else
      T = T_0;
    done = (rho < rho_min);
  }
  usrmyclose(outM);
  cout << std::scientific;
  cout << std::setprecision(6);
  cout << "radius " << z << "  Mass " << M << '\n';
  cout << "density " << rho << " temp " << T << " sound speed "
       << sqrt(ss2) << '\n';
  return(true);
}
//end usrcore.c

//file usrshockpos.c
/*
+ usrshockpos - Calculate shock position.
     Description :
        Returns  1D  shock  position   by  searching  from  the  right
        boundary.
*/
float usrshockpos(bool verbose, int dir, float nx, float qa[], float pa[]){
  int id, ilev;
  bool found = false;
  float x, dx, xc[MAXD], rho, xlmin, xrmax, pos, orho;
  struct cell *cptr;
  for (id = 1; id <= nd; id++)
    xc[id-1] = 0.0;
  dx = lev[nlevs-1]->dxbar[dir-1];
  if (nd == 1)
    dir = 1;
  else{
    if (dir == 1)
      xc[1] = nx;
    else
      xc[0] = nx;
  }
  x = xr[dir-1] - 0.5*dx;    
  xc[dir-1] = x;
  usrqis(xc, qa);
  xlmin = big;
  xrmax = -big;
  for (ilev = 0; ilev < nlevs; ilev++){
    xlmin = min(xlmin, lev[ilev]->xl[dir-1]);
    xrmax = max(xrmax, lev[ilev]->xr[dir-1]);		
  }
  x = xrmax - 0.5*dx;
  orho = qa[iqd];
  while (x > xlmin){
    xc[dir-1] = x;     
    cptr = find_fcell(ilev, xc);
    if (cptr != NULL){
      if (!cptr->halo){
	dx = lev[ilev]->dxbar[dir-1];
	rho = cptr->qa[iqd];	
	if (rho > 1.1*qa[iqd]){
	  if (rho > orho)
	    orho = rho;
	  else{
	    usrpq(cptr->qa, pa);
	    found = true;
	    break;
	  }
	}
      }
    }
    x -= dx;
  }
  if (found)
    pos = x;
  else
    pos = -big;
  if (verbose && mype == 0){
    cout.setf(ios::scientific);
    cout << setprecision(4);
    cout << "dir " << dir << " nx " << nx << " pos " << pos
      << " pg " << pa[iqe] << " D " <<
      pa[iqu0+dir]/(1.0 - qa[iqd]/pa[iqd]) <<'\n';    
  }
  return(global_max(pos));
}
//end usrshockpos.c

//file usrsonpos.c
/*
+ usrsonpos - Calculate sonic point position.
     Description :
        Returns 1D  sonic point position  by searching from  the right
        boundary. xs is shock position, detv is detonation velocity;
*/
float usrsonpos(float xs, float detv){
  int id, ilev, dir, nx;
  bool found = false;
  float x, dx, rho, xc[MAXQ], qa[MAXQ], pa[MAXQ], xlmin,
    xrmax, pos;
  struct cell *cptr;
  for (id = 1; id <= nq; id++)
    xc[id-1] = 0.0;
  dx = lev[nlevs-1]->dxbar[0];  
  xlmin = big;
  xrmax = -big;
  for (ilev = 0; ilev < nlevs; ilev++){
    xlmin = min(xlmin, lev[ilev]->xl[0]);
    xrmax = max(xrmax, lev[ilev]->xr[0]);		
  }
  x = xs - 10.0*dx;
  while (x > xlmin){
    xc[0] = x;     
    cptr = find_fcell(ilev, xc);
    if (cptr != NULL){
      if (!cptr->halo){
	dx = lev[ilev]->dxbar[0];
	usrpq(cptr->qa, pa);
	if (mabs((pa[iqu0+1] - detv))/usrss(pa) >= 1.0){
	  found = true;
	  break;
	}
      }
    }
    x -= dx;
  }
  if (found)
    pos = x;
  else
    pos = -big;
  return(global_max(pos));
}
//end usrsonpos.c

//file usrpfcell.c
/*
+ usrpfcell - Prints cell solution.
     Description :
        Prints the solution at the cell cptr).
*/
void usrpfcell(cell *cptr)
{
  int i, imax;
  float pa[MAXQ], ydot[MAXSPECIES+1], en, maxen;
  bool old_isothermal;  
  usrpq(cptr->qa, pa);
  old_isothermal = isothermal;
  isothermal = false;
  usrreact(pa, ydot);
  isothermal = old_isothermal;  
  cout << setprecision(8);
  maxen = -1.0;
  imax = -1;
  for (i = 0; i < nq; i++){
    cout << "i" << setw(3) << i << " qa " <<  setw(15) <<
      cptr->qa[i] <<  setw(3) << " pa " <<  setw(15) << pa[i];
    if ((i > iqal0) && (i < iqal0+nal)){
      en = - ydot[i-iqal0-1]*Qi[i-iqal0-1];
      if (fabs(en) > maxen){
	maxen = fabs(en);
	imax = i-iqal0-1;
      }
      cout << " xa " << setw(14) << pa[i]*Ai[i-iqal0-1] << setw(3) <<
	" en " << setw(15) << en << " " << ionam[i-iqal0-1];
    }
    else if (i == iqd)
      cout << " density";
    else if (i == iqu0+1)
      cout << " x velocity";
    else if ((nd == 2) && (i == iqu0+2))
      cout << " y velocity";    
    else if (i == iqe)
      cout << " pressure";
    else if (i == iqal0+nal)    
      cout << " React_prog ";
    else if (i == iqte)
      cout << " Temp";
    cout << '\n';
  }
  cout << "Abs max energy rate " << maxen << " for " <<
    ionam[imax] <<'\n';
  en = 0.0;
  for (i = 0; i < niso; i++)
    en -= pa[iqal0+1+i]*Qi[i];
  cout << "Energy generation rate " << ydot[niso]*usrcv(pa) << '\n';;
  cout<< "binding energy " << en << '\n'; 
  return;
}
//end usrpfcell.c

//file usrinitintp.c
/*
+ usrinitintp - Reads initial denity and temperature.
     Description :
        This is called by the command ATMINIT and reads the density
        and temperature in an atmosphere created by the ATMOS
        command. The files are of the form filename_ad/iso_rho/T.xq.
*/
bool usrinitintp(char filename[])
{
  int i, j, count;
  float x, tv;
  char name[SZNAM], tail[SZNAM], str1[SZNAM], str2[SZNAM], str3[SZNAM];
  FILE *infile;
  strcpy(name, filename);  
  strcat(name,"_rho.xq");      
  if (!(infile=fopen(name, "r")))
    return(false);
  safe_fscanf(infile, "%d %d", &i, &j);
  safe_fscanf(infile, "%s", str1);
  safe_fscanf(infile, "%s %s", str2, str3);
  count = 0;
  /* Count lines */
  while (!(safe_fscanf(infile,"%lf %lf", &x, &tv) == EOF)){
    count++;
  }
  /* Do rho */
  intpts = count;
  rewind(infile);
  safe_fscanf(infile, "%d %d", &i, &j);
  safe_fscanf(infile, "%s", str1);
  safe_fscanf(infile, "%s %s", str2, str3);
  for (i = 0; i < count; i++){
    safe_fscanf(infile,"%lf %lf", &x, &tv);
    if (i == 1)
      intdx = x;
    if (i == count - 1)
      intxmax = x;    
    intrho[i] = tv;
  }
  /* Do T */
  usrmyclose(infile);  
  infile = NULL;
  strcpy(name, filename);  
  strcat(name,"_T.xq");
  if (!(infile=fopen(name, "r")))
    return(false);
  safe_fscanf(infile, "%d %d", &i, &j);
  safe_fscanf(infile, "%s", str1);
  safe_fscanf(infile, "%s %s", str2, str3);
  for (i = 0; i < count; i++){
    safe_fscanf(infile,"%lf %lf", &x, &tv);
    inttemp[i] = tv;
  }  
  usrmyclose(infile);
  infile = NULL;
  return(true);
}    
//end usrinitintp.c

// file intp.c
/*
+ intp - interpolates on an array
     Description :
        Locates  array position  in a  table and  interpolates between
        variable value bracketing values in 1D
  Inputs
	- quan 	- variable value
	- qmin 	- lower limit of array
	- qmax 	- upper limit of array
	- nq	- number of array indices (intervals = nq-1)
	- dq 	- array step value
  Outputs
	- i	- lower array index bracketing the variable value
	- fraci	- fractional value between bracketing indices

  Addition: CJW 13/1/15
*/
void intp(float quan, float qmin, float qmax, float nq, float dq,
  int &i, float &fraci)
{
  if ( quan <= qmin+0.00001 ){
    if(debug)
      cout << "! Interpolation value ("<<quan<<") near lower limit. Using lower limit: " << qmin << '\n';
    i = 0;
    fraci = 0;
    return;
  }else if ( quan >= qmax-0.00001 ){
    if(debug)
      cout << "! Interpolation value ("<<quan<<") near upper limit. Using upper limit: " << qmax << '\n';
    i = nq-2;
    fraci = 1.0;
    return;
  }else{
    i = (int)((quan - qmin)/dq);
    fraci =  ((quan - qmin)-(i*dq))/dq;
  }
  return;
}
//end intp.c

//file usrgsl.c
/*
+ usrgsl - Stiff integrator.
     Description :
        Integrates the  solution over time  dt starting with  pa[]. At
        the end,  pa[] is  updated.  Uses the  gsl.  Returns  false on
        failure.
*/
bool usrgsl(cell *cptr, float pa[], float dt)
{
  int i, nvar, nsubl, mintimestep, maxtimestep;
  size_t dim;
  float tl = 0.0,  epsabs = gslatol, epsrel = gslrtol, y[MAXSPECIES+1],
    ydot[MAXSPECIES+1], srate, dtl;
  struct gsl_param *param_ptr;
  gsl_odeiv2_driver *driver;  
  param_ptr = new gsl_param;
  /* Number of eqns */
  nvar = niso + 1;
  param_ptr->nvar = nvar;
  dim = (size_t)nvar;
  /* Check for large rates */
  usrreact(pa, ydot);
  srate = -1.0;
  /* Isotope rate */      
  for (i = 0; i < niso; i++)
    srate = max(srate, fabs(Ai[i]*ydot[i]));		// ydot*A = dX/dt rate of mass fraction
  /* Energy rate */
  srate = max(srate, fabs(ydot[niso]/pa[iqte]));
  dtl = dt;
  if (dt*srate < cns1)
    /* Rate not large */
    nsubl = nsub;
  else{
    /* Large rate */
    if (cptr->refined) {
      // Set source to zero
      dtl = 0.0;
      return true;		// added by BM
    }
    else if (curlev < solev){
      // Reduce evolution time for coarse grids
      dtl = max(0.0, cns2/srate);
    }
    // Set local nsub 
    nsubl = nsub + int(dtl*srate/cns1);  
  }
  /* check for min/max allowed time step */
  if (gslhmax > 0.0)
    if (dtl/nsubl > gslhmax)
      dtl = gslhmax*nsubl;
  if (gslhmin > 0.0) {
    if (dtl/nsubl < gslhmin)
      dtl = gslhmin*nsubl;
    if (dtl > dt) {
      dtl = dt;	// cannot go greater than dt
    }
  }
  /* Put these results in cptr */
  if (cptr != NULL){
    cptr->ilev = curlev;
    cptr->dtrat = dtl/dt;
    cptr->nsub = nsubl;
  }
  /* Abundances */
  for (i = 0; i < niso; i++)
    y[i] = pa[iqal0+1+i];
  /* Complete solution */
  for (i = 0; i < nq; i++)
    param_ptr->pa[i] = pa[i];
  /* Temperature (in units of 1.0e9 K) */
  y[niso] = 1.0e-9*pa[iqte];
  /* Define ode system  */
  gsl_odeiv2_system sys = {usrgslreact, usrjac, dim, param_ptr};
  /* Allocate ode driver */
  if (strcmp(ode, "RK2") == 0)
    /* Explicit rk2 integrator */    
    driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk2,
      dtl/nsubl, epsabs, epsrel);
  else if (strcmp(ode, "RK2IMP") == 0)
    /* Implicit rk2 integrator */
    driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk2imp,
      dtl/nsubl, epsabs, epsrel);
  else
    /* Implicit MSBDF integrator */
    driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf,
      dtl/nsubl, epsabs, epsrel);
  /* check for minimum time step requirements, gslhmin cannot be > dt */
  if (gslhmin > 0.0) {
    if (gslhmin < dtl/nsubl)
      mintimestep = gsl_odeiv2_driver_set_hmin(driver, gslhmin);
    else
      mintimestep = gsl_odeiv2_driver_set_hmin(driver, dtl/nsubl);
  }
  if (gslhmax > 0.0) {
    maxtimestep = gsl_odeiv2_driver_set_hmax(driver, gslhmax);
  }
  /* Integrate from tl to tl+dt */
  int status = gsl_odeiv2_driver_apply(driver, &tl, tl+dtl, y);
  if (status != GSL_SUCCESS){
    printf ("error in ode integrator, return value=%d\n", status);
    return(false);
  }
  gsl_odeiv2_driver_free(driver);
  /* Delete param_ptr */
  delete param_ptr;
  /* Update pa */
  for (i = 0; i < niso; i++)
    pa[iqal0+1+i] = y[i];
  /* Rescale temperature */
  pa[iqte] = 1.0e9*y[niso];
  return(true);
}
//end usrgsl.c

//file usrgslreact.c
/*
+ usrgslreact - gsl wrapper for reaction rates
     Description :
        Right hand side  function for gsl ode  integrator.  It returns
        the reaction sources  in ydot using y as  input and param_ptr
        to get pa[].
*/
int usrgslreact(float t, const float y[], float ydot[],
  void *param)
{
  (void)(t); /* avoid unused parameter warning */
  int i;
  float pa[MAXQ];
  struct gsl_param *param_ptr;
  param_ptr = (struct gsl_param*)param;
  /* Zero ydot for safety */
  for (i = 0; i < niso + 1; i++)
    ydot[i] = 0.0;
  /* Fill up to iqal0 from param_ptr->pa */
  for (i = 0; i <= iqal0; i++)
    pa[i] = param_ptr->pa[i];
  /* Use y for abundances and temperature */
  for (i = 0; i < niso; i++)
    pa[iqal0+1+i] = y[i];
  if (hydrostatic)
    /* Fix temperature */
    pa[iqte] = param_ptr->pa[iqte];
  else
  /* Rescale temperature */    
    pa[iqte] = 1.0e9*y[niso];
  /* Fix up solution for intermediate steps */
  usrfxp(pa);
  /* Get abundance and temperature source terms */
  usrreact(pa, ydot);
  /* Temperature in units of 1.0e9 K */
  ydot[niso] *= 1.0e-9;
  return GSL_SUCCESS;
}
//end usrgslreact.c

//file usrjac.c
/*
+ usrjac - Jacobian 
     Description :
        Returns  the  Jacobian  in  dfdy[]  given  the  abundances  in
        y[]. The reaction   rates, RR are in the  param structure.  The
        Jacobian   is   evaluated    by   finite   differences   using
        usrgslreact.x */
int usrjac (float t, const float y[], float *dfdy, float dfdt[],
  void *param)
{
  bool negy;
  int nvar, i, j;
  size_t dim;
  //  float rel = 3.162278e-8, ax = 1.0;
    //float rel = 1e-6, ax = 1.0;
    float rel = 1e-8, ax = 1.0;
  float yt[MAXSPECIES+1], h, fminus[MAXSPECIES+1], fplus[MAXSPECIES+1];
  struct gsl_param *param_ptr;
  param_ptr = (struct gsl_param*)param;
  jaccall++;
  /* No of variables */
  nvar = param_ptr->nvar;
  dim = (size_t)nvar;
  /* Zero dfdt for safety */
  for (i = 0; i < nvar; i++)
    dfdt[i] = 0.0;
  /* Set up gsl matrix */
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, dim, dim);
  gsl_matrix *m_ptr = &dfdy_mat.matrix;
  /* Copy of y */
  for (i = 0; i < nvar; i++)
    yt[i] = y[i];
  /* Central difference for Jacobian */
  for (j = 0; j < nvar; j++){
    /* Increment */
    h = rel*max(y[j], ax);
    /* Right difference */
    yt[j]  = y[j] + h;
    usrgslreact(t, yt, fplus, param);
    /* Left difference */
    yt[j]  = y[j] - h;
    negy = yt[j] <= 0.0;
    if (!negy)
      usrgslreact(t, yt, fminus, param);
    /* Restore yt */
    yt[j] = y[j];
    /* Set jth column of Jacobian */
    for (i = 0; i < nvar; i++){
      if (negy)
	m_ptr->data[i*nvar+j] = fplus[i]/h;
      // gsl_matrix_set (m_ptr, i, j, fplus[i]/h);
      else
	m_ptr->data[i*nvar+j] = 0.5*(fplus[i] - fminus[i])/h;	
      // gsl_matrix_set (m_ptr, i, j, 0.5*(fplus[i] - fminus[i])/h);
    }
  }
  return GSL_SUCCESS;;
}
//end usrjac.c



//file usrpref.c
/*
+ usrpref - primitive reference solution.
     Description:
        Provides  a   primitive  reference  solution   from  primitive
        solution in mcell mcptr.
*/
void usrpref(mcell *mcptr, float pref[])
{
  int i, id;
  float rho, ss;
  for (i = 0; i < nq; i++)
    pref[i] = mabs(mcptr->qa[i]);
  rho = mcptr->qa[iqd];
  ss = sqrt(mcptr->qa[iqe]/rho);
  for (id = 1; id <= nd; id++)
    pref[iqu0+id] = mabs(mcptr->qa[iqu0+id]) + ss;
  if (keps){
    pref[iqk] = 0.01*pref[iqe]/rho;
    pref[iqeps] = pref[iqk];
  }
  for (i = iqal0 + 1; i <= iqal0 + nal; i++)
    pref[i] = 1.0;
  if (grav)
    pref[iqpot] = big;
  return;
}
//end usrpref.c






//file modsol.c
/*
+ modsol - Find detonation speed
       Description :
       Use shooting to find detonation  speed by bisection if bisect =
       true, Otherwise  it just  calculations solution  for detonation
       speed D. Returns true if it succeeds.       
*/
bool modsol(struct cell *pcptr, bool verbose, bool bisect)
{
  int i, nstep, nstepl, nsubl, cnt;
  float pr[MAXQ], pa[MAXQ], pao[MAXQ], srate, dt, omach, mach, imach,
    cnsl, ss, u, phi0, phi, ydot[MAXSPECIES+1], xc[MAXD], dU_r, dU_l;
  char scmd[SZNAM];
  struct cell *cptr;
  cptr = new cell;
  while (D_r - D_l > par[0]*D){
    /* Get upstream solution */
    for (i = 0; i < nd; i++)
      xc[i] = 0.0;
    xc[0] = 1.0;
    usrqis(xc, cptr->qa);
    usrpq(cptr->qa, pr);    
    /* Set initial abundances */
    for (i = 0; i < niso; i++)
      pa[iqal0+1+i] = pr[iqal0+1+i];
    usrshock(pr, pr[iqd], pr[iqte], D, pa[iqd], pa[iqte], pa[iqu0+1],
      pa[iqe]);
    //*(int*)0 = 0;
    /* Copy of solution */
    for (i = 0; i < nq; i++)
      pao[i] = pa[i];
    omach = pao[iqu0+1]/usrss(pao); // old mach no;    
    pa[iqal0+nal] = 0.0; // x
    t = 0.0;
    nstep = 0;
    nstep = 0;
    cnsl = cns; // initial source courant no
    modsrc(pao, ydot);
    ss = usrss(pao);
    u = pao[iqu0+1];
    phi0 = 1.0e9*ydot[nal]*(ss*ss - u*u)/u;
    omach = u/ss; // old mach no;
    imach = omach; // initial mach no
    if (mon){
      usrclosexq(monmod, 1);
      strcpy(scmd, "rm ");
      strcat(scmd, monpref);
      strcat(scmd, "*.xq");
      sys_cmd(scmd);
      usropenxq(monmod, 1);
      pr[iqu0+1] = D;
      usrqp(pr, pcptr->qa);
      usrmon(monmod);
      usrqp(pa, pcptr->qa);
      usrmon(monmod);      
    }
    /* Integrate solution */
    dt = 1e-16/u;
    while (nstep < int(par[1])){
      /* Timestep */
      modsrc(pa, ydot);
      //      usrreact(pa, ydot);
      srate = -1.0;
      /* Isotope rate */      
      for (i = 0; i < niso; i++)
	srate = max(srate, fabs(Ai[i]*ydot[i]));
      /* Energy rate */
      srate = max(srate, fabs(1e9*ydot[niso]/pa[iqte]));
      //      srate = max(srate, fabs(ydot[niso]/pa[iqte]));
      //      dt = cnsl/srate;
      if (dt*srate < cns1)
	/* Rate not large */
	nsubl = nsub;
      else
	/* Set local nsub */
	nsubl = nsub + int(dt*srate/cns1);    
      /* Put nsubl in cptr */
      pcptr->nsub = nsubl;
      cout << std::scientific;
      if ((nstep == 0) && verbose){
	mach = pa[iqu0+1]/usrss(pa);	
        cout << "initial mach " << mach << " dt " << dt << " dx " <<
	  dt*pa[iqu0+1] << " phi " << phi0 << " cnsl " << cnsl << '\n';
      }
      cnt = 0;

      while (!modgsl(nsubl, pa, dt)){
	/* modgsl failure - reduce step */
	dt *= 0.5;		// BM
	if (verbose)
	  cout << "modgsl failure reduced dt " << dt << '\n';	
	cnt++;
	if (cnt > 100){
  	  if (verbose)
	    cout << "modgsl failure after 100 tries " << dt << '\n';	
	  cout << " x " << pa[iqal0+nal] << " dt " << dt << " phi " <<
	    phi <<" nstep " << nstep << " cnsl " << cnsl <<'\n';
	  delete cptr;
	  return(false);
	}
	/* Restore previous solution */
	for (i = 0; i < nq; i++)
	  pa[i] = pao[i]; 
      }
      t += dt;      
      /* Fix up solution */
      usrqp(pa, cptr->qa);
      usrfxq(cptr);
      usrpq(cptr->qa, pa);
      if (mon){
	usrqp(pa, pcptr->qa);
	usrmon(monmod);
      }      
      /* Check for sonic point */
      ss = usrss(pa);
      u = pa[iqu0+1];
      phi = 1.0e9*ydot[nal]*(ss*ss - u*u)/u;
      mach = u/ss;
      /* Check for error when the mach number goes to zero! BM
       TODO this likely needs improvement, but helps continue the calculations */
      if (mach < 0.1*imach) {
        cout << "Mach approaching zero, so resetting." << endl;
        // Try new D, closer to further bound
        dU_r = fabs(D_r - D);
        dU_l = fabs(D - D_l);
        if (dU_r > dU_l)
          D = 0.5*(D_r + D);		// Do not reset, just perturb D
        else
          D = 0.5*(D_l + D);		// Do not reset, just perturb D
        break;
      }
      nstep++;
      if ((mach - omach) > 0.1*(1 - mach))
        dt *= par[4];
      else
	dt *= par[6];
      //      usrreact(pa, ydot);
      //      if ((phi < 0) && (fabs(mach - imach) > 0.1)){
      if (mach - imach > 0.1){
	if (mach > par[2]){
	  /* Underdriven */
	  if (bisect){
	    D_l = D;
	    D = 0.5*(D_l + D_r);
	  }	  
	  if (verbose)
	    cout << "under \n";
	  break;
	}
	else if ((phi > 0) || (mach - omach < par[3]*mach)){
	  /* Overdriven */
	  if (bisect){
	    D_r = D;
	    D = 0.5*(D_l + D_r);
	  }
	  if (verbose)	  
	    cout << "over \n";
	  break;
	}
      }
      for (i = 0; i < nq; i++)
	pao[i] = pa[i];      
      omach = mach;
    }
    if (verbose){
      cout << std::setprecision(10);      
      cout << "mach " << mach << " omach " << omach <<
	" x " << pa[iqal0+nal] << " dt " << dt << '\n';
      cout <<  "phi " << phi <<" nstep " << nstep
	   << " cnsl " << cnsl <<'\n';
      cout << "D_l " << D_l << " D_r " << D_r << " D " << D << " lson "
	   << pa[iqal0+nal] << '\n';
    }
    if (nstep >= int(par[1])){
      delete cptr;
      if (verbose)
	cout << "modsol no solution in " << int(par[1]) << " steps \n";	

      return(false);
    }
    else if (!bisect)
      break; 
  }
  /* Output info  */
  cout << "curv D lson \n";  
  cout << setprecision(6) << curv << " " << setprecision(10) << D <<
    " " << pa[iqal0+nal] << '\n';
  cout << std::setprecision(6);
  delete cptr;
  return(true); 
}
//end modsol.c
 
//file modgsl.c 
/*
+ modgsl - stiff integrator
   Description:
    Integrates the solution  over time dt starting  with pa[].  At
    the end, pa[] is updated. Uses gsl. Returns false on failure.
*/
bool modgsl(int nsubl, float pa[],float dt)
{
  int i;
  float tl= 0, epsabs = 1.0e-8, epsrel = 1.0e-8, y[MAXSPECIES+1],
    ydot[MAXSPECIES+1], srate, ss;    
  struct gsl_param* param_ptr;
  gsl_odeiv2_driver* driver;
  param_ptr = new(gsl_param);
  /* Number of eqns */
  param_ptr->nvar = niso + 3;
  reactcall = 0;
  jaccall = 0;
  /* Abundances */
  for (i = 0; i < niso; i++)
    y[i] = pa[iqal0 + 1 + i];
  /* Complete solution */
  for (i = 0; i < MAXQ; i++)
    param_ptr->pa[i] = pa[i];
  /* Temperature - in units of 10^9 */
  y[niso] = 1e-9*pa[iqte];
  /* Density - in units of 10^9 */
  y[niso+1] = 1e-9*pa[iqd];
  /* Velocity - in units of 10^9 */
  y[niso+2] = 1e-9*pa[iqu0+1];
  /* Define ode system  */
  size_t dim;
  dim = param_ptr->nvar; // density, temp and velocity
  gsl_odeiv2_system sys ={modgsldiff,modjac,dim,param_ptr};
  /* Allocate ode driver */
  if (strcmp(ode,"RK2") == 0)
    driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rk2,
      dt / nsubl,epsabs,epsrel);
  else if (strcmp(ode,"RK2IMP") == 0)
    /* Implicit rk2 integrator */
    driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rk2imp,
      dt / nsubl,epsabs,epsrel);
  else if (strcmp(ode,"RK4IMP") == 0)
    /* Implicit rk2 integrator */
    driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rk4imp,
      dt / nsubl,epsabs,epsrel);
  else if (strcmp(ode, "BSIMP") == 0)
    /* Implicit Bulirsch-Stoer integrator */
    driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_bsimp,
      dt / nsubl,epsabs,epsrel);
  else if (strcmp(ode,"MSADAMS") == 0)
    /* A variable-coefficient linear multistep Adams integrator */
    driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_msadams,
      dt / nsubl,epsabs,epsrel);
  else
    /* Implicit MSBDF integrator */
    driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_msbdf,
      dt / nsubl,epsabs,epsrel);
  /* turn off internal error handling - BM
      - instead gsl will return error code to MG
       this avoids seg faults inside gsl */
  gsl_set_error_handler_off();
  /* Integrate from tl to tl+dt */
  int status = gsl_odeiv2_driver_apply(driver,&tl,tl + dt,y);
  if (status != GSL_SUCCESS) {
        printf("error in ode integrator, return value=%d\n",status);
    return(false);
  }
  gsl_odeiv2_driver_free(driver);
  /* Delete param_ptr */
  delete param_ptr;
  /* Update pa */
  for (i = 0; i < niso; i++)
    pa[iqal0 + 1 + i] = y[i];
  /* Temperature - rescale */
  pa[iqte] = 1e9*y[niso];
  /* Density - rescale */
  pa[iqd] = 1e9*y[niso+1];
  /* Velocity - rescale */
  pa[iqu0+1] = 1e9*y[niso+2];
  /* Distance from shock */
  pa[iqal0+nal] += pa[iqu0+1]*dt;
  return(true);
}
//end modgsl.c

//file modgsldiff.c
/*
+ modgsldiff - gsl wrapper for derivatives
   Description:
    Right  hand   side  function  for  gsl   ode  integrator.  Returns
    derivatives  in  ydot  using  y  as input  and  param_ptr  to  get
    parameter pa.
*/
int modgsldiff(float t,const float y[],float ydot[],
  void* param)
{
  (void)(t); /* avoid unused parameter warning */
  int i;
  float pa[MAXQ];
  float tv,c;
  struct gsl_param* param_ptr;
  param_ptr = (struct gsl_param*)param;
  // BM added from Malcolm's code 
  for (i = 0; i < niso + 2; ++i) {
    bool notOK = std::isnan(y[i]) || !std::isfinite(y[i]);
    //std::cout << y[i] << std::endl;
    if (notOK) {
      std::cout << "Nan in gsldiff, return failure" << std::endl;
      return -4;
    }
  }
  /* Zero ydot for safety */
  for (i = 0; i < param_ptr->nvar; i++)
    ydot[i] = 0.0;
  /* Fill up to iqal0 from param_ptr->pa */
  for (i = 0; i <= iqal0; i++)
    pa[i] = param_ptr->pa[i];
  /* Use y for abundances and temperature */
  for (i = 0; i < niso; i++)
    pa[iqal0 + 1 + i] = y[i];
  /* Temperature - scale up */
  pa[iqte] = 1e9*y[niso];
  /* Density - scale up */
  pa[iqd] = 1e9*y[niso+1];
  /* Velocity = scale up */
  pa[iqu0+1] = 1e9*y[niso+2];
  /* Temperature range check to prevent failures in usrhelm*/
  double temp = pa[iqte];
  if (temp < 0.0 || (log10(temp) <= tlo) || (log10(temp) >= thi)) {
    //    std::cout << "Temp range fail, log(T): " << log10(temp) << std::endl;
    return GSL_EBADFUNC;
  }
  // Added from Malcolm's code BM
  if (std::isnan(temp)) {
    //std::cout << "Sonic in gsldiff, return failure" << std::endl;
    return -2;
  }
  /* Apply fix for intermediate steps */
  modfix(pa);
  /* Source terms */
  modsrc(pa, ydot);
  // BM added checks from Malcolm's code
  /* MSBDF solver can get stuck at singular points
    If we detect that we have gone sonic, return an error code 
    and automatically we will reduce the time-step 
    */
    /* Check for NaNs */
  for (i = 0; i < niso + 2; ++i) {
    bool notOK = std::isnan(ydot[i]) || !std::isfinite(ydot[i]);
    //std::cout << ydot[i] << std::endl;
    if (notOK) {
      std::cout << "Nan in gsldiff, return failure" << std::endl;
      return -5;
    }
  }
  // check for supersonic flow
  tv = pa[iqu0+1];
  c = usrss(pa);
  tv = tv*tv - c*c;
  if (tv > 0.0) {
    //std::cout << "Sonic in gsldiff, return failure" << std::endl;
    return -3;
  }

  return GSL_SUCCESS;
}
//end modgsldiff.c

//file modjac.c
/*
+ modjac - gsl wrapper for jacobian
      Description :
        Returns  the  Jacobian  in  dfdy[]  given  the  parameters  in
        y[]. The reaction  rates, RR are in the  param structure.  The
        Jacobian   is   evaluated    by   finite   differences   using
        usrgslreact.x
*/
int modjac(float t,const float y[],float* dfdy,float dfdt[],
  void* param)
{
  bool negy;
  int nvar,i,j;
  //float rel = 3.162278e-8,ax = 1.0;
  float rel = 1e-10,ax = 1.0;
  float yt[MAXSPECIES+1],h,fminus[MAXSPECIES+1],fplus[MAXSPECIES+1];
  struct gsl_param* param_ptr;
  param_ptr = (struct gsl_param*)param;
  jaccall++;
  /* No of variables */
  nvar = param_ptr->nvar;
  /* Zero dfdt for safety */
  for (i = 0; i < nvar; i++)
    dfdt[i] = 0.0;
  /* Set up gsl matrix */
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,nvar,nvar);
  gsl_matrix* m_ptr = &dfdy_mat.matrix;
  /* Copy of y */
  for (i = 0; i < nvar; i++)
    yt[i] = y[i];
  /* Central difference for Jacobian */
  for (j = 0; j < nvar; j++) {
    /* Increment */
    h = rel * std::max(y[j],ax);
    /* Right difference */
    yt[j] = yt[j] + h;
    modgsldiff(t,yt,fplus,param);
    yt[j] = y[j];
    /* Left difference */
    yt[j] = yt[j] - h;
    negy = yt[j] <= 0.0;
    if (!negy)
      modgsldiff(t,yt,fminus,param);
    yt[j] = y[j];
    /* Set jth column of Jacobian */
    for (i = 0; i < nvar; i++) {
      if (negy)
        m_ptr->data[i * nvar + j] = fplus[i] / h;
      //	gsl_matrix_set (m_ptr, i, j, fplus[i]/h);
      else
        m_ptr->data[i * nvar + j] = 0.5 * (fplus[i] - fminus[i]) / h;
      //	gsl_matrix_set (m_ptr, i, j, 0.5*(fplus[i] - fminus[i])/h);
    }
  }
  return GSL_SUCCESS;;
}
//end modjac.c

//file modsrc.c
/*
+ modsrc - Source terms
      Description :
        Calculate  source  terms for steady solution.
*/
void modsrc(float pa[],float ydot[])
{
  int i;
  float rho,temp,u,Abar,pg,tv,dpdt,dedt,dpdd,phi,drhodt,dTdt,dudt,kterm,
    dedy,dpdy,he[10];
  rho = pa[iqd];
  temp = pa[iqte];
  u = pa[iqu0+1];
  Abar = usrAbar(pa);
  pg = 0.0;
  usrhelm(temp,rho,usrYE(pa),he);
  pg += he[0];
  dpdd = he[3];
  dpdt = he[4];
  dedt = he[6];
  usrion(temp,rho,Abar,he);
  pg += he[0];  
  dpdd += he[3];
  dpdt += he[4];
  dedt += he[6];
  usrrad(temp,rho,he);
  pg += he[0];  
  dpdd += he[3];
  dpdt += he[4];
  dedt += he[6];
  /* Set pressure in pa */
  pa[iqe] = pg;  
  usrreact(pa, ydot);
  phi = 0.0;
  dedy = 1.5*k_b*N_A*temp;
  dpdy = dedy*rho/1.5;
  /* Note as Q[i] are negative, have to add dqdy term */
  for (i = 0; i < niso; i++)
    phi += (dedy + Qi[i])*ydot[i];
  phi *= dpdt/dedt;
  for (i = 0; i < niso; i++)
    phi -= dpdy*ydot[i];
  /* Add term to phi for curvature */
  kterm = curv*rho*u*u*(D - u);
  phi += kterm;
  //  phi = min(phi, 0.0);
  tv = usrss(pa);
  /* Density time derivative */  
  drhodt = phi/(tv*tv - u*u);
  /* Velocity time derivative */
  dudt = (-u*drhodt - kterm/u)/rho;
  /* Temperature time derivative */ 
  dTdt = (u*u - dpdd)*drhodt;
  for (i = 0; i < niso; i++)
    dTdt -= dpdy*ydot[i];
  /* Add curvature term for temperature derivative */
  dTdt += kterm;
  dTdt /= dpdt;  
  /* Update ydot array */
  ydot[niso] = 1e-9*dTdt;
  ydot[niso+1] = 1e-9*drhodt;
  ydot[niso+2] = 1e-9*dudt;
}
//end modsrc.c

//file modfix.c
/*
+ modfix - Fix dependent variables
       Description :
        Set dependent physical variables from other solution variables
*/
void modfix(float pa[MAXQ])
{
  pa[iqd] = max(1.0e-3,pa[iqd]);
  pa[iqte] = max(1e5,pa[iqte]);
  /* Abundances */
    usrfxp(pa);
  /* Pressure from thermo state */
  float Abar,Ye;
  Abar = usrAbar(pa);
  Ye = usrYE(pa);
  pa[iqe] = usrpg(pa[iqte],pa[iqd],Abar,Ye);
}
//end modfix.c



// A universal safe fscanf wrapper that handles strings, numbers, and any layout BM
// Updated safe fscanf that returns a success status for loops BM
int safe_fscanf(FILE *stream, const char *format, ...) {
    va_list args;
    va_start(args, format);
    
    int items_read = vfscanf(stream, format, args);
    
    va_end(args);

    // 1. If it cleanly hits the End of File, return 0 to stop while loops safely
    if (items_read == EOF) {
        return 0; 
    }

    // 2. If it read nothing or failed to match the format, it's a critical corruption error
    if (items_read < 0) {
        fprintf(stderr, "CRITICAL ERROR: safe_fscanf failed to read expected data!\n");
        // Put your parallel abort here, e.g., MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1); 
    }

    // 3. Success! Return 1 to keep loops going
    return 1;
}
