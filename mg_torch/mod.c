#undef MAIN
#undef USER
#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "const.h"
#include "struct.h"
#include "funcs.h"
#include "model.h"
#include "user.h"
/* GSL headers*/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
//#define _LOGNANCHECKS;
extern void usrpq(float qa[], float pa[]);
extern void usrqp(float pa[], float qa[]);
extern bool usrqis(float xc[], float qa[]);
extern void usrfxq(cell *cptr);
extern void usrfxp(float pa[]);
extern float usrAbar(float pa[]);
extern float usrYE(float pa[]);
extern float usrpg(float temp, float rho, float Abar, float Ye);
extern float usrss(float pa[]);
extern void usrshock(float pa[], float rho_r, float T_r, float v_r,
  float &rho_l, float &T_l, float &v_l, float &p_l);
extern void usrhelm(float temp, float rho, float Ye, float he[]);
extern void usrreact(float pa[], float ydot[]);
extern void usrion(float temp, float rho, float Abar, float he[]);
extern void usrrad(float temp, float rho, float he[]);
extern void usrmon(char monmod[]);
extern bool usropenxq(char monmod[], int mode);
extern void usrclosexq(char monmod[], int mode);
extern bool usrrewxq(char monmod[]);
extern void sys_cmd(char cmd[]);
bool modsol(struct cell *cptr, bool verbose);
bool modgsl(int nsub, float pa[],float dt);
int modgsldiff(float t,const float y[],float ydot[],
  void* param);
int modjac(float t,const float y[],float* dfdy,float dfdt[],
  void* param);
void modsrc(float pa[],float ydot[]);
void modfix(float pa[]);
/* GSL structure  */
struct gsl_param {
  float pa[MAXQ];
  int nvar;
};
/* Offsets etc */
int iqdn,iqx,iqt, idrho, idte;
float mass_flux = 1.0; // Default for hydrostatic runs to behave
float mom_flux = 1.0; // Momentum flux
bool log_shk_rels = false;
/* Physical constants */
const float k_b = 1.380649e-16, N_A = 6.0221409e23;
/* Range and spacings for EOS */
const float   tlo = 4.0, dlo = -10.0, thi = 11.0, dhi = 11.0;
const float  tstp  = (thi - tlo)/float(JMAX-1),
  dstp  = (dhi - dlo)/float(IMAX-1);
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
    cnsl, ss, u, phi0, phi, ydot[MAXSPECIES+1], xc[MAXD];
  char scmd[SZNAM];
  struct cell *cptr;
  cptr = new cell;
  while (D_r - D_l > par[0]*D){
    /* Get upstream solution */
    for (i = 0; i < nd; i++)
      xc[i] = 0.0;
    usrqis(xc, cptr->qa);
    usrpq(cptr->qa, pr);    
    /* Set initial abundances */
    for (i = 0; i < niso; i++)
      pa[iqal0+1+i] = pr[iqal0+1+i];
    usrshock(pr, pr[iqd], pr[iqte], D, pa[iqd], pa[iqte], pa[iqu0+1],
      pa[iqe]);
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
	dt *= 0.5;
	if (verbose)
	  cout << "modgsl failure reduced dt " << dt << '\n';	
	cnt++;
	if (cnt > 100){
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
  /* Integrate from tl to tl+dt */
  int status = gsl_odeiv2_driver_apply(driver,&tl,tl + dt,y);
  if (status != GSL_SUCCESS) {
    //    printf("error in ode integrator, return value=%d\n",status);
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
  /* Apply fix for intermediate steps */
  modfix(pa);  
  /* Source terms */
  modsrc(pa, ydot);
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

