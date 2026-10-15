#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_roots.h>
#include <time.h>

//#define double float


#define IMAX 211
#define JMAX 71
/* Helmholtz free energy and derivatives */ 
double f[JMAX+1][IMAX+1], fd[JMAX+1][IMAX+1], ft[JMAX+1][IMAX+1],
  fdd[JMAX+1][IMAX+1], ftt[JMAX+1][IMAX+1], fdt[JMAX+1][IMAX+1],
  fddt[JMAX+1][IMAX+1], fdtt[JMAX+1][IMAX+1], fddtt[JMAX+1][IMAX+1];

/* Range and spacings */
const double   tlo = 4.0, dlo = -10.0, thi = 11.0, dhi = 11.0;
const double  tstp  = (thi - tlo)/double(JMAX-1),
  dstp  = (dhi - dlo)/double(IMAX-1);
/* Density and temperature points */
double d[IMAX+1], t[JMAX+1];
/* Physical constants */
const double k_b = 1.380649e-16, N_A = 6.0221409e23,
  stefan = 7.565914e-15;

inline double min(double a, double b)
{
 return (a < b ? a : b);
}

inline double max(double a, double b)
{
 return (a > b ? a : b);
}

int main();
void helm(double temp, double den, double ye, double he[]);
void ion(double temp, double den, double abar, double he[]);
void rad(double temp, double den, double he[]);
void init_helm();
double psi0(double z);
double dpsi0(double z);
double ddpsi0(double z);
double psi1(double z);
double dpsi1(double z);
double ddpsi1(double z);
double psi2(double z);
double dpsi2(double z);
double ddpsi2(double z);
double herm5(int i, int j, double w0t, double w1t, double w2t, 
 double w0mt, double w1mt, double w2mt, double w0d, double w1d,
 double w2d, double w0md, double w1md, double w2md);
void nadyo(double temp, double den, double ye, double &pres, double &ein,
 double &sound);

void bisectHelmP(double pres, double den, double ye, double &temp, double &ein, double &sound);
void bisectHelmE(double ein, double den, double ye, double &temp, double &pres, double &sound);
void bisectNadyoP(double pres, double den, double ye, double &temp, double &ein, double &sound);
void bisectNadyoE(double ein, double den, double ye, double &temp, double &pres, double &sound);

double iterateHelmP(double logtemp, void * param);

struct EOSparams {
  double pres, den, ye, ein, sound;
};



int main()
{
  double temp1, temp, den, ye, abar, he[10], pres, ein, einh, sound, presr, Tl, Th, seconds;
  int status;
  int iter = 0, max_iter = 100;

  clock_t t1,t2;		// structured variables for time
  struct stat buffer;
  long n1=0, n2=0;



  /* Initialize */
  init_helm();
  ye = 1.0;
  temp = 1.0e4;
  pres = 10e13;
  abar = 14.0;
  //while(temp > 0.0){
  temp = 1e4;
  den = 5962050.56;
  ye = 0.5;
  while(temp > 0.0){
    std::cout << "Input temperature and density and ye" << '\n';
    std::cin >> temp;
    std::cin >> den;
    std::cin >> ye;

    std::cout << std::scientific;
    std::cout << std::setprecision(6);
    if (temp > 0.0){
      temp1=temp;
      t1 = clock();		// get initial clock value
      helm(temp, den, ye, he);
      t2 = clock();		// get next clock value
      presr = he[0];
      einh = he[1];
      std::cout << '\n' << "Timmes Method (Table lookup)" << '\n';
      std::cout << "Pres from T and den:" << '\n';
      std::cout << "  temp = " << temp << " pres = " << presr << '\n';
      std::cout << "  ein = " << einh << '\n';
      std::cout << "  sound = " << he[9] << '\n';
      seconds = (double(t2 - t1)/double(CLOCKS_PER_SEC));	// convert to seconds
      std::cout << "    that took " <<  seconds << " seconds" << '\n';

      // Other pressure information, disabled for now
      ion(temp, den, abar, he);
      rad(temp, den, he);

      //n1 = times(&t1);		// get initial clock value
      t1 = clock();		// get initial clock value

      struct EOSparams *p;
      p = new EOSparams;
      p->pres = presr;
      p->den = den;
      p->ye = ye;

      gsl_root_fsolver* solver;
      //gsl_function_fdf fwrapper;
      //solver=gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);

      gsl_function fwrapper;
      solver=gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
      //solver=gsl_root_fsolver_alloc(gsl_root_fsolver_falsepos);
      //solver=gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
      fwrapper.function = iterateHelmP;
      fwrapper.params = p;
      //Tl = 1e7;
      //Th = 1e9;
      Tl = log(1e4);
      Th = log(1e11);

      gsl_root_fsolver_set(solver, &fwrapper, Tl, Th);
      //gsl_root_fdfsolver_set (s, &FDF, x);
      iter = 0;
      do {
        iter++;
        status = gsl_root_fsolver_iterate(solver);
        temp = gsl_root_fsolver_root(solver);
	Tl = gsl_root_fsolver_x_lower(solver);
	Th = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(Tl, Th, 1E-12, 1E-5);
      } while (iter < max_iter && status < 0);
      if (iter == max_iter) std::cout << "WARNING:  MAX ITER REACHED" << '\n';

      temp = exp(temp);
      ein = p->ein;
      sound = p->sound;
      gsl_root_fsolver_free(solver);
      delete p;

      //n2 = times(&t2);		// get next clock value
      t2 = clock();		// get next clock value

      std::cout << "Temp from P and den (GSL):" << '\n';
      std::cout << "  temp = " << temp << '\n';
      std::cout << "  ein = " << ein << '\n';
      std::cout << "  sound = " << sound << '\n';
      seconds = (double(t2 - t1)/double(CLOCKS_PER_SEC));	// convert to seconds
      std::cout << "    that took " <<  seconds << " seconds" << '\n';



      t1 = clock();		// get initial clock value
      bisectHelmP(presr, den, ye, temp, ein, sound);
      t2 = clock();		// get next clock value
      std::cout << "Temp from P and den:" << '\n';
      std::cout << "  temp = " << temp << '\n';
      std::cout << "  ein = " << ein << '\n';
      std::cout << "  sound = " << sound << '\n';
      seconds = (double(t2 - t1)/double(CLOCKS_PER_SEC));	// convert to seconds
      std::cout << "    that took " <<  seconds << " seconds" << '\n';



      t1 = clock();		// get initial clock value
      bisectHelmE(einh, den, ye, temp, pres, sound);
      t2 = clock();		// get next clock value
      std::cout << "T and P from ein and den:" << '\n';
      std::cout << "  temp = " << temp << '\n';
      std::cout << "  pres = " << pres << '\n';
      std::cout << "  sound = " << sound << '\n';
      seconds = (double(t2 - t1)/double(CLOCKS_PER_SEC));	// convert to seconds
      std::cout << "    that took " <<  seconds << " seconds" << '\n';



      std::cout << '\n' << "Nadyozhin approx:" << '\n';
      t1 = clock();		// get initial clock value
      nadyo(temp1, den, ye, pres, einh, sound);
      t2 = clock();		// get next clock value
      std::cout << "Pres from T and den:" << '\n';
      std::cout << "  pres = " << pres << '\n';
      std::cout << "  ein = " << einh << '\n';
      std::cout << "  sound = " << sound << '\n';
      seconds = (double(t2 - t1)/double(CLOCKS_PER_SEC));	// convert to seconds
      std::cout << "    that took " <<  seconds << " seconds" << '\n';


      t1 = clock();		// get initial clock value
      bisectNadyoP(pres, den, ye, temp, ein, sound);
      t2 = clock();		// get next clock value
      std::cout << "Temp from P and den:" << '\n';
      std::cout << "  temp = " << temp << '\n';
      std::cout << "  ein = " << ein << '\n';
      std::cout << "  sound = " << sound << '\n';
      seconds = (double(t2 - t1)/double(CLOCKS_PER_SEC));	// convert to seconds
      std::cout << "    that took " <<  seconds << " seconds" << '\n';

      t1 = clock();		// get initial clock value
      bisectNadyoE(einh, den, ye, temp, pres, sound);
      t2 = clock();		// get next clock value
      std::cout << "T and P from ein and den:" << '\n';
      std::cout << "  temp = " << temp << '\n';
      std::cout << "  pres = " << pres << '\n';
      std::cout << "  sound = " << sound << '\n';
      seconds = (double(t2 - t1)/double(CLOCKS_PER_SEC));	// convert to seconds
      std::cout << "    that took " <<  seconds << " seconds" << '\n';
      std::cout << '\n';

      //temp = temp*1.1;
    }
  }

}

/*
This routine performs a  thermodynamically consistent interpolation in
a tabular electron-positron equation  of state using biquintic hermite
basis functions.  Based on F.X.  Timmes, D.F. Swesty, APjS,  126, 501,
2000. Data is stored in helm_table.dat

input :
temp = temperature(in K)
den = density (in g cm~3)
ye = electrons per baryon = zbar/abar

helm_data.h contains the  table of helmholtz free energy  and eight of
its partial derivatives : f,  df_d, df_t, df_dd, df_tt, df_dt, df_ddt,
df_dtt and df_ddtt

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

  Range:  1e8 < T <1e10
          1e-3 < d < 1e11


  */
void helm(double temp, double den, double ye, double he[])
{
  int i, j, iat, jat;
  double free, df_d, df_t, df_dd, df_tt, df_dt, dt, dt2, dt2i, dd,
  dd2, ddi, dd2i, xt, xd, mxt, mxd, si0t, si1t, si2t, si0mt, si1mt,
  si2mt, si0d, si1d, si2d, si0md, si1md, si2md, dsi0t, dsi1t, dsi2t,
  dsi0mt, dsi1mt, dsi2mt, dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md, 
  ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, ddsi0d, ddsi1d,
  ddsi2d, ddsi0md, ddsi1md, ddsi2md, z, w0t, w1t, w2t, w0mt, w1mt, w2mt,
  w0d, w1d, w2d, w0md, w1md, w2md, din, ptot, dpdt, dpdd, stot, dsdt,
    dsdd, etot, dedt, dedd, dedd_p, dedp_d;

  /* f and its derivatives 
  double f[JMAX][IMAX], fd[JMAX][IMAX], ft[JMAX][IMAX], fdd[JMAX][IMAX],
  ftt[JMAX][IMAX], fdt[JMAX][IMAX], fddt[JMAX][IMAX], 
  fdtt[JMAX][IMAX],fddtt[JMAX][IMAX]; */

  /* Find the table locations, entering with ye * density */

  jat = int((log10(temp) - tlo)/tstp) + 1;
  din = ye*den;
  iat = int((log10(din) - dlo)/dstp) + 1;
  /* Various differences */
  dt = t[jat+1] - t[jat];
  dt2 = dt*dt;
  dd = d[iat+1] - d[iat];
  dd2 = dd*dd;
  xt = max((temp - t[jat])/dt, 0.0);
  xd = max((din - d[iat])/dd,0.0);
  mxt = 1.0 - xt;
  mxd = 1.0 - xd;

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
  dpdd = ye*(din*din*df_dd + 2.0*din*df_d);
  stot = -df_t*ye;
  dsdt = -df_tt*ye;
  dsdd = -df_dt*ye*ye;

  etot = ye*free + temp*stot;
  dedt = temp*dsdt;
  dedd = ye*ye*df_d + temp*dsdd;
  he[0] = ptot;
  he[1] = etot;
  he[2] = stot;
  he[3] = dpdd;
  he[4] = dpdt;
  he[5] = dedd;
  he[6] = dedt;
  he[7] = dsdd;
  he[8] = dsdt;

  /* (de/dd)_p */
  dedd_p = dedd - dedt*dpdd/dpdt;
  /* (de/dp)_d */
  dedp_d = dedt/dpdt;
  he[9] = sqrt((ptot/(den*den) - dedd_p)/dedp_d);

  /* Flash sound speed */
  double pres, presi, chit, chid, x7, af;
  pres = ptot;
  presi = 1.0e0/pres;
  chit  = temp*presi*dpdt;
  chid  = dpdd*den*presi;
  x7 = pres*chit/(temp*dedd*den);
  af = sqrt((chit*x7 + chid)*pres/den);
  //std::cout << std::scientific;
  //std::cout << std::setprecision(6);
  //std::cout << "Helm alt sound " << sqrt(dpdd - dpdt*dsdd/dsdt) << '\n';
  //std::cout << "Flash sound " << af << '\n';
  //std::cout << "Helmholtz: pressure " << he[0] << " energy " << he[1] << 
  //  " sound " << he[9] <<'\n';
  return;
}
/*
Calculates ion thermodynamics

input :
temp = temperature(in K)
den = density (in g cm~3)
abar = average atomic mass of nucleons

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
void ion(double temp, double den, double abar, double he[])
{
  double tv1, tv2, pres, etot, dpdt, dpdd, dedt, dedd, stot, dsdt, dsdd;
  tv1 = k_b*N_A/abar;
  tv2 = tv1*den;
  pres = tv2*temp;
  dpdt = tv2;
  dpdd = tv1*temp;

  etot = 1.5*tv1*temp;
  dedt = 1.5*tv1;   
  dedd = 0.0;

  stot = 0.0; /* Warning do not use */
  dsdt = 1.5*tv1/temp;
  dsdd = - tv1/den;
  he[0] = pres;
  he[1] = etot;
  he[2] = stot;
  he[3] = dpdd;
  he[4] = dpdt;
  he[5] = dedd;
  he[6] = dedt;
  he[7] = dsdd;
  he[8] = dsdt;
  std::cout << std::scientific;
  std::cout << std::setprecision(6);
  std::cout << "Ion      : pressure " << he[0] << " energy " << he[1] << 
    " sound " << sqrt(dpdd - dpdt*dsdd/dsdt) <<'\n';
  return;
} 
/*
Calculates radiation thermodynamics

input :
temp = temperature(in K)
den = density (in g cm~3)

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
void rad(double temp, double den, double he[])
{
  double tv1, tv2, pres, etot, dpdt, dpdd, dedt, dedd, stot, dsdt, dsdd;
  tv1 = temp*temp;
  tv1 = tv1*tv1;
  tv2 = stefan*tv1;
  pres = tv2/3.0;
  dpdt = 4.0*pres/temp;
  dpdd = 0.0;

  etot = tv2/den;
  dedt = 4.0*etot/temp;   
  dedd = -etot/den;

  stot = 4.0*etot/(3.0*temp);
  dsdt = 3.0*stot/temp;
  dsdd = - stot/den;
  he[0] = pres;
  he[1] = etot;
  he[2] = stot;
  he[3] = dpdd;
  he[4] = dpdt;
  he[5] = dedd;
  he[6] = dedt;
  he[7] = dsdd;
  he[8] = dsdt;
  std::cout << std::scientific;
  std::cout << std::setprecision(6);
  std::cout << "Rad      : pressure " << he[0] << " energy " << he[1] << 
    " sound " << sqrt(dpdd - dpdt*dsdd/dsdt) <<'\n';
  return;
} 
void init_helm()
{
  FILE *input;
  int i, j;
  double tf, tfd, tft, tfdd, tftt, tfdt, tfddt, tfdtt, tfddtt;
  input = fopen("helm_table.dat","r");
  for (j = 1; j <= JMAX; j++){
    for (i = 1; i <= IMAX; i++){
      fscanf(input,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &tf, &tfd, 
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
 fclose(input);

 /* Density and temperature arrays */
 for (i = 1; i <= IMAX; i++)
   d[i] = powf(10.0, dlo + (i-1)*dstp);
 for (j = 1; j <= JMAX; j++)
   t[j] = powf(10.0, tlo + (j-1)*tstp);

}
  /* quintic hermite basis statement functions*/
  /* psi0 and its derivatives */
double psi0(double z)
{
  return(z*z*z*(z*(-6.0*z + 15.0) - 10.0) + 1.0);
}
double dpsi0(double z)
{
  return(z*z*(z*(-30.0*z + 60.0) - 30.0));
}
double ddpsi0(double z)
{
  return(z*(z*(-120.0*z + 180.0) - 60.0));
}
  /* psi1 and its derivatives */
double psi1(double z)
{
  return(z*(z*z*(z*(-3.0*z + 8.0) - 6.0) + 1.0));
}
double dpsi1(double z)
{
  return(z*z*(z*(-15.0*z + 32.0) - 18.0) + 1.0);
}
double ddpsi1(double z)
{
  return(z*(z*(-60.0*z + 96.0) - 36.0));
}
  /* psi2 and its derivatives */
double psi2(double z)
{
  return(0.5*z*z*(z*(z*(-z + 3.0) - 3.0) + 1.0));
}
double dpsi2(double z)
{
  return(0.5*z*(z*(z*(-5.0*z + 12.0) - 9.0) + 2.0));
}
double ddpsi2(double z)
{
  return(0.5*(z*(z*(-20.0*z + 36.0) - 18.0) + 2.0));
}
  /* bicubic hermite polynomial statement function */
double herm5(int i, int j, double w0t, double w1t, double w2t, 
  double w0mt, double w1mt, double w2mt, double w0d, double w1d,
  double w2d, double w0md, double w1md, double w2md)
{
  double tv;

  tv  = f[j][i]*w0d*w0t + f[j][i+1]*w0md*w0t
    + f[j+1][i]*w0d*w0mt + f[j+1][i+1]*w0md*w0mt
    + ft[j][i]*w0d*w1t + ft[j][i+1]*w0md*w1t
    + ft[j+1][i]*w0d*w1mt + ft[j+1][i+1]*w0md*w1mt
    + ftt[j][i]*w0d*w2t + ftt[j][i+1]*w0md*w2t
    + ftt[j+1][i]*w0d*w2mt + ftt[j+1][i+1]*w0md*w2mt
    + fd[j][i]*w1d*w0t + fd[j][i+1]*w1md*w0t
    + fd[j+1][i]*w1d*w0mt + fd[j+1][i+1]*w1md*w0mt
    + fdd[j][i]*w2d*w0t + fdd[j][i+1]*w2md*w0t
    + fdd[j+1][i]*w2d*w0mt + fdd[j+1][i+1]*w2md*w0mt
    + fdt[j][i]*w1d*w1t + fdt[j][i+1]*w1md*w1t
    + fdt[j+1][i]*w1d*w1mt + fdt[j+1][i+1]* w1md*w1mt
    + fddt[j][i]*w2d*w1t + fddt[j][i+1]*w2md*w1t
    + fddt[j+1][i]*w2d*w1mt + fddt[j+1][i+1]*w2md*w1mt
    + fdtt[j][i]*w1d*w2t + fdtt[j][i+1]*w1md*w2t
    + fdtt[j+1][i]*w1d*w2mt + fdtt[j+1][i+1]*w1md*w2mt
    + fddtt[j][i]*w2d*w2t + fddtt[j][i+1]*w2md*w2t
    + fddtt[j+1][i]*w2d*w2mt + fddtt[j+1][i+1]*w2md*w2mt;
  return(tv);
}
/* Nadyozhin approx (apjs106) */
void nadyo(double temp, double den, double ye, double &pres, double &ein,
  double &sound)
{
  const double Q = 4.82744e-8, B = 2.68192e13, a1 = 45.0/16.0,
    a2 = 2865.0/512.0, a3 = 15.0/4.0, a4 = 165.0/32, a5 = 5.0/4.0,
    a6 = 9.0/16.0, a7 = 15.0/16.0, a8 = 955.0/512.0, a9 = 21.0/4.0, 
    a10 = 405.0/32;
  double al, tv1, tv2, tv3;
  /* 
  al = 5.929896509e9/temp;
  tv1 = exp(-2.0*al);
  tv2 = k_b*N_A*ye*temp;
  tv3 = den*ye;
  pres = den*tv2*(1.0 + Q*tv3*powf(al, 1.5)*(1.0 - (a1
    - a2/al)/al) + B*tv1*(1.0 + (a3 + a4/al)/al)/
    (al*al*al*tv3*tv3)); 
  ein = 1.5*tv2*(1.0 + a5*(1.0 - (1.0 - a6/al)/al)/al
    + Q*tv3*powf(al, 1.5)*(1.0 - (a7 + a8/al)/al) + 
    2.0*B/(3.0*powf(tv3/al, 2.0))*tv1*
    (1.0 +(a9 + a10/al)/al));
  */
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15,
    t17, t18, t19, t20, t21, t22, t23, t24, t25, t27, t29, t31, t33,
    t38, t39, t40, t45, t49, t63;
  al = 5.929896509e9/temp;
  t4 = pow(al,1.5);
  t5 = 1/al;
  t7 = al*al;
  t8 = 1/t7;
  t13 = den*den;
  t19 = exp(-2.0*al);
  /* Pressure */
  pres = k_b*temp*den*N_A*(1.0+Q*den*t4*(1.0-a1*t5+a2*t8)+
    B/t13/t7/al*t19*(1.0+a3*t5+a4*t8));


  /* Internal Energy */
  t2 = al;
  t3 = 1/t2;
  t11 = pow(t2,1.5);
  t18 = den*den;
  t21 = t2*t2;
  t24 = exp(-2.0*t2);
  ein = 1.5*k_b*temp*N_A*(1.0+a5*(0.1E1-(1.0-a6*t3)*t3)*t3+
    Q*den*t11*(1.0-(a7+a8*t3)*t3)+
    2.0*B/3.0/t18/t21*t24*(1.0+(a9+a10*t3)*t3));
  double  dpdd, dpdt, dedd, dedt, dedd_p, dedp_d;
  /* (dp/dd)_T */
  t1 = k_b*temp;
  t3 = pow(al,1.5);
  t4 = 1/al;
  t6 = al*al;
  t7 = 1/t6;
  t9 = 1.0-a1*t4+a2*t7;
  t12 = den*den;
  t18 = exp(-2.0*al);
  t23 = 1/t6/al*t18*(1.0+a3*t4+a4*t7);
  dpdd = t1*N_A*(1.0+Q*den*t3*t9+B/t12*t23)+t1*den*N_A*(Q*t3*t9
    -2.0*B/t12/den*t23);
  /* (dp/dT)_d */
  t2 = Q*den;
  t3 = al;
  t4 = pow(t3,1.5);
  t5 = 1/t3;
  t7 = t3*t3;
  t8 = 1/t7;
  t10 = 1.0-a1*t5+a2*t8;
  t13 = den*den;
  t15 = B/t13;
  t17 = 1/t7/t3;
  t19 = exp(-2.0*t3);
  t20 = t17*t19;
  t23 = 1.0+a3*t5+a4*t8;
  t31 = pow(t3,0.5);
  t33 = -al/temp;
  t45 = t7*t7;
  t49 = t19*t23*t33;
  dpdt = k_b*den*N_A*(1.0+t2*t4*t10+t15*t20*t23)+
    k_b*temp*den*N_A*(0.15E1*t2*t31*t10*t33+t2*t4*(a1*t8*t33-
    2.0*a2*t17*t33)-3.0*t15/t45*t49-2.0*t15*t17*t49+t15*t20*(-a3*t8*t33
    -2.0*a4*t17*t33));
  /* (de/dd)_T */
  t2 = al;
  t3 = pow(t2,1.5);
  t5 = 1/t2;
  t11 = den*den;
  t15 = t2*t2;
  t18 = exp(-2.0*t2);
  dedd = 1.5*k_b*temp*N_A*(Q*t3*(1.0-(a7+a8*t5)*t5)-4.0/3.0*B/t11/
    den/t15*t18*(1.0+(a9+a10*t5)*t5));
  /* (de/dT)_d */
  t2 = al;
  t3 = 1/t2;
  t5 = 1.0-a6*t3;
  t8 = a5*(1.0-t5*t3);
  t10 = Q*den;
  t11 = pow(t2,1.5);
  t13 = a7+a8*t3;
  t15 = 1.0-t13*t3;
  t18 = den*den;
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
  sound = sqrt((pres/(den*den) - dedd_p)/dedp_d);
      //std::cout << "Nadyozhin: temp " << temp << " den " << den << " ye " << ye << " pressure " << pres << " energy " << ein <<
      //  " sound " << sound <<'\n';


  return;
}



double iterateHelmP(double logtemp, void * param) {
  double f, he[10], tol;
  struct EOSparams* p = (EOSparams*)param;
  helm(exp(logtemp), p->den, p->ye, he);
  p->ein = he[1];
  p->sound = he[9];
  f = (he[0] - p->pres)/p->pres;
  return f;
}




void bisectHelmP(double pres, double den, double ye, double &temp, double &ein, double &sound) {
  double a, b, x, y, f, he[10], tol;
  tol = 1E-15;

  // upper and lower limits for T
  a = log(1e4);
  b = log(1e11);

  do {
    x = 0.5*(a+b);
    helm(double(exp(a)), den, ye, he);
    f = (he[0] - pres)/pres;
    temp = exp(x);
    helm(temp, den, ye, he);
    y = (he[0] - pres)/pres;
    if (f*y < 0.0)
      b = x;
    else
      a = x;
  } while (fabs(y) > tol);

  ein = he[1];
  sound = he[9];
  return;
}



void bisectHelmE(double ein, double den, double ye, double &temp, double &pres, double &sound) {
  double a, b, x, y, f, he[10], tol;
  tol = 1E-12;

  // upper and lower limits for T
  a = log(1e4);
  b = log(1e11);

  do {
    x = 0.5*(a+b);
    helm(double(exp(a)), den, ye, he);
    f = (he[1] - ein)/ein;
    temp = exp(x);
    helm(temp, den, ye, he);
    y = (he[1] - ein)/ein;
    if (f*y < 0.0)
      b = x;
    else
      a = x;
  } while (fabs(y) > tol);

  pres = he[0];
  sound = he[9];
  return;
}


void bisectNadyoP(double pres, double den, double ye, double &temp, double &ein, double &sound) {
  double a, b, x, y, f, pres2, tol;
  tol = 1E-12;

  // upper and lower limits for T
  a = log(1e4);
  b = log(1e11);

  do {
    x = 0.5*(a+b);
    //helm(double(exp(a)), den, ye, he);
    nadyo(double(exp(a)), den, ye, pres2, ein, sound);
    f = (pres2 - pres)/pres;
    temp = exp(x);
    //helm(temp, den, ye, he);
    nadyo(temp, den, ye, pres2, ein, sound);
    y = (pres2 - pres)/pres;
    if (f*y < 0.0)
      b = x;
    else
      a = x;
  } while (fabs(y) > tol);

  return;
}






void bisectNadyoE(double ein, double den, double ye, double &temp, double &pres, double &sound) {
  double a, b, x, y, f, einn, tol;
  tol = 1E-12;

  // upper and lower limits for T
  a = log(1e5);
  b = log(1e11);

  do {
    x = 0.5*(a+b);
    //helm(double(exp(a)), den, ye, he);
    nadyo(double(exp(a)), den, ye, pres, einn, sound);
    f = (einn - ein)/ein;
    temp = exp(x);
    //helm(temp, den, ye, he);
    nadyo(temp, den, ye, pres, einn, sound);
    y = (einn - ein)/ein;
    if (f*y < 0.0)
      b = x;
    else
      a = x;
  } while (fabs(y) > tol);
  return;
}



