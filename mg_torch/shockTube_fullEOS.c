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
#include <fstream>

/* Use GSL for root finding */
#include <gsl/gsl_roots.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
 

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

void GSLbisectHelmP7(double den, double ye, double &temp, double he[]);
double GSLiterateHelmP(double logtemp, void * param);


void GSLbisectHelmS7(double &den, double ye, double &temp, double he[]);
double GSLiterateHelmS(double logden, void * param);


void GSLbisectHelmP7d(double &den, double ye, double temp, double he[]);
double GSLiterateHelmPd(double logden, void * param);

int GSLiterateHelmSTD(const gsl_vector *x, void *param, gsl_vector *ret);
void GSLiterHelmS7(double &den, double ye, double &temp, double he[]);

int dGSLiterateHelmSTD(const gsl_vector *x, void *param, gsl_matrix *ret);
int fdGSLiterateHelmSTD(const gsl_vector *x, void *param, gsl_vector *ret_vec, gsl_matrix *ret_mat);

void eosTotal(double temp, double den, double ye, double he[], double abar);
/* shock tube problem functions */

void getRiemannState(double dxdt, double u_star,
  double cl_star, double cr_star, double pl_star, double pr_star, double dl_star, double dr_star, double Vl, double Vr,
  double dl, double ul, double pl, double cl, double sl, double dr, double ur, double pr, double cr, double sr,
  double &d, double &u, double &p);



/* Data structure */
struct EOSparams {
  double ye, den, temp;
  double he[10];
};



int main()
{
  double temp1, temp, den, ye, abar, he[10], hel[10], her[10], pres, ein, einh, sound, presr, Th, seconds;
  int status;
  int iter = 0, max_iter = 100;

  clock_t tm1,tm2;		// structured variables for time
  double cpu_time;
  struct stat buffer;
  long n1=0, n2=0;

  /* Initialize */
  init_helm();
  ye = 0.5;


  // start the clock
  tm1 = clock();




  /* This program will now solve the classic shock tube problem
     to get analytical soln for the hydrodynamics
  Note EOS gives:
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
 */

  double ul, Tl, dl, pl, cl, sl, hl, ur, Tr, dr, pr, cr, sr, hr, cvl, cvr, chidl, chidr, chiTl, chiTr, Rl, Rr, cpl, cpr, sigma, z, sel, ser, ustar, ushockl, ushockr;
  double Jplus, Jminus, plstar, dlstar, clstar, prstar, drstarm, crstar, F, dpldu, dprdu, w1, w2, Mshock, hlstar, h2lstar, hrstar, h2rstar, dhdu, dlstarp, plstarp, drstarp, prstarp, w1p, w2p, condition, ushockp;
  double tol, ucheck, dh, hlstarp, hrstarp, h2lstarp, h2rstarp, drstar, dp, tlstar, trstar, trstarp, tlstarp;

  /* ICs */
  ul = 0.0;
  Tl = 1.0e9;
  dl = 1.0e7;
  ur = 0.0;
  Tr = 1.0e8;
  dr = 1.0e6;

  /* Now get initial pressures, sound speeds, enthalpy, and entropu */
  //helm(Tl, dl, ye, hel);	// Equation of state TODO add other components
  eosTotal(Tl, dl, ye, hel, 13.7);
  pl = hel[0];
  cl = hel[9];
  sl = hel[2];
  hl = hel[1]+(pl/dl);


  //helm(Tr, dr, ye, her);	// Equation of state TODO add other components
  eosTotal(Tr, dr, ye, her, 13.7);
  pr = her[0];
  cr = her[9];
  sr = her[2];
  hr = her[1]+(pr/dr);

  /* Provide guess for ustar, based on Gottlieb and Groth but adapted for general EOS */
  cvl = hel[6];
  cvr = her[6];
  chidl = hel[3]*dl/pl;
  chidr = her[3]*dr/pr;
  chiTl = hel[4]*Tl/pl;
  chiTr = her[4]*Tr/pr;
  Rl = chiTl*chiTl*pl/(dl*Tl*chidl);
  Rr = chiTr*chiTr*pr/(dr*Tr*chidr);
  cpl = Rl+cvl;
  cpr = Rr+cvr;
  if (pl >= pr)
    sigma = cpl/cvl;
  else
    sigma = cpr/cvr;
  z = (cr*cvr*Rl/cl*cvl*Rr)*pow(pl/pr,sigma);
  //sel = ul + 2*cvl*cl/Rl;		// this is the Gottlieb and Groth estimate
  //ser = ur - 2*cvr*cr/Rr;
  sel = ul + cl;			// for SN1A, u+c and u-c give better estimate of ustar
  ser = ur - cr;
  ustar = (z*sel+ser)/(1+z);		// for now estimate as avaerage of wave speeds



  ushockl = 0.0;
  ushockr = 0.0;
  bool done=false;
  int count = 0;
  tol = 2E-3;
  do {
    std::cout << "ustar (guess) = " << ustar << '\n';
    if (count > 100000) {
      std::cout << "error in riemann convergence. try increasing atol.\n";
    }
    /* left moving rarefaction */
    if (ustar > ul) {
	std::cout << " left moving rarefaction... " << '\n';
      // let Fl=0
      Jplus = ul;	// take as constant
      dp = -1.0E-5*dl*cl*(ustar-ul);	// Obtain dp based on difference in velocity
      plstar = pl;
      dlstar = dl;
      clstar = cl;
      F = 0.0;
      do {
	F += dp/(dlstar*clstar);
	ucheck = Jplus-F;
	// update variables
	plstar += dp;
	/* TODO to generalize, get T from S, and use T and P to get rho. get c from T.*/
	hel[0] = plstar;
	hel[2] = sl;
	GSLiterHelmS7(dlstar, ye, tlstar, hel);
        clstar = hel[9];
      } while (ucheck<ustar);
      dpldu = -dlstar*clstar;	// derivative of p w.r.t. u
    }
    /* left moving shock */
    else {
 	std::cout << " left moving shock... " << '\n';
      // initial guess for ushock (-ve)
      //ushockl = 0.75*ustar + 0.25*ul;	// this is Browne estimate based on rho2=5rho1, condition ushock < ustar
      ushockl = 0.5*ustar;		// For some reason this is a better buess in this case
      bool done2=false;
      int count2 = 0;
      do {
        if (count2 > 100000) {
          std::cout << "error in shock convergence. try increasing atol.\n";
        }
        // transform to shock frame
        w1 = ul-ushockl;
        w2 = ustar - ushockl;
        Mshock = w1/cl;
        dlstar = dl*w1/w2;
        plstar = pl+(dl*w1*w1)-(dlstar*w2*w2);
	/* TODO to generalize just solve P/r=RT */
	hel[0] = plstar;
        GSLbisectHelmP7(dlstar, ye, tlstar, hel);	// gets T from P and rho
        clstar = hel[9];
        hlstar = hel[1]+(plstar/dlstar);
        h2lstar = hl+(0.5*w1*w1) - (0.5*w2*w2);
        dh = hlstar - h2lstar;	
        // perturb ushock
        ushockp = 1.05*ushockl;	// perturb the ushock value  by 5% to get dh/du
        w1p = ul-ushockp;
        w2p = ustar - ushockp;
        dlstarp = dl*w1p/w2p;
        plstarp = pl+(dl*w1p*w1p)-(dlstarp*w2p*w2p);
	hel[0] = plstarp;
        GSLbisectHelmP7(dlstarp, ye, tlstarp, hel);	// TODO gets T from P and rho
        hlstarp = hel[1]+(plstarp/dlstarp);
        h2lstarp = hl+(0.5*w1p*w1p) - (0.5*w2p*w2p);
        dhdu = (hlstarp-hlstar - h2lstarp + h2lstar)/(ushockp-ushockl);
	// check if h2lstar = hlstar, if not new guess for ushock
        condition = fabs(1.0-(hlstar/h2lstar));
        if (condition < tol) {
          done2 = true;
        }
        /* new guess for ushock */
        else {
          ushockl = ushockl - dh/dhdu;
        }
        count2++;
      } while (!done2);
      dpldu = 2*cl*Mshock*Mshock*Mshock/(1+Mshock*Mshock);	// derivative of p w.r.t. u
    }
    /* right moving rarefaction */
    if (ustar < ur) {
	std::cout << " right moving rarefaction... " << '\n';
      // let Fr=0
      Jminus = ur;	// constant
      dp = 1.0E-5*dr*cr*(ustar-ur);		// Obtain dp based on difference in velocity
      prstar = pr;
      drstar = dr;
      crstar = cr;
      F = 0.0;
      do {
	F += dp/(drstar*crstar);
	ucheck = Jminus+F;
	// update variables
	prstar += dp;
	/* TODO to generalize, get T from S, and use T and P to get rho. get c from T.*/
	her[0] = prstar;
	her[2] = sr;
	GSLiterHelmS7(drstar, ye, trstar, her);
        crstar = her[9];
      } while (ucheck>ustar);
      dprdu = -drstar*crstar;	// derivative of p w.r.t. u
    }
    /* right moving shock */
    else {
	std::cout << " right moving shock... " << '\n';
      // initial guess for ushock (+ve)
      //ushockr = 1.25*ustar - 0.25*ur;	// Browne estimate based on rho2=5rho1 but condition ushockr > ustar...
      ushockr = 2.0*ustar;		// For some reason this is a better buess in this case
      bool done2=false;
      int count2 = 0;
      do {
        if (count2 > 100000) {
          std::cout << "error in shock convergence. try increasing atol.\n";
        }
        // transform to shock frame
        w1 = ushockr-ur;
        w2 = ushockr-ustar;
        Mshock = w1/cr;
        drstar = dr*w1/w2;
        prstar = pr+(dr*w1*w1)-(drstar*w2*w2);
	/* TODO to generalize just solve P/r=RT */
	her[0] = prstar;
        GSLbisectHelmP7(drstar, ye, trstar, her);	// gets T from P and rho
        crstar = her[9];
        hrstar = her[1]+(prstar/drstar);
        h2rstar = hr+(0.5*w1*w1) - (0.5*w2*w2);
        dh = hrstar - h2rstar;	
        // perturb ushock by 5% to get dh/du
        ushockp = 1.05*ushockr;	// perturb the ushock value 
        w1p = ushockp-ur;
        w2p = ushockp- ustar;
        drstarp = dr*w1p/w2p;
        prstarp = pr+(dr*w1p*w1p)-(drstarp*w2p*w2p);
	her[0] = prstarp;
        GSLbisectHelmP7(drstarp, ye, trstarp, her);	// TODO gets T from P and rho
        hrstarp = her[1]+(prstarp/drstarp);
        h2rstarp = hr+(0.5*w1p*w1p) - (0.5*w2p*w2p);
        dhdu = (hrstarp-hrstar - h2rstarp + h2rstar)/(ushockp-ushockr);
	// check if h2lstar = hlstar, if not new guess for ushock
        condition = fabs(1.0-(hrstar/h2rstar));
	//std::cout << "  condition = " << condition << '\n';
        if (condition < tol) {
          done2 = true;
        }
        /* new guess for ushock */
        else {
          ushockr = ushockr - dh/dhdu;
        }
        count2++;
      } while (!done2);
      dprdu = 2*cr*Mshock*Mshock*Mshock/(1+Mshock*Mshock);	// derivative of p w.r.t. u

    }

    /* check for convergence */
    condition = fabs(1.0-(plstar/prstar));
    if (condition < tol) {
       done = true;
    }
    /* new guess for ustar */
    else {
      ustar = ustar - ((plstar-prstar)/(dpldu-dprdu));		// new guess for ustar based on difference in p, and their derivatives
    }
    count++;
  } while (!done);

  tm2 = clock();

  cpu_time = ((double) (tm2 - tm1)) / CLOCKS_PER_SEC;
  std::cout << "ustar = " << ustar << '\n';
  std::cout << "  that took = " << cpu_time << "s" << '\n';



  /* now get profile, reconstruct the fan if necessary

 loop through dxdt....
  t = 0.0008 in test case
  xl = -500000  xr = 500000
  dxdt_max = +- 625000000

 */
  double dxdt = 0.0;
  double time = 0.0008;
  double bd = 500000;		// x location of +/- boundaries
  double npts = 200;		// number of points in reconstruction
  double delta = 2*bd/(npts*time);	// this is the change in dxdt from one data point to the next
  double x, d, u, p;

  /* create file */
  std::ofstream cellout;
  cellout.open("shocktube_results.dat",std::ios::out); 
  cellout << "#x rho u p\n";
  cellout.close();





  /* first get flux at x= 0 */
  tm1 = clock();
  getRiemannState(dxdt, ustar, clstar, crstar, plstar, prstar, dlstar, drstar, ushockl, ushockr,
  dl, ul, pl, cl, sl, dr, ur, pr, cr, sr, d, u, p); 	// returns d u p at dxdt
  tm2 = clock();
  cpu_time = ((double) (tm2 - tm1)) / CLOCKS_PER_SEC;
  std::cout << "at cell interface:" << '\n';

  std::cout << "  u = " << u << '\n';
  std::cout << "  d = " << d << '\n';
  std::cout << "  p = " << p << '\n';
  std::cout << "  Vr = " << ushockr << '\n';
  std::cout << "  xshock = " << ushockr*time << '\n';
  std::cout << "  Mshock = " << ushockr/cr << '\n';
  std::cout << "  that took = " << cpu_time << "s" << '\n';

  /* Now loop throuhg dxdt to reconstruct entire profile */
  dxdt = -bd/(time);
  std::cout << "reconstructing profile..." << '\n';

  tm1 = clock();
  do {
      x = dxdt*time;
      getRiemannState(dxdt, ustar, clstar, crstar, plstar, prstar, dlstar, drstar, ushockl, ushockr, dl, ul, pl, cl, sl, dr, ur, pr, cr, sr, d, u, p);

      // output to file
      std::ofstream cellout;
      cellout.open("shocktube_results.dat",std::ios::app); 
      cellout << x << " " << d << " " << u << " " << p << "\n";
      cellout.close();

      dxdt += delta;
  }
  while (dxdt < (bd/time));
  tm2 = clock();
  cpu_time = ((double) (tm2 - tm1)) / CLOCKS_PER_SEC;

  std::cout << "shocktube_results.dat created" << '\n';
  std::cout << "  that took = " << cpu_time << "s" << '\n';

  return 0;
}










/* This function returns the density velocity and pressure at a given location of dxdt
Vl Vr is the velocity of right left shocks (these can be zero for exp. fan) */
void getRiemannState(double dxdt, double u_star,
  double cl_star, double cr_star, double pl_star, double pr_star, double dl_star, double dr_star, double Vl, double Vr,
  double dl, double ul, double pl, double cl, double sl, double dr, double ur, double pr, double cr, double sr,
  double &d, double &u, double &p) {

    double cminus, cplus, cfan, pfan, dfan, tfan, ufan, Jplus, Jminus, dp, F;
    double he[10];
    double ye=0.5;

    /* if (x,t) is on the left of the contact surface... */
    if (u_star > dxdt) {
      /* leftward moving rarefaction */
      if (u_star>ul) {
        // check left u_star-cl_star
         cminus = u_star-cl_star;
        /* if inside expansion fan */
        if (cminus > dxdt) {
          /* check for pre-expansion state */
          if ((ul-cl) > dxdt) {
            u = ul;
            p = pl;
            d = dl;
          }
          else {
	    /* inegrate in fan until C- curves match TODO make it match new procedure */
            // let Fl=0, int until u-c = dxdt
            Jplus = ul;	// take as constant
            dp = -1.0E-5*dl*cl*(u_star-ul);
            pfan = pl;
            dfan = dl;
            cfan = cl;
            F = 0.0;
            do {
	      F += dp/(dfan*cfan);
	      ufan = Jplus-F;
	      // update variables
	      pfan += dp;
	      /* TODO to generalize, get T D from S, and use T and P to get rho. get c from T.*/
	      he[0] = pfan;
	      he[2] = sl;
	      GSLiterHelmS7(dfan, ye, tfan, he);
              cfan = he[9];
              cminus = ufan-cfan;
            } while (cminus<dxdt);

	    u = ufan;
	    //cfan = u-dxdt;
  	    p = pfan;
            d = dfan;
          }
        }
        /* otherwise solution is known */
        else {
          u = u_star;
          p = pl_star;
          d = dl_star;
        }
      }
      /* leftward moving shock */
      else {
        /* pre-shock state */
        if (dxdt < Vl) {
          u = ul;
          p = pl;
          d = dl;
        }
        /* post shock state */
        else {
          u = u_star;
          p = pl_star;
          d = dl_star;
        }
      }
    }
    /* if (x,t) is on the right of the contact surface... */
    else {
      /* rightward moving rarefaction */
      if (u_star<ur) {
        // check right u_star+cl_star
        cplus = u_star+cr_star;
        /* if inside expansion fan */
        if (cplus < dxdt) {
          /* check for pre-expansion state */
          if ((ur+cr) < dxdt) {
            u = ur;
            p = pr;
            d = dr;
          }
          else {					// TODO


            Jminus = ur;	// constant
            dp = 1.0E-5*dr*cr*(u_star-ur);
            pfan = pr;
            dfan = dr;
            cfan = cr;
            F = 0.0;
            do {
	      F += dp/(dfan*cfan);
	      ufan = Jminus+F;
	      // update variables
	      pfan += dp;
	      /* TODO to generalize, get T from S, and use T and P to get rho. get c from T.*/
	      he[0] = pfan;
	      he[2] = sr;
	      GSLiterHelmS7(dfan, ye, tfan, he);
              cfan = he[9];
              cplus = ufan+cfan;
            } while (cminus>dxdt);
	    u = ufan;
	    //cfan = dxdt-u;
  	    p = pfan;
            d = dfan;
          }
        }
        /* otherwise solution is known */
        else {
          u = u_star;
          p = pr_star;
          d = dr_star;
        }
      }
      /* rightward moving shock */
      else {
        /* pre-shock state */
        if (dxdt > Vr) {
          u = ur;
          p = pr;
          d = dr;
        }
        /* post shock state */
        else {
          u = u_star;
          p = pr_star;
          d = dr_star;
        }
      }
    }
    return;
}


























/////////////////////////////////////////////////////////////////////////////////////////
////////// EOS STUFF /////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


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
  double a, b, c, d, sioncon, amu, h, pi, deni, tempi, tv1, tv2, pres, etot, dpdt, dpdd, dedt, dedd, stot, dsdt, dsdd;
  tv1 = k_b*N_A/abar;
  tv2 = tv1*den;
  pres = tv2*temp;
  dpdt = tv2;
  dpdd = tv1*temp;

  etot = 1.5*tv1*temp;
  dedt = 1.5*tv1;
  dedd = 0.0;

  stot = 0.0; /* Warning do not use */
  //Why is this 0? Should use Sackur-Tetrode eq!
  h = 6.6260689633e-27;
  amu = 1.66053878283e-24;
  pi = 3.14159265358979323846264338327950;
  deni = 1.0/den;
  tempi = 1.0/temp;
  sioncon = (2.0*pi*amu*k_b)/(h*h);
  a = abar*abar*sqrt(abar) * deni/N_A;
  b = temp*sioncon;
  c = a*b*sqrt(b);
  d = log(c);
  stot = (pres*deni + etot)*tempi + tv1*d;

  dsdt = (dpdt*deni + dedt)*tempi - (pres*deni + etot) * tempi*tempi + 1.5*tv1/temp;
  dsdd = (dpdd*deni - pres*deni*deni + dedd)*tempi - tv1/den;
  //dsdt = 1.5*tv1/temp;
  //dsdd = - tv1/den;
  he[0] = pres;
  he[1] = etot;
  he[2] = stot;
  he[3] = dpdd;
  he[4] = dpdt;
  he[5] = dedd;
  he[6] = dedt;
  he[7] = dsdd;
  he[8] = dsdt;
  he[9] = 0; // Added by Erol, to prevent errors in eosTotal()
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
  //std::cout << "Rad      : pressure " << he[0] << " energy " << he[1] << 
  //  " sound " << sqrt(dpdd - dpdt*dsdd/dsdt) <<'\n';
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

/* this function returns the temperature, ye, internal energy, and sound speed given pressure, density, and Y values, using bisection method */
void GSLbisectHelmP7(double den, double ye, double &temp, double he[]) {
  double Tl, Th;
  int status;
  int iter = 0, max_iter = 100;
      struct EOSparams *p;
      p = new EOSparams;
      p->den = den;
      p->ye = 0.5;	// 0.5 for all iso7 species
      p->he[0] = he[0];


      gsl_root_fsolver* solver;
      gsl_function fwrapper;
      solver=gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
      //solver=gsl_root_fsolver_alloc(gsl_root_fsolver_falsepos);
      //solver=gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
      fwrapper.function = GSLiterateHelmP;
      fwrapper.params = p;
      Tl = log(1e5);
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
				status = gsl_root_test_interval(Tl, Th, 1E-6, 1E-5);
      } while (iter < max_iter && status < 0);
      if (iter == max_iter) std::cout << "WARNING:  MAX ITER REACHED" << '\n';

      temp = exp(temp);
      he[1] = p->he[1];
      he[2] = p->he[2];
      he[3] = p->he[3];
      he[4] = p->he[4];
      he[5] = p->he[5];
      he[6] = p->he[6];
      he[7] = p->he[7];
      he[8] = p->he[8];
      he[9] = p->he[9];


      gsl_root_fsolver_free(solver);
      delete p;

      return;
}


double GSLiterateHelmP(double logtemp, void * param) {
  double f, he[10], tol, ptot, etot, stot, dpddtot, dpdttot, deddtot, dedttot, dsddtot, dsdttot, sound, abar;
  struct EOSparams* p = (EOSparams*)param;


  // TODO need abar, just set to something, see sod problem
  abar = 13.7;

  ptot = 0.0;
  etot = 0.0;
  stot = 0.0;
  dpddtot = 0.0;
  dpdttot = 0.0;      
  dsddtot = 0.0;
  dsdttot = 0.0;
  deddtot = 0.0;
  dedttot = 0.0;
  helm(exp(logtemp), p->den, p->ye, he);
  ptot += he[0];
  etot += he[1];
  stot += he[2];
  dpddtot += he[3];
  dpdttot += he[4];   
  deddtot += he[5]; 
  dedttot += he[6];  
  dsddtot += he[7];
  dsdttot += he[8];      
  ion(exp(logtemp), p->den, abar, he);
  ptot += he[0];
  etot += he[1];
  stot += he[2];
  dpddtot += he[3];
  dpdttot += he[4];      
  deddtot += he[5];
  dedttot += he[6];
  dsddtot += he[7];
  dsdttot += he[8];            
  rad(exp(logtemp), p->den, he);
  ptot += he[0];
  etot += he[1];
  stot += he[2];
  dpddtot += he[3];
  dpdttot += he[4];      
  deddtot += he[5];
  dedttot += he[6];
  dsddtot += he[7];
  dsdttot += he[8];
  sound = sqrt(dpddtot - dpdttot*dsddtot/dsdttot);

  p->he[1] = etot;
  p->he[2] = stot;	// check
  p->he[3] = dpddtot;
  p->he[4] = dpdttot;
  p->he[5] = deddtot;
  p->he[6] = deddtot;
  p->he[7] = dsddtot;
  p->he[8] = dsdttot;
  p->he[9] = sound;


/*
  he[0] = pres;
  he[1] = etot;
  he[2] = stot;
  he[3] = dpdd;
  he[4] = dpdt;
  he[5] = dedd;
  he[6] = dedt;
  he[7] = dsdd;
  he[8] = dsdt;
*/






  f = (ptot - p->he[0])/p->he[0];
  return f;
}




/* this function returns the temperature and density given ye, pressure, and entropy, using bisection method */
void GSLbisectHelmS7(double &den, double ye, double &temp, double he[]) {
  double Tl, Th;
  int status;
  int iter = 0, max_iter = 100;
      struct EOSparams *p;
      p = new EOSparams;
      p->ye = 0.5;	// 0.5 for all iso7 species
      p->he[0] = he[0];	// pres
      p->he[2] = he[2];	// entropy


      gsl_root_fsolver* solver;
      gsl_function fwrapper;
      solver=gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
      //solver=gsl_root_fsolver_alloc(gsl_root_fsolver_falsepos);
      //solver=gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
      fwrapper.function = GSLiterateHelmS;
      fwrapper.params = p;
      Tl = log(1e5);
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
				status = gsl_root_test_interval(Tl, Th, 1E-6, 1E-5);
      } while (iter < max_iter && status < 0);
      if (iter == max_iter) std::cout << "WARNING:  MAX ITER REACHED" << '\n';

      temp = exp(temp);

      he[1] = p->he[1];
      he[3] = p->he[3];
      he[4] = p->he[4];
      he[5] = p->he[5];
      he[6] = p->he[6];
      he[7] = p->he[7];
      he[8] = p->he[8];
      he[9] = p->he[9];
      den = p->den;
      gsl_root_fsolver_free(solver);
      delete p;


      return;
}


 

double GSLiterateHelmS(double logtemp, void * param) {
  double f, den, he[10], tol, ptot, etot, dpddtot, dpdttot, dsddtot, dsdttot, temp;
  struct EOSparams* p = (EOSparams*)param;

  ptot = 0.0;
  etot = 0.0;
  dpddtot = 0.0;
  dpdttot = 0.0;      
  dsddtot = 0.0;
  dsdttot = 0.0;

  temp = exp(logtemp);
  he[0] = p->he[0];

  GSLbisectHelmP7d(den, p->ye, temp, he);


  p->he[1] = he[1];
  p->he[3] = he[3];
  p->he[4] = he[4];
  p->he[5] = he[5];
  p->he[6] = he[6];
  p->he[7] = he[7];
  p->he[8] = he[8];
  p->he[9] = he[9];
  p->den = den;
  p->temp = temp;
  //f = (he[0] - p->he[0])/p->he[0];
  f = (he[2] - p->he[2])/p->he[2];		// we now get correct P, given D and T.  check for correct entropy.  Iterate D.



  return f;
}










/* this function returns the temperature, ye, internal energy, and sound speed given pressure, density, and Y values, using bisection method */
void GSLbisectHelmP7d(double &den, double ye, double temp, double he[]) {
  double Dl, Dh;
  int status;
  int iter = 0, max_iter = 100;
      struct EOSparams *p;
      p = new EOSparams;
      p->temp = temp;
      p->ye = 0.5;	// 0.5 for all iso7 species
      p->he[0] = he[0];


      gsl_root_fsolver* solver;
      gsl_function fwrapper;
      solver=gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
      //solver=gsl_root_fsolver_alloc(gsl_root_fsolver_falsepos);
      //solver=gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
      fwrapper.function = GSLiterateHelmPd;
      fwrapper.params = p;
      Dl = log(1e-3);
      Dh = log(1e11);

      gsl_root_fsolver_set(solver, &fwrapper, Dl, Dh);
      //gsl_root_fdfsolver_set (s, &FDF, x);
      iter = 0;
      do {
        iter++;
        status = gsl_root_fsolver_iterate(solver);
        den = gsl_root_fsolver_root(solver);
				Dl = gsl_root_fsolver_x_lower(solver);
				Dh = gsl_root_fsolver_x_upper(solver);
				status = gsl_root_test_interval(Dl, Dh, 1E-6, 1E-5);
      } while (iter < max_iter && status < 0);
      if (iter == max_iter) std::cout << "WARNING:  MAX ITER REACHED" << '\n';

      den = exp(den);
      he[1] = p->he[1];
      he[2] = p->he[2];
      he[3] = p->he[3];
      he[4] = p->he[4];
      he[5] = p->he[5];
      he[6] = p->he[6];
      he[7] = p->he[7];
      he[8] = p->he[8];
      he[9] = p->he[9];


      gsl_root_fsolver_free(solver);
      delete p;

      return;
}



double GSLiterateHelmPd(double logden, void * param) {
  double f, he[10], tol, ptot, etot, stot, dpddtot, dpdttot, dsddtot, dsdttot, deddtot, dedttot, sound, abar;
  struct EOSparams* p = (EOSparams*)param;

  abar = 13.7;

  ptot = 0.0;
  etot = 0.0;
  stot = 0.0;
  dpddtot = 0.0;
  dpdttot = 0.0;      
  dsddtot = 0.0;
  dsdttot = 0.0;
  deddtot = 0.0;
  dedttot = 0.0;
  helm(p->temp, exp(logden), p->ye, he);
  ptot += he[0];
  etot += he[1];
  stot += he[2];
  dpddtot += he[3];
  dpdttot += he[4];   
  deddtot += he[5]; 
  dedttot += he[6];  
  dsddtot += he[7];
  dsdttot += he[8];      
  ion(p->temp, exp(logden), abar, he);
  ptot += he[0];
  etot += he[1];
  stot += he[2];
  dpddtot += he[3];
  dpdttot += he[4];   
  deddtot += he[5]; 
  dedttot += he[6];  
  dsddtot += he[7];
  dsdttot += he[8];      
  rad(p->temp, exp(logden), he);
  ptot += he[0];
  etot += he[1];
  stot += he[2];
  dpddtot += he[3];
  dpdttot += he[4];   
  deddtot += he[5]; 
  dedttot += he[6];  
  dsddtot += he[7];
  dsdttot += he[8];      
  sound = sqrt(dpddtot - dpdttot*dsddtot/dsdttot);

  p->he[1] = etot;
  p->he[2] = stot;	// check
  p->he[3] = dpddtot;
  p->he[4] = dpdttot;
  p->he[5] = deddtot;
  p->he[6] = deddtot;
  p->he[7] = dsddtot;
  p->he[8] = dsdttot;
  p->he[9] = sound;

  f = (ptot - p->he[0])/p->he[0];
  return f;
}















/* this function returns the temperature and density given ye, pressure, and entropy, using bisection method */
void GSLiterHelmS7(double &den, double ye, double &temp, double he[]) {
  double Tl, Th;
  int status;
  int iter = 0, max_iter = 100;
      struct EOSparams *p;
      p = new EOSparams;
      p->ye = 0.5;	// 0.5 for all iso7 species
      p->he[0] = he[0];	// pres
      p->he[2] = he[2];	// entropy


    gsl_multiroot_function f;
    gsl_multiroot_fsolver *workspace_f;
    gsl_multiroot_function_fdf fdf;
    gsl_multiroot_fdfsolver *workspace;
    gsl_vector *x;
   /* allocate a, x, b */
    x = gsl_vector_alloc(2);

    /* Define Solver
    workspace_f = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 2);
    f.f = &GSLiterateHelmSTD;
    f.n = 2;
    f.params = p; */

    /* Define Solver */
    workspace = gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_hybridsj, 2);
    fdf.f = &GSLiterateHelmSTD;
    fdf.df = &dGSLiterateHelmSTD;
    fdf.fdf = &fdGSLiterateHelmSTD;
    fdf.n = 2;
    fdf.params = p;



   /* set initial value */
        gsl_vector_set(x, 0, 1.0e8);
        gsl_vector_set(x, 1, den);


    /* set solver */
    //gsl_multiroot_fsolver_set(workspace_f, &f, x);
    gsl_multiroot_fdfsolver_set(workspace, &fdf, x);

 /* main loop */
    for(iter = 0; iter < max_iter; iter++)
    {

        //std::cout << "starting iter = " << iter << '\n';
        //status = gsl_multiroot_fsolver_iterate(workspace_f);
	status = gsl_multiroot_fdfsolver_iterate(workspace);
        // std::cout << "status = " << status << '\n';

        /*printf("%d iter: ",iter);
        for(int i = 0; i < 2; i++)
            printf("%10.3e ", gsl_vector_get(workspace->x, i));
        printf("\n");*/
 
	//status =  gsl_multimin_test_gradient(x, 1E-6)
	status = gsl_multiroot_test_residual(workspace->f, 1E-4);
         //std::cout << "status2 = " << status << '\n' << '\n';
	if (status > -1)
	  break;


        if((status == GSL_EBADFUNC) || (status == GSL_ENOPROG))
        {
            printf("Status: %s\n", gsl_strerror(status));
            break;
        }
    }


//solution
 temp = gsl_vector_get(workspace->x, 0);
 den = gsl_vector_get(workspace->x, 1);


      he[1] = p->he[1];
      he[3] = p->he[3];
      he[4] = p->he[4];
      he[5] = p->he[5];
      he[6] = p->he[6];
      he[7] = p->he[7];
      he[8] = p->he[8];
      he[9] = p->he[9];


      //gsl_multiroot_fsolver_free(workspace_f);
      gsl_multiroot_fdfsolver_free(workspace);
      delete p;


      return;
}



int GSLiterateHelmSTD(const gsl_vector *x, void *param, gsl_vector *ret)
{
  double f, den, he[10], tol, ptot, etot, dpddtot, dpdttot, dsddtot, dsdttot, temp;
  struct EOSparams* p = (EOSparams*)param;

/*  ptot = 0.0;
  etot = 0.0;
  dpddtot = 0.0;
  dpdttot = 0.0;      
  dsddtot = 0.0;
  dsdttot = 0.0;*/
  
    temp = gsl_vector_get(x, 0);
    den = gsl_vector_get(x, 1);

   // std::cout << "  trying T = " << temp << " den = " << den << '\n';
   //helm(temp, den, p->ye, he);
   eosTotal(temp, den, p->ye, he, 13.7);
 
  p->he[1] = he[1];
  p->he[3] = he[3];
  p->he[4] = he[4];
  p->he[5] = he[5];
  p->he[6] = he[6];
  p->he[7] = he[7];
  p->he[8] = he[8];
  p->he[9] = he[9];




    gsl_vector_set(ret, 0, (he[0] - p->he[0])/p->he[0]);
    gsl_vector_set(ret, 1, (he[2] - p->he[2])/p->he[2]);

   //std::cout << "     err_p = " << (he[0] - p->he[0])/p->he[0] << '\n';
   //std::cout << "     err_s = " << (he[2] - p->he[2])/p->he[2] << '\n';


 
    return GSL_SUCCESS;
}




int dGSLiterateHelmSTD(const gsl_vector *x, void *param, gsl_matrix *ret)
{
  double f, den, he[10], tol, ptot, etot, dpddtot, dpdttot, dsddtot, dsdttot, temp;
  struct EOSparams* p = (EOSparams*)param;

/*  ptot = 0.0;
  etot = 0.0;
  dpddtot = 0.0;
  dpdttot = 0.0;      
  dsddtot = 0.0;
  dsdttot = 0.0;*/
 
    temp = gsl_vector_get(x, 0);
    den = gsl_vector_get(x, 1);
    //helm(temp, den, p->ye, he);
    eosTotal(temp, den, p->ye, he, 13.7);

    gsl_matrix_set(ret, 0, 0, he[4]/p->he[0]);	// dP/dT div by P
    gsl_matrix_set(ret, 0, 1, he[3]/p->he[0]);	// dP/dD div by P
 
    gsl_matrix_set(ret, 1, 0, he[8]/p->he[2]);	// dS/dT div by S
    gsl_matrix_set(ret, 1, 1, he[7]/p->he[2]);	 // dP/dD div by S
 
    return GSL_SUCCESS;
}




/* simultaneous eval of f and df */
int fdGSLiterateHelmSTD(const gsl_vector *x, void *param, gsl_vector *ret_vec, gsl_matrix *ret_mat)
{
    GSLiterateHelmSTD(x, param, ret_vec);
    dGSLiterateHelmSTD(x, param, ret_mat);
 

    return GSL_SUCCESS;
}

// function for the total EOS, includes helm(), rad and ion contributions
void eosTotal(double temp, double den, double ye, double he[], double abar){
  // void helm(double temp, double den, double ye, double he[]);
  // void ion(double temp, double den, double abar, double he[]);
  // void rad(double temp, double den, double he[]);
  double he_temp[10];
  for(int i=0;i<9;i++)
      he[i]=0.0;
  helm(temp, den, ye, he);
  ion(temp, den, abar, he_temp);
  for(int i=0;i<9;i++)
      he[i] += he_temp[i];
  rad(temp, den, he_temp);
  for(int i=0;i<9;i++)
      he[i] += he_temp[i];
  he[9] = he[3] - he[4] * he[7] / he[8];
  he[9] = sqrt(he[9]);

  return;
}


