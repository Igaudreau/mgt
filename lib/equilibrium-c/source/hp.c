/*
C library for equilibrium chemistry calculations
    References:
        "Computer Program for Calculation of Complex Equilibrium Compositions and Applications"
        Nasa Reference Publication 1311, October 1995
        Sanford Gordon and Bonnie J. McBride

        "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"
        NASA/TP - 2002-211556, September 2002
        Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "thermo.h"
#include "linalg.h"
#include "common.h"
#include "hp.h"
#include "pt.h"




static double temperature_guess(int nsp, double ht, double pt, double n, double* ns, double* lewis){
    /*
    Guess a first iteration temperature assuming constant Cp from 298 K and target p and h
    Inputs:
        nsp   : Number of species
        ht    : Target mixture specific enthalpy (J/kg)
        pt    : Target pressure (Pa)
        n     : Mixture mixture weight [nsp]
        ns    : Intiial composition mixture weights [nsp]
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]

    Output:
        T : Temperature Guess (K)
    */
    int s;
    double* lp;
    double cp,h0,Cps298,H0s298,nss,T;

    cp = 0.0;
    h0 = 0.0;
    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        Cps298 = compute_Cp0_R(298.15, lp)*Ru;
        H0s298 = compute_H0_RT(298.15, lp)*Ru*298.15;

        nss = ns[s];
        cp += nss*Cps298;
        h0 += nss*H0s298;
    }
    T = 298.15+(ht - h0)/cp;

    //printf("s0: %e\n", s0);
    T = fmin(fmax(T, 200.0),20000.0); // Limit in case of bad initial state.


    return T;
}

int solve_hp(double pt,double ht,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
               double* X1, double* Teq, int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed pressure and specific entropy
    Inputs:
        pt    : target pressre (Pa)
        ht    : target specific enthalpy (J/kg)
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1  : Equilibrium Mole Fraction [nsp]  
        Teq: Equilibrium Temperature 
    */


    double *H_RTs, *Cp_Rs, *ns, *bi0; // Dynamic arrays
    double n,M1,T,error;
    int s,k,errorcode;
    H_RTs = (double*) malloc(sizeof(double)*nsp);     // Species Molar Enthalpy 
    Cp_Rs = (double*) malloc(sizeof(double)*nsp);     // Species Specific Heat @ Const. pressure 
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    double *lp;
    double h, cp, f, df;
    double Thigh, Tlow;
    const double atol=1e-6;

    // For now call PT solver and iterate using Newton iteration
    // guess a temperature
    composition_guess(a, M, X0, nsp, nel, ns, &n, bi0);
    T = temperature_guess(nsp, ht, pt, n, ns, lewis);

    // first apply bisection to get close to the solution
    Thigh = 9000.0;
    Tlow = 200.0;
    for (k=0; k<20; k++){
        solve_pt(pt, T, X0, nsp, nel, lewis, M, a, X1, verbose);		// equibrate at this T and p

        // calculate specific h and cp (dh/dT)
        composition_guess(a, M, X1, nsp, nel, ns, &n, bi0);		// obtains specific ns (in mol/kg)
        h = 0.0;
        cp = 0.0;
        for (s=0; s<nsp; s++){
            lp = lewis + 9*3*s;
            H_RTs[s] = compute_H0_RT(T, lp);		// Divided by RuT
            Cp_Rs[s] = compute_Cp0_R(T, lp);		// Divided by Ru
            h    += ns[s]*H_RTs[s];
            cp   += ns[s]*Cp_Rs[s];
        }
        h *= Ru*T;	// in J/kg
        cp *= Ru;	// in J/kgK
/*        printf("AFTER BISECTION\n");
        printf("T = %f \n",T);
        printf("h = %f \n",h);
        printf("cp = %f \n",cp);*/

        f = h - ht;		// iterate until target f reached

        if (f > 0.0) {
          // h, T is too high
          Thigh = T;
          T = T-0.5*(T-Tlow);
        }
        else {
          // h, T is too low
          Tlow = T;
          T = T+0.5*(Thigh-T);
        }
    }

    // compare h to ht and per form the Newton's method here.
    errorcode=0;
    error=1e99;
    f = 0.0;
    df = 1.0;
    // Begin Iterations
    for (k=0; k<attempts; k++){
        T = T - 0.1*(f/df);		// for some reason need to dampen the Newton iertation to ensure stability...
        solve_pt(pt, T, X0, nsp, nel, lewis, M, a, X1, verbose);		// equibrate at this T and p

        // calculate specific h and cp (dh/dT)
        composition_guess(a, M, X1, nsp, nel, ns, &n, bi0);		// obtains specific ns (in mol/kg)
        h = 0.0;
        cp = 0.0;
        for (s=0; s<nsp; s++){
            lp = lewis + 9*3*s;
            H_RTs[s] = compute_H0_RT(T, lp);		// Divided by RuT
            Cp_Rs[s] = compute_Cp0_R(T, lp);		// Divided by Ru
            h    += ns[s]*H_RTs[s];
            cp   += ns[s]*Cp_Rs[s];
        }
        h *= Ru*T;	// in J/kg
        cp *= Ru;	// in J/kgK
        f = h - ht;		// iterate until target f reached
        df = cp;

/*        printf("AFTER NEWTON\n");
        printf("T = %f \n",T);
        printf("h = %f \n",h);
        printf("cp = %f \n",cp);
        printf("f = %e \n",f);
        printf("df = %e \n",df);
        printf("atol = %e \n",atol);*/

        error = fabs(f);	// check for absolute convergence


        if (verbose>0){
            printf("iter %2d: [%f]",k,T);
            for (s=0; s<nsp; s++) printf(" %f",ns[s]);
            printf("  (%e)\n", error);
        }

        // Exit loop if convergence is achieved
        if (error<atol) break;

        // Exit loop if nans appearing in dlnns
        if (isnan(error)) {
            printf("Solver nan'd, exiting!\n");
            errorcode=1;
            break;
        }

        if (k>=attempts-1) {
            printf("Solver not converged, exiting!\n");
            errorcode=1;
            break;
        }
    }

    if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, error);
    if ((verbose>0)&&(errorcode>0)) printf("Convergence failure in %d iter, error: %e\n", k, error);


    // Compute output composition
    n = 0.0;
    for (s=0; s<nsp; s++) n += ns[s];
    M1 = 1.0/n;
    for (s=0; s<nsp; s++) X1[s] = M1*ns[s];
    *Teq = T;

    free(H_RTs);
    free(Cp_Rs);
    free(ns);
    free(bi0);
    return errorcode;
}

#ifdef TEST
int main(){
    printf("Called main in hp.c! Bad!\n");
    return 0;
}
#endif
