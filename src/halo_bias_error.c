#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"


/* This is the bias factor for halos. Mass (in h^{-1} M_sol) is the input.
 * This is not the scale-dependent bias factor, which can be calculated later.
 * 
 * The parameterization of b(M) is taken from Sheth Mo & Tormen (2001 MNRAS 232 1)
 * but the actual parameters have been replaced by the values of Appendix A in 
 * Tinker, Weinberg, Zheng, & Zehavi. 2005 Apj (submitted)
 */

double bias(double mass)
{
  static int prev_cosmo=-1;
  static long IDUM=-3333;
  double rm,sig,k,neff,b,logk,khi,klo,phi,plo,nu,psp,x,fac,fac1,fac2;
  static int flag=0,na=4;
  static double pnorm,*xa,*ya,*za;
  int i;

  /* Original SMT parameters */
  double a=0.707,bb=0.5,c=0.6;

  /* Tinker et al parameters */
  a=0.707;
  bb=0.35;
  c=0.8;

  /* Fitting to Mike Warren's simulations. */
  bb=0.28; 

  /* Use the globel parameters. */
  a=BIAS_A; bb=BIAS_B; c=BIAS_C;

  if(!flag)
    {
      xa=dvector(1,na);
      ya=dvector(1,na);
      za=dvector(1,na);
      flag=1;
      xa[1] = log(1/2.51);
      xa[2] = log(1/1.49);
      xa[3] = log(1/0.905);
      xa[4] = log(1/0.501);
    }

  if(prev_cosmo != RESET_COSMOLOGY)
    {
      prev_cosmo=RESET_COSMOLOGY;
      for(i=1;i<=4;++i) 
	ya[i] = gasdev(&IDUM);
      spline(xa,ya,na,1.0E+30,1.0E+30,za);
    }


  /* First normalize the power spectrum
   */
  pnorm=SIGMA_8/sigmac(8.0);
  rm=pow(3.0*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  sig=pnorm*sigmac(rm);

  nu=DELTA_CRIT/sig;
  b=1+1.0/(sqrt(a)*DELTA_CRIT)*(sqrt(a)*a*nu*nu + sqrt(a)*bb*pow(a*nu*nu,1-c) - 
				(pow(a*nu*nu,c)/(pow(a*nu*nu,c)+bb*(1-c)*(1-c/2.))));

  splint(xa,ya,za,na,log(1/sig),&fac1);
  fac2 = 0.0158*pow(sig,-2.7) + 0.00251*pow(sig,2.2);
  fac = 1 + fac1*fac2;
  if(fac<0.01)fac=0.01;
  b*=fac;

  return(b);


  /* This is the Seljak & Warren (2004) bias.
   * There's an alpha_s in the correction term, which here I set to zero. (RUNNING.)
   * See App. A of Tinker et al for a comparison of this bias with above bias.
   */
  x=mass/MSTAR;
  b=(0.53+0.39*pow(x,0.45)+0.13/(40*x+1) + 5.0e-4*pow(x,1.5));
  b=b+log10(x)*(0.4*(OMEGA_M-0.3+SPECTRAL_INDX-1)+0.3*(SIGMA_8-0.9+HUBBLE-0.7)+0.8*0.0);
  return(b);

  /* This is the old Mo & White (1996) formula
   */
  return(1+DELTA_CRIT/sig/sig-1/DELTA_CRIT);
 

}

/* Just like the halo mass function, we'll set this up such that
 * you just need to interpolate.
 *
 * Now this has been changed to calculate the spatial-scale dependence
 * of the bias factor. If you don't want the scale dep. b, then just
 * input r<0.
 *
 * If the global flag LINEAR_PSP==0, uses the scale dependence calculated for
 * for halo bias relative to the non-linear matter \xi_m(r):
 *
 * f^2(r) = (1.0+xi*1.17)^1.49/(1.0+xi*0.69)^2.09  --> b(r) = b0*f(r)
 *
 * For LINEAR_PSP==1, use scale dependence determined for the linear P(k):
 *
 * f(r) = 1 + exp[-(r/A)^0.7] --> where A is a parameter that we're gonna have to 
 *                                determine in more detail, but now A \approx 1
 */
double bias_interp(double m, double r)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2;
  int i,n=100;
  double dm,max=16.3,min=9,a,b,m1,m2,dm1,xi;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting bias for %f %f\n",OMEGA_M,SIGMA_8);
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      dm=(double)(max-min)/n;
      for(i=1;i<=n;++i)
	{
	  x[i]=pow(10.0,min+i*dm);
	  y[i]=log(bias(x[i]));
	  if(DELTA_HALO!=200)
	    {
	      m1 = halo_mass_conversion2(x[i]*1.001,halo_c200(x[i]*1.001),DELTA_HALO,200.0);
	      m2 = halo_mass_conversion2(x[i]*0.999,halo_c200(x[i]*0.999),DELTA_HALO,200.0);
	      dm1 = (x[i]*1.001 - x[i]*0.999)/(m1-m2);
	      y[i] = y[i] + log(dm1);
	      x[i]=log(halo_mass_conversion2(x[i],halo_c200(x[i]),DELTA_HALO,200.0));
	    }
	  else
	    x[i]=log(x[i]);
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;

    }

  m=log(m);
  splint(x,y,y2,n,m,&a);
  if(r<0 || r>10)return(exp(a));

  if(LINEAR_PSP)
    {
      b=exp(-pow(r,0.7)) + 1;
    }
  else
    {
      xi=xi_interp(r);
      b=pow(1.0+xi*1.17,1.49*0.5)*pow(1.0+xi*0.69,-2.09*0.5);
    }
  a=exp(a)*b;
  return(a);

}


/* This is the integrand which qromo or qtrap would call
 * to calculate the large-scale galaxy bias.
 * The integral is a number-weighted average of the halo
 * bias function, integrated over the halo mass function.
 */
double func_galaxy_bias(double m)
{
  return(dndM_interp(m)*N_avg(m)*bias_interp(m,-1.));
}
