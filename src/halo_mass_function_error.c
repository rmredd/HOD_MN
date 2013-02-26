#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* This is a calculation of the differential halo mass function of
 * Jenkins et al 2001 (MNRAS 321, 372). 
 *
 * Also included is the mass function from my own analysis, which is 
 * a variant of the Sheth-Tormen mass function (leaving all 3 parameters
 * as free parameters, with a 2-parameter high-mass cutoff of the form
 * suggested by Reed et al.
 *
 */

double halo_mass_function(double mass)
{
  double sig,logm,a,slo,shi,rm,rlo,rhi,mlo,mhi,dsdM,n,nuprime,nufnu,p,A;
  static int flag=0, prev_cosmo=-1,na=4;
  static double pnorm,*xa,*ya,*za;

  double fac1,fac2,fac;
  static long IDUM=-555;
  int i;
  static double
    a1 = 0.325277,
    a2 = 0.492785,
    a3 = 0.310289,
    a4 = 1.317104,
    a5 = 2.425681;
    
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
  logm=log10(mass);

  mlo=0.99*mass;
  mhi=1.01*mass;
  rlo=pow(3.0*mlo/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  rhi=pow(3.0*mhi/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  slo=pnorm*sigmac(rlo);
  shi=pnorm*sigmac(rhi);
  dsdM=(shi-slo)/(mhi-mlo);

  if(BEST_FIT)
    {
      n = a1*sqrt(2*a2/PI)*(1+pow(sig*sig/(a2*DELTA_CRIT*DELTA_CRIT),a3))
	*DELTA_CRIT/sig*exp(-a2*DELTA_CRIT*DELTA_CRIT/(2*sig*sig))
	*exp(-a4/sig/pow(fabs(cosh(2*sig)),a5));
      
      splint(xa,ya,za,na,log(1/sig),&fac1);
      fac2 = 0.00562*pow(sig,-2.9) + 0.00158*pow(sig,2.2);
      fac = 1 + fac1*fac2;
      if(fac<0.01)fac=0.01;
      n = n*fac;
      /* printf("%e %e %e %e\n",sig,n,n/fac,fac); */
      n = -n*OMEGA_M*RHO_CRIT/mass/sig*dsdM;
      return(n);
    }

  /* Jenkins et al.
   */
  a=-JENKINS_A*OMEGA_M*RHO_CRIT/mass/sig;
  n=a*dsdM*exp(-pow(fabs(JENKINS_B-log(sig)),JENKINS_C));
  return(n);
  
}


/* It may be a bit costly to run the above function every time you need
 * dn/dM, so here we put the values into an array and then interpolate.
 */
double dndM_interp(double m)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2;
  int i,n=1000;
  double dm,max=16.7,min=9,a,m1,m2,dm1;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting mass function for %f %f\n",OMEGA_M,SIGMA_8);
      fflush(stdout);

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
	  y[i]=log(halo_mass_function(x[i]));
	  if(DELTA_HALO!=200)
	    {
	      m1 = halo_mass_conversion2(x[i]*1.001,halo_c200(x[i]*1.001),DELTA_HALO,200.0);
	      m2 = halo_mass_conversion2(x[i]*0.999,halo_c200(x[i]*0.999),DELTA_HALO,200.0);
	      dm1 = (x[i]*1.001 - x[i]*0.999)/(m1-m2);
	      y[i] = y[i] + log10(dm1);
	      x[i]=log10(halo_mass_conversion2(x[i],halo_c200(x[i]),DELTA_HALO,200.0));
	    }
	  else
	    x[i]=log(x[i]);
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;
    }

  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(exp(a));

}




