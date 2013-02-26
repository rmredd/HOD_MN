#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "header.h"

double *mu_x,*mu_y,*mu_z;
int mu_n;

/* Equation (10) of vdB et al, assuming alpha = 1 (standard NFW)
 */
double func_jeans1(double x)
{
  return(x/(1+x)/(1+x));
}

/* Equation (11) of vdB et al.
 */
double func_jeans2(double x)
{
  double y;
  splint(mu_x,mu_y,mu_z,mu_n,x,&y);
  return(y/(x*x*x*(1+x)*(1+x)));
}
double jeans_dispersion(double mass, double rx, double v[])
{
  int i,j,k,n,nc,nr;
  double rmin,rmax,dlogr;
  static double cmin,cmax,dlogc;
  double rvir,rs,s1,sig,mu1,vcirc,cvir;

  static long IDUM2 = -4555;

  static double fac=-1;
  static int prev_cosmo=-1;

  if(fac<0)
    {
      mu_n = nc = 1000;
      cmin = 0.1;
      cmax = 100.0;
      dlogc = (log(cmax) - log(cmin))/(nc-1);
      
      mu_x = dvector(1,nc);
      mu_y = dvector(1,nc);
      mu_z = dvector(1,nc);
      
      for(i=1;i<=nc;++i)
	{
	  mu_x[i] = exp((i-1)*dlogc)*cmin;
	  mu_y[i] = qromo(func_jeans1,0,mu_x[i],midpnt);
	}
      spline(mu_x,mu_y,nc,1.0E+30,1.0E+30,mu_z);
      fac=sqrt(4.499E-48)*3.09E19;
    }

  cvir = halo_concentration(mass);
  rvir = pow(3*mass/(4*PI*DELTA_HALO*OMEGA_M*RHO_CRIT),THIRD);
  rs = rvir/cvir;
  rx = rx/rs;
  vcirc = fac*fac*mass/rvir;
  splint(mu_x,mu_y,mu_z,mu_n,cvir,&mu1);

  s1 = qromo(func_jeans2,rx,cmax,midpnt);
  sig = sqrt(cvir*vcirc/mu1*rx*(1+rx)*(1+rx)*s1);

  if(isnan(sig))
    sig = sqrt(vcirc/2);

  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sig*VBIAS;

  return sig;

  nr = 100;
  rmin = rvir*0.01;
  rmax = rvir*0.99;
  dlogr = (log(rmax) - log(rmin))/(nr-1);
  rs = rvir/cvir;

  vcirc = fac*fac*mass/rvir;
  splint(mu_x,mu_y,mu_z,mu_n,cvir,&mu1);

  for(i=0;i<nr;++i)
    {
      rx = exp(i*dlogr)*rmin/rs;
      s1 = qromo(func_jeans2,rx,cmax,midpnt);
      sig = cvir/mu1*rx*(1+rx)*(1+rx)*s1*2;
      printf("%f %f\n",rx*rs/rvir,sqrt(sig));
    }
  exit(0);

}
