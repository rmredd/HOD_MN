#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header.h"

double m_g3,r_g3,*x_g3,*y_g3,*z_g3;
double n_g3;
int flag_g3=0;

double func_ng(double m);
double func_ng2(double m);

/* The restr
 */

double restricted_number_density(double r)
{
  static int flag=1;
  static double *x,*y,*y2;
  int i,n=50,j;
  double mlimit,dlogm,logm,mmin,sum=0,t0,t1,s1,s2,s3,m,r1,r2,ng2,rlim,rmin,rmax;

  if(flag)
    {
      n_g3 = n;
      x_g3=dvector(1,n);
      y_g3=dvector(1,n);
      z_g3=dvector(1,n);
      flag=0;
    }

  /* Reset the static variables in this function.
   */
  func_ng2(-1);

  r_g3=r;
  ng2=GALAXY_DENSITY*GALAXY_DENSITY;

  /* Calculate the maximum allowable halo mass, which had 
   * rvir = r_g3 - rvir(M_low).
   */
  r1=pow(3.*HOD.M_low/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
  rlim = r_g3 - r1;
  mlimit=log(4./3.*DELTA_HALO*RHO_CRIT*PI*rlim*rlim*rlim*OMEGA_M);
  if(mlimit>log(HOD.M_max))mlimit=log(HOD.M_max);
  mmin=log(HOD.M_low);

  if(HOD.color==2)
    {
      dlogm=(mlimit-mmin)/(n-1);
      m = mmin;
      for(i=1;i<=n;++i)
	{
	  if(N_avg(exp(m))>0)break;
	  m += dlogm;
	}
      mmin = m;
      r1=pow(3.*exp(mmin)/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
      rlim = r_g3 - r1;
      mlimit=log(4./3.*DELTA_HALO*RHO_CRIT*PI*rlim*rlim*rlim*OMEGA_M);
    }

  if(EXCLUSION==2) {
    dlogm=(mlimit-mmin)/(n-1);
    x_g3[1] = mmin;
    y_g3[1] = qromo(func_galaxy_density,mmin,mmin+dlogm,midpnt);
    for(i=2;i<=n;++i) 
      {
	x_g3[i] = i*dlogm+mmin;
	y_g3[i] = y_g3[i-1] + qromo(func_galaxy_density,(i-1)*dlogm+mmin,mmin+i*dlogm,midpnt);
      }    
    spline(x_g3,y_g3,n,1.0E+30,1.0E+30,z_g3);
    s1 = qromo(func_ng2,mmin,mlimit,midpnt);
    return(sqrt(s1));
  }

  /* Calculate the double integral at specified masses.
   */
  dlogm=(mlimit-mmin)/(n-1);
  for(i=1;i<=n;++i)
    {
      logm=(i-0.5)*dlogm+mmin;
      m_g3=exp(logm);
      r2 = pow(3*m_g3/(4*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
      if(EXCLUSION==3) {
	if(ellipsoidal_exclusion_probability(r1/r2,r_g3/(r1+r2))==0)break; }
      else {
	if(r1+r2>r_g3)break; }
      s1=qtrap(func_ng,mmin,mlimit,1.0E-4);
      sum+=s1*m_g3*dlogm;
      if(s1==0)break;
      if(sum>=ng2)break;
    }
  return sqrt(sum);
}

double func_ng2(double m)
{
  static double fac2=-1,fac1=-1;
  double s1,rv1,n,N,m1,mx;

  if(m<0)
    {
      fac1=fac2=-1;
      return(0);
    }

  m1=exp(m);
  if(fac2<0)
    fac2=pow(3.0/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
  if(fac1<0)
    fac1=4./3.*PI*RHO_CRIT*DELTA_HALO*OMEGA_M;

  rv1 = r_g3 - pow(m1,1.0/3.0)*fac2;
  rv1 = rv1;
  mx = fac1*rv1*rv1*rv1;

  n=dndM_interp(m1);
  N=N_avg(m1);
  splint(x_g3,y_g3,z_g3,n_g3,log(mx),&s1);
  return(n*N*s1*m1);
}


double func_ng(double m)
{
  static double fac2=-1;
  double s1,rv1,rv2,exfac=1,n,N;

  m=exp(m);
  if(fac2<0)
    fac2=pow(3.0/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rv1=pow(m_g3,1.0/3.0)*fac2;
  rv2=pow(m,1.0/3.0)*fac2;

  if(EXCLUSION==3) 
    {
      if(0.5*(rv1+rv2)>r_g3)return(0);
      if(1.5*(rv1+rv2)>r_g3)exfac=ellipsoidal_exclusion_probability(rv2/rv1,r_g3/(rv2+rv1));
    }
  else
    {
      if(rv1+rv2>r_g3)return(0);
    }

  n=dndM_interp(m)*dndM_interp(m_g3);
  N=N_avg(m)*N_avg(m_g3);
  return(exfac*n*N*m);
}

/* This is the probability that two halos do not overlap, given their
 * radii and separation. Of course, for spherical halos P(x) is a step function
 * at x = (r1+r2)/r_sep  = 1, but for ellipsoidal halos there is a chance 
 * that they could be closer. In detail, P(x) changes depending on the mass
 * ratio of the halos, but using tabulated values does not appear to make
 * significant difference in the results for xi_2h(r). The function below is
 * a fit to Monte Carlo results for a halos with a distribution of axis ratios
 * which is lognormal in e_b = (1-b/a) and e_c = (1-c/a) with dispersions of 0.2
 * mean <b/a>=0.9 and <c/a>=0.8 (pretty reasonable values).
 */
double ellipsoidal_exclusion_probability(double rv, double r)
{
  static int flag=0,nr=101,nratio=31;
  static double **xprob,*rad,*ratio,rhi,rlo,mhi,mlo,dr,dm;
  float x1,x2,x3;
  int i,j,im,ir;
  FILE *fp;

  if(rv<1)rv=1.0/rv;
  
  r=(r-0.8)/0.29;
  if(r>1)return(1.0);
  if(r<0)return(0.0);
  return(3*r*r-2*r*r*r);
}  
