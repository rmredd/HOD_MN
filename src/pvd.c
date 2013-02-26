#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header.h"

double two_halo_pvd(double r, double *mu, double *vt);
double one_halo_pvd(double r);
double r_g12;
double pvd1h_func1(double m);
double check_2halo_pairs(double r, double *n);

void pairwise_velocity_dispersion()
{
  int i,j,k,nr=30;
  double r,rmin,rmax,dlogr,pvd,pvd1,pvd2,pvd3,x1,x2,x3,mu,pvdt2,pvdt;
  FILE *fp;
  char fname[1000];
  
  sprintf(fname,"%s.PVD",Task.root_filename);
  fp = fopen(fname,"w");

  OUTPUT=2;

  rmin = 0.05;
  rmax = 30.0;
  dlogr = (log(rmax)-log(rmin))/(nr-1);

  /* check_2halo_pairs(2.0,&x3); */

  for(i=1;i<=nr;++i)
    {
      r = exp(dlogr*(i-1))*rmin;
      x1 = one_halo_real_space(r);
      x2 = 1 + two_halo_real_space(r);
      pvd1 = one_halo_pvd(r);
      if(r>R_MIN_2HALO)
	{
	  pvd2 = two_halo_pvd(r,&mu,&pvdt2);
	  /*
	  pvd3 = check_2halo_pairs(r,&x3);
	  printf("BOO %f %f %e %e\n",pvd2,sqrt(pvd3),x2,x3);
	  pvd2 = sqrt(pvd2*pvd2 + pvd3*x3)/(1+x3);
	  x2 = x2+x3;
	  */
	}
      else 
	pvd2 = 0;
      pvd = sqrt((pvd1*pvd1*x1 + pvd2*pvd2*x2)/(x1+x2));
      pvdt = sqrt((pvd1*pvd1*x1 + pvdt2*pvdt2*x2)/(x1+x2));
      x3 = mu*x2/(x1+x2);
      fprintf(fp,"%f %f %f %f %f %f %f %e %e\n",r,pvd,pvd1,pvd2,pvdt,pvdt2,x3,x1,x2);
      fflush(fp);
    }
  fclose(fp);
}

double two_halo_pvd(double r, double *mu, double *vt)
{
  double v,dv=1,p1=0,p2=0,p=0,p0=0;
  int i,n=4000;

  for(i=-n;i<n;++i)
    {
      v = (i+0.5)*dv;
      p = galaxy_prob_vz(v,r,0.0);
      if(isnan(p))continue;
      p1 += v*p*dv;
      p2 += v*v*p*dv;
      p0 += p*dv;
    }
  if(p0==0) {*vt = 0; *mu=0; return 0;}
  p1 = p1/p0;
  *vt = sqrt(p2/p0 - p1*p1);
  
  p0 = p1 = p2 = 0;
  for(i=-n;i<n;++i)
    {
      v = (i+0.5)*dv;
      p = galaxy_prob_vz(v,r,PI/2*0.99);
      if(isnan(p))continue;
      p1 += v*p*dv;
      p2 += v*v*p*dv;
      p0 += p*dv;
      //printf("%d %f %f %e %e\n",i,r,v,p,p0);
    }
  *mu = p1/p0;
  return(sqrt(p2/p0-p1/p0*p1/p0));
}

double check_2halo_pairs(double r, double *n)
{
  double s1,s2,s3,func_pvd1(),func_pvd2(),func_pvd3();

  s1 = qromo(func_pvd1,log(HOD.M_low),log(HOD.M_low*10),midpnt);
  s2 = qromo(func_pvd2,log(1e+14),log(HOD.M_max),midpnt);
  s3 = qromo(func_pvd3,log(1e+14),log(HOD.M_max),midpnt);
  printf("%e %e %e\n",s1*s2,s3/s2,(1+two_halo_real_space(r))*GALAXY_DENSITY*GALAXY_DENSITY);
  *n = s1*s2/(GALAXY_DENSITY*GALAXY_DENSITY);
  return(s3/s2);
}

double func_pvd1(double m)
{
  m=exp(m);
  return(N_avg(m)*dndM_interp(m)*m);
}
double func_pvd2(double m)
{
  m=exp(m);
  return(N_sat(m)*dndM_interp(m)*m);
}
double func_pvd3(double m)
{
  double sigv_sat;
  static double fac=-1;
  m=exp(m);
  if(fac<0)
    fac=sqrt(4.499E-48/2.0)*
      pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sigv_sat=VBIAS*fac*pow(m,0.3333333);  
  return(N_sat(m)*dndM_interp(m)*m*(sigv_sat*sigv_sat + 4.0e4));
}
    
double one_halo_pvd(double r)
{
  double rhi, mlo, fac, s1;
  
  r_g12 = r;
  rhi=(1.9*pow(3*HOD.M_max/(4*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0));
  if(r>rhi)return 0;

  fac=1.0/(2*PI*r_g12*r_g12*GALAXY_DENSITY*GALAXY_DENSITY);
  
  mlo = 4./3.*PI*RHO_CRIT*DELTA_HALO*OMEGA_M*pow(r*.5,3.0);
  if(mlo<HOD.M_low)
    mlo = HOD.M_low;
  
  s1=fac*qromo(pvd1h_func1,log(mlo),log(HOD.M_max),midpnt)/one_halo_real_space(r);
  
  return sqrt(s1); 
}


double pvd1h_func1(double m)
{
  double N,n,fac2,rvir,f_ss,f_cs,cvir,x,rfof,ncen,nsat,sigv_sat,sigv_cen;
  static double fac=-1;
  static int prev_cosmo=0;

  if(fac<0 || RESET_COSMOLOGY!=prev_cosmo)
    {
      prev_cosmo=RESET_COSMOLOGY;
      fac=sqrt(4.499E-48/2.0)*
	pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
    }

  m=exp(m);
  sigv_sat=VBIAS*fac*pow(m,0.33333333333)*ROOT2;
  sigv_cen=(fac*pow(m,0.33333333333))*sqrt(VBIAS_C*VBIAS_C + VBIAS*VBIAS);
  cvir=halo_concentration(m)*CVIR_FAC;
  n=dndM_interp(m);  
  nsat=N_sat(m);
  ncen=N_cen(m);
  
  rvir=2*pow(3.0*m/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  /* Break up the contribution of pairs into
   * central-satellite (cs) and satellite-satellite (ss) pairs.
   */
  f_ss=dFdx_ss(r_g12/rvir,cvir)*moment_ss(m)*0.5*sigv_sat*sigv_sat;
  f_cs=dFdx_cs(r_g12/rvir,cvir)*nsat*ncen*sigv_cen*sigv_cen;
  x=n*(f_ss+f_cs)/rvir*m;

  return(x);

}
