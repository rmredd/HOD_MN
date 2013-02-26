#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"


double func_cvir1(double a);
double func_cvir2(double a);
double cvir_2pop(double mass, float K, double delta, double *mdelta);
double c200c_maccio(double mass, int iwmap, double delta, double *mdelta);

double collapse_redshift(double z);
double cvir_pnorm_g1,
  cvir_sigma_g1;
double cvir_g1, cvir_anow;

/* This calculates and tabulates the halo concentrations
 * as a function of halo mass. Uses the "Bullock model", 
 * described in a little more detail below.
 */

double halo_concentration(double m)
{
  static int flag=1,n,prev_cosmo=0;
  static double *x,*y,*y2;
  int i;
  float x1,x2,cfac;
  double a,dm,x3,x4,m1;
  FILE *fp;
  char fname[1000];

  if(flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      MSTAR = mstar();
      if(OUTPUT)
	fprintf(stdout,"Calc cvir with DELTA_HALO= %f\n",DELTA_HALO);
      prev_cosmo=RESET_COSMOLOGY;
      n=50;
      if(flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      cvir_pnorm_g1=SIGMA_8/sigmac(8.0);

      flag=0;
      
      dm=(log(HOD.M_max)-log(1.0E8))/(n-1);
      for(i=1;i<=n;++i)
	{
	  // below is one of our models
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i]=cvir_2pop(x[i],0,DELTA_HALO,&m1);
	  //printf("CCC %e %e\n",m1,y[i]);
	  //fflush(stdout);
	  x[i] = log(m1);
	  y[i] = log(y[i]);
	  continue;

	  // below is for the Maccio model
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i]=c200c_maccio(x[i],0,DELTA_HALO,&m1);
	  x[i] = log(m1);
	  y[i] = log(y[i]);
	  continue;


	  //below is biullock cvir model
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i]=cvir_model(x[i]);
	  x[i]=log(halo_mass_conversion(x[i],&y[i],DELTA_HALO));
	  y[i]=log(y[i]);
	  //y[i] = log(10);
	}
      spline(x,y,n,1.0E+30,1.0E+30,y2);
    }
  cfac = 0.13*log10(m/1.0e12) + 0.125;
  cfac = 1;
  if(m>80*MSTAR) m=80*MSTAR;
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(exp(a)*cfac);
  
}


/* This is the cvir model, which should reproduce the results
 * of cvir3.f, which can be obtained off James Bullock's web site.
 *
 * Basic idea; For a halo of mass M, find the collapse redshift
 * of that halo be finding the redshift at which rms mass variance
 * sigma(Mf) = DELTA_CRIT (using linear growth factor), where Mf is some
 * fraction of the redshift zero halo mass.
 *
 * cvir = k*(1+z_coll(M*f))
 *
 * Model params:
 *  - k = 3.12    - fitting parameter (taken from cvir3.f)
 *  - f = 0.001   - fraction of halo mass that collapsed at z_coll
 */

double cvir_model(double mass)
{
  double cvir,z_coll,zbrent(),rad;
  double k=3.12;
  double f=0.001;
 
  rad=pow(3.0*f*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  cvir_sigma_g1=cvir_pnorm_g1*sigmac(rad);

  if(collapse_redshift(0.0)<0)return(k);
  z_coll=zbrent(collapse_redshift,0.0,200.0,1.0E-3);
  cvir=k*(1+z_coll);
  return(cvir);
}

/* This is the input function to find where sig(M)*D(z)-DELTA_CRIT = 0
 * which is the redshift at which a halo of mass M collapsed.
 */
double collapse_redshift(double z)
{
  double D;
  D=growthfactor(z);
  return(cvir_sigma_g1*D-DELTA_CRIT);
}

/* Some quantities are specific to an overdensity of 200 (i.e., the Jenkins mass
 * function and the halo bias relation of Tinker et al. )
 * Specfically, these are actually for FOF 0.2 halos, for which 200 is the
 * current best approximation. (Although Warren et al quotes 250, which is the most recent
 * value.)
 *
 * Therefore, the halo concentrations for Delta=200 halos need to be easily 
 * accessible for halo mass conversion of these quantities. The values are 
 * tabulates here, as opposed to the above routine which tabulates the 
 * concentrations for a user-specified overdensity.
 */

double halo_c200(double m)
{
  static int flag=1,n,prev_cosmo=0;
  static double *x,*y,*y2;
  int i;
  float x1,x2;
  double a,dm,x3,x4,m1;
  FILE *fp;
  char fname[1000];

  if(flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(OUTPUT)
	fprintf(stdout,"Calc cvir with DELTA_HALO= %f\n",200.0);
      prev_cosmo=RESET_COSMOLOGY;
      n=50;
      if(flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      cvir_pnorm_g1=SIGMA_8/sigmac(8.0);

      flag=0;
      
      dm=(log(HOD.M_max)-log(1.0E8))/(n-1);
      for(i=1;i<=n;++i)
	{
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i]=cvir_2pop(x[i],0,200.0,&m1);
	  x[i] = log(m1);
	  y[i] = log(y[i]);
	  continue;

	  //below is bullock
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i]=cvir_model(x[i]);
	  x[i]=log(halo_mass_conversion(x[i],&y[i],200.0));
	  y[i]=log(y[i]);
	}
    }
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(exp(a));
  
}


double func_cvir1(double a)
{
  double c;
  c = a;
  return exp(a)/RT2PI/0.34*exp(-(c - cvir_g1 + 0.34*0.34/2)*(c - cvir_g1 + 0.34*0.34/2)/(2*0.34*0.34));
}
double func_cvir2(double a)
{
  double c;
  c = a;
  return 1.0/RT2PI/0.34*exp(-(c - cvir_g1 + 0.34*0.34/2)*(c - cvir_g1 + 0.34*0.34/2)/(2*0.34*0.34));
}

/* This is a more robust version of andrey's cvir model
 */
double cvir_2pop(double mass, float K, double delta, double *mdelta)
{
  double f = 0.3;
  double rad, z_coll, cvir, cdelta, dc, rdelta, a_coll, sig, x, a_now, c_min, x2,
    mlo, mhi, rlo, rhi, slo, shi, dlogsdlogM;


  K = 4.8; // for RELAXED halos
  K = 4.4; // for ALL halos

  cvir_pnorm_g1=SIGMA_8/sigmac(8.0);
  rad=pow(3.0*f*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  sig=cvir_pnorm_g1*sigmac(rad);

  mlo=0.99*mass;
  mhi=1.01*mass;
  rlo=pow(3.0*mlo/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  rhi=pow(3.0*mhi/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  slo=cvir_pnorm_g1*sigmac(rlo);
  shi=cvir_pnorm_g1*sigmac(rhi);
  dlogsdlogM=log(shi/slo)/log(mhi/mlo);
  
  a_coll = 0.6*pow(1.69/sig,0.75);
  a_coll = 1.0/(sig*sqrt(fabs(dlogsdlogM)/0.185));
  z_coll = 1/a_coll - 1;
  cvir_anow = a_now = 1/(REDSHIFT+1);

  
  cvir_g1 =  log(K*(1+z_coll)/(1+REDSHIFT));
  //c_min = log(K*1.2*pow(1+REDSHIFT,-0.3));
  c_min = log(K);

  cvir = qromo(func_cvir1,c_min,10,midpnt);
  sig = 0.34;
  x =erf(fabs(c_min - cvir_g1 + sig*sig/2)/(2*sig*sig));
  x2 = qromo(func_cvir2,c_min,10,midpnt);
  
  // for relaxed halos
  //cvir = cvir/x2;

  //printf("FRAC %f %e %f\n",REDSHIFT,mass,x2);

  // for all halos
  x=OMEGA_M*pow(1+REDSHIFT,3.0)/(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
  cvir = cvir + (1-x2)*K*1.15*pow(1+REDSHIFT,-0.3);

  //cvir = cvir + (1-x2)*K*1.2*pow(1+REDSHIFT,-0.3*(0.3/OMEGA_M)); // for RELAXED halos

  //cvir = cvir + (1-x2)*K*0.85; // this is for ALL halos

  // standard result
  //cvir = K*(1+z_coll)/(1+REDSHIFT);

  x=x-1;
  dc=(18*PI*PI+82*x-39*x*x)/(1+x);
  *mdelta = halo_mass_conversion2(mass, cvir, dc, delta);
  rad = pow(3.0*mass/(4.0*PI*RHO_CRIT*OMEGA_M*dc),THIRD);
  rdelta = pow(3.0*(*mdelta)/(4.0*PI*delta*OMEGA_M*RHO_CRIT),1.0/3.0);
  cdelta = cvir*rdelta/rad;
  return cdelta;

}

/* In the top part of this function, I'm implementing equation (10) in
 *  Maccio et al, which is c200c for WMAP5 cosmology. I then convert
 *  to the given overdensity and return that
 */
double c200c_maccio(double mass, int iwmap, double delta, double *mdelta)
{
  double K, f = 0.01;
  double rad, z_coll, c200, cdelta, dc, rdelta;

  c200 = 0.830 - 0.098*log10(mass/1.0e12);
  c200 = pow(10.0,c200);

  // get the value of 200crit
  dc = 200/OMEGA_M/pow(1+REDSHIFT,3.0)*(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
  
  // convert the mass to the desired overdensity
  *mdelta = halo_mass_conversion2(mass, c200, dc, delta);

  // scale the concentration to the desired overdensity
  rad = pow(3.0*mass/(4.0*PI*200*RHO_CRIT/pow(1+REDSHIFT,3.0)*(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M))),1.0/3.0);
  rdelta = pow(3.0*(*mdelta)/(4.0*PI*delta*OMEGA_M*RHO_CRIT),1.0/3.0);
  cdelta = c200*rdelta/rad;
  return cdelta;


  switch (iwmap) {
  case 1: K = 4.0; break;
  case 2: K = 3.9; break;
  case 3: K = 3.8; break;
  case -1: K = 3.8; break;
  case -2: K = 3.6; break;
  case -3: K = 3.4; break;
  default: endrun("No model specified for Maccio c200\n");
  }

  cvir_pnorm_g1=SIGMA_8/sigmac(8.0);
  rad=pow(3.0*f*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  cvir_sigma_g1=cvir_pnorm_g1*sigmac(rad);

  //printf("%e %f %f %e\n",mass,rad,cvir_sigma_g1,collapse_redshift(0.0));

  if(collapse_redshift(-0.99)<0)return(K);
  z_coll=zbrent(collapse_redshift,-0.99,200.0,1.0E-3);

  c200 =  K*pow(OMEGA_M*pow(1+z_coll,3.0)+(1-OMEGA_M),THIRD)/pow(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M),THIRD);
  //printf("%e %f %f %f %f %e\n",mass,c200,z_coll,REDSHIFT,OMEGA_M,pow(OMEGA_M*pow(1+z_coll,3.0)+(1-OMEGA_M),3.0) );

  dc = 200/OMEGA_M/pow(1+REDSHIFT,3.0)*(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M));
  
  *mdelta = halo_mass_conversion2(mass, c200, dc, delta);
  rad = pow(3.0*mass/(4.0*PI*200*RHO_CRIT/pow(1+REDSHIFT,3.0)*(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M))),1.0/3.0);
  rdelta = pow(3.0*(*mdelta)/(4.0*PI*delta*OMEGA_M*RHO_CRIT),1.0/3.0);
  cdelta = c200*rdelta/rad;

  // printf("ZC %e %f\n",mass,z_coll);

  return cdelta;
}

