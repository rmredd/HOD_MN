#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header.h"

double NFW_func(double x);
double NFW_density3(double r);
double NFW_galaxy_density(double r);
double NFW_density2(double r);
double NFW_galaxy_density2(double r);
double I_integral(double t);
double I_integral2(double t);
double mean_NFW_dispersion(double m);
double satellite_dispersion(double r, double rs, double rvir);

double RS;

double velocity_dispersion_profile(double sig, double cvir, double rvir, double r)
{
  double sigr,rs,mass,s_nfw,sigx,c_dm,x,Navg,x1;
  double xmin,xmax,dlogx,fac,a,dlogr,rho0,sum;
  int i,j,k,n2=1000;

  I_integral(1.0);
  I_integral2(1.0);
  
  
  c_dm=cvir/CVIR_FAC;
  mass = 4./3.*PI*rvir*rvir*rvir*RHO_CRIT*OMEGA_M*DELTA_HALO;
  s_nfw = mean_NFW_dispersion(mass);
  RS = rs = rvir/c_dm;
  x = r/rs*CVIR_FAC;
  sigr = sig/s_nfw*ROOT2*sig*sqrt(c_dm/NFW_func(c_dm)/CVIR_FAC*(x)*(1+x)*(1+x)*
				  I_integral2(x));	  
  /* return(sigr); */
  
  CVIR_FAC=1;
  fac=sqrt(4.499E-48/2)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  dlogr = (log(1.0e16)-log(1.0e12))/(99);
  mass = 1.0e14;
  for(i=0;i<100;++i)
    {
      /* mass=exp(dlogr*i)*1.0e12; */
      CVIR_FAC = (i+20)/100.0;
      sigx=mean_NFW_dispersion((double)mass);
      sigx=fac*pow(mass,1.0/3.0);
      c_dm = halo_concentration(mass);
      rvir = pow(3*mass/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
      RS = rs = rvir/c_dm;
      sigr=0;
      Navg = N_avg(mass);
      rho0 = N_avg(mass)/qromo(NFW_galaxy_density2,1.0E-6,rvir,midpnt);
      for(j=1;j<=100;++j)
	{
	  r=j/100.0*rvir;
	  x = r/rs;
	  sigr += fac*pow(mass,1.0/3.0)*ROOT2*sqrt(c_dm/NFW_func(c_dm)/CVIR_FAC*(x)*(1+x)*(1+x)*
						   I_integral2((double)(x)))*
	    4*PI*r*r*rvir/100*NFW_galaxy_density(r)*rho0/Navg;	  
	}
      printf("%f %e %e %e %f\n",log10(mass),sigr,sigx,sigr/sigx,CVIR_FAC);
    }
  exit(0);

  /*
  s_nfw = mean_NFW_dispersion((double)mass);
  sigr = sig/s_nfw*ROOT2*sig*
    sqrt(cvir/NFW_func(cvir)*(r/rs)*(1+r/rs)*(1+r/rs)*I_integral((double)(r/rs)));

  sigr = sig/s_nfw*ROOT2*sig*
    sqrt(rvir/NFW_func(cvir)*satellite_dispersion(r,rs,rvir));
  return(sigr);
  */

  mass=1.0E14;
  cvir=halo_concentration((double)mass)*CVIR_FAC;
  c_dm=cvir/CVIR_FAC;
  rvir=pow(3.0*mass/(4*PI*DELTA_HALO*OMEGA_M*RHO_CRIT),1.0/3.0);
  fac=sqrt(4.499E-48/2)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sig=fac*pow(mass,1.0/3.0);


  RS = rs = rvir/c_dm;
  Navg = N_avg(mass);
  rho0 = N_avg(mass)/qromo(NFW_galaxy_density2,1.0E-6,rvir,midpnt);
  fprintf(stderr,"%f %e %f\n",N_avg(mass),rho0,rs);

  dlogr=(log(c_dm)-log(0.01))/n2;
  if(dlogr>1)exit(0);
  sum=0;
  x1=0;
  s_nfw = mean_NFW_dispersion((double)mass);
  
  for(j=1;j<=n2;++j)
    {
      r=exp((j-0.5)*dlogr)*0.01*RS;
      x=r/rs*CVIR_FAC;

      sigr = sig/s_nfw*ROOT2*sig*
	sqrt(rvir/NFW_func(c_dm)*satellite_dispersion(r,rs,rvir));
      sigx = sig/s_nfw*ROOT2*sig*sqrt(cvir/NFW_func(cvir)*(r/rs)*(1+r/rs)*(1+r/rs)*
				      I_integral((double)(r/rs)));	  
      sigr = sig/s_nfw*ROOT2*sig*sqrt(c_dm/NFW_func(c_dm)/CVIR_FAC*(x)*(1+x)*(1+x)*
				      I_integral2((double)(x)));	  
      
      sum += NFW_galaxy_density(r)*rho0*sigr*r*r*r*dlogr;
      x1 += NFW_galaxy_density(r)*rho0*r*r*r*dlogr;
      printf("%f %f %f %f %f %f %f %e %f\n",r,sigr,sigx,sum*4*PI/Navg,sig,s_nfw,x1*4*PI/Navg,
	     NFW_galaxy_density(r)*rho0,RS);
    }
  exit(0);
  
  return(sigr);
}

double satellite_dispersion(double r, double rs, double rvir)
{
  double x1,x2,x;
  double func_jeans_sat();

  x=CVIR_FAC*r/rs;
  /*x1=x*(1+x)*(1+x)*qromo(func_jeans_sat,(double)(r),(double)(rvir),midpnt);*/
  x1=x*(1+x)*(1+x)*qromo(func_jeans_sat,(double)(r),1.0E+30,midinf);
  return(x1);
}

double func_jeans_sat(double r)
{
  double x;
  x=r/RS;
  return(NFW_func(x)/(r*r*CVIR_FAC*x*(1+CVIR_FAC*x)*(1+CVIR_FAC*x)));
}

double mean_NFW_dispersion(double m)
{
  static int flag=0,n=100;
  static double *x,*y,*y2;
  double xmin,xmax,dlogx,fac,a,r,rs,dlogr,cvir,rvir,rho0,sigr,sum,mass,sig;
  int i,j,k,n2=1000;

  if(flag)goto SKIP_INIT;
  flag=1;

  xmin=log(1.0e9);
  xmax=log(1.0e16);
  dlogx=(xmax-xmin)/(n-1);

  x=dvector(1,n);
  y=dvector(1,n);
  y2=dvector(1,n);

  fac=sqrt(4.499E-48/2)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;

  for(i=1;i<=n;++i)
    {
      x[i]=(i-1)*dlogx+xmin;
      mass=exp((i-1)*dlogx+xmin);

      cvir=halo_concentration((double)mass);
      rvir=pow(3.0*mass/(4*PI*DELTA_HALO*OMEGA_M*RHO_CRIT),1.0/3.0);
      sig=fac*pow(mass,1.0/3.0);

      RS = rs = rvir/cvir;
      rho0 = mass/qromo(NFW_density2,1.0E-6,rvir,midpnt);

      dlogr=(log(cvir)-log(0.01))/n2;
      if(dlogr>1)exit(0);
      sum=0;
      
      for(j=1;j<=n2;++j)
	{
	  r=exp((j-0.5)*dlogr)*0.01*RS;
	  sigr = ROOT2*sig*sqrt(cvir/NFW_func(cvir)*(r/rs)*(1+r/rs)*(1+r/rs)*
				I_integral((double)(r/rs)));	  
	  sum += NFW_density3(r)*rho0*sigr*r*r*r*dlogr;
	}
      y[i]=sum*4*PI/mass;
    }
  spline(x,y,n,1.0E+30,1.0E+30,y2);
 SKIP_INIT:

  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(a);
}

double I_integral(double t)
{
  static int flag=0,n=1000;
  static double *x,*y,*y2;
  double xmin,xmax,dlogx,func_I(),a;
  int i;

  if(flag)goto SKIP_INIT;

  flag=1;
  
  xmin=-10.0;
  xmax=10.0;
  dlogx=(xmax-xmin)/(n-1);

  x=dvector(1,n);
  y=dvector(1,n);
  y2=dvector(1,n);

  for(i=1;i<=n;++i)
    {
      x[i]=exp((i-1)*dlogx+xmin);
      y[i]=qromo(func_I,x[i],1.0E+30,midinf);
    }
  spline(x,y,n,1.0E+30,1.0E+30,y2);

 SKIP_INIT:

  splint(x,y,y2,n,t,&a);
  return(a);
}

double func_I(double t)
{
  double f;
  f=log(1+t)-t/(1+t);
  return(f/(t*t*t*(1+t)*(1+t)));
}

double I_integral2(double t)
{
  static int flag=0,n=1000;
  static double *x,*y,*y2;
  double xmin,xmax,dlogx,func_I2(),a;
  int i;

  if(flag)goto SKIP_INIT;

  flag=1;
  
  xmin=-10.0;
  xmax=10.0;
  dlogx=(xmax-xmin)/(n-1);

  x=dvector(1,n);
  y=dvector(1,n);
  y2=dvector(1,n);

  for(i=1;i<=n;++i)
    {
      x[i]=exp((i-1)*dlogx+xmin);
      y[i]=qromo(func_I,x[i],1.0E+30,midinf);
    }
  spline(x,y,n,1.0E+30,1.0E+30,y2);

 SKIP_INIT:

  splint(x,y,y2,n,t,&a);
  return(a);
}

double func_I2(double t)
{
  double f,t1;
  t1=t/CVIR_FAC;
  f=log(1+t1)-t1/(1+t1);
  return(f/(t*t*t*(1+t)*(1+t)));
}

double NFW_func(double x)
{
  return(log(1+x)-x/(1+x));
}

double NFW_density3(double r)
{
  return(1.0/(r/RS*(1+r/RS)*(1+r/RS)));
}

double NFW_galaxy_density(double r)
{
  return(1.0/(r/RS*CVIR_FAC*(1+r/RS*CVIR_FAC)*(1+r/RS*CVIR_FAC)));
}

double NFW_density2(double r)
{
  return(4*PI*r*r/(r/RS*(1+r/RS)*(1+r/RS)));
}

double NFW_galaxy_density2(double r)
{
  return(4*PI*r*r/(r/RS*CVIR_FAC*(1+r/RS*CVIR_FAC)*(1+r/RS*CVIR_FAC)));
}

