#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* This file contains routines to fit a red/blue color-selected sample.
 * It does a simultaneous fit to the red, blue, and all samples. It basically
 * acts as a wrapper around the normal wp_minimization routines.
 */
void initial_color_values(double *a, double **pp, double *yy);
void wp_input(void);
void wp_color_input(void);
double chi2_wp_color(double *a);
double dd_func_red_fraction(double x);

/* external functions
 */
void mcmc_color_minimization(void);

void fit_color_samples()
{
  int n,niter,i,j;
  double *a,**pp,*yy,FTOL=1.0E-3,chi2min,s1,dlogm,m;
  FILE *fp;
  char aa[1000];

  mcmc_color_minimization();

  fprintf(stderr,"\n\nCHI2 MINIMIZATION OF W_P(R_P) COLOR  DATA..........\n");
  fprintf(stderr,    "--------------------------------------------\n\n");

  HOD.blue_fraction = 0.5857; /* <- Millenium fraction (sloan);; SDSS fraction -> 0.565; */
  HOD.blue_fraction = 0.6555; /* <- Millenium fraction (B-V>0.8) */
  HOD.blue_fraction = 0.565; /* SDSS -19,-20 */
  HOD.blue_fraction = 0.492; /* SDSS -20,-21 */
  HOD.blue_fraction = 0.379; /* SDSS -21 */

  wp_color.ON = 1;

  if(POWELL)
    FTOL=1.0E-3;
  else
    FTOL=1.0E-5;

  for(n=0,i=1;i<=7;++i)
    {
      n+=HOD.free[i];
      if(!OUTPUT)continue;
      printf("wp_min> free[%i] = %d\n",i,HOD.free[i]);
    }
  /* The parameters that govern the blue fraction aren't 
   * listed in the HOD.free array, so add them in.
   * NB: There are four parameters for these two functions,
   * one of which is fit by the number densities.
   */
  n+=3;

  if(OUTPUT)printf("wp_min> Number of free parameters: %d\n",n);

  wp_color_input();

  wp.ncf=n;
  a=dvector(1,n);
  if(POWELL)
    pp=dmatrix(1,n,1,n);
  else
    pp=dmatrix(1,n+1,1,n);
  yy=dvector(1,n+1);

  initial_color_values(a,pp,yy);

  if(POWELL) 
    {
      if(OUTPUT)printf("wp_min> starting powell.\n");
      powell(a,pp,n,FTOL,&niter,&chi2min,chi2_wp_color);
      chi2min = chi2_wp_color(a);
    }
  else
    {
      if(OUTPUT)printf("wp_min> starting amoeba.\n");
      amoeba(pp,yy,n,FTOL,chi2_wp_color,&niter);
      for(i=1;i<=n;++i)a[i]=pp[1][i];
      chi2min = chi2_wp_color(a);
    }	

  s1=qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt);
  GALAXY_BIAS=s1/GALAXY_DENSITY;

  printf("POWELL %e %e %e ",chi2min,HOD.M_min,HOD.fblue0_cen);
  if(HOD.pdfc==7)
    printf("%e ",HOD.M_cen_max);
  for(i=1;i<=n;++i)printf("%e ",a[i]);
  printf(" %f\n",GALAXY_BIAS);

  /* Output the fit and the HOD curve.
   */
  sprintf(aa,"%s.fit",Task.root_filename);
  fp=fopen(aa,"w");
  fprintf(fp,"%e %e %e ",chi2min,HOD.M_min,HOD.fblue0_cen);
  if(HOD.pdfc==7)
    fprintf(fp,"%e ",HOD.M_cen_max);
  for(i=1;i<=n;++i)fprintf(fp,"%e ",a[i]);
  fprintf(fp," %f\n",GALAXY_BIAS);
  fclose(fp);

  sprintf(aa,"%s.HOD",Task.root_filename);
  fp=fopen(aa,"w");
  dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
  for(i=1;i<=100;++i)
    {
      m=exp((i-1)*dlogm)*HOD.M_low;
      fprintf(fp,"%e %e %e %e\n",m,N_cen(m),N_sat(m),N_avg(m));
    }
  fclose(fp);

}

double chi2_wp_color(double *a)
{
  double chi2,nsat_blue,ncen_blue,dlogm,m,func_blue_fraction(),temp1,temp2,fsat_blue,fsat_red,fsat_all;
  int i,j;
  static int iter=0;

  COVAR = 0;

  HOD.color = 0;
  GALAXY_DENSITY2 = GALAXY_DENSITY = wp_color.ngal_full;
  wp.np = wp_color.n_full;

  for(i=wp.ncf-2;i<=wp.ncf;++i)
    if(a[i]<0)return(1.0e+7);

  for(i=1;i<=wp.np;++i)
    {
      wp.r[i] = wp_color.r_full[i];
      wp.x[i] = wp_color.x_full[i];
      wp.e[i] = wp_color.e_full[i];
      for(j=1;j<=wp.np;++j)
	  wp.covar[i][j] = wp_color.covar_full[i][j];
    }

  chi2 = chi2_wp(a);

  HOD.M_low0 = set_low_mass();

  fsat_all = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;

  /*
  dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
  for(i=1;i<=100;++i)
    {
      m=exp((i-1)*dlogm)*HOD.M_low;
      printf("HODFULL %e %f %f %f\n",m,N_cen(m)+N_sat(m),N_cen(m),N_sat(m));
    }
  printf("GALDEN FULL %e\n",qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt));
  sprintf(Task.root_filename,"gal_full");
  populate_simulation();
  */


  if(!iter)
    for(i=1;i<=wp.np;++i)
      for(j=1;j<=wp.np;++j)
	wp_color.covar_full[i][j] = wp.covar[i][j];

  HOD.fblue0_sat = a[wp.ncf - 2];
  HOD.sigma_fblue_sat = a[wp.ncf - 1];
  HOD.sigma_fblue_cen = a[wp.ncf - 0];
  HOD.fblue0_cen = 1.0;
  
  HOD.color = 1;
  HOD.fblue0_cen = pow(10.0,zbrent(func_blue_fraction,-5.0,1.0,1.0E-5));

  if(OUTPUT)
    fprintf(stdout,"old fblue0_cen= %e %e\n",HOD.fblue0_cen,HOD.M_min);

  if(DENSITY_DEPENDENCE) 
    {
      HOD.color = 2;
      temp1 = HOD.M_min;
      temp2 = HOD.M_low;
      populate_simulation();
      HOD.M_min = temp1;
      HOD.M_low = temp2;
      HOD.fblue0_cen = pow(10.0,zbrent(dd_func_red_fraction,-5.0,1.0,1.0E-5));
      HOD.color = 1;
      if(OUTPUT)
	fprintf(stdout,"new fblue0_cen= %e %e\n",HOD.fblue0_cen,
		(1. - HOD.M_min_fac*(1.0 - HOD.fblue0_cen)));
      //      populate_simulation();
      //exit(0);
    }
  else
    HOD.fblue0_cen = pow(10.0,zbrent(func_blue_fraction,-5.0,1.0,1.0E-5));

  if(ERROR_FLAG)
    {
      ERROR_FLAG=0;
      return(1.0e7);
    }


  /* This is only if you have square functions for the central occupation
   */
  /*
  nsat_blue = qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt);
  ncen_blue = qromo(func_central_density,log(HOD.M_low),log(HOD.M_cen_max),midpnt);
  HOD.fblue0_cen = (wp_color.ngal_blue - nsat_blue)/ncen_blue;
  */

  GALAXY_DENSITY2 = GALAXY_DENSITY = wp_color.ngal_blue;
  wp.np = wp_color.n_blue;

  for(i=1;i<=wp.np;++i)
    {
      wp.r[i] = wp_color.r_blue[i];
      wp.x[i] = wp_color.x_blue[i];
      wp.e[i] = wp_color.e_blue[i];
      for(j=1;j<=wp.np;++j)
	  wp.covar[i][j] = wp_color.covar_blue[i][j];
    }

  chi2 += chi2_wp(a);
  fsat_blue = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;

  /*
  dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
  for(i=1;i<=100;++i)
    {
      m=exp((i-1)*dlogm)*HOD.M_low;
      printf("HODBLUE %e %f %f %f\n",m,N_cen(m)+N_sat(m),N_cen(m),N_sat(m));
    }
  printf("GALDEN BLUE %e\n",qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt));

  sprintf(Task.root_filename,"gal_blue");
  populate_simulation();
  */

  if(!iter)
    for(i=1;i<=wp.np;++i)
      for(j=1;j<=wp.np;++j)
	wp_color.covar_blue[i][j] = wp.covar[i][j];

  HOD.color = 2;
  GALAXY_DENSITY2 = GALAXY_DENSITY = wp_color.ngal_red;
  wp.np = wp_color.n_red;

  for(i=1;i<=wp.np;++i)
    {
      wp.r[i] = wp_color.r_red[i];
      wp.x[i] = wp_color.x_red[i];
      wp.e[i] = wp_color.e_red[i];
      for(j=1;j<=wp.np;++j)
	  wp.covar[i][j] = wp_color.covar_red[i][j];
    }

  chi2 += chi2_wp(a);
  fsat_red = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;

  /*
  dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
  for(i=1;i<=100;++i)
    {
      m=exp((i-1)*dlogm)*HOD.M_low;
      printf("HODRED %e %f %f %f\n",m,N_cen(m)+N_sat(m),N_cen(m),N_sat(m));
    }
  printf("GALDEN RED %e\n",qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt));

  sprintf(Task.root_filename,"gal_red");
  populate_simulation();

  exit(0);
  */
  if(!iter)
    for(i=1;i<=wp.np;++i)
      for(j=1;j<=wp.np;++j)
	wp_color.covar_red[i][j] = wp.covar[i][j];
  
  iter++;

  wp.fsat_all = fsat_all;
  wp.fsat_red = fsat_red;
  wp.fsat_blue = fsat_blue;

  printf("COLOR_PARAMS %d %e %e %e %e %e %e %e %e\n",iter,chi2,HOD.M_min,HOD.M1,HOD.alpha,HOD.fblue0_cen,HOD.fblue0_sat,HOD.sigma_fblue_cen,HOD.sigma_fblue_sat);
  printf("COLOR_ITER %d %e %e %f %f %f\n",iter,chi2,HOD.fblue0_cen,fsat_all,fsat_blue,fsat_red);
  fflush(stdout);
  return(chi2);

}

double func_blue_fraction(double x)
{
  double n;
  HOD.fblue0_cen=pow(10.0,x);
  n = qtrap(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),1.0E-4);
  //printf("%f %e %e %e %e\n",HOD.fblue0_cen,n,wp_color.ngal_blue,central_blue_fraction(HOD.M_min),HOD.M_low); 
  return n - wp_color.ngal_blue;
}

void initial_color_values(double *a, double **pp, double *yy)
{
  static int flag=0;
  int i,j;
  double d[100];

  COVAR=0;
  
  if(!flag) {
    i=0;j=0;
    if(HOD.free[++i])a[++j]=HOD.M_min;
    if(HOD.free[++i])a[++j]=HOD.M1;
    if(HOD.free[++i])a[++j]=HOD.alpha;
    if(HOD.free[++i])a[++j]=HOD.M_cut;
    if(HOD.free[++i])a[++j]=HOD.sigma_logM;
    if(HOD.free[++i])a[++j]=CVIR_FAC;
    if(HOD.pdfc>=7){
      if(HOD.free[++i])a[++j]=HOD.M_cen_max; }
    else {
      if(HOD.free[++i])a[++j]=HOD.MaxCen; }

    a[++j] = HOD.fblue0_sat;
    a[++j] = HOD.sigma_fblue_sat;
    a[++j] = HOD.sigma_fblue_cen;

    printf("INITIAL VALUES: ");
    for(i=1;i<=wp.ncf;++i)printf("%e ",a[i]);
    printf("\n");
  }
  flag++;

  /* Make the starting stepsize 10% of the initial values.
   */
  for(i=1;i<=wp.ncf;++i)
    d[i]=a[i]*0.25/flag;


  if(POWELL)
    {
      for(i=1;i<=wp.ncf;++i)
	{
	  for(j=1;j<=wp.ncf;++j)
	    {
	      pp[i][j]=0;
	      if(i==j)pp[i][j]+=d[j];
	    }
	}
    }
  else
    {
      for(j=1;j<=wp.ncf;++j)
	pp[1][j]=a[j];
      yy[1]=chi2_wp_color(a);
    
      for(i=1;i<=wp.ncf;++i)
	{
	  a[i]+=d[i];
	  if(i>1)a[i-1]-=d[i-1];
	  for(j=1;j<=wp.ncf;++j)
	    yy[i+1]=chi2_wp_color(a);	  
	  pp[i+1][j]=a[j];
	}
      a[wp.ncf]-=d[wp.ncf];
    }
}

void wp_color_input()
{
  float x1,x2,x3;
  FILE *fp;
  int i,j,n=11;
  char a[1000],b[1000];


  /* The first call to wp_input() will be the standard full sample.
   */
  wp_input();

  /* Put the data into the wp_color arrays
   */
  wp_color.r_full = dvector(1,wp.np);
  wp_color.x_full = dvector(1,wp.np);
  wp_color.e_full = dvector(1,wp.np);
  wp_color.covar_full = dmatrix(1,wp.np,1,wp.np);
  wp_color.n_full = wp.np;
  if(wp.np>n)wp_color.n_full = n;
  wp_color.ngal_full = GALAXY_DENSITY;

  for(i=1;i<=wp.np;++i)
    {
      wp_color.r_full[i] = wp.r[i];
      wp_color.x_full[i] = wp.x[i];
      wp_color.e_full[i] = wp.e[i];
      for(j=1;j<=wp.np;++j)
	wp_color.covar_full[i][j] = wp.covar[i][j];
    }
  free_dvector(wp.r,1,wp.np);
  free_dvector(wp.x,1,wp.np);
  free_dvector(wp.e,1,wp.np);
  free_dmatrix(wp.covar,1,wp.np,1,wp.np);

  /* The second call is the blue sample.
   */
  sprintf(a,"%s",wp.fname_wp);
  sprintf(b,"%s",wp.fname_covar);

  sprintf(wp.fname_wp,"%s_blue",a);
  sprintf(wp.fname_covar,"%s_blue",b);
  wp_input();

  wp_color.r_blue = dvector(1,wp.np);
  wp_color.x_blue = dvector(1,wp.np);
  wp_color.e_blue = dvector(1,wp.np);
  wp_color.covar_blue = dmatrix(1,wp.np,1,wp.np);
  wp_color.n_blue = wp.np;
  if(wp.np>n)wp_color.n_blue = n;
  wp_color.ngal_blue = GALAXY_DENSITY*HOD.blue_fraction;

  for(i=1;i<=wp.np;++i)
    {
      wp_color.r_blue[i] = wp.r[i];
      wp_color.x_blue[i] = wp.x[i];
      wp_color.e_blue[i] = wp.e[i];
      for(j=1;j<=wp.np;++j)
	wp_color.covar_blue[i][j] = wp.covar[i][j];
    }
  free_dvector(wp.r,1,wp.np);
  free_dvector(wp.x,1,wp.np);
  free_dvector(wp.e,1,wp.np);
  free_dmatrix(wp.covar,1,wp.np,1,wp.np);

  /* The third call is the red sample.
   */
  sprintf(wp.fname_wp,"%s_red",a);
  sprintf(wp.fname_covar,"%s_red",b);
  wp_input();

  wp_color.r_red = dvector(1,wp.np);
  wp_color.x_red = dvector(1,wp.np);
  wp_color.e_red = dvector(1,wp.np);
  wp_color.covar_red = dmatrix(1,wp.np,1,wp.np);
  wp_color.n_red = wp.np;
  if(wp.np>n)wp_color.n_red = n;
  wp_color.ngal_red = GALAXY_DENSITY*(1-HOD.blue_fraction);

  for(i=1;i<=wp.np;++i)
    {
      wp_color.r_red[i] = wp.r[i];
      wp_color.x_red[i] = wp.x[i];
      wp_color.e_red[i] = wp.e[i];
      for(j=1;j<=wp.np;++j)
	wp_color.covar_red[i][j] = wp.covar[i][j];
    }
  free_dvector(wp.r,1,wp.np);
  free_dvector(wp.x,1,wp.np);
  free_dvector(wp.e,1,wp.np);
  free_dmatrix(wp.covar,1,wp.np,1,wp.np);

  /* Open permanent arrays that are larger than any individual array in
   * the color sequence.
   */
  wp.r = dvector(1,15);
  wp.x = dvector(1,15);
  wp.e = dvector(1,15);
  wp.covar = dmatrix(1,15,1,15);
}

/*
double N_cen_blue(double m)
{
  if(m<HOD.M_low)return(0);
  x = (log10(m) - log10(HOD.M_min_blue))/HOD.sigma_logM_blue;
  return(HOD.MaxCen_blue*exp(-x*x/2));    
}
double N_cen_red(double m)
{
  if(m<HOD.M_low)return(0);
  x = (log10(m) - log10(HOD.M_min_red))/HOD.sigma_logM_red;
  return(HOD.MaxCen_red*exp(-x*x/2));    
}
*/
