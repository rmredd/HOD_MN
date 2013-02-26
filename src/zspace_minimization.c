#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif


#include "header.h"

/* This is a series of routines to calculate the projected correlation
 * function from the real-space one (HOD calculation) and then chi2
 * minimize the parameters based on the SDSS data.
 */

/* Local functions.
 */
double func_chi2(double *a);
void initial_zspace_values(double *a, double **pp, double *yy);

/* External functions.
 */
void wp_input(void);

  
/********************************************************************
 * Below is the actual minimization of the wp + z-space data.
 *
 * HOD.free[] is a vector which holds 1/0 as to whether or not a parameter is going
 * to be help constant during the chi^2 minimization. 1==vary 0==constant
 *
 *
    1)  if(HOD.free[++i])opar[++j]=(HOD.M_min);
    2)  if(HOD.free[++i])opar[++j]=(HOD.M1);
    3)  if(HOD.free[++i])opar[++j]=HOD.alpha;
    4)  if(HOD.free[++i])opar[++j]=(HOD.M_cut);
    5)  if(HOD.free[++i])opar[++j]=HOD.sigma_logM;
    6)  if(HOD.free[++i])opar[++j]=CVIR_FAC;
    7)  if(HOD.free[++i])opar[++j]=HOD.MaxCen;
    8)  if(HOD.free[++i])opar[++j]=OMEGA_M;
    9)  if(HOD.free[++i])opar[++j]=SIGMA_8;
    10) if(HOD.free[++i])opar[++j]=VBIAS;
    11) if(HOD.free[++i])opar[++j]=VBIAS_C;

 */

void zspace_minimization(char *fname)
{
  int n,niter,i,j;
  double *a,**pp,*yy,FTOL=1.0E-3,chi2min,qromo(),midpnt(),s1;

  fprintf(stderr,"\n\nCHI2 MINIMIZATION OF W_P(R_P)+xi(s,p) DATA..........\n");
  fprintf(stderr,    "-----------------------------------------------------\n\n");

  OUTPUT=0;
  wp_input();
  Work.imodel=2;

  /* Find the number of free parameters in the minimization
   * for the real-space correlation function.
   */
  for(n=0,i=1;i<=100;++i)
    {
      n+=HOD.free[i];
      if(OUTPUT)
	printf("zspace_min> free[%i] = %d\n",i,HOD.free[i]);
    }
  wp.ncf=n;

  a=dvector(1,n);
  if(POWELL)
    pp=dmatrix(1,n,1,n);
  else
    pp=dmatrix(1,n+1,1,n);
  yy=dvector(1,n+1);

  initial_zspace_values(a,pp,yy);

  if(POWELL) 
    {
      powell(a,pp,n,FTOL,&niter,&chi2min,func_chi2);
      chi2min = func_chi2(a);
    }
  else
    amoeba(pp,yy,n,FTOL,func_chi2,&niter);

  s1=qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt);
  GALAXY_BIAS=s1/GALAXY_DENSITY;

  if(!ThisTask) {
    printf("ZSPACEMIN %e %e ",chi2min,HOD.M_min);
    for(i=1;i<=n;++i)printf("%e ",a[i]);
    printf(" %f\n",GALAXY_BIAS);
  }

  output_parameter_file(fname);
    
}

double func_chi2(double *a)
{
  static double *opar;
  static int iter=0,flag=1;
  double x1,x2;
  int i,j;

  if(flag)
    {
      flag=0;
      opar=dvector(1,100);
      i=0;j=0;
      if(HOD.free[++i])opar[++j]=(HOD.M_min);
      if(HOD.free[++i])opar[++j]=(HOD.M1);
      if(HOD.free[++i])opar[++j]=HOD.alpha;
      if(HOD.free[++i])opar[++j]=(HOD.M_cut);
      if(HOD.free[++i])opar[++j]=HOD.sigma_logM;
      if(HOD.free[++i])opar[++j]=CVIR_FAC;
      if(HOD.free[++i])opar[++j]=HOD.MaxCen;
      if(HOD.free[++i])opar[++j]=OMEGA_M;
      if(HOD.free[++i])opar[++j]=SIGMA_8;
      if(HOD.free[++i])opar[++j]=VBIAS;
      if(HOD.free[++i])opar[++j]=VBIAS_C;
    }
  RESET_COSMOLOGY++;
  j=0;
  for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
  i=N_HOD_PARAMS;
  if(HOD.free[++i])OMEGA_M = a[++j];
  if(HOD.free[++i])SIGMA_8 = a[++j];
  if(HOD.free[++i])VBIAS   = a[++j];
  if(HOD.free[++i])VBIAS_C = a[++j];
  if(HOD.free[++i])GAMMA   = a[++j];
  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];

 
  if(OMEGA_M<=0 || SIGMA_8<=0 || VBIAS<=0)return(1.0e7);

  x1=chi2_wp(a);
  if(x1>9.0E6)
    {
	printf("ETRAP: %d chi2 for wp = %e\n",iter+1,x1);
      return(x1);
    }
  x2=chi2_zspace(a);

  if(!ThisTask) {
    printf("TRY %d ",++iter);
    for(i=1;i<=wp.ncf;++i)
      printf("%.4e ",a[i]);
    printf("%e %e %e\n",x1,x2,x1+x2);
    fflush(stdout);
  }

  return(x1+x2);
}


void initial_zspace_values(double *a, double **pp, double *yy)
{
  int i,j=0,i1,j1;
  double x1,x2,d[100],start_dev[100];
  long IDUM = -556;

  i=0;j=0;
  if(HOD.free[++i]){ a[++j]=(HOD.M_min);start_dev[j]=0.001; }
  if(HOD.free[++i]){ a[++j]=(HOD.M1);start_dev[j]=0.01; } //.0005
  if(HOD.free[++i]){ a[++j]=HOD.alpha;start_dev[j]=0.03; } //.005
  if(HOD.free[++i]){ a[++j]=(HOD.M_cut);start_dev[j]=0.01; } //.001
  if(HOD.free[++i]){ a[++j]=(HOD.sigma_logM);start_dev[j]=0.01; }
  if(HOD.free[++i]){ a[++j]=CVIR_FAC;start_dev[j]=0.02; }
  if(HOD.pdfc==7) {
    if(HOD.free[++i])a[++j]=(HOD.M_cen_max); start_dev[j]=0.001; }
  else {
    if(HOD.free[++i])a[++j]=HOD.MaxCen; start_dev[j]=0.02; }
  if(HOD.free[++i]){ a[++j]=(HOD.M_sat_break);start_dev[j]=0.001; }
  if(HOD.free[++i]){ a[++j]=HOD.alpha1;start_dev[j]=0.02; }

      if(HOD.free[++i])a[++j]=OMEGA_M;
      if(HOD.free[++i])a[++j]=SIGMA_8;
      if(HOD.free[++i])a[++j]=VBIAS;
      if(HOD.free[++i])a[++j]=VBIAS_C;
      if(HOD.free[++i])a[++j]=GAMMA;
      if(HOD.free[++i])a[++j]=SPECTRAL_INDX;

  /*
  i=0;j=0;
  if(HOD.free[++i])a[++j]=HOD.M_min;
  if(HOD.free[++i])a[++j]=HOD.M1;
  if(HOD.free[++i])a[++j]=HOD.alpha;
  if(HOD.free[++i])a[++j]=HOD.M_cut;
  if(HOD.free[++i])a[++j]=HOD.sigma_logM;
  if(HOD.free[++i])a[++j]=CVIR_FAC;
  if(HOD.free[++i])a[++j]=HOD.MaxCen;
  if(HOD.free[++i])a[++j]=OMEGA_M;
  if(HOD.free[++i])a[++j]=SIGMA_8;
  if(HOD.free[++i])a[++j]=VBIAS;
  if(HOD.free[++i])a[++j]=VBIAS_C;
  */

  if(!ThisTask){
    printf("INITIAL VALUES: ");
    for(i=1;i<=wp.ncf;++i)printf("%e ",a[i]);
    printf("\n");
  }
  for(i=1;i<=wp.ncf;++i)
    if(a[i]!=0)
      d[i]=a[i]*0.10;
    else
      d[i]=0.2;

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
      yy[1]=func_chi2(a);
    
      for(i=1;i<=wp.ncf;++i)
	{
	  a[i]+=d[i];
	  if(i>1)a[i-1]-=d[i-1];
	  yy[i+1]=func_chi2(a);	  
	  for(j=1;j<=wp.ncf;++j)
	    pp[i+1][j]=a[j];
	}
      a[wp.ncf]-=d[wp.ncf];
    }
}
