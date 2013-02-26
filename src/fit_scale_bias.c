#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

void initial_scale_bias_parameters(double **pp, double *yy);
double chi2_scale_bias(double *a);
char *BIAS_FILENAME;
double LARGE_SCALE_BIAS;

void fit_scale_bias(int argc, char **argv)
{
  FILE *fp;
  int n,i,j,k,niter;

  double FTOL = 1.e-3;
  double **pp,*yy;

  BIAS_FILENAME = argv[3];
  LARGE_SCALE_BIAS = atof(argv[4]);

  wp.ncf = 4;

  pp=dmatrix(1,wp.ncf+1,1,wp.ncf);
  yy=dvector(1,wp.ncf+1);
  initial_scale_bias_parameters(pp,yy);
  amoeba(pp,yy,wp.ncf,FTOL,chi2_scale_bias,&niter);

}

void initial_scale_bias_parameters(double **pp, double *yy)
{
  int i,j;
  double d[100],a[100];

  a[1] = 1.17;
  a[2] = 1.79;
  a[3] = 1.01;
  a[4] = 2.03;
  
  for(i=1;i<=wp.ncf;++i)
    d[i]=a[i]*0.2;

  for(j=1;j<=wp.ncf;++j)
    pp[1][j]=a[j];
  yy[1]=chi2_scale_bias(a);

  for(i=1;i<=wp.ncf;++i)
    {
      a[i]+=d[i];
      if(i>1)a[i-1]-=d[i-1];
      yy[i+1]=chi2_scale_bias(a);
	
      for(j=1;j<=wp.ncf;++j)
	pp[i+1][j]=a[j];
    }
  a[wp.ncf]-=d[wp.ncf];

}


double chi2_scale_bias(double *a)
{
  static FILE *fp;
  static int n, flag = 1;
  static float *r,*x,*e;

  int i;
  double xi,chi2=0,b;
  char aa[1000];

  if(flag)
    {
      fp = openfile(BIAS_FILENAME);
      n = filesize(fp);

      r = vector(1,n);
      x = vector(1,n);
      e = vector(1,n);

      for(i=1;i<=n;++i) {
	fscanf(fp,"%f %f %f",&r[i],&x[i],&e[i]);
	x[i] = x[i]/LARGE_SCALE_BIAS; //bias_interp(7.70e10,-1);
	e[i] = 0.05*x[i];
	fgets(aa,1000,fp); }
      flag = 0;
    }

  for(i=1;i<=n;++i)
    {
      if(r[i]>20)continue;
      xi = xi_interp(r[i]);
      b = pow(1+a[1]*xi,a[2])/pow(1+a[3]*xi,a[4]);
      chi2 += (b - x[i])*(b - x[i])/e[i]/e[i];
    }
  if(isnan(chi2))chi2 = 1.0E7;
  printf("ITER %e %f %f %f %f\n",chi2,a[1],a[2],a[3],a[4]);
  return(chi2);
}
