#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

int g31_number;

double poisson_prob(int n, double nave);

double func1_m2n(double m);
double func2_m2n(double m);

void mass2number()
{
  int i,j,k,n;
  double s1,s2,mbar;
  char fname[1000];
  FILE *fp;

  sprintf(fname,"%s.M2N",Task.root_filename);
  fp = fopen(fname,"w");

  for(i=2;i<=10;++i)
    {
      g31_number = i;
      s1 = qromo(func1_m2n,log(HOD.M_low),log(HOD.M_max),midpnt);
      s2 = qromo(func2_m2n,log(HOD.M_low),log(HOD.M_max),midpnt);
      mbar = s1/s2;
      fprintf(fp,"%d %e %e\n",i,mbar,mbar/i);
    }

  for(i=20;i<=150;i+=10)
    {
      g31_number = i;
      s1 = qromo(func1_m2n,log(HOD.M_low),log(HOD.M_max),midpnt);
      s2 = qromo(func2_m2n,log(HOD.M_low),log(HOD.M_max),midpnt);
      mbar = s1/s2;
      fprintf(fp,"%d %e %e\n",i,mbar,mbar/i);
    }
  fclose(fp);
}

double func1_m2n(double m)
{
  double ncen,nsat,x;
  m = exp(m);
  ncen = N_cen(m);
  nsat = N_sat(m);
  x = poisson_prob(g31_number-1,nsat);
  //  printf("%e %e %e %e %e\n",m,ncen,nsat,poisson_prob(g31_number-1,nsat));
  if(isnan(x))x = 0;
  return ncen*x*m*m*dndM_interp(m);
}
double func2_m2n(double m)
{
  double ncen,nsat,x;
  m = exp(m);
  ncen = N_cen(m);
  nsat = N_sat(m);
  x = poisson_prob(g31_number-1,nsat);
  if(isnan(x))x = 0;
  return ncen*x*m*dndM_interp(m);
}
