#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

double func_ncen(double m)
{
  m=exp(m);
  return(N_cen(m));
}
double func_alpha(double a)
{
  double x;
  HOD.alpha = a;
  x = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
  return(GALAXY_DENSITY - x);
}
double func_m1(double m)
{
  double x;
  HOD.M1=exp(m);
  x = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
  return(GALAXY_DENSITY - x);
}  

int neg_binomial()
{
  int i,j,k;
  double fsat,alpha,r;

  double func_ncen(),func_m1();

  fsat = 1 - qromo(func_ncen,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  
  HOD.M_cut/=pow(1.2589,30.0);

  for(i=1;i<=60;++i)
    {
      HOD.M_cut = HOD.M_cut*1.2589;
      /* alpha = zbrent(func_alpha,0.5,3.0,1.0E-3); */
      alpha = zbrent(func_m1,log(1.0e11),log(1.0e15),1.0E-3);
      printf("NB%d %e %e\n",i,HOD.M_cut,exp(alpha));

      RESET_FLAG_1H = RESET_FLAG_2H = 1;

      for(j=-10;j<=17;++j)
	{
	  r = pow(10.0,j/10.0);
	  printf("XI%d %f %e\n",i,r,one_halo_real_space(r)+two_halo_real_space(r));
	}
    }
  exit(0);
  return 0;
}
