
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

double nbody_xi_from_file(double r);

void vpf()
{
  int i;
  double vpf, nbar, xibar, r, func_xibar();

  for(i=1;i<=13;++i)
    {
      r=i/1.0;
      nbar = 4./3.*PI*r*r*r*GALAXY_DENSITY;
      xibar = qromo(func_xibar,log(0.01),log(r),midpnt)*3.0/(r*r*r);
      /* xibar = qtrap(func_xibar,log(0.01),log(r),1.0E-5)*3/(r*r*r); */
      vpf = pow(1+nbar*xibar,-1/xibar);
      printf("%.1f %e %e %e %e\n",r,vpf,nbar,xibar,
	     one_halo_real_space(r)+two_halo_real_space(r));
    }
  exit(0);
}

double func_xibar(double r)
{
  r=exp(r);
  return(1/r/r*r*r*r);
  return(r*r*r*(one_halo_real_space(r)+two_halo_real_space(r)));
  return(nbody_xi_from_file(r)*r*r*r);

  r=exp(r);
  return(3*r*r*r*(one_halo_real_space(r)+two_halo_real_space(r)));
  return(3*r*r*r*(nbody_xi_from_file(r)));
}

double nbody_xi_from_file(double r)
{
  static FILE *fp;
  static int n=0;
  static double *rr,*xx,*yy;
  char aa[1000];
  float x1,x2,x3,x4;
  double a;
  int i;

  if(!n)
    {
      fp=openfile("xi.B7");
      n = filesize(fp) - 1;
      
      rr = dvector(1,n);
      xx = dvector(1,n);
      yy = dvector(1,n);
      
      fgets(aa,1000,fp);
      
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f %f %f",&x1,&x2,&x3,&x4);
	  fgets(aa,1000,fp);
	  rr[i] = x3;
	  xx[i] = x4;
	}
      spline(rr,xx,n,1.0E+30,1.0E+30,yy);
    }
  splint(rr,xx,yy,n,r,&a);
  return(a);

}
