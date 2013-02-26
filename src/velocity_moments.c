#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header.h"

double halo_velocity_moments(double r, double m1, double m2, double *vrad, 
			  double *vtan, double dlogr);

void output_velocity_moments(int imass)
{
  double r,dlogr,rmin=0.36,rmax=45.0,vrad[5],vtan[5],m,npair;
  int i,n=20;

  m = RHO_CRIT*pow(RESOLUTION,3.0)*OMEGA_M*pow(2.0,imass+0.5)*100;
  dlogr = (log(rmax)-log(rmin))/(n-1);
  for(i=0;i<n;++i)
    {
      r = exp(i*dlogr)*rmin;
      npair = halo_velocity_moments(r,m,m,vrad,vtan,dlogr);
      if(!isnan(vrad[4]) && !isnan(vtan[4]))
	printf("%f %f %f %f %f %f %f %e\n",r,vrad[1],vrad[2],vrad[3],vrad[4],vtan[2],vtan[4],npair);
    }
  exit(0);
}

double halo_velocity_moments(double r, double m1, double m2, double *vrad, 
			  double *vtan, double dlogr)
{
  static int flag=0;
  static double *prob;
  double binsize=20,maxvel,nbins=401,*cf,wgal[3],sgal[4],phi;
  double v,dv=1,p1=0,p2=0,p=0,p0=0,p3=0,p4=0,nhalo;
  int i,n=4000,j;

  nhalo = dndM_interp(m1)*dndM_interp(m2)*pow(BOX_SIZE,3.0)*log(2)*log(2)*m1*m2/2*
    (1+xi_interp(r)*bias_interp(m1,r)*bias_interp(m2,r))*4*PI*r*r*r*dlogr;

  if(!flag++)
    prob = dvector(1,nbins);

  dv = binsize;
  n = nbins;
  maxvel = binsize*floor(nbins/2);
  wgal[0] = wgal[1] = wgal[2] = 0;
  sgal[0] = sgal[1] = sgal[2] = sgal[3] = 0;
  
  vdelta_v4(m1,m2,r,0.0,binsize,-maxvel,nbins,cf,
	    cf,cf,prob,wgal,sgal);

  for(j=i=1;i<=n;++i)
    if(prob[i]>prob[j])j=i;
  

  for(i=1;i<=n;++i)
    {
      v=-maxvel+(i-1)*binsize;
      p = prob[i];
      if(p*dv<1/nhalo)p=0;
      /* printf("PROB %f %e\n",v,p); */
      p0 += p*dv;
      p2 += v*v*p*dv;
      p4 += v*v*v*v*p*dv;
    }
  p2 = p2/p0;
  vtan[2] = sqrt(p2);
  p4 = p4/p0*pow(p2,-2.0)-3.0;
  vtan[4] = p4;
  
  vdelta_v4(m1,m2,r,PI/2.0,binsize,-maxvel,nbins,cf,
	    cf,cf,prob,wgal,sgal);

  p0 = p1 = p2 = p3 = p4 = 0;
  for(i=1;i<=n;++i)
    {
      v=-maxvel+(i-1)*binsize;
      p = prob[i]; 
      if(p*dv<1/nhalo)p=0;
      p0 += p*dv;
      p1 += v*p*dv;
      p2 += v*v*p*dv;
    }
  p1 = p1/p0;
  p2 = sqrt(p2/p0-p1*p1);
  for(i=1;i<=n;++i)
    {
      v=-maxvel+(i-1)*binsize;
      p = prob[i];
      if(p*dv<1/nhalo)p=0;
      p3 += pow(v-p1,3.0)*p*dv;
      p4 += pow(v-p1,4.0)*p*dv;
    }
  p3 = (p3/p0)*pow(p2,-3.0);
  p4 = (p4/p0)*pow(p2,-4.0)-3.0; 
  vrad[1] = p1;
  vrad[2] = p2;
  vrad[3] = p3;
  vrad[4] = p4;
  return(nhalo);
}


