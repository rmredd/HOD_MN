#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"

int main(int argc, char **argv)
{
  int i,j,k,n;
  float *x,*y,*z,*vx,*vy,*vz,*mass,*rvir,*vmax;
  FILE *fp;
  float dx,dy,dz,r,xpos,ypos,zpos;
 
  fp = openfile(argv[1]);
  n = filesize(fp);
  
  x = vector(1,n);
  y = vector(1,n);
  z = vector(1,n);
  vx = vector(1,n);
  vy = vector(1,n);
  vz = vector(1,n);

  mass = vector(1,n);
  rvir = vector(1,n);
  vmax = vector(1,n);

  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%d %f %f %f %f %f %f %f %f %f",
	     &j,&mass[i],&rvir[i],&vmax[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i]);
    }
  for(i=1;i<=n;++i)
    {
      if(i%1000==0)fprintf(stderr,"%d/%d\n",i,n);
      xpos = x[i];
      ypos = y[i];
      zpos = z[i];
      
      for(j=i+1;j<=n;++j)
	{
	  if(mass[i]<0)continue;
	  dx = xpos - x[j];
	  dy = ypos - y[j];
	  dz = zpos - z[j];
	  r = sqrt(dx*dx + dy*dy + dz*dz);
	  if(r<=rvir[j] || r<=rvir[i]) 
	    {
	      if(vmax[i]<vmax[j])mass[i]=-1;
	      else mass[j]=-1;
	    }
	}
    }

  k=0;
  for(i=1;i<=n;++i)
    {
      if(mass[i]<0)continue;
      printf("%7d %.5e %9.3f %9.3f %9.3f %9.3f %9.3f %10.3f %10.3f %10.3f\n",
	     ++k,mass[i],rvir[i],vmax[i],x[i],
	     y[i],z[i],vx[i],vy[i],vz[i]);
    }


}
