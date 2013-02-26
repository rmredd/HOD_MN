#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifdef PARALLEL
#include <mpi.h>
#endif
#include "header.h"

double integrated_bin(double xlo, double ylo, double dx1, double dy1, int n)
{
  int i,j;
  double dx,dy,xsum=0,vol=0,x1,x2,rs,rp;

  dx=dx1/n;
  dy=dy1/n;

  for(i=1;i<=n;++i)
    for(j=1;j<=n;++j)
      {
	rs=xlo+(i-0.5)*dx;
	rp=ylo+(j-0.5)*dy;
	x2=two_halo(rs,rp);
	x1=one_halo(rs,rp);
	xsum+=2*dy*rs*dx*(x1+x2);
      }
  dx=xlo+dx1;
  dy=ylo+dy1;
  xsum/=((dx*dx-xlo*xlo)*dy1);
  return(xsum);

}
