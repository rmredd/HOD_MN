#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "header.h"

void linlin_bins()
{
  FILE *fp,*fp2;
  char fname[100];
  int nsize_rpi,nsize_rsig,i,j;
  double rp,rs,r,dx1,dx2,dx3,dx4,dx5,rslo,rshi,drs,binsize_rpi=2;
  float x1;


  nsize_rpi=15;
  nsize_rsig=15;
  binsize_rpi=2.0;

  sprintf(fname,"%s.linlin",Task.root_filename);
  fp=fopen(fname,"w");

  rshi=0;
  for(i=1;i<=nsize_rsig;++i)
    {
      rslo=rshi;
      rshi=rslo+binsize_rpi;
      drs=rshi-rslo;
      rs=0.5*(rshi+rslo);

      for(j=1;j<=nsize_rpi;++j)
	{
	  rp=(j-0.5)*binsize_rpi;
	  r=sqrt(rs*rs+rp*rp);
	  if(!KAISER)
	    {
	      dx4=two_halo(rs,rp);
	      dx1=one_halo(rs,rp);
	    }
	  else
	    {
	      dx1=0;
	      dx4=0;
	    }
	  /*
	  dx5=integrated_bin(rslo,(j-1)*binsize_rpi,drs,binsize_rpi,10);
	  */
	  dx5=xi2d_interp(rslo,(j-1)*binsize_rpi,rslo+drs,j*binsize_rpi);
	  fprintf(fp,"%e %e %e %e %e\n",rs,rp,dx5,dx1,dx4);
	  fprintf(stdout,"%e %e %e %e %e\n",rs,rp,dx5,dx1,dx4);
	  fflush(stdout);
	  fflush(fp);
	}      
    }
  fclose(fp);
  

}
