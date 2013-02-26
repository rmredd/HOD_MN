#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "header.h"

/* This will calculate the radial separation (pi) at which the value of \xi(sigma,pi) falls
 * by a factor of 2 relative to the value at \xi(0,pi) (the plots demonstrate that xi 
 * reaches a horizontal asymptote at small pi).
 *
 * I guess I'll do two versions of this:
 *  1) assume that xi(sigma,pi) has already been calculated and just inpterpolate.
 *  2) calculate xi(sigma,pi) as needed.
 *
 * Input:
 *  r_sigma--> the transverse separation at which to do the measure.
 * 
 * Output:
 *  r_xi/2(sigma)--> in Mpc/h
 */

double RHALF_FACTOR = 0.5;

void calc_rhalf(double r[], double rhalf[], int nr)
{
  double *recv;
  int i,istep=1,istart=0;
  char fname[100];
  FILE *fp2;

#ifdef PARALLEL
  istep = NTask;
  istart = ThisTask;
#endif

  recv=dvector(0,nr-1);

  for(i=0;i<nr;++i)
    rhalf[i]=0;

  for(i=istart;i<nr;i+=istep)
    rhalf[i]=small_scale_measure(r[i]);

#ifdef PARALLEL
  MPI_Allreduce(rhalf,recv,nr,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(i=0;i<nr;++i)
    rhalf[i]=recv[i];
#endif

  free_dvector(recv,0,nr-1);

}

double small_scale_measure(double rs)
{
  int niter=0,i;
  double rpi,xi_0,x1a,x2a,x3a,x1b,x2b,x3b,rlo=0.2,rhi=10.0,xi,ratio,tol=0.001,error,
    xa[3],ya[3],a,b;


  xa[0]=0.05;
  xa[1]=0.12;
  xa[2]=0.19;

  /* To calculate xi_0, do it at three small pi and average.
   */
  ya[0]=x1a=one_halo(rs,0.05);
  ya[1]=x2a=one_halo(rs,0.12);
  ya[2]=x3a=one_halo(rs,0.19);

  /* The first two-halo call results in NANs. Need to check this out.
   */
  x1b=two_halo(rs,0.05);
  ya[0]+=x1b=two_halo(rs,0.05);
  ya[1]+=x2b=two_halo(rs,0.12);
  ya[2]+=x3b=two_halo(rs,0.19);

  /*
  printf("Calculating r_xi/2(sigma=%.3f):\n",rs);
  printf(" xi(%.3f,0.10)= %e + %e = %e\n",rs,x1a,x1b,x1a+x1b);
  printf(" xi(%.3f,0.11)= %e + %e = %e\n",rs,x2a,x2b,x2a+x2b);
  printf(" xi(%.3f,0.12)= %e + %e = %e\n",rs,x3a,x3b,x3a+x3b);
  */

  for(i=0;i<3;++i)
    {
      xa[i]=log(xa[i]);
      ya[i]=log(ya[i]);
    }
  least_squares(xa,ya,3,&a,&b);
  xi_0=exp(a+b*log(0.1));


  /*xi_0=(x1a+x1b+x2a+x2b+x3a+x3b)/3.0;*/
  if(OUTPUT)
    {
      printf(" xi(%.3f,pi->0)= %e\n",rs,xi_0);
      fflush(stdout);
    }
  if(isnan(xi_0))return(0);

  /* Check to make sure that the upper limit (in rpi) actually has \xi
   * less than xi(0)/2
   */
  xi=one_halo(rs,rhi)+two_halo(rs,rhi);
  while(xi>xi_0*RHALF_FACTOR)
    {
      /*printf("ssm> xi(%.2f)= %.2e > %.2e, resetting rhi to %.2f\n",rhi,xi,xi_0*0.5,2*rhi);*/
      rhi*=2;
      xi=one_halo(rs,rhi)+two_halo(rs,rhi);
    }

  rpi=0.5*(rlo+rhi);
  xi=one_halo(rs,rpi)+two_halo(rs,rpi);
  ratio=xi/xi_0;
  error=fabs(ratio-RHALF_FACTOR)/RHALF_FACTOR;
  while(error>tol)
    {
      niter++;
      if(niter>30)
	{
	  printf("Too many iterations in small_scale_measure.\n");
	  printf("Exiting at %f accuracy\n",error);
	  break;
	}
      if(ratio<RHALF_FACTOR)
	rhi=rpi;
      else
	rlo=rpi;
      rpi=0.5*(rlo+rhi);
      xi=one_halo(rs,rpi)+two_halo(rs,rpi);
      ratio=xi/xi_0;
      error=fabs(ratio-RHALF_FACTOR)/RHALF_FACTOR;
      if(OUTPUT)
	printf("SM [%d] rpi= %f xi= %e ratio= %f error= %f\n",niter,rpi,xi,ratio,error);
    }
  if(OUTPUT)
    printf("R_XI/2(%.3f)= %f [%d iterations]\n",rs,rpi,niter);
  return(rpi);

}

