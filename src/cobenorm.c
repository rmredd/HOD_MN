#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"


/* This is the one that Andrew Zentner gave me:
 * taken from Bunn, Liddle & White (and appears to 
 * give different answers from the Bunn & White functions.
 */
double cobenorm(double Om)
{
  double Ol, n, f, pwr, g, r=0, d0, 
    cobe_norm;

  n = SPECTRAL_INDX;

  Ol = 1 - Om;

  f= 0.750 - 0.130*(Ol*Ol);

  pwr = -0.8 - 0.05*log(Om);

  g = 1.0-0.180*(1.0-n)*Ol - 0.03*r*Ol;

  d0 = 1.91e-5*(pow(7.0,((n-1.0)/2.0)));

  cobe_norm = d0*(exp(1.010*(1.0-n))/
		  sqrt(1.0+r*f))*(pow(Om,pwr))*g;

  return cobe_norm;
}

/* 

      Ol = 1.d0-Om

      f=0.75d0-0.13d0*(Ol**2.d0)

      pwr = -0.8d0 - 0.05d0*log(Om)

      g = 1.d0-0.18d0*(1.d0-n)*Ol - 0.03d0*r*Ol

      d0 = 1.91d-5*(7.d0**((n-1.d0)/2.d0))

      cobe_norm = d0*(exp(1.01d0*(1.d0-n))/
     & dsqrt(1.d0+r*f))*(Om**pwr)*g

*/

/* The function below was supplied by Risa Wechsler;
 * I set it up to always assume flat universe.
 */
double cobenorm_risa(double omega_m)
/* Return the Bunn & White (1997) fit for delta_H */
/* Given lambda, omega_m, qtensors, and tilt */
/* Open model with tensors is from Hu & White */
{
  //  cout<<omega_m<<" "<<lambda<<endl;
  double n,
    //  omega_m,
    lambda;
  int qtensors = 0;

  /* n = tilt-1; */

  n = SPECTRAL_INDX;
  /* omega_m = OMEGA_M; */
  lambda = 1-omega_m;


  if (fabs(omega_m+lambda-1.0)<1e-5) {	/* Flat universe */
    if (qtensors)
      return 1.94e-5*pow(omega_m, -0.785-0.05*log(omega_m))*
	exp(n+1.97*n*n);
    else
      return 1.94e-5*pow(omega_m, -0.785-0.05*log(omega_m))*
	exp(-0.95*n-0.169*n*n);
   } else if (fabs(lambda)<1e-5) {	/* No lambda */
    if (qtensors)
      return 1.95e-5*pow(omega_m,-0.35-0.19*log(omega_m)-0.15*n)*
	exp(+1.02*n+1.7*n*n);
     else return 1.95e-5*pow(omega_m, -0.35-0.19*log(omega_m)-0.17*n)*
	   exp(-n-0.14*n*n);
   } else return 1e-5*(2.422-1.166*exp(omega_m)+0.800*exp(lambda)
		      +3.780*omega_m-2.267*omega_m*exp(lambda)+0.487*SQR(omega_m)+
		      0.561*lambda+3.329*lambda*exp(omega_m)-8.568*omega_m*lambda+
		      1.080*SQR(lambda));
}


/* Read in the cosmo parameters from a completed chain and calculate the chi^2 
 * of the matter power spectrum
 * with respect to the COBE normalization.
 */
double cobe_prior(double omega_m)
{
  double pk_model,pk_cobe,chi2;
  static int niter=0;
  pk_model = linear_power_spectrum(0.0023/7)/pow(transfnc(0.0023/7),2.0);
  pk_cobe = cobenorm(omega_m);
  pk_cobe*=pk_cobe;
  chi2 = (pk_cobe - pk_model)*(pk_cobe - pk_model)/(0.07*pk_cobe*0.07*pk_cobe);
  niter++;
  printf("COBE%d %e %e %e\n",niter,pk_cobe,pk_model,chi2);

  return chi2;
}


/* Read in the cosmo parameters from a completed chain and calculate the chi^2 
 * of the matter power spectrum
 * with respect to the COBE normalization.
 */
void cobe_prior_from_file(char *filename)
{
  FILE *fp;
  int n,i,j,k,i1,i2,ip=0;
  float xx[10];
  char aa[4];

  double pk_model,pk_cobe,chi2,pnorm,p2,xk;
  double pka[20][41];

  fp = openfile(filename);
  n = filesize(fp);

  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%4s %d %d",aa,&i1,&i2);
      for(j=0;j<10;++j)fscanf(fp,"%f",&xx[j]);
      SIGMA_8 = xx[6];
      SPECTRAL_INDX = xx[7];
      OMEGA_TEMP = xx[5];

      RESET_COSMOLOGY++;
      pk_model = linear_power_spectrum(0.0023/7)/pow(transfnc(0.0023/7),2.0);
      pk_cobe = cobenorm(OMEGA_TEMP)*cobenorm(OMEGA_TEMP);
      chi2 = (pk_cobe - pk_model)*(pk_cobe - pk_model)/(0.07*pk_cobe*0.07*pk_cobe);

      printf("%d %e %e %e %f %f %f\n",
	     i1,pk_model,pk_cobe,chi2,SIGMA_8,SPECTRAL_INDX,OMEGA_TEMP);
      fflush(stdout);

      if(i%(n/20)==0 && ip<20) {
	for(k=0,j=-40;j<=0;++j)
	  {
	    xk = pow(10.0,j/10.0);
	    pka[ip][k] = linear_power_spectrum(xk)/(xk*xk*xk)*2*PI*PI;
	    k++;
	  }
	ip++;
      }
    }

  for(i=0;i<=41;++i)
    {
      printf("BOO %f ",(i-40)/10.0);
      for(j=0;j<20;++j)
	printf("%e ",pka[j][i]);
      printf("\n");
    }

  exit(0);
}
