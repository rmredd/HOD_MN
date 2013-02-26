#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* It is not elegant to do things in this manner, but the code was
 * written with auto-correlation functions in mind first, and so copy+edit for
 * a second HOD was the most efficient method of modifying the code for the
 * cross-correlation function.
 *
 * NB-- No blue/red galaxies for the second HOD.
 */

/* This file holds several useful functions for the HOD, such as the number of
 * central and satellite halos, their second moment.
 *
 * This also intializes the HOD parameters.
 */

/* internal routines.
 */
double fixfunc1_2(double m);
double fixfunc2_2(double m);
double func_mlow2(double m);
double func_mhi2(double m);

/*****************
 * DIFFERENT HOD PARAMETERIZATIONS: (1) Central galaxies
 *
 * This integer HOD.pdfc controls the parameterization of the central occupation.
 * (I realize 'pdf' is not the right acronym, but it's a holdover from old code.)
 *
 * For soft central cutoffs, the PDF is the nearest integer distribution, AKA
 * Bernoulli distribution.
 *
 * 0 = No galaxies, just halos. (always returns 1)
 * 1 = Hard cutoff at M_min. Unity above M_min. (DEFAULT)
 * 2 = Soft cutoff of the form 0.5*(1+erf((log10(m) - log10(HOD.M_min))/HOD.sigma_logM)
 * 3 = Soft cutoff of the form N_cen = exp(-HOD.M_min/m)
 * 4 = Hard cutoff, but N_cen -> 0 for M>M_min (i.e., blue galaxies) 
 *      MaxCen*exp(-(lgM-lgM_min)**2/2/(sigma_logM)**2)
 *      NB! -> The value of N_cen(M_min) can be < 1. (It will== MaxCen)
 * 5 = Hard Cutoff, Step function, but N_cen != 1, but some number < 1.
 * 6 = Same as 4, but symmetric Gaussian around M_min (This is for mag bins.)
 * 7 = A sqaure function. N_cen=1 for M_min <= m <= M_cen_max (for mag bin data)
 * 8 = A sqaure function (like case 7) but with Gaussian cutoffs on either edge
 *      instead of sharp cutoffs.
 *
 * 9 = Magnitude bin model, but the upper mass cutoff is defined by the lower
 *      mass cutoff of the next-highest bin. (see function in ml_minimization.c)
 *
 * M_low is a parameter created to keep the integrals from having to go to M=0 for
 * soft central cutoffs. N_cen(mlow) = 1.0E-3;
 * If the value of sigma_logM gets above 1 or so, then M_low can be rediculously 
 * small, so I have placed a lower limit on M_low of 1.0E+7 Msol/h
 */
double N_cen2(double m)
{
  double x,f=1;

  switch(HOD2.pdfc) {
  case 0:
    return 1;
  case 1: 
    if(m<HOD2.M_min)return(0);
    return(f*1.0);
    break;
  case 2:
    return(f*0.5*(1+erf((log10(m) - log10(HOD2.M_min))/HOD2.sigma_logM)));
    break;
  case 3:
    return(f*exp(-HOD2.M_min/m));
    break;
  case 4:
    if(m<HOD2.M_min)return(0);
    x = (log10(m) - log10(HOD2.M_min))/HOD2.sigma_logM;
    return(f*HOD2.MaxCen*exp(-x*x/2));
  case 5:
    if(m<HOD2.M_min)return(0);
    return(f*HOD2.MaxCen);
  case 6:
    x = (log10(m) - log10(HOD2.M_min))/HOD2.sigma_logM;
    return(f*HOD2.MaxCen*exp(-x*x/2));    
  case 7:
    if(m>=HOD2.M_min && m<=HOD2.M_cen_max)return(f);
    return(0);
  case 8:
    if(m>=HOD2.M_min && m<=HOD2.M_cen_max)return(f);
    if(m<HOD2.M_low)return(0);
    if(m<HOD2.M_min)
      x = (log10(m) - log10(HOD2.M_min))/HOD2.sigma_logM;
    else
      x = (log10(m) - log10(HOD2.M_cen_max))/HOD2.sigma_logM;
    return(f*exp(-x*x/2));    
  case 9:
    return(f*N_cen_i(m,HOD2.i_wp));
  default:
    endrun("Illegal value of HOD2.pdfc.");
  }
  return 0;
}

/*****************
 * DIFFERENT HOD PARAMETERIZATIONS: (1) Satellite galaxies
 *
 * This integer HOD2.pdfs controls the parameterization of the satellite occupation.
 *
 * 0 = halos only, no galaxies (always returns zero)
 * 1 = power law: N_sat = pow(m/HOD2.M1,HOD2.alpha) [cut off at M_low]
 * 2 = power law with soft cutoff: N_sat = pow((m-HOD2.M_cut)/HOD2.M1,alpha) 
 * 3 = power law with exp cutoff: N_sat = pow(m/HOD2.M1,alpha)*exp(-M_cut/(m-M_min))
 * 4 = broken power law with exp cutoff (alpha changes at some high mass value)
 * 5 = broken power law (m-mcut)/m1 parameterization
 * 
 * The PDF is always assumed to be Poisson for satellite galaxies.
 */
double N_sat2(double m)
{
  double m1,f=1;

  switch(HOD2.pdfs) {
  case 0:
    return 0;
  case 1: 
    if(m<HOD2.M_min)return(0);
    return(f*pow(m/HOD2.M1,HOD2.alpha));
    break;
  case 2:
    if(m<HOD2.M_low || m<HOD2.M_cut)return(0);
    return(f*pow((m-HOD2.M_cut)/HOD2.M1,HOD2.alpha));
    break;
  case 3:
    if(m<HOD2.M_min)return(0);
    return(f*exp(-HOD2.M_cut/(m-HOD2.M_min))*pow(m/HOD2.M1,HOD2.alpha));
    break;
  case 4:
    if(m<HOD2.M_min)return(0);
    if(m<HOD2.M_sat_break)
      return(f*exp(-HOD2.M_cut/(m-HOD2.M_min))*pow(m/HOD2.M1,HOD2.alpha));
    else
      {
	m1 = exp(log(HOD2.M_sat_break) - HOD2.alpha/HOD2.alpha1*log(HOD2.M_sat_break/HOD2.M1));
	/*
	m1 = exp(log(HOD2.M_sat_break) - 1/HOD2.alpha1*
		 (HOD2.alpha*log(HOD2.M_sat_break/HOD2.M1) - HOD2.M_cut/(HOD2.M_sat_break-HOD2.M_min)));
	*/
	return(f*exp(-HOD2.M_cut/(m-HOD2.M_min))*pow(m/m1,HOD2.alpha1));
      }
    break;
  case 5:
    if(m<HOD2.M_low || m<HOD2.M_cut)return(0);
    if(m<HOD2.M_sat_break)
      return(f*pow((m-HOD2.M_cut)/HOD2.M1,HOD2.alpha));
    else
      {
	m1 = HOD2.M_sat_break*pow((HOD2.M_sat_break-HOD2.M_cut)/HOD2.M1,-HOD2.alpha/HOD2.alpha1);
	return(f*pow(m/m1,HOD2.alpha1));
      }
    break;
  default:
    endrun("Illegal value of HOD2.pdfs.");
  }
  return 0;
}


/* If what is needed is the total number of galaxies in a halo.
 */
double N_avg2(double m)
{
  return(N_cen2(m)+N_sat2(m));
}


/* This is the <(N_sat-1)(N_sat)> moment, which is
 * for the number of pairs of satellite galaxies.
 * For a Poisson distribution, <N(N-1)> = N^2
 */
double moment_ss2(double m)
{
  double n;
  n=N_sat2(m);
  return(n*n);
}

/* This is a function to set the HOD parameters until 
 * I get some input code or batch file set up.
 * 
 */
void set_HOD2_params()
{
  int i,j=1;
  double m,error=1.0,tol=1.0E-4,prev,s1,mlo;


  if(HOD2.pdfc == 2 || HOD2.pdfc == 3 || HOD2.pdfc == 6 || HOD2.pdfc == 8 || HOD2.pdfc == 9)
    SOFT_CENTRAL_CUTOFF=1;

  /* Error trap both the galaxy density and M_min both left unspecified.
   */
  if(HOD2.M_min<=0 && GALAXY_DENSITY2<=0 && HOD2.free[0]==0)
    endrun("ERROR: Must specify either M_min or GALAXY_DENSITY2");

  /* If the user has specified M_min and M1, calculate the galaxy density.
   */
  if(HOD2.M_min>0 && HOD2.M1>0)
    {
      HOD2.M_low = -1;
      if(SOFT_CENTRAL_CUTOFF)
	HOD2.M_low = exp(zbrent(func_mlow2,log(HOD2.M_min*1.0E-6),log(HOD2.M_min*1.1),1.0E-5));
      else
	HOD2.M_low = HOD2.M_min;
      GALAXY_DENSITY2=qromo(func_galaxy_density2,log(HOD2.M_low),log(HOD2.M_max),midpnt);
      if(OUTPUT) {
	fprintf(stdout,"M_low= %e\n",HOD2.M_low);
	fprintf(stdout,"ng= %e\n",GALAXY_DENSITY2); }
    }      

  /* If M_min<=0 then use the specified galaxy density to calculate M_min
   */
  if(HOD2.M_min<=0)
    {
      if(HOD2.pdfc==7 && HOD2.pdfc==8)
	HOD2.M_min=pow(10.0,zbrent(fixfunc1_2,8.0,log10(HOD2.M_cen_max*0.99),1.0E-5));
      else
	HOD2.M_min=pow(10.0,zbrent(fixfunc1_2,8.0,14.8,1.0E-5));
	
      HOD2.M_low = -1;
      if(SOFT_CENTRAL_CUTOFF)
	HOD2.M_low = exp(zbrent(func_mlow2,log(HOD2.M_min)-5*HOD2.sigma_logM*2.3,
			       log(HOD2.M_min),1.0E-5));
      else
	HOD2.M_low = HOD2.M_min;
      if(HOD2.M_low<1.0E7)HOD2.M_low=1.0E+7;
      if(OUTPUT) {
	fprintf(stdout,"M_min %e [ng= %e]\n",HOD2.M_min,GALAXY_DENSITY2);
	fprintf(stdout,"M_low= %e\n",HOD2.M_low); }
    }

  /* If M1<=0 then use the specified galaxy density to calculate M1
   */
  if(HOD2.M1<=0)
    {
      HOD2.M_low = -1;
      if(SOFT_CENTRAL_CUTOFF)
	HOD2.M_low = exp(zbrent(func_mlow2,log(HOD2.M_min)-5*HOD2.sigma_logM*2.3,
			       log(HOD2.M_min*1.1),1.0E-5));
      else
	HOD2.M_low = HOD2.M_min;
      if(HOD2.M_low<1.0E7)HOD2.M_low=1.0E+7;
      HOD2.M1=pow(10.0,zbrent(fixfunc2_2,log10(HOD2.M_low),15.8,1.0E-5));
      if(OUTPUT) {
	fprintf(stdout,"M1 %e [ng= %e]\n",HOD2.M1,GALAXY_DENSITY2);
	fprintf(stdout,"M_min = %e M_low= %e\n",HOD2.M_min,HOD2.M_low); }
    }

  HOD2.M_hi = set_high_central_mass2();
  return;

}

/* If soft central cutoff, then put a lower limit on mass integrals
 * that begins at the mass where N_cen=0.001.
 */
double func_mlow2(double m)
{
  /* Have a check in case the passed mass is equal to M_min, but the value of
   * N_cen is < 0.001 (which might be the case for very small sigma_logM)
   */
  if(fabs(exp(m)-HOD2.M_min)<0.001*HOD2.M_min)
    if(N_cen2(exp(m))<0.001)return(0);
  /* fprintf(stderr,"MLO %e %e %e %e %e\n",exp(m),N_cen(exp(m)),HOD2.M_min,HOD2.M_cen_max,HOD2.M_low); */
  return(N_cen2(exp(m))-0.001);
}

/* It is straightforward to calculate what M_low
 * should be given the other parameters on N_cen.
 */
double set_low_mass2()
{
  double m;

  if(!SOFT_CENTRAL_CUTOFF)return(HOD2.M_min);
  switch(HOD2.pdfc){
  case 8:
  case 6:
    m = log10(HOD2.M_min) - sqrt(-2*HOD2.sigma_logM*HOD2.sigma_logM*log(0.001));
    m = pow(10.0,m);
    return(m);
  default:
    m = exp(zbrent(func_mlow2,log(HOD2.M_min)-5*HOD2.sigma_logM*2.3,
		   log(HOD2.M_min),1.0E-5));
    return(m);
  }
  return(0);
}

/* If modeling magnitude bin samples, then there will be a high mass
 * scale for central galaxies as well as a low mass scale. This finds
 * the mass at which N_cen(m)=0.001, where m>M_min.
 */
double set_high_central_mass2()
{
  double m,n;

  if(HOD2.pdfc==7)
    return(HOD2.M_cen_max);
  if(!(HOD2.pdfc==6 || HOD2.pdfc==8 || HOD2.pdfc==9))
    return(HOD2.M_max);

  m = HOD2.M_min;
  n = N_cen(m);

  while(n>0.001)
    {
      m*=2;
      n = N_cen(m);
      if(m>HOD2.M_max)return(HOD2.M_max);
    }
  m = exp(zbrent(func_mhi2,log(m/2),log(m),1.0E-5));
  return(m);
}

double func_mhi2(double m)
{
  m=exp(m);
  return(N_cen2(m)-0.001);
}


/* This is a copy of the above function that can be called from any routine.
 * (But this one integrates over dlogm
 */
double func_galaxy_density2(double m)
{
  double n1,n2,m0;

  m=exp(m);
  n1=dndM_interp(m);
  n2=N_avg2(m);
  return(n1*n2*m);
}

/* This is the equation for zbrent to solve. What value
 * of M_min gives the correct galaxy density?
 * For HODs that have sharp central cutoffs, use M_min as the lower
 * limit of the integration. For soft cutoffs, first find M_low.
 */
double fixfunc1_2(double m)
{
  double n,mlo;

  HOD2.M_min=m=pow(10.0,m);
  HOD2.M_low=0;
  if(SOFT_CENTRAL_CUTOFF)
    mlo = (zbrent(func_mlow2,log(HOD2.M_min)-5*HOD2.sigma_logM*2.3,log(HOD2.M_min),1.0E-5));
  else
    mlo = log(HOD2.M_min);
  if(exp(mlo)<1.0E7)mlo=log(1.0E+7);
  n=qromo(func_galaxy_density,mlo,log(HOD2.M_max),midpnt);
  /* fprintf(stderr,"MMIN %e %e %e %e\n",m,exp(mlo),n,GALAXY_DENSITY2); */
  return(GALAXY_DENSITY2-n);
}

/* This function is sent to zbrent to determine what M1 is based
 * on the number density. Both M_min and M_low have already been specified.
 */
double fixfunc2_2(double m)
{
  double n;
  HOD2.M1=m=pow(10.0,m);  
  n=qromo(func_galaxy_density,log(HOD2.M_low),log(HOD2.M_max),midpnt);
  return(GALAXY_DENSITY2-n);
}


/* This function is to be passed to qromo to integrate the number density
 * of satellite galaxies.
 */
double func_satfrac2(double m)
{
  m=exp(m);
  return(N_sat2(m)*dndM_interp(m)*m);
}

/* This function is to be passed to qromo to integrate the number density
 * of satellite galaxies.
 */
double func_satellite_density2(double m)
{
  m=exp(m);
  return(N_sat2(m)*dndM_interp(m)*m);
}

/* This function is to be passed to qromo to integrate the number density
 * of satellite galaxies.
 */
double func_central_density2(double m)
{
  m=exp(m);
  return(N_cen2(m)*dndM_interp(m)*m);
}
