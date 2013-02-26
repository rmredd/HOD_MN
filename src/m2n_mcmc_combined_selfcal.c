#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"


/* UPDATES: 
 *
 * Thu Dec 3 17:21:30 PST 2009 -- beginning getting rid of
 * FAKE DATA revisions. Added individual redshifts for each clustering
 * sample.
 * Tue Dec  8 13:28:18 CST 2009 -- Using new M/N extractions: changing the mass (and thus radius) 
 * when getting the number of galaxies in the virial radius.
 */

// parameters of Eduardo's P(M|N200)
/*
#define B_rozo -0.12596
#define alpha_rozo 1.149
*
#define B_rozo -1.1
#define alpha_rozo 0.74
#define sig_rozo 0.35
*/

/* Here we're SELF-CALIBRATING the scale-dependent bias. We're keeping the internal priors
 * but removing the scale-dep bias, obviously.
 */

/* External functions from wp_minimization.c
 */
void wp_input(void);
double poisson_prob(int n, double nave);

/* Internal functions.
 */
double m2n_chi2_wp_wrapper(double *a, int icall);
double m2n_initialize(double *a, double **cov1, double *avg1, double *start_dev);
void m2n_input(void);
void put_parameters_in_place(double *a, int icall);
double chi2_m2n(int icall);
double chi2_number_profiles(void);
int parameter_out_of_range(double *a);
void put_parameters_in_place(double *a, int icall);
double internal_priors(void);

double prob_ngals_ntrue(double m, int ngals, int ntrue);
double prob_ntrue_mass(double m, double ntrue);
double m2n_func1(double m);
double m2n_func1a(double m);
double m2n_func2(double m);
double m2n_func3(double m);
double comoving_volume(double zlo, double zhi);
double efunc1(double z);
void input_maxbcg_counts(void);

void test_pdf(void);

/* Variables.
 */
int USE_IWEIGHT;
double H2OFZ;

/* Rozo variables.
 */
//double B_rozo = -1.02, B_rozo0 = -1.02, B_rozo_err = 0.15;
//double alpha_rozo = 0.74, alpha_rozo0 = 0.74, alpha_rozo_err = 0.1;
//double sig_rozo = 0.35, sig_rozo0 = 0.35, sig_rozo_err = 0.1;

//Fri Mar  6 16:24:51 CST 2009
// confused about above. appears I've been using P(N|M) above? instead of P(M|N) which i wanted
//double B_rozo = 0.95, B_rozo0 = 0.95, B_rozo_err = 0.12;
//double alpha_rozo = 1.06, alpha_rozo0 = 1.06, alpha_rozo_err = 0.11;
//double sig_rozo = 0.45, sig_rozo0 = 0.45, sig_rozo_err = 0.10;

//Fri Mar  6 17:13:35 CST 2009
// NO! i need to be using P(N|M)-- but above seems wrong; using new params from Eduardo w/ WMAP5 priors
// also this is <lnN|M> so no need to -sig^2/2
double B_rozo = 1.09, B_rozo0 = 1.09, B_rozo_err = 0.09; //Mpivot = 2.06e13 Msol
double alpha_rozo = 0.75, alpha_rozo0 = 0.75, alpha_rozo_err = 0.024;
double sig_rozo = 0.35, sig_rozo0 = 0.35, sig_rozo_err = 0.07;
double PIVOT_MASS;

int USE_VARIABLE_ROZO = 1;
int COBENORM = 0;
int USE_MONTE_CARLO_SYSTEMATICS = 2; // this is for evolution of HOD, LF, weak lensing errors
double SYSTEMATIC_OFFSET;

int USE_INTERNAL_PRIORS=1;


//NB NB NB:: making these REALLY small
double bias_amp_err=0.06,
  mf_amp_err=0.05,
  scalebias_amp_err=0.15;
  

/******************************************************************
 *
 * HOD.free[] also controls which variables will be held constant/vary
 * during MCMC minimization. Since this routine will also so z-space
 * minimization if requested, indices>6 are cosmological.
 *
 *  i     variable
 * ---    --------
 * [1] ->  M_min
 * [2] ->  M1
 * [3] ->  alpha
 * [4] ->  M_cut
 * [5] ->  sigmaM
 * [6] ->  CVIR_FAC
 * [7] ->  MaxCen (or M_cen_max)
 * [8] ->  M_sat_break
 * [9] ->  alpha1
 *
 * [10]->  OMEGA_M
 * [11]->  SIGMA_8
 * [12]->  VBIAS
 * [13]->  VBIAS_C
 * [14]->  GAMMA
 * [15]->  SPECTRAL_INDX
 * [16]->  HUBBLE PARAMETER [used for P(k)]
 *
 * [0] -> The galaxy_density will be considered data with errors on it,
 *         and therefore no variable will be fixed by the galaxy density.
 * 
 */
void m2n_mcmc()
{
  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*avg1,chi2,chi2prev,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,**chain,*start_dev,*eval_prev;
  int n,i,j,k,nrot,niter=0,count=0,imax_chain=100000,NSTEP=50,NSTEP_MAX=10000,convergence=0;
  long IDUM=-555;

  int *pcheck,pcnt,ptot=20,firstflag=1,*iweight,total_weight;
  double t0,tprev,temp,chi2a,chi2b,chi2wp,s8original,chi2cmb;

  double delta_halo_rhoc = 200;


  M2N.scalebias_selfcal=1;
  M2N.bias_pow1 = 1.49/2.0;
  M2N.bias_amp1 = 1.17;


  //set up pivot mass: 2.06e13*HUBBLE fo M200mean--> change to M200crit
  PIVOT_MASS = halo_mass_conversion2(2.06e13*HUBBLE,halo_c200(2.06e13*HUBBLE),200.0,507.0);
  PIVOT_MASS = 2.06e13*HUBBLE;

  // initialize use of SYSTEMATIC ERROR
  USE_ERRORS = 1;
  M2N.IDUM = -5555;
  
  // fix the value of OEMGA_M fot T(k)
  M2N.fix_omegam = 0;
  M2N.constant_omegam = 0.27;

  opar=dvector(1,100);

  MCMC=2;
  OUTPUT = 0;
  pcheck=calloc(ptot,sizeof(int));

  /* Since we're at constant halo overdensity wrt RHO_CRIT,
   * set the overdensity of the halo
   */
  H2OFZ = 0.27*pow(1.25,3.0)+0.73;
  DELTA_HALO = delta_halo_rhoc/0.27*(0.27*1.25*1.25*1.25 + (1-0.27))/(1.25*1.25*1.25);

  //FAKE_DATA
  //DELTA_HALO = 507;

  /* read in the wp data for a single luminosity threshold sample.
   */
  // wp_input();

  /* read in the M2N data.
   */
  m2n_input();

  Work.imodel=2;
  Work.chi2=1;


  srand48(32498793);

  /* Find the number of free parameters in the minimization
   * for the real-space correlation function.
   */
  // hardwire: 4+5+5 HOD params, 4 cosmo
  n = 18;

  //FAKE_DATA removing the cosmology
  //n = 14;

  if(USE_ERRORS)
    n+=4;
  if(USE_VARIABLE_ROZO)
    n+=3;
  if(USE_MONTE_CARLO_SYSTEMATICS==2)
    n+=1;
  wp.ncf=n;

  if(OUTPUT)
    printf("mcmc_min> %d  free parameters\n",n);

  a=dvector(1,n);
  start_dev=dvector(1,n);
  aprev=dvector(1,n);
  atemp=dvector(1,n);
  cov1=dmatrix(1,n,1,n);
  avg1=dvector(1,n);

  tmp=dmatrix(1,n,1,n);
  tmp1=dmatrix(1,n,1,1);
  evect=dmatrix(1,n,1,n);
  eval=dvector(1,n);
  eval_prev=dvector(1,n);

  chain=dmatrix(1,imax_chain,1,n);
  iweight = ivector(1,imax_chain);
  for(i=1;i<=imax_chain;++i)
    iweight[i] = 0;

  IDUM=IDUM_MCMC;

  chi2prev=m2n_initialize(a,cov1,avg1,start_dev);
  niter++;
  for(i=1;i<=n;++i)
    {
      aprev[i] = a[i];
      chain[1][i] = a[i];
    }

  pcnt=0;
  pcheck[pcnt]=1;

  stepfac=1;

  //FAKE_DATA
  //stepfac = 0.1;

  while(niter<NSTEP)
    {
      pcnt++;
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  stepfac = stepfac*pow(0.9,5-j);
	  if(!ThisTask)printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      stepfac=0.5;
      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];

      /* Since we're at constant halo overdensity wrt RHO_CRIT,
       * set the overdensity of the halo
       */
      //DELTA_HALO = delta_halo_rhoc/OMEGA_M*(OMEGA_M*1.25*1.25*1.25 + (1-OMEGA_M))/(1.25*1.25*1.25);

      // put cosmo+other globals in place
      put_parameters_in_place(a,0);
      
      // Check to see if any of our parameters are out of range.
      if(parameter_out_of_range(a)){
	printf("RANGE 0 %d\n",parameter_out_of_range(a));
	for(i=1;i<=-n;++i)printf("RANGE %d %e\n",i,a[i]);
	continue; }

      // put HOD params in place, one at a time.
      s8original = SIGMA_8;
      chi2wp = chi2a = 0;
      for(i=1;i<=3;++i)
	{
	  //reset cosmology for z=0.1
	  if(i==1) { SIGMA_8 = s8original*growthfactor(0.068); REDSHIFT=0.068; }
	  if(i==2) { SIGMA_8 = s8original*growthfactor(0.104); REDSHIFT=0.104; }
	  if(i==3) { SIGMA_8 = s8original*growthfactor(0.126); REDSHIFT=0.126; }

	  //FAKE_DATA
	  //SIGMA_8 = s8original*growthfactor(0.25);

	  RESET_COSMOLOGY++;

	  put_parameters_in_place(a,i);
	  chi2wp += m2n_chi2_wp_wrapper(a,i);
	  
	  // reset cosmology for z=0.25
	  REDSHIFT = 0.25;
	  SIGMA_8 = s8original*growthfactor(0.25);
	  RESET_COSMOLOGY++;
	  chi2a = chi2_m2n(i);
	}
      chi2 = chi2a+chi2wp;
      if(USE_INTERNAL_PRIORS)
	chi2 = chi2 + internal_priors();
      if(COBENORM)
	{
	  SIGMA_8 = s8original;
	  RESET_COSMOLOGY++;
	  chi2 = chi2 + cobe_prior(OMEGA_M);
	}

      if(!ThisTask){
	printf("TRY %d ",++count);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	printf("%e\n",chi2);fflush(stdout);
      }
      pcheck[pcnt]=1;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	{
	  /* This for loop puts the prev element in the chain is
	   * the current trial point is rejected.
	   */
	  if(USE_IWEIGHT)
	    iweight[niter+1]++;
	  pcheck[pcnt]=0;
	  continue;
	  
	}

      niter++;
      iweight[niter]++;

      for(i=1;i<=n;++i)
	chain[niter][i]=a[i];
      for(i=1;i<=n;++i)
	avg1[i] += a[i];
      for(i=1;i<=n;++i)
	aprev[i] = a[i];
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cov1[i][j] += a[i]*a[j];
      chi2prev=chi2;

      if(!ThisTask){
	printf("ACCEPT %d %d ",niter,count);
	for(i=1;i<=n;++i)
	  printf("%e ",a[i]);
	printf("%e %e %e\n",chi2,chi2a,chi2wp);fflush(stdout);
	printf("HSTATS %d %e %e %e %e\n",niter,HOD.M_min,number_weighted_halo_mass(),
	       number_weighted_central_mass(),
	       qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);

	fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	printf("FSAT %d %e %e %e %e\n",niter,fsat,HOD.M_min,HOD.sigma_logM,qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
      }

    }

  stepfac=1.6/sqrt(n);
  pcnt=-1;
  t0 = second();

  NSTEP = niter;

  while(niter<imax_chain)
    {
      pcnt++;
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  stepfac=1.6/sqrt(n);
	  if(!ThisTask)printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      stepfac=1.6/sqrt(n)*1.2;
      //FAKE_DATA
      //stepfac=1.9/sqrt(n)*1.2;
      //FAKE_DATA L384
      //stepfac=2.5/sqrt(n)*1.2;

      // set convergence to 8000
      if(niter==8000) convergence=1;

      if(convergence)goto SKIP_MATRIX;

      for(j=1;j<=n;++j)
	{
	  avg1[j]=0;
	  for(k=1;k<=n;++k)
	    cov1[j][k]=0;
	}
      total_weight = 0;
      for(i=1;i<=niter;++i)
	{
	  for(j=1;j<=n;++j)
	    {
	      avg1[j]+=chain[i][j]*iweight[i];
	      for(k=1;k<=n;++k)
		cov1[j][k]+=chain[i][j]*chain[i][k]*iweight[i];
	    }
	  total_weight+=iweight[i];
	}

      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  tmp[i][j] = cov1[i][j]/total_weight - avg1[i]*avg1[j]/(total_weight*total_weight);

      jacobi(tmp,n,eval,evect,&nrot);
      gaussj(evect,n,tmp1,1);

    SKIP_MATRIX:
      for(i=1;i<=n;++i)
	atemp[i] = gasdev(&IDUM)*sqrt(eval[i])*stepfac;

      for(i=1;i<=n;++i)
	for(a[i]=0,j=1;j<=n;++j)
	  a[i] += atemp[j]*evect[j][i];

      for(i=1;i<=n;++i) 
	a[i] += aprev[i];

      /* Since we're at constant halo overdensity wrt RHO_CRIT,
       * set the overdensity of the halo
       */
      //DELTA_HALO = delta_halo_rhoc/OMEGA_M*(OMEGA_M*1.25*1.25*1.25 + (1-OMEGA_M))/(1.25*1.25*1.25);

      // put cosmo+other globals in place
      put_parameters_in_place(a,0);
      
      // Check to see if any of our parameters are out of range.
      if(parameter_out_of_range(a)){
	printf("RANGE 0 %d\n",parameter_out_of_range(a));
	for(i=1;i<=-n;++i)printf("RANGE %d %e\n",i,a[i]);
	continue; }

      // put HOD params in place, one at a time.
      s8original = SIGMA_8;
      chi2wp = chi2a = 0;
      for(i=1;i<=3;++i)
	{
	  //reset cosmology for z=0.1
	  if(i==1) { SIGMA_8 = s8original*growthfactor(0.068); REDSHIFT=0.068; }
	  if(i==2) { SIGMA_8 = s8original*growthfactor(0.104); REDSHIFT=0.104; }
	  if(i==3) { SIGMA_8 = s8original*growthfactor(0.126); REDSHIFT=0.126; }

	  //FAKE_DATA
	  //SIGMA_8 = s8original*growthfactor(0.25);

	  RESET_COSMOLOGY++;

	  put_parameters_in_place(a,i);
	  chi2wp += m2n_chi2_wp_wrapper(a,i);
	  
	  // reset cosmology for z=0.25
	  REDSHIFT = 0.25;
	  SIGMA_8 = s8original*growthfactor(0.25);
	  RESET_COSMOLOGY++;
	  chi2a = chi2_m2n(i);
	}
      chi2 = chi2a+chi2wp;
      if(USE_INTERNAL_PRIORS)
	chi2 = chi2 + internal_priors();
      if(COBENORM)
	{
	  SIGMA_8 = s8original;
	  RESET_COSMOLOGY++;
	  chi2 = chi2 + cobe_prior(OMEGA_M);
	}

      tprev = t0;
      t0 = second();
      ++count;
      if(!ThisTask) {
	printf("TRY %d ",count);
	for(i=1;i<=n;++i)
	  printf("%.4e ",a[i]);
	if(RESTART==2) {
	  printf("%e %e %.2f\n",chi2,chi2/(1+exp(-count/100.0)),
		 timediff(tprev,t0));fflush(stdout); }
	else {
	  printf("%e %.2f\n",chi2,
		 timediff(tprev,t0));fflush(stdout); }
      }

      pcheck[pcnt]=0;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	{
	  if(USE_IWEIGHT)
	    iweight[niter+1]++;
	  continue;
	}
      pcheck[pcnt]=1;

      //      if(NSTEP<NSTEP_MAX)NSTEP++;
      niter++;
      if(!convergence)NSTEP = niter;
      iweight[niter]++;

      if(niter%NSTEP_MAX==0 && !convergence && niter>NSTEP_MAX)
	{
	  convergence = 1;
	  for(i=1;i<=n;++i)
	    {
	      x1=fabs(eval[i]-eval_prev[i])/eval_prev[i];
	      if(x1>0.01)convergence = 0;
	      printf("CONVERGENCE CHECK %d %d %e %e %e\n",niter/NSTEP_MAX,i,x1,eval[i],eval_prev[i]);
	    }
	  for(i=1;i<=n;++i)
	    eval_prev[i] = eval[i];
	  convergence = 0;

	  if(convergence)
	    printf("CONVERGENCE ACCOMPLISHED %d %d \n",niter,count);	    
	}
      if(niter==NSTEP_MAX)
	{
	  for(i=1;i<=n;++i)
	    eval_prev[i] = eval[i];
	}


      for(i=1;i<=n;++i)
	chain[niter][i]=a[i];
      for(i=1;i<=n;++i)
	avg1[i] += a[i];
      for(i=1;i<=n;++i)
	aprev[i] = a[i];
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cov1[i][j] += a[i]*a[j];
      chi2prev=chi2;

      if(!ThisTask) {
	printf("ACCEPT %d %d ",niter,count);
	for(i=1;i<=n;++i)
	  printf("%e ",a[i]);
	printf("%e %e %e\n",chi2,chi2a,chi2wp);fflush(stdout);
	
      }

    }
}


/* This is a routine where i can hand-set some boundaries for the 
 * cosmological parameters (and such)
 */
int parameter_out_of_range(double *a)
{
  int i;

  if(OMEGA_M<0.1 || OMEGA_M>0.5) return 1;
  if(SIGMA_8<0.45 || SIGMA_8>1.2) return 2;
  if(SPECTRAL_INDX<0.86 || SPECTRAL_INDX>1.06) return 3;
  if(HUBBLE < 0.55 || HUBBLE > 0.85) return 4;

  // check cvir_fac for all 3 HODs
  if(a[4]<0.2 || a[4]>2) return 5;
  if(a[9]<0.2 || a[9]>2) return 6;
  if(a[14]<0.2 || a[14]>2) return 7;

  if(USE_ERRORS)
    {
      if(M2N.bias_amp > 1.30 || M2N.bias_amp < 0.7) return 8;
      if(M2N.mf_amp > 1.15 || M2N.mf_amp < 0.85) return 9;
      // if(M2N.scalebias_amp > 1.45 || M2N.scalebias_amp < 0.55) return 10;
    }
  if(USE_VARIABLE_ROZO)
    {
      if(B_rozo > B_rozo0+3*B_rozo_err || B_rozo < B_rozo0-3*B_rozo_err)return 11;
      if(alpha_rozo > alpha_rozo0+3*alpha_rozo_err || alpha_rozo < alpha_rozo0-3*alpha_rozo_err)return 12;
      if(sig_rozo > sig_rozo0+3*sig_rozo_err || sig_rozo < sig_rozo0-3*sig_rozo_err)return 13;
    }
  if(USE_MONTE_CARLO_SYSTEMATICS==2)
    {
      if(fabs(SYSTEMATIC_OFFSET-1)>0.127*3)return 14;
    }
  return 0;
}

double m2n_chi2_wp_wrapper(double *a, int icall)
{
  static long IDUM4=-2342;
  static int flag=1, niter=0;
  static double *b;
  double chi2;
  int i,j;

  if(flag)
    {
      b=dvector(1,100);
      flag=0;
    }

  /* check to make sure that none of the HOD parameters
   * are NEGATIVE. (with exception of sigma_logM)
   */
  for(j=0,i=1;i<=N_HOD_PARAMS;++i) {
    if(a[++j]<=0 && j!=8 && j!=13) { printf("NEG %d %d %e\n",i,j,a[j]); return(1.0E7); }
  }

  // if ICALL==1, then get the SYSTEMATIC_OFFSET (if required)
  // use 3-sigma clipping
  if(icall==1 && USE_MONTE_CARLO_SYSTEMATICS==1)
    {
      SYSTEMATIC_OFFSET = gasdev(&IDUM4)*0.127 + 1;
      while((fabs(SYSTEMATIC_OFFSET-1)-3*0.127)>0)
	SYSTEMATIC_OFFSET = gasdev(&IDUM4)*0.127 + 1;

      // monte carlo the DM systematics
      if(!USE_ERRORS)
	{
	  M2N.bias_amp = gasdev(&IDUM4)*bias_amp_err + 1;
	  while((fabs(M2N.bias_amp-1)-3*bias_amp_err)>0)
	    M2N.bias_amp = gasdev(&IDUM4)*bias_amp_err + 1;

	  M2N.mf_amp = gasdev(&IDUM4)*mf_amp_err + 1;
	  while((fabs(M2N.mf_amp-1)-3*mf_amp_err)>0)
	    M2N.mf_amp = gasdev(&IDUM4)*mf_amp_err + 1;

	  //M2N.scalebias_amp = gasdev(&IDUM4)*scalebias_amp_err + 1;
	  //while((fabs(M2N.scalebias_amp-1)-3*scalebias_amp_err)>0)
	  // M2N.scalebias_amp = gasdev(&IDUM4)*scalebias_amp_err + 1;
	  printf("DARKMATTERSYS%d %f %f %f\n",M2N.bias_amp, M2N.mf_amp, M2N.scalebias_amp);
	}

      printf("SYSOFFSET%d %f\n",niter,SYSTEMATIC_OFFSET);
    }
  
  j = 0;
  if(icall==1) i = 0;
  if(icall==2) i = 4;
  if(icall==3) i = 9;

  b[++j] = pow(10.0,a[++i]); // M1
  b[++j] = a[++i]; // alpha
  b[++j] = pow(10.0,a[++i]); // Mcut
  if(icall>1)  b[++j] = pow(10.0,a[++i]); // sigma_logm
  b[++j] = a[++i]; // cvir_fac

  // put the data in place.
  wp.np = M2N.wpn[icall];
  for(i=1;i<=wp.np;++i)
    {
      wp.r[i] = M2N.wpr[icall][i];
      wp.x[i] = M2N.wpx[icall][i];
      wp.e[i] = M2N.wpe[icall][i];
    }

  //if FAKE_DATA, the comment this part out
  for(i=1;i<=wp.np;++i)
    for(j=1;j<=wp.np;++j)
      wp.covar[i][j] = M2N.wpcovar[icall][i][j];
  

  for(i=0;i<=100;++i)
    HOD.free[i] = 0;
  HOD.free[2] = HOD.free[3] = HOD.free[4] = 1;
  if(icall>1) HOD.free[5] = 1;
  else HOD.sigma_logM = 0.2;

  switch(icall) {
  case 1: GALAXY_DENSITY = 0.011134; break;
  case 2: GALAXY_DENSITY = 0.0031634; break;
  case 3: GALAXY_DENSITY = 0.0011486; break;
  }

  chi2 = chi2_wp(b);

  printf("M1MMIN%d %d %f %f\n",niter,icall,satellite_mass_scale()/HOD.M_min,HOD.M1/HOD.M_min);

  for(i=1;i<=wp.np;++i)
    printf("WP%d %d %f %e %e %e %e\n",niter,icall,wp.r[i],wp.x[i],M2N.wpmodel[i],wp.e[i],chi2);
  if(icall==3)niter++;
 
  //FAKE_DATA
  //chi2/=10;
 
  return chi2;
}

double m2n_initialize(double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0;
  double x1,x2,omega_m,chi2a,chi2,chi2wp,s8original,chi2cmb=0;
  long IDUM = -556;

  printf("Top of m2n_initialize\n");

  //set up N_HOD_PARAMS
  N_HOD_PARAMS = 14;

  //FAKE_DATA
  //FAKE_DATA
  //FAKE_DATA L384
  i = 0;
  a[++i] = log10(8.280946e+12); //M1 
  start_dev[i] = 0.001;
  a[++i] = 9.434173e-01;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(3.410614e+12); //Mcut
  start_dev[i] = 0.001;
  a[++i] = 1.0;           //CVIR_FAC
  start_dev[i] = 0.02;

  a[++i] = log10(4.079740e+13); //M1
  start_dev[i] = 0.001;
  a[++i] = 1.082929e+00;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(1.036308e+12); //Mcut
  start_dev[i] = 0.001;
  a[++i] =log10(0.3);    //sigma_logm
  start_dev[i] = 0.01;
  a[++i] = 1.0;           //CVIR_FAC
  start_dev[i] = 0.02;

  a[++i] = log10(1.080903e+14); //M1
  start_dev[i] = 0.001;
  a[++i] = 1.071983e+00;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(3.821450e+11); //Mcut
  start_dev[i] = 0.001;
  a[++i] =  log10(0.4);    //sigma_logm
  start_dev[i] = 0.01;
  a[++i] = 1.0;           //CVIR_FAC
  start_dev[i] = 0.02;


  //FAKE_DATA - MILLENNIUM
  i = 0;
  a[++i] = log10(6.313680e+12); //M1 
  start_dev[i] = 0.001;
  a[++i] = 1.10;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(2.740872e+11); //Mcut
  start_dev[i] = 0.001;
  a[++i] = 1;           //CVIR_FAC
  start_dev[i] = 0.02;

  a[++i] = log10( 2.028518e+13); //M1
  start_dev[i] = 0.001;
  a[++i] = 1.16;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(3.e8); //Mcut
  start_dev[i] = 0.001;
  a[++i] = -0.672;    //sigma_logm
  start_dev[i] = 0.01;
  a[++i] = 1;           //CVIR_FAC
  start_dev[i] = 0.02;

  a[++i] = log10(5.696855e+13); //M1
  start_dev[i] = 0.001;
  a[++i] = 1.207;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(3.726604e+12); //Mcut
  start_dev[i] = 0.001;
  a[++i] =  -0.4128;    //sigma_logm
  start_dev[i] = 0.01;
  a[++i] = 1;           //CVIR_FAC
  start_dev[i] = 0.02;


  //FAKE_DATA
  
  a[++i] = 0.25;//0.25; //omega_m
  start_dev[i] = 0.01;
  a[++i] = 0.9; //sigma8
  start_dev[i] = 0.01;
  a[++i] = 1.0; //ns
  start_dev[i] = 0.02;
  a[++i] = 0.72; //hubble
  start_dev[i] = 0.02;
  

  //PRODUCTION RUNS-- REAL DATA

  i = 0;
  //since we pretty much know what were doing here, let's just hard-wire things
  a[++i] = log10(6.2e12); //M1 
  start_dev[i] = 0.001;
  a[++i] = 0.95;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(3.7e11); //Mcut
  start_dev[i] = 0.001;
  a[++i] = 0.5;           //CVIR_FAC
  start_dev[i] = 0.02;

  a[++i] = log10(2.1e13); //M1
  start_dev[i] = 0.001;
  a[++i] = 0.95;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(3.7e12); //Mcut
  start_dev[i] = 0.001;
  a[++i] = log10(0.2);    //sigma_logm
  start_dev[i] = 0.01;
  a[++i] = 0.5;           //CVIR_FAC
  start_dev[i] = 0.02;

  a[++i] = log10(6.1e13); //M1
  start_dev[i] = 0.001;
  a[++i] = 0.95;          //alpha
  start_dev[i] = 0.03;
  a[++i] = log10(6.7e12); //Mcut
  start_dev[i] = 0.001;
  a[++i] = log10(0.4);    //sigma_logm
  start_dev[i] = 0.01;
  a[++i] = 0.5;           //CVIR_FAC
  start_dev[i] = 0.02;

  a[++i] = 0.26; //omega_m
  start_dev[i] = 0.01;
  a[++i] = 0.84; //sigma8
  start_dev[i] = 0.01;
  a[++i] = 0.96; //ns
  start_dev[i] = 0.02;
  a[++i] = 0.72; //hubble
  start_dev[i] = 0.02;


  if(USE_ERRORS)
    {
      a[++i] = 1.0;
      start_dev[i] = bias_amp_err/5;
      a[++i] = 1.0;
      start_dev[i] = mf_amp_err/5;
      a[++i] = M2N.bias_pow1;
      start_dev[i] = a[i]*0.01;
      a[++i] = M2N.bias_amp1;
      start_dev[i] = a[i]*0.01;
    }
  if(USE_VARIABLE_ROZO)
    {
      a[++i] = alpha_rozo0;
      start_dev[i] = 0.01;
      a[++i] = B_rozo0;
      start_dev[i] = 0.01;
      a[++i] = sig_rozo0;
      start_dev[i] = 0.01;
    }
  if(USE_MONTE_CARLO_SYSTEMATICS==2)
    {
      a[++i] = 1.0;
      start_dev[i] = 0.03;
    }

  a[1] = 1.278522e+01;
  a[2] = 9.067063e-01;
  a[3] = 1.263600e+01;
  a[4] = 1.670926e+00;
  a[5] = 1.336860e+01;
  a[6] = 9.364447e-01;
  a[7] = 1.262506e+01;
  a[8] = -4.339940e-01;
  a[9] = 9.381532e-01;
  a[10] = 1.384847e+01;
  a[11] = 9.717490e-01;
  a[12] = 1.065310e+01;
  a[13] = -1.469427e-01;
  a[14] = 1.773086e+00;
  a[15] = 2.794929e-01;
  a[16] = 7.964727e-01;
  a[17] = 9.697585e-01;
  a[18] = 7.325640e-01;
  a[19] = 1.050308e+00;
  a[20] = 1.060286e+00;
  // a[21] = 6.261735e-01; // these two are the scale bias now.
  //a[22] = 6.261735e-01;
  a[23] = 7.583900e-01;
  a[24] = 1.103295e+00;
  a[25] = 3.106082e-01;
  a[26] = 1.000134e+00; 


  a[19] = a[20] = 1;

  if(!ThisTask)
    {
      printf("INITIAL VALUES: ");
      for(i=1;i<=wp.ncf;++i)printf("%e ",a[i]);
      printf("\n");
    }

  for(i=1;i<=wp.ncf;++i)
    {
      avg1[i]=a[i];
      for(j=1;j<=wp.ncf;++j)
	cov1[i][j]=a[i]*a[j];
    }

  // put cosmo+other globals in place
  put_parameters_in_place(a,0);
  
  // put HOD params in place, one at a time.
  // use the individual median redshifts for each
  // clustering sample:
  s8original = SIGMA_8;
  chi2a = chi2wp = 0;
  for(i=1;i<=3;++i)
    {
      //reset cosmology for z=0.1
      if(i==1) { SIGMA_8 = s8original*growthfactor(0.068); REDSHIFT=0.068; }
      if(i==2) { SIGMA_8 = s8original*growthfactor(0.104); REDSHIFT=0.104; }
      if(i==3) { SIGMA_8 = s8original*growthfactor(0.126); REDSHIFT=0.126; }

      //FAKE_DATA
      //SIGMA_8 = s8original*growthfactor(0.25);

      RESET_COSMOLOGY++;
      
      put_parameters_in_place(a,i);
      chi2wp += m2n_chi2_wp_wrapper(a,i);
      

      // reset cosmology for z=0.25
      REDSHIFT = 0.25;
      SIGMA_8 = s8original*growthfactor(0.25);
      printf("S8 %f %f\n",s8original,SIGMA_8);
      RESET_COSMOLOGY++;
      chi2a = chi2_m2n(i); // only third call matters
    }
  if(COBENORM)
    {
      SIGMA_8 = s8original;
      RESET_COSMOLOGY++;
      chi2 = chi2 + cobe_prior(OMEGA_M);
    }
  chi2 = chi2a+chi2wp+chi2cmb;
  if(USE_INTERNAL_PRIORS)
    chi2 = chi2 + internal_priors();

  if(!ThisTask) {
    printf("TRY 0 ");
    for(i=1;i<=wp.ncf;++i)
      printf("%.4e ",a[i]);
    printf("%e\n",chi2wp+chi2a);fflush(stdout);
    printf("INITIAL CHI2: %e\n",chi2);
    fflush(stdout);
  }
  return(chi2);
}

/* get the input data.
 * --------------------
 *
 * M2N data. (first off)
 *
 * NB-- the radius is 200rho_crit(z) = 200rho_crit(0)*H^2(z), then
 *   converted to comoving by multiply by 1+z (z=0.25 for this sample).
 */

void m2n_input()
{
  int i,ifile,j,i1;
  FILE *fp;
  char aa[1000];
  double mag, m2n_new;
  float x1; // dummy for mass error
  float *ngals, *err_m, *err_n;
  double **tmp,**tmp2;

  float correction_factor[6] = { 1.334, 1.18, 1.141, 1.098, 1.074, 1.011 };
  float correction_error[6] = { 0.057, 0.038, 0.037, 0.0097, 0.004, 0.003 };
  double err;

  for(ifile=1;ifile<=3;++ifile)
    {

      //FAKE_DATA
      if(ifile==2)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/L384_MOCKS/DATA/M2N_20.5.dat"); 
      if(ifile==1)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/L384_MOCKS/DATA/M2N_19.5.dat"); 
      if(ifile==3)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/L384_MOCKS/DATA/M2N_21.0.dat"); 

      //FAKE_DATA
      if(ifile==2)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/MILLENNIUM_DATA/M2N_20.5.dat"); 
      if(ifile==1)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/MILLENNIUM_DATA/M2N_19.5.dat"); 
      if(ifile==3)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/MILLENNIUM_DATA/M2N_21.0.dat"); 

      //PRODUCTION RUNS
      if(ifile==2)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/DENSITY_PROFILES/worktemp/M2N_z0threshold_20.5_ngals10.dat"); 
      if(ifile==1)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/DENSITY_PROFILES/worktemp/M2N_z0threshold_19.5_ngals10.dat"); 
      if(ifile==3)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/DENSITY_PROFILES/worktemp/M2N_z0threshold_21.0_ngals10.dat"); 

      fp = openfile(M2N.m2n_filename);
      M2N.ndata_i[ifile] = filesize(fp);
      M2N.ndata += M2N.ndata_i[ifile];
      fclose(fp);
    }

  M2N.mass = vector(1,M2N.ndata);
  M2N.model = dvector(1,M2N.ndata);
  M2N.radius = vector(1,M2N.ndata);
  M2N.m2n = vector(1,M2N.ndata);
  M2N.err = vector(1,M2N.ndata);
  M2N.Ngals_lo = ivector(1,M2N.ndata);
  M2N.Ngals_hi = ivector(1,M2N.ndata);
  
  M2N.model_mass = dvector(1,M2N.ndata);
  M2N.model_m2n = dvector(1,M2N.ndata);
  
  err_m = vector(1,M2N.ndata);
  err_n = vector(1,M2N.ndata);
  ngals = vector(1,M2N.ndata);

  i = 0;
  for(ifile=1;ifile<=3;++ifile)
    {
      //FAKE_DATA
      if(ifile==2)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/L384_MOCKS/DATA/M2N_20.5.dat"); 
      if(ifile==1)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/L384_MOCKS/DATA/M2N_19.5.dat"); 
      if(ifile==3)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/L384_MOCKS/DATA/M2N_21.0.dat"); 

      //FAKE_DATA
      if(ifile==2)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/MILLENNIUM_DATA/M2N_20.5.dat"); 
      if(ifile==1)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/MILLENNIUM_DATA/M2N_19.5.dat"); 
      if(ifile==3)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/MILLENNIUM_DATA/M2N_21.0.dat"); 

      //PRODUCTION RUNS -- REAL DATA
      if(ifile==2)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/DENSITY_PROFILES/worktemp/M2N_z0threshold_20.5_ngals10.dat"); 
      if(ifile==1)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/DENSITY_PROFILES/worktemp/M2N_z0threshold_19.5_ngals10.dat"); 
      if(ifile==3)
	sprintf(M2N.m2n_filename,"/data/tinker/SDSS/DENSITY_PROFILES/worktemp/M2N_z0threshold_21.0_ngals10.dat"); 


      fp = openfile(M2N.m2n_filename);
      for(j=1;j<=M2N.ndata_i[ifile];++j) 
	{
	  ++i;
	  fscanf(fp,"%f %f %f %f %d %d %f %f %f",&M2N.mass[i], &err_m[i], &M2N.m2n[i], &M2N.err[i], 
		 &M2N.Ngals_lo[i], &M2N.Ngals_hi[i], &M2N.Ngals_mean[i], &ngals[i], &err_n[i]);
	  fgets(aa,1000,fp);

	  //printf("M2N_CORR: %d %e %e ",i,M2N.m2n[i], M2N.err[i]);
	  
	  // get the error from miscentering
	  err = M2N.m2n[i]*correction_error[(i-1)%6];
	  M2N.err[i] = sqrt(M2N.err[i]*M2N.err[i] + err*err);
	  
	  // add in the miscentering correction
	  M2N.m2n[i] *= correction_factor[(i-1)%6];
	  
	  printf("%e %e %f\n",M2N.m2n[i], M2N.err[i], correction_factor[(i-1)%6]);

	  //printf("FRAC %d %e %e\n",i,M2N.err[i]/M2N.m2n[i],err_m[i]/M2N.mass[i]);

	  //change the HOD values by evolution z=0.25 to z=0.1
	  //FAKE_DATA
	  //m2n_new = M2N.m2n[i]*(1 - 0.1332 + 0.06398*log10(ngals[i]));

	  // PRODUCTION RUNS-- no evolution on HOD (add as error later)
	  m2n_new = M2N.m2n[i];

	  //add 10% error for evolution of luminosity function
	  // M2N.err[i] = m2n_new*sqrt(pow(M2N.err[i]/M2N.m2n[i],2.0) + 0.05*0.05);

	  // PRODCTION RUNS-- add 18% error to MASSES
	  // no-- it's done in the M/N extraction
	  //M2N.m2n[i] = m2n_new*1.18;

	  //M2N.radius[i] = pow(3*M2N.mass[i]/(4*PI*200*RHO_CRIT*H2OFZ),THIRD)*1.25;
	  M2N.radius[i] = pow(3*M2N.mass[i]/(4*PI*RHO_CRIT*DELTA_HALO*OMEGA_M),THIRD);

	  //FAKE_DATA (but no longer used anyway)
	  //M2N.radius[i] = pow(3*M2N.mass[i]/(4*PI*507.0*RHO_CRIT*OMEGA_M),THIRD);
	}
      fprintf(stderr,"Done reading [%d] lines from [%s]\n",M2N.ndata_i[ifile],M2N.m2n_filename);
      fclose(fp);
    }

  // Construct the covariance matrix using total differential error
  M2N.covar = dmatrix(1,M2N.ndata,1,M2N.ndata);
  
  for(i=1;i<=M2N.ndata;++i)
    for(j=1;j<=M2N.ndata;++j)
      M2N.covar[i][j] = 0;

  //diagonals
  for(i=1;i<=M2N.ndata;++i)
    M2N.covar[i][i] = M2N.err[i]*M2N.err[i];

  //off-diagonals, 1->2 and 1->3
  for(i=1;i<=M2N.ndata_i[1];++i)
    {
      j=i+M2N.ndata_i[1];
      M2N.covar[i][j] = err_m[i]*err_m[i]/(ngals[i]*ngals[j]);
      j=i+M2N.ndata_i[1]+M2N.ndata_i[2];
      if(j<=M2N.ndata) //-21 has only 5 points
	M2N.covar[i][j] = err_m[i]*err_m[i]/(ngals[i]*ngals[j]);
    }
  //off-diagonals, 2->3
  for(i1=1;i1<=M2N.ndata_i[2];++i1)
    {
      i=i1+M2N.ndata_i[1];
      j=i1+M2N.ndata_i[1]+M2N.ndata_i[2];
      if(j<=M2N.ndata) //-21 has only 5 points
	M2N.covar[i][j] = err_m[i]*err_m[i]/(ngals[i]*ngals[j]);
    }

  //make matrix symmetric
  
  for(i=1;i<=M2N.ndata;++i)
    for(j=1;j<=M2N.ndata;++j)
      {
	if(i<=j)continue;
	M2N.covar[i][j] = M2N.covar[j][i];
      }
  
  //print out (temp)
  /*
  for(i=1;i<=M2N.ndata;++i)
    for(j=1;j<=M2N.ndata;++j)
      printf("COVAR %d %d %e\n",i,j,M2N.covar[i][j]);
  exit(0);
  */

  //invert
  tmp=dmatrix(1,M2N.ndata,1,1);
  tmp2=dmatrix(1,M2N.ndata,1,M2N.ndata);
  for(i=1;i<=M2N.ndata;++i)
    for(j=1;j<=M2N.ndata;++j)
      tmp2[i][j]=M2N.covar[i][j];
  gaussj(tmp2,M2N.ndata,tmp,1);
  for(i=1;i<=M2N.ndata;++i)
    for(j=1;j<=M2N.ndata;++j)
      M2N.covar[i][j]=tmp2[i][j];
  free_dmatrix(tmp,1,M2N.ndata,1,1);
  free_dmatrix(tmp2,1,M2N.ndata,1,M2N.ndata);

  // Now get the correlation functions.
  for(ifile=1;ifile<=3;++ifile)
    {
      if(ifile==1) mag = 19.5;
      if(ifile==2) mag = 20.5;
      if(ifile==3) mag = 21.0;

      //FAKE_DATA
      COVAR = 0;
      sprintf(wp.fname_wp,"/data/tinker/SDSS/L384_MOCKS/DATA/wp_%.1f.dat",mag);
      sprintf(wp.fname_wp,"/data/tinker/SDSS/MILLENNIUM_DATA/wp_%.1f.dat",mag);

      // PRODUCTION RUNS
      COVAR = 1;
      sprintf(wp.fname_wp,"/data/tinker/SDSS/DR7_DATA/wp_%.1f.dat",mag);
      sprintf(wp.fname_covar,"/data/tinker/SDSS/DR7_DATA/wp_covar_%.1f.dat",mag);
      
      // TRUNCATED PRODUCTION RUNS
      COVAR = 1;
      sprintf(wp.fname_wp,"/data/tinker/SDSS/DR7_DATA/truncated_data/wp_%.1f.dat",mag);
      sprintf(wp.fname_covar,"/data/tinker/SDSS/DR7_DATA/truncated_data/wp_covar_%.1f.dat",mag);
      

      wp_input();

      M2N.wpn[ifile] = wp.np;
      M2N.wpr[ifile] = dvector(1,wp.np);
      M2N.wpx[ifile] = dvector(1,wp.np);
      M2N.wpe[ifile] = dvector(1,wp.np);
      M2N.wpcovar[ifile] = dmatrix(1,wp.np,1,wp.np);

      // put the data in place
      for(i=1;i<=wp.np;++i)
	{
	  M2N.wpr[ifile][i] = wp.r[i];
	  M2N.wpx[ifile][i] = wp.x[i];
	  M2N.wpe[ifile][i] = wp.e[i];
	}

      // invert the covariance matrix
      // if FAKE_DATA, then comment out this

      tmp=dmatrix(1,wp.np,1,1);
      tmp2=dmatrix(1,wp.np,1,wp.np);
      for(i=1;i<=wp.np;++i)
	for(j=1;j<=wp.np;++j)
	  tmp2[i][j]=wp.covar[i][j];
      gaussj(tmp2,wp.np,tmp,1);
      for(i=1;i<=wp.np;++i)
	for(j=1;j<=wp.np;++j)
	  M2N.wpcovar[ifile][i][j]=tmp2[i][j];
      free_dmatrix(tmp,1,wp.np,1,1);
      free_dmatrix(tmp2,1,wp.np,1,wp.np);

      // keep the allocations for the last wpinput call
      if(ifile<3)
	{
	  free_dvector(wp.r,1,wp.np);
	  free_dvector(wp.x,1,wp.np);
	  free_dvector(wp.e,1,wp.np);
	  // if FAKE_DATA, comment this line out
	  free_dmatrix(wp.covar,1,wp.np,1,wp.np);
	}
    }


}


/* Calculate the chi2 for the M2N data.
 * ----------------------------------------------------------------------------
 *
 * For each Ngals, calculate the mean halo mass and the mean number of galaxies.
 *
 *  M_bar(N_gals) = \int dM dn/dM P(N_gals|Mass) Mass
 *     divided by  \int dM dn/dM P(N_gals|Mass) 
 *
 *  N_bar(N_gals) = \int dM dn/dM P(N_gals|Mass) N_true(Mass)
 *     divided by  \int dM dn/dM P(N_gals|Mass) 
 *
 * Where P(N_true|Mass) is the Poisson distribution for satellites + nearest int for centrals.
 *
 * Where P(N_gals|N_true) is something we have yet to determine.
 *
 * For bins in Ngals that are wider than one, do a weighted sum between all the 
 * values of Ngals. 
 *
 * How many systems do you find at a fixed Ngals?
 *
 *   n(N_gals) = \int dM dn/dM P(N_gals|Mass)
 *
 */
double chi2_m2n(int icall)
{
  int i,j,k,n, ifirst;
  double mbar, nbar, nsys, m2n, chi2, x1, vol, x2, x3;
  static int iter=0;

  chi2 = 0;

  // set the first value by the magnitude of the sample
  if(icall==1)ifirst = 1;
  if(icall==2)ifirst = 1 + M2N.ndata_i[1];
  if(icall==3)ifirst = 1 + M2N.ndata_i[1] + M2N.ndata_i[2];
  for(M2N.current_bin=1,i=ifirst;i<=ifirst+M2N.ndata_i[icall]-1;++i,++M2N.current_bin)
    {
      mbar = nbar = nsys = x1 = 0;
      for(j=M2N.Ngals_lo[i];j<=M2N.Ngals_hi[i];++j)
	{
	  M2N.current_Ngals = j;
	  mbar += qromo(m2n_func1,log(HOD.M_min),log(HOD.M_max),midpnt);
	  nbar += qromo(m2n_func2,log(HOD.M_min),log(HOD.M_max),midpnt);
	  nsys += qromo(m2n_func3,log(HOD.M_min),log(HOD.M_max),midpnt);
	  x1 += qromo(m2n_func1a,log(HOD.M_min),log(HOD.M_max),midpnt);
	}
      m2n = mbar/nbar;
      if(USE_MONTE_CARLO_SYSTEMATICS)
	m2n*=SYSTEMATIC_OFFSET;

      //FAKE_DATA
      //m2n = M2N.mass[i]/N_sat(M2N.mass[i]);

      M2N.model_m2n[i] = m2n;
      M2N.model_mass[i] = mbar/nsys;
      x1 /= nsys;
      chi2 += (M2N.m2n[i] - m2n)*(M2N.m2n[i] - m2n)/(M2N.err[i]*M2N.err[i]);
      M2N.model[i] = m2n;
      x2 = (exp(0.95)*pow(M2N.Ngals_mean[i]/40.,1.06)*1e14*HUBBLE);
      x3 = M2N.model_mass[i]/M2N.mass[i];
      //convert that mass to 500crit from 200crit.
      /*
      DELTA_HALO*=2.5;
      RESET_COSMOLOGY++;
      x2 = halo_mass_conversion2(x2,halo_concentration(x2),DELTA_HALO,507.0);
      x3 = M2N.model_mass[i]/x2;
      DELTA_HALO/=2.5;
      RESET_COSMOLOGY++;
      */
      printf("CHIM2N%d %d %e %e %e %e %e %e %e %e\n",iter,i,M2N.m2n[i],M2N.err[i],m2n,M2N.mass[i],M2N.model_mass[i],chi2,x2,x3);
    }
  //exit(0);

  // if we're at the third call, then do the total covariance matrix
  if(icall<3)return 0;

  chi2 = 0;
  for(i=1;i<=M2N.ndata;++i)
    for(j=1;j<=M2N.ndata;++j)
	chi2+=(M2N.model[i]-M2N.m2n[i])*(M2N.model[j]-M2N.m2n[j])*M2N.covar[i][j];
  printf("CHI2M2Nx%d %e\n",iter,chi2);

  iter++;
  return chi2;
}


/* THis is the integrand of the mean mass integral:
 * dM dn/dM P(N_gals|M) Mass
 *
 * The PDF of P(M|N_gals) is a lognormal with 
 * mean = <M|Ngals> = exp(B) pow(N_gals/40,alpha)
 * sigma = 0.48
 *
 * with B = -0.12
 * with alpha = 1.15
 */
double m2n_func1(double m)
{
  int i;
  double x, logm, mu, sig, c, rvir, rs, rhos, mtrue, logj;

  logm = m;
  m = exp(m);
  i = M2N.current_bin;
  logj = log(M2N.current_Ngals);

  sig = sig_rozo;
  mu = B_rozo - sig*sig/2 + alpha_rozo*(logm-log(1.0e14*HUBBLE)) + log(40.0);
  x = exp(-(mu - logj)*(mu - logj)/(2*sig*sig))/(RT2PI*sig)/m;

  mu = log(exp(B_rozo)*pow(M2N.current_Ngals/40.,alpha_rozo)*1e14*HUBBLE)- sig*sig/2;
  x = exp(-(mu-logm)*(mu-logm)/(2*sig*sig))/(RT2PI*sig)/m;

  mu = B_rozo + alpha_rozo*(logm-log(PIVOT_MASS));
  x = exp(-(mu - logj)*(mu - logj)/(2*sig*sig))/(RT2PI*sig)/M2N.current_Ngals;

  // extrapolate the HOD profile out/in to the mean radius of the bin
  c = halo_concentration(m);
  rvir = pow(3*m/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD);
  rs = rvir/c;
  rhos = m/(4*PI*rvir*rvir*rvir*HK_func(rs/rvir));
  mtrue = 4*PI*rhos*pow(M2N.radius[i],3.0)*HK_func(rs/M2N.radius[i]);

  //FAKE_FATA
  //mtrue = m;

  // printf("%e %e %e %e %e\n",m,mtrue,x,mu);

  return mtrue*m*dndM_interp(m)*x;
}

double m2n_func1a(double m)
{
  int i;
  double x, logm, mu, sig, c, rvir, rs, rhos, mtrue, logj;

  logm = m;
  m = exp(m);
  i = M2N.current_bin;  
  logj = log(M2N.current_Ngals);


  sig = sig_rozo;
  mu = B_rozo - sig*sig/2 + alpha_rozo*(logm-log(1.0e14*HUBBLE)) + log(40.0);
  x = exp(-(mu - logj)*(mu - logj)/(2*sig*sig))/(RT2PI*sig)/m;

  return m*m*dndM_interp(m)*x;
}


/* THis is the integrand of the mean number integral:
 * dM dn/dM P(N_gals|M) N_true(M)
 *
 * NB! need to do a correction for the fact that the radius for N
 * isn't R200 exactly for all masses-- it's R200 for the mean mass in the bin.
 *
 */
double m2n_func2(double m)
{
  int i;
  double x, logm, mu, sig, c, rvir, rs, rhos, nsat, logj;

  logm = m;
  m = exp(m);
  i = M2N.current_bin;
  logj = log(M2N.current_Ngals);

  //mu = exp(-0.12)*pow(M2N.current_Ngals/40.0,1.15);
  //mu = B_rozo + alpha_rozo*log(M2N.current_Ngals/40.0) + log(4e14*HUBBLE);
  sig = sig_rozo;
  mu = B_rozo - sig*sig/2 + alpha_rozo*(logm-log(1.0e14*HUBBLE)) + log(40.0);
  x = exp(-(mu - logj)*(mu - logj)/(2*sig*sig))/(RT2PI*sig)/m;

  mu = log(exp(B_rozo)*pow(M2N.current_Ngals/40.,alpha_rozo)*1e14*HUBBLE)- sig*sig/2;
  x = exp(-(mu-logm)*(mu-logm)/(2*sig*sig))/(RT2PI*sig)/m;

  mu = B_rozo + alpha_rozo*(logm-log(PIVOT_MASS));
  x = exp(-(mu - logj)*(mu - logj)/(2*sig*sig))/(RT2PI*sig)/M2N.current_Ngals;

  //FAKE_DATA
  //return m*N_sat(m)*dndM_interp(m)*x;

  // extrapolate the HOD profile out/in to the mean radius of the bin
  c = CVIR_FAC*halo_concentration(m);
  rvir = pow(3*m/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD);
  rs = rvir/c;
  rhos = N_sat(m)/(4*PI*rvir*rvir*rvir*HK_func(rs/rvir));
  nsat = 4*PI*rhos*pow(M2N.radius[i],3.0)*HK_func(rs/M2N.radius[i]);
  //printf("%e %e %e %e\n",m,rvir,M2N.radius[i],nsat/N_sat(m));

  return m*nsat*dndM_interp(m)*x;
}

/* THis is the integrand of the number of systems
 * dM dn/dM P(N_gals|M)
 *
 */
double m2n_func3(double m)
{
  int i;
  double x, logm, mu, sig, logj;

  logm = m;
  m = exp(m);
  logj = log(M2N.current_Ngals);

  //mu = exp(-0.12)*pow(M2N.current_Ngals/40.0,1.15);
  sig = sig_rozo;
  mu = B_rozo - sig*sig/2 + alpha_rozo*(logm-log(1.0e14*HUBBLE)) + log(40.0);
  x = exp(-(mu - logj)*(mu - logj)/(2*sig*sig))/(RT2PI*sig)/m;

  mu = log(exp(B_rozo)*pow(M2N.current_Ngals/40.,alpha_rozo)*1e14*HUBBLE)- sig*sig/2;
  x = exp(-(mu-logm)*(mu-logm)/(2*sig*sig))/(RT2PI*sig)/m;
  //printf("%e %e %e %f %f\n",m,dndM_interp(m),x,mu,logm);

  mu = B_rozo + alpha_rozo*(logm-log(PIVOT_MASS));
  x = exp(-(mu - logj)*(mu - logj)/(2*sig*sig))/(RT2PI*sig)/M2N.current_Ngals;


  return m*dndM_interp(m)*x;
}




/*
 * Takes the vector of variables and places the values into the proper globals
 * - the "icall" parameter tells which threshold we're using. icall==0 means non-HOD parameters.
 */
void put_parameters_in_place(double *a, int icall)
{
  static int flag=1;
  static double s8original;
  int i;
  
  if(flag)
    {
      s8original = SIGMA_8;
      flag = 0;
    }

  if(icall==0) 
    {    
      i = 14;

      //FAKE_DATA 
      // ??? what do i do with this? I need to keep this, obviously, so why is it labeled FAKE?
      OMEGA_M = a[++i];
      SIGMA_8 = a[++i];
      SPECTRAL_INDX = a[++i];
      HUBBLE = a[++i];
      RESET_COSMOLOGY++;

      //FAKE_DATA!!!! -- constant cosmology
      //SIGMA_8 = s8original;


      if(USE_ERRORS)
	{
	  M2N.bias_amp = a[++i];
	  M2N.mf_amp = a[++i];
	  //M2N.scalebias_amp = a[++i];
	  M2N.bias_pow1 = a[++i];
	  M2N.bias_amp1 = a[++i];
	}
      if(USE_VARIABLE_ROZO)
	{
	  alpha_rozo0 = a[++i];
	  B_rozo0 = a[++i];
	  sig_rozo0 = a[++i];
	}
      if(USE_MONTE_CARLO_SYSTEMATICS==2)
	{
	  SYSTEMATIC_OFFSET = a[++i];
	}
      return;
    }

  if(icall==1) //-19.5
    {
      i = 0;
      HOD.M1 = a[++i];
      HOD.alpha = a[++i];
      HOD.M_cut = a[++i];
      CVIR_FAC = a[++i];
      return;
    }
  if(icall==2) //-20.5
    {
      i = 4;
      HOD.M1 = a[++i];
      HOD.alpha = a[++i];
      HOD.M_cut = a[++i];
      HOD.sigma_logM = a[++i];
      CVIR_FAC = a[++i];
      return;
    }
  if(icall==1) //-19.5
    {
      i = 9;
      HOD.M1 = a[++i];
      HOD.alpha = a[++i];
      HOD.M_cut = a[++i];
      HOD.sigma_logM = a[++i];
      CVIR_FAC = a[++i];
      return;
    }
}


double internal_priors()
{
  static int niter=0;
  double   chifac = 0;

  chifac += ((M2N.mf_amp-1.0)*(M2N.mf_amp-1.0)/
	     (mf_amp_err*mf_amp_err));
  chifac += ((M2N.bias_amp-1.0)*(M2N.bias_amp-1.0)/
	     (bias_amp_err*bias_amp_err));
  //comment out for selfcalibration
  // chifac += ((M2N.scalebias_amp-1.0)*(M2N.scalebias_amp-1.0)/
  //	     (scalebias_amp_err*scalebias_amp_err));

  chifac += pow((alpha_rozo-alpha_rozo0)/alpha_rozo_err,2.0);
  chifac += pow((B_rozo-B_rozo0)/B_rozo_err,2.0);
  chifac += pow((sig_rozo-sig_rozo0)/sig_rozo_err,2.0);

  printf("CHIFAC%d %e\n",++niter,chifac);

  return chifac;
}
