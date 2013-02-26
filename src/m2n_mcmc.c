#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

// parameters of Eduardo's P(M|N200)
/*
#define B_rozo -0.12596
#define alpha_rozo 1.149
*
#define B_rozo -1.1
#define alpha_rozo 0.74
#define sig_rozo 0.35
*/

/* External functions from wp_minimization.c
 */
void wp_input(void);
double poisson_prob(int n, double nave);

/* Internal functions.
 */
double m2n_chi2_wp_wrapper(double *a);
double m2n_initialize(double *a, double **cov1, double *avg1, double *start_dev);
void m2n_input(void);
double chi2_m2n(void);
double chi2_number_profiles(void);
int parameter_out_of_range(double *a);

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
double B_rozo = -1.02, B_rozo0 = -1.02, B_rozo_err = 0.15;
double alpha_rozo = 0.74, alpha_rozo0 = 0.74, alpha_rozo_err = 0.1;
double sig_rozo = 0.35, sig_rozo0 = 0.35, sig_rozo_err = 0.1;
int USE_VARIABLE_ROZO = 1;

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
  double t0,tprev,temp,chi2a,chi2b;

  double delta_halo_rhoc = 200;

  int icvir;

  //test_pdf();

  // initialize use of SYSTEMATIC ERROR
  USE_ERRORS = 1;
  M2N.IDUM = -5555;
  
  // fix the value of OEMGA_M fot T(k)
  M2N.fix_omegam = 0;
  M2N.constant_omegam = 0.27;

  opar=dvector(1,100);

  MCMC=Task.MCMC;

  pcheck=calloc(ptot,sizeof(int));

  /* Since we're at constant halo overdensity wrt RHO_CRIT,
   * set the overdensity of the halo
   */
  H2OFZ = 0.27*pow(1.25,3.0)+0.73;
  DELTA_HALO = delta_halo_rhoc/OMEGA_M*(OMEGA_M*1.25*1.25*1.25 + (1-OMEGA_M))/(1.25*1.25*1.25);

  /* read in the wp data for a single luminosity threshold sample.
   */
  wp_input();

  /* read in the M2N data.
   */
  m2n_input();
  //input_maxbcg_counts();

  Work.imodel=2;
  Work.chi2=1;


  srand48(32498793);

  /* Find the number of free parameters in the minimization
   * for the real-space correlation function.
   */
  for(n=0,i=1;i<100;++i)
    {
      n+=HOD.free[i];
      /* if(i>N_HOD_PARAMS && HOD.free[i])MCMC=3;*/
      if(OUTPUT)
	printf("mcmc_min> free[%i] = %d\n",i,HOD.free[i]);
    }
  if(USE_ERRORS)
    n+=3;
  if(USE_VARIABLE_ROZO)
    n+=3;
  wp.ncf=n;

  /* Find out which free parameter is for CVIR_FAC
   */
  j=0;
  if(HOD.free[6])
    for(i=0;i<6;++i)
      if(HOD.free[i])j++;
  icvir=j+1;

  if(HOD.free[0])
    {
      wp.ngal = GALAXY_DENSITY;
      wp.ngal_err = 0.1*wp.ngal;
      FIX_PARAM = 0;
    }

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
      stepfac=0.7;
      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];

      
      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
	  i=N_HOD_PARAMS;
	  if(HOD.free[++i])OMEGA_M         = a[++j];
	  if(HOD.free[++i])SIGMA_8         = a[++j];
	  if(HOD.free[++i])VBIAS           = a[++j];
	  if(HOD.free[++i])VBIAS_C         = a[++j];
	  if(HOD.free[++i])GAMMA           = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];
	  if(HOD.free[++i])HUBBLE          = a[++j];
	}
      if(VBIAS_C<0)continue;

      if(USE_ERRORS)
	{
	  M2N.mf_amp = a[++j];
	  M2N.bias_amp = a[++j];
	  M2N.scalebias_amp = a[++j];
	}
      if(USE_VARIABLE_ROZO)
	{
	  alpha_rozo = a[++j];
	  B_rozo = a[++j];
	  sig_rozo = a[++j];
	}

      /* Hard-wire CVIR variation
       */
      if(HOD.free[6])
	CVIR_FAC = a[icvir];

      /* Since we're at constant halo overdensity wrt RHO_CRIT,
       * set the overdensity of the halo
       */
      DELTA_HALO = delta_halo_rhoc/OMEGA_M*(OMEGA_M*1.25*1.25*1.25 + (1-OMEGA_M))/(1.25*1.25*1.25);

      // Check to see if any of our parameters are out of range.
      if(parameter_out_of_range(a)){ muh(1); continue; }

 
      /* Draw random value of cvir from prior.
       */
      /* if(CVIR_FAC<0.3 || CVIR_FAC>1.2)continue; */
      /* CVIR_FAC = 0.9*drand48()+0.3;  */
      /* GAMMA = gasdev(&IDUM)*0.02 + 0.15; */


      chi2=m2n_chi2_wp_wrapper(a);
      // reset cosmology for z=0.25
      SIGMA_8 = SIGMA_8*growthfactor(0.25);
      RESET_COSMOLOGY++;
      if(MCMC>1 && chi2<1.0E7)chi2+= chi2a = chi2_m2n();
      if(MCMC>1 && chi2<1.0E7)chi2+= chi2b = chi2_number_profiles();

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

	  /* For the initialization, don't use this: we need
	   * separate elements for estimating the covariance matrix.
	   */
	  /*
	  for(i=1;i<=n;++i)
	    a[i] = aprev[i];
	  chi2 = chi2prev;
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
	printf("%e %e %e %e\n",chi2,chi2a,chi2b,chi2-chi2a-chi2b);fflush(stdout);
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

      /* We seem to be having a problem with this.
       * So, broadcast the model params from the root processor.
       */
#ifdef PARALLEL      
      MPI_Bcast(&a[1],n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);
#endif

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
	  i=N_HOD_PARAMS;
	  if(HOD.free[++i])OMEGA_M         = a[++j];
	  if(HOD.free[++i])SIGMA_8         = a[++j];
	  if(HOD.free[++i])VBIAS           = a[++j];
	  if(HOD.free[++i])VBIAS_C         = a[++j];
	  if(HOD.free[++i])GAMMA           = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];
	  if(HOD.free[++i])HUBBLE          = a[++j];
	}
      if(VBIAS_C<0)continue;

      if(USE_ERRORS)
	{
	  M2N.mf_amp = a[++j];
	  M2N.bias_amp = a[++j];
	  M2N.scalebias_amp = a[++j];
	}
      if(USE_VARIABLE_ROZO)
	{
	  alpha_rozo = a[++j];
	  B_rozo = a[++j];
	  sig_rozo = a[++j];
	}
      
      /* Hard-wire CVIR variation
       */
      if(HOD.free[6])
	CVIR_FAC = a[icvir];

      /* Since we're at constant halo overdensity wrt RHO_CRIT,
       * set the overdensity of the halo
       */
      DELTA_HALO = delta_halo_rhoc/OMEGA_M*(OMEGA_M*1.25*1.25*1.25 + (1-OMEGA_M))/(1.25*1.25*1.25);

      // Check to see if any of our parameters are out of range.
      if(parameter_out_of_range(a))continue;


      /* Draw random value of cvir from prior.
       */
      /* CVIR_FAC = a[n]; */
      /* if(CVIR_FAC<0.3 || CVIR_FAC>1.2)continue; */
      /* CVIR_FAC = 0.7*drand48()+0.3; */
      /* GAMMA = gasdev(&IDUM)*0.02 + 0.15; */
      //      printf("GAMMA %d %f %f\n",count+1,GAMMA,CVIR_FAC);

      chi2=m2n_chi2_wp_wrapper(a);
      // reset cosmology for z=0.25
      SIGMA_8 = SIGMA_8*growthfactor(0.25);
      RESET_COSMOLOGY++;
      if(MCMC>1 && chi2<1.0E7)chi2 += chi2a = chi2_m2n();
      if(MCMC>1 && chi2<1.0E7)chi2 += chi2b = chi2_number_profiles();

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
      if(0) {
	printf("CPU%02d %d ",ThisTask,count);
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
	  /*
	  for(i=1;i<=n;++i)
	    a[i] = aprev[i];
	  chi2 = chi2prev;
	  */
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
	printf("%e %e %e %e\n",chi2,chi2a,chi2b,chi2-chi2a-chi2b);fflush(stdout);
	
	if(MCMC==1)
	  {
	    printf("HSTATS %d %e %e %e %e\n",niter,HOD.M_min,number_weighted_halo_mass(),
		   number_weighted_central_mass(),
		   qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
	    
	    fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	    printf("FSAT %d %e %e %e %e\n",niter,fsat,HOD.M_min,HOD.sigma_logM,qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
	  }
      }

    }
}


/* This is a routine where i can hand-set some boundaries for the 
 * cosmological parameters (and such)
 */
int parameter_out_of_range(double *a)
{
  int i;
  i=N_HOD_PARAMS;
  if(HOD.free[++i])//OMEGA_M         = a[++j];
    if(OMEGA_M<0.15 || OMEGA_M>0.4) return 1;
  if(HOD.free[++i])//SIGMA_8         = a[++j];
    if(SIGMA_8<0.55 || SIGMA_8>1.1) return 1;
  if(HOD.free[++i])//VBIAS           = a[++j];
    ;
  if(HOD.free[++i])//VBIAS_C         = a[++j];
    ;
  if(HOD.free[++i])//GAMMA           = a[++j];
    ;
  if(HOD.free[++i])//SPECTRAL_INDX   = a[++j];
    if(SPECTRAL_INDX<0.88 || SPECTRAL_INDX>1.05) return 1;
  if(HOD.free[++i])//HUBBLE          = a[++j];
    if(HUBBLE < 0.6 || HUBBLE > 0.8) return 1;
  if(HOD.free[6])//CVIR_FAC
    if(CVIR_FAC<0.2 || CVIR_FAC>2) return 1;
  if(USE_ERRORS)
    {
      if(M2N.bias_amp > 1.10 || M2N.bias_amp < 0.9) return 1;
      if(M2N.mf_amp > 1.10 || M2N.mf_amp < 0.9) return 1;
      if(M2N.scalebias_amp > 1.30 || M2N.scalebias_amp < 0.7) return 1;
    }
  if(USE_VARIABLE_ROZO)
    {
      if(B_rozo > B_rozo0+2*B_rozo_err || B_rozo < B_rozo0-2*B_rozo_err)return 1;
      if(alpha_rozo > alpha_rozo0+2*alpha_rozo_err || alpha_rozo < alpha_rozo0-2*alpha_rozo_err)return 1;
      if(sig_rozo > sig_rozo0+2*sig_rozo_err || sig_rozo < sig_rozo0-2*sig_rozo_err)return 1;
    }
  return 0;
}

double m2n_chi2_wp_wrapper(double *a)
{
  static int flag=1;
  static double *b;
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
    if(HOD.free[i] && i!=5) { 
      if(a[++j]<=0) { printf("NEG %d %d %e\n",i,j,a[j]); return(1.0E7); } }
    if(HOD.free[i] && i==5) {
      ++j; }
  }
  
  /* check to make sure that none of the cosmological 
   * parameters are negative, either.
   */
  for(i=N_HOD_PARAMS+1;i<100;++i)
    if(HOD.free[i])
      if(a[++j]<=0) { printf("NEG %d %d %e\n",i,j,a[j]); return(1.0E7); } 

  i=0;j=0;
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M_min */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M1 */
  if(HOD.free[++i]){j++;b[j]=a[j];}           /* alpha */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M_cut */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* sigma_logM */
  if(HOD.free[++i]){j++;b[j]=a[j];}           /* cvir_fac */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* MaxCen */
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);} /* M_sat_break */
  if(HOD.free[++i]){j++;b[j]=a[j];}           /* alpha1 */

  return(chi2_wp(b));
}

double m2n_initialize(double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0;
  double x1,x2,omega_m;
  long IDUM = -556;

  omega_m = 1;
  //if(MCMC>1)
  //omega_m = OMEGA_M;

  i=0;j=0;
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_min/omega_m);start_dev[j]=0.001; }
  if(HOD.free[++i]){ a[++j]=log10(HOD.M1/omega_m);start_dev[j]=0.001; } //.0005
  if(HOD.free[++i]){ a[++j]=HOD.alpha;start_dev[j]=0.03; } //.005
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_cut/omega_m);start_dev[j]=0.01; } //.001
  if(HOD.free[++i]){ a[++j]=log10(HOD.sigma_logM);start_dev[j]=0.01; }
  if(HOD.free[++i]){ a[++j]=CVIR_FAC;start_dev[j]=0.02; }
  if(HOD.pdfc==7) {
    if(HOD.free[++i])a[++j]=log10(HOD.M_cen_max/omega_m); start_dev[j]=0.001; }
  else {
    if(HOD.free[++i])a[++j]=HOD.MaxCen; start_dev[j]=0.02; }
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_sat_break/omega_m);start_dev[j]=0.001; }
  if(HOD.free[++i]){ a[++j]=HOD.alpha1;start_dev[j]=0.02; }

  if(MCMC>1)
    {
      if(HOD.free[++i])a[++j]=OMEGA_M;
      if(HOD.free[++i])a[++j]=SIGMA_8;
      if(HOD.free[++i])a[++j]=VBIAS;
      if(HOD.free[++i])a[++j]=VBIAS_C;
      if(HOD.free[++i])a[++j]=GAMMA;
       if(HOD.free[++i])a[++j]=SPECTRAL_INDX;
       if(HOD.free[++i])a[++j]=HUBBLE;
    }
  if(USE_ERRORS)
    {
      a[++j]=1;
      a[++j]=1;
      a[++j]=1;
    }
  if(USE_VARIABLE_ROZO)
    {
      a[++j] = alpha_rozo0;
      a[++j] = B_rozo0;
      a[++j] = sig_rozo0;
    }

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

  if(MCMC>1)
    {
      RESET_COSMOLOGY++;
      j=0;
      for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
      i=N_HOD_PARAMS;
      if(HOD.free[++i]){ OMEGA_M = a[++j]; start_dev[j] = 0.01; } 
      if(HOD.free[++i]){ SIGMA_8 = a[++j]; start_dev[j] = 0.01; } 
      if(HOD.free[++i]){ VBIAS   = a[++j]; start_dev[j] = 0.01; } 
      if(HOD.free[++i]){ VBIAS_C = a[++j]; start_dev[j] = 0.02; } 
      if(HOD.free[++i]){ GAMMA   = a[++j]; start_dev[j] = 0.015; } 
      if(HOD.free[++i]){ SPECTRAL_INDX    = a[++j]; start_dev[j] = 0.02; }
      if(HOD.free[++i]){ HUBBLE    = a[++j]; start_dev[j] = 0.02; }
    }

  if(USE_ERRORS)
    {
      M2N.bias_amp = a[++j]; start_dev[j] = 0.01;
      M2N.mf_amp = a[++j]; start_dev[j] = 0.01;
      M2N.scalebias_amp = a[++j]; start_dev[j] = 0.02;
    }
  if(USE_VARIABLE_ROZO)
    {
      alpha_rozo = a[++j]; start_dev[j] = 0.01;
      B_rozo = a[++j]; start_dev[j] = 0.01; 
      sig_rozo = a[++j]; start_dev[j] = 0.01;
    }

  x1=m2n_chi2_wp_wrapper(a);
  // reset cosmology for z=0.25
  SIGMA_8 = SIGMA_8*growthfactor(0.25);
  RESET_COSMOLOGY++;
  if(MCMC>1 && x1<1.0E7)x1+=chi2_m2n();
  if(MCMC>1 && x1<1.0E7)x1+=chi2_number_profiles();
  
  if(!ThisTask) {
    printf("TRY 0 ");
    for(i=1;i<=wp.ncf;++i)
      printf("%.4e ",a[i]);
    printf("%e\n",x1+x2);fflush(stdout);
    printf("INITIAL CHI2: %e\n",x1);
    fflush(stdout);
  }
  return(x1);
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
  int i;
  FILE *fp;
  char aa[1000];
  float x1; // dummy for mass error

  sprintf(M2N.m2n_filename,"/home/tinker/SDSS/DENSITY_PROFILES/M2N_erin_ngals10.threshold_0.dat"); //-19.5
  //sprintf(M2N.m2n_filename,"/home/tinker/SDSS/DENSITY_PROFILES/M2N_erin_ngals10.threshold_3.dat"); //-20.5
  //sprintf(M2N.m2n_filename,"/home/tinker/SDSS/DENSITY_PROFILES/M2N_erin_ngals10.threshold_5.dat"); //-21.0

  sprintf(M2N.m2n_filename,"/Users/tinker/cosmo/MN_RATIO/DENSITY_PROFILES/M2N_z0threshold_19.5_ngals10.dat"); 

  fp = openfile(M2N.m2n_filename);
  M2N.ndata = filesize(fp);

  M2N.mass = vector(1,M2N.ndata);
  M2N.radius = vector(1,M2N.ndata);
  M2N.m2n = vector(1,M2N.ndata);
  M2N.err = vector(1,M2N.ndata);
  M2N.Ngals_lo = ivector(1,M2N.ndata);
  M2N.Ngals_hi = ivector(1,M2N.ndata);

  M2N.model_mass = dvector(1,M2N.ndata);
  M2N.model_m2n = dvector(1,M2N.ndata);

  for(i=1;i<=M2N.ndata;++i) 
    {
      fscanf(fp,"%f %f %f %f %d %d",&M2N.mass[i], &x1, &M2N.m2n[i], &M2N.err[i], 
	     &M2N.Ngals_lo[i], &M2N.Ngals_hi[i]);
      fgets(aa,1000,fp);
      M2N.radius[i] = pow(3*M2N.mass[i]/(4*PI*200*RHO_CRIT*H2OFZ),THIRD)*1.25;
    }
  

  fprintf(stderr,"Done reading [%d] lines from [%s]\n",M2N.ndata,M2N.m2n_filename);
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
double chi2_m2n()
{
  int i,j,k,n;
  double mbar, nbar, nsys, m2n, chi2, x1, vol;
  static int iter=0;

  chi2 = 0;

  printf("%f %f %f\n",alpha_rozo, B_rozo, sig_rozo);

  // calculate the number densities of the clusters for the current cosmology
  /*
  vol = comoving_volume(0.1,0.3);
  for(i=1;i<=M2N.counts_nbins;++i)
    M2N.ndens_N200[i] = M2N.counts_N200[i]/vol;
  */
  for(i=1;i<=M2N.ndata;++i)
    {
      mbar = nbar = nsys = x1 = 0;
      M2N.current_bin = i;
      for(j=M2N.Ngals_lo[i];j<=M2N.Ngals_hi[i];++j)
	{
	  M2N.current_Ngals = j;
	  mbar += qromo(m2n_func1,log(HOD.M_min),log(HOD.M_max),midpnt);
	  nbar += qromo(m2n_func2,log(HOD.M_min),log(HOD.M_max),midpnt);
	  nsys += qromo(m2n_func3,log(HOD.M_min),log(HOD.M_max),midpnt);
	  x1 += qromo(m2n_func1a,log(HOD.M_min),log(HOD.M_max),midpnt);
	}
      m2n = mbar/nbar;
      M2N.model_m2n[i] = m2n;
      M2N.model_mass[i] = mbar/nsys;
      x1 /= nsys;
      chi2 += (M2N.m2n[i] - m2n)*(M2N.m2n[i] - m2n)/(M2N.err[i]*M2N.err[i]);
      
      if(OUTPUT)
	printf("CHIM2N%d %d %e %e %e %e %e %e %e %e\n",iter,i,M2N.m2n[i],M2N.err[i],m2n,M2N.mass[i],M2N.model_mass[i],chi2,x1,exp(-0.12)*pow(M2N.Ngals_lo[i]/40.,1.15)*4e14*HUBBLE);
      //exit(0);
    }
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

  // extrapolate the HOD profile out/in to the mean radius of the bin
  c = halo_concentration(m);
  rvir = pow(3*m/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD);
  rs = rvir/c;
  rhos = m/(4*PI*rvir*rvir*rvir*HK_func(rs/rvir));
  mtrue = 4*PI*rhos*pow(M2N.radius[i],3.0)*HK_func(rs/M2N.radius[i]);

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
  //printf("%e %e %e %f %f\n",m,dndM_interp(m),x,mu,logm);

  return m*dndM_interp(m)*x;
}




double chi2_number_profiles()
{
  double chi2, c, rvir, rhos, rho, rr, m, rs; 
  float x1;
  int i, j, ibin[6] = { 6, 7, 8, 9, 10, 11 };
  char fname[1000];
  FILE *fp;

  static int flag = 1, nbins, nrad=14, iter=0;
  static float **ngal, **rad, **err;

  return 0;

  if(flag)
    {
      nbins = M2N.ndata;
      ngal = matrix(1,nbins,1,nrad);
      rad = matrix(1,nbins,1,nrad);
      err = matrix(1,nbins,1,nrad);
      flag = 0;

      for(i=1;i<=nbins;++i)
	{
	  sprintf(fname,"/home/tinker/SDSS/DENSITY_PROFILES/ngals3D_bin%02d.dat",ibin[i-1]);
	  fp = openfile(fname);
	  for(j=1;j<=nrad;++j)
	    {
	      fscanf(fp,"%f %f %f %f",&rad[i][j],&ngal[i][j],&err[i][j],&x1);
	      rad[i][j] /= 1.25; // put in comoving units
	      ngal[i][j]/=(M2N.mass[i]/M2N.m2n[i]); //normalize to unity
	      err[i][j]/=(M2N.mass[i]/M2N.m2n[i]); //normalize to unity
	    }
	  fprintf(stderr,"Done reading [%d] lines from [%s]\n",nrad,fname);
	}
    }

  for(i=1;i<=M2N.ndata;++i)
    {
      x1 = 0;
      m = M2N.model_mass[i];
      c = CVIR_FAC*halo_concentration(m);
      rvir = pow(3*m/(4*PI*RHO_CRIT*OMEGA_M*DELTA_HALO),THIRD);
      rs = rvir/c;
      rhos = 1/(4*PI*rvir*rvir*rvir*HK_func(1/c)); // normalize to unity
      for(j=1;j<=nrad;++j)
	{
	  if(rad[i][j]>M2N.radius[i])break;
	  rr = rad[i][j]/rs;
	  rho = rhos/(rr*(1+rr)*(1+rr));
	  printf("PROFILE%d %e %e %e %e\n",i,rad[i][j],ngal[i][j],err[i][j],rho);
	  x1 += (ngal[i][j]-rho)*(ngal[i][j]-rho)/(err[i][j]*err[i][j]);
	}
      printf("CHIPROF%d %d %e %f %f %e %e\n",iter,i,x1, M2N.radius[i],rvir,
	     M2N.mass[i],M2N.model_mass[i]);

      chi2 += x1;
    }
  exit(0);
  ++iter;
  return chi2;
}

double meanN_givenM(double m)
{
  static int flag = 1, n=1000;
  static double *mass, *nbar, *yy;
  int i,j,k;
  double logm, dlogm, mlo, mhi, x, pn, ptot, mu, a, sig;

  if(flag)
    {
      flag = 0;
      mass = dvector(1,n);
      nbar = dvector(1,n);
      yy = dvector(1,n);

      mlo = 1.0E12;
      mhi = 1.0E16;
      dlogm = log(mhi/mlo)/(n-1);
      for(i=1;i<=n;++i)
	{
	  mass[i] = exp(dlogm*(i-1))*mlo;
	  logm = log(mass[i]);
	  
	  pn = ptot = 0;
	  sig = sig_rozo;
	  for(j=1;j<=220;++j)
	    {
	      mu = B_rozo - sig*sig/2 + alpha_rozo*log(i/40.0) + log(4e14*HUBBLE);      
	      x = exp(-(mu - logm)*(mu - logm)/(2*sig*sig))/(RT2PI*sig)/m;
	      ptot += x; 
	      pn += x*i;
	    }
	  nbar[i] = log(pn/ptot);
	  mass[i] = logm;
	}
      spline(mass,nbar,n,1.0E+30,1.0E+30,yy);
    }
  splint(mass,nbar,yy,n,log(m),&a);
  return exp(a);
}

/* testing to see if i can back out the P(N|M) from P(M|N)
 */
void test_pdf()
{
  int i,j,k;
  double m, sig, mu, x, logm, ptot=0, xx[200], Nbar, pn = 0;

  m = 2.0e13;
  logm = log(m);
  for(i=1;i<=100;++i)
    {
      sig = sig_rozo;
      mu = B_rozo - sig*sig/2 + alpha_rozo*log(i/40.0) + log(4e14*HUBBLE);      
      x = exp(-(mu - logm)*(mu - logm)/(2*sig*sig))/(RT2PI*sig)/m;
      ptot += x; 
      pn += x*i;
      xx[i] = x;
    }
  mu = exp(B_rozo + log(4e14*HUBBLE))/1.6;    
  Nbar = log(pow(m/mu,1/alpha_rozo)*40);
  Nbar = log(pn/ptot);
  fprintf(stderr,"%f %f\n",exp(Nbar),pn/ptot);
  sig = sig_rozo/alpha_rozo;
  mu = Nbar - sig*sig/2;
  for(i=1;i<=100;++i)
    {
      x = exp(-(mu - log(i))*(mu - log(i))/(2*sig*sig))/(RT2PI*sig)/i;
      printf("PDF %d %e %e\n",i,xx[i]/ptot,x);
    }
  exit(0);
}

// 41266 degrees in the sky
void input_maxbcg_counts()
{
  FILE *fp;
  char fname[100];
  int i,j,n,nf;

  sprintf(fname,"/home/tinker/SDSS/MAXBCG_ML/CATALOG/public_counts.dat");
  sprintf(fname,"/Users/tinker/cosmo/MN_RATIO/MAXBCG_ML/CATALOG/public_counts.dat");
  fp = openfile(fname);
  nf = filesize(fp);
  n = 220;
  M2N.counts_N200 = ivector(1,n);
  M2N.ndens_N200 = dvector(1,n);
  for(i=1;i<=n;++i)
    M2N.counts_N200[i] = 0;
  muh(0);
  for(i=1;i<=nf;++i)
    fscanf(fp,"%d %d",&j,&M2N.counts_N200[j]);
  fclose(fp);
  M2N.counts_nbins = n;
  fprintf(stderr,"Done reading [%d] lines from [%s]\n",nf,fname);
}

double comoving_volume(double zlo, double zhi)
{
  return 7500./41266.*pow(qromo(efunc1,zlo,zhi,midpnt)*c_on_H0,3.0)*4./3.*PI;
}
double efunc1(double z)
{
  return pow(OMEGA_M*(1+z)*(1+z)*(1+z) + (1-OMEGA_M),-0.5);
}
