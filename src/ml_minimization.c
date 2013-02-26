#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

#define MAGSOL 5.33

/* Data structures to hold multiple correlation functions
 * and their corresponding covariance matrices.
 */

struct WP_ARRAY {
  double *x;
  double *r;
  double *e;
  double *eigen;
  double **covar;

  double mag1;
  double mag2;
  double mean_mag;
  double ngal;
  double engal;
  double ngal_hod;
  int np;

  double M_min;
  double sigma_logM;
  double M1;
  double M_low;
  double M_cut;
  double mmin_hi;

} wp_array[15];

struct gao_parameters {
  double f0;
  double s0;
} GAO;

double **covar_array[15];

/* Local globals for qromo
 */
double *mass_g9,
  *lum_g9,
  *sig_g9,
  *y_lum_g9,
  *y_sig_g9,
  GLOBAL_SIGMA_LOGM = 0.15,
  curlum_g9;
int n_g9,
  COBE_NORM=1,
  VARY_TF=1,
  RANDOMIZE_NDENS=0,
  VARY_NGAL=0,
  VARY_GAMMA=0,
  VARY_CVIR=0,
  VARY_ALPHA=0,
  VARY_ALPHA1=0;


/* External function for density dependence.
 */
void dd_hod_functions(float *hmass, float *hdensity, int nhalo);


/* Local functions.
 */
void ml_input_data(int);
void switch_wp_data(int n);
double schechter(double mag, double mstar, double phi0,double alpha);
double func_integrate_schechter(double m);
double func_integrate_lum_schechter(double m);
double ml_initialize(int nwp, int nhod, double *a, double **cov1, double *avg1, double *start_dev);

double chi2_ml_function(int count);
void mass2light_functions(int n_wp, double *a);
double total_halo_light(double mass, int n_wp, double *a, double *ss);
double func_ml_lum(double m);
double func_ml_norm(double m);
double chi2_CNOC_data(int count, double *ml_cnoc, double *lum_cnoc);

double chi2_ml_wrapper(double *a);

void ml_function_test(int n_wp, double *a);
void ml_function_test_monte_carlo(int n_wp, double *a);
double poisson_prob(int n, double nave);
void dispersion_test(int n_wp, double *a, double mass);
double random_luminosity(double m1, double m2);
void calc_multiplicity_function(int n_wp, int n_start);
double chi2_multiplicity_function(int count);
void analyze_chain(void);

double chi2_number_density(int);

void output_wp_fits(int n_wp, double *a);

void read_chain_from_file(double **ax, int nx, int np, int *astart);

void PVD_one_halo(int ngal, float ga[1000][6], int nwp, double mass);

/* Local Functions for simulation population.
 */
void ml_mock_catalog(int n_wp, double *a);
double random_2df_magnitude_error(void);
void ml_mock_catalog_gao(int n_wp, double *a);
void ml_mock_gadget(int n_wp, double *a);
double N_cen_gao(double m, int ii, double *mmin);
double N_cen_hiden(double m);
void exponential_velocity(float v[]);
double non_gaussian_velocity(void);

double func_one_sat(double m);

/* External functions.
 */
double chi2_wp_wrapper(double *a);
void wp_input(void);
double func_BW2(double m); /* <--- calcs satellite galaxy density */

void sigma_logM_compare();

void ml_minimization()
{
  int i,j,k,n_wp=9,n_hod=0,n1,n2,n3;

  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*a1,*avg1,chi2,chi2prev,*start_dev,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,chi2array[15],
    xk,**chain,x,x2,x3,*eval_prev;
  int n,nrot,niter=0,count=0,imax_chain=30000,NSTEP = 50,NSTEP_MAX=10000,convergence=0,
    n_hod_params;
  long IDUM=-555;

  double chi2_ngal, chi2_ml, chi2_multi, chi2_wp, chi2_pk;

  int *pcheck,pcnt,ptot=20,nancheck,astart[10];

  FILE *fp;
  char fname[1000];
  int covar_outflag=1;

  fprintf(stderr,"\n\nCHI2 MINIMIZATION OF M/L DATA..........\n");
  fprintf(stderr,    "--------------------------------------------\n\n");

  /* If HOD.free[0]==-1, then randomize the number densities by 
   * 1-sigma away from the true value.
   */
  if(HOD.free[0]==-1)
    {
      HOD.free[0]=0;
      RANDOMIZE_NDENS=1;
    }

  if(wp.n_wp==4)NSTEP = 250;

  ml_input_data(wp.n_wp);

  Work.imodel=2;
  Work.chi2=1;
  HOD.pdfc = 9;

  /* In Zehavi et al, this is 40 Mpc/h,
   * but Norberg et al use 50.
   */
  wp.pi_max=50.0;

  srand48(32498793);

  astart[0] = 1;

  n_wp = wp.n_wp;

  /* If HOD.free[1]==2, then the galaxy densities
   * will vary within their errors, but we'll still
   * calculate M_min from the given HOD params and ngal.
   */
  if(HOD.free[1]==2)
    {
      VARY_NGAL=1;
      HOD.free[1]=0;
    }

  /* Number of magnitude bins:
   * default is 9, which means half-mag bins.
   * alternative is 4, which is full-mag bins.
   */

  if(wp.n_wp==9 && !HOD.free[1])
    {
      wp.ncf = 3;
      for(i=1;i<=6;++i)
	astart[i] = astart[i-1] + 2;
      for(i=7;i<=8;++i)
	astart[i] = astart[i-1] + 3;
    }
  if(wp.n_wp==9 && HOD.free[1])
    {
      wp.ncf = 4;
      for(i=1;i<=6;++i)
	astart[i] = astart[i-1] + 3;
      for(i=7;i<=8;++i)
	astart[i] = astart[i-1] + 4;
    }

  wp.ncf=0;
  for(i=1;i<=N_HOD_PARAMS;++i)
    if(HOD.free[i] && i!=2 && i!=6)wp.ncf++;
  printf("mcmc_min> %d free HOD params for each w_p.\n",wp.ncf);
  n_hod_params = wp.ncf;

  if(wp.n_wp==4)
    {
      for(i=1;i<wp.n_wp;++i)
	astart[i] = astart[i-1] + wp.ncf;
    }

  n = astart[wp.n_wp-1]+wp.ncf-1;
  printf("mcmc_min> %d  HOD  parameters\n",n);

  if(HOD.free[3]) {
    VARY_ALPHA=1;
    HOD.free[3] = 0;
    n++;
  }
  if(HOD.free[9]) {
    VARY_ALPHA1=1;
    HOD.free[9] = 0;
    n++;
  }
  if(HOD.free[6]) {
    VARY_CVIR=1;
    HOD.free[6] = 0;
    n++;
  }
  if(HOD.free[12]==2) {
    VARY_GAMMA=1;
    HOD.free[12]=0;
  }
  
  for(i=N_HOD_PARAMS+1;i<100;++i)
    if(HOD.free[i]) { n++; MCMC++; }

  printf("mcmc_min> %d  free parameters (MCMC=%d)\n",n,MCMC);

  pcheck=calloc(ptot,sizeof(int));

  a=dvector(1,n);
  start_dev=dvector(1,n);
  a1=dvector(1,n_hod);
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

  chi2prev=ml_initialize(wp.n_wp,n_hod,a,cov1,avg1,start_dev);

  if(IDUM_MCMC==-1)
    {
      OUTPUT=2;
      mass2light_functions(n_wp,a);
    }
  if(IDUM_MCMC==-2)
    {
      ml_mock_catalog(n_wp,a);
      exit(0);
    }
  if(IDUM_MCMC==-3)
    {
      calc_multiplicity_function(n_wp,2);
      calc_multiplicity_function(n_wp,5);
      exit(0);
    }
  if(IDUM_MCMC==-4)
    {
      output_wp_fits(n_wp,a);
      exit(0);
    }
  if(IDUM_MCMC==-5)
    ml_mock_catalog_gao(n_wp,a);
  if(IDUM_MCMC==-6)
    ml_function_test(n_wp,a);
  if(IDUM_MCMC==-7)
    ml_function_test_monte_carlo(n_wp,a);
  if(IDUM_MCMC==-8)
    analyze_chain();
  if(IDUM_MCMC==-10)
    ml_mock_gadget(n_wp,a);
  if(IDUM_MCMC==-11)
    sigma_logM_compare();

  niter++;
  for(i=1;i<=n;++i)
    aprev[i] = a[i];

  for(i=1;i<=n;++i)
    chain[niter][i]=a[i];

  IDUM=IDUM_MCMC;

  pcnt=0;
  pcheck[pcnt]=1;

  if(RESTART)
    {
      VARY_ALPHA=1;
      read_chain_from_file(chain,NSTEP_MAX,n,astart);
      NSTEP = NSTEP_MAX;
      niter = NSTEP_MAX;
    }


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
	  stepfac = 1;
	}
      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  i = N_HOD_PARAMS;
	  j = astart[wp.n_wp-1]+n_hod_params-1+VARY_ALPHA+VARY_CVIR+VARY_ALPHA1;
	  k = j;
	  if(HOD.free[++i])OMEGA_TEMP = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX = a[++j];

	  if(VARY_GAMMA)
	    GAMMA = gasdev(&IDUM)*0.02 + 0.15;
	  if(VARY_TF)
	    OMEGA_M = OMEGA_TEMP;
	}
      printf("COSMO %f %f %f\n",OMEGA_TEMP,SIGMA_8,SPECTRAL_INDX);

      chi2=0;

      /* Restrict CVIR_FAC (and maybe GAMMA)
       */
      i = astart[wp.n_wp-1]+n_hod_params+VARY_ALPHA;
      if(VARY_CVIR)
	if(a[i]<0.3 || a[i]>1.2)continue;
      if(SPECTRAL_INDX>1.5 || SPECTRAL_INDX<0.5)continue;
      /*
      if(HOD.free[10])
	if(GAMMA<0.11 || GAMMA>0.19)continue;
      */

      ++count;
      for(i=n_wp-1;i>=0;--i)
	{
	  n_hod = n_hod_params;
	  HOD.free[5] = 1;
	  if(i<6 && wp.n_wp==9)
	    { n_hod = n_hod_params-1; HOD.free[5] = 0; HOD.sigma_logM = GLOBAL_SIGMA_LOGM; }
	  for(k=0,j=astart[i];j<=astart[i]+n_hod-1;++j)
	    {
	      a1[++k] = a[j];
	    }
	  if(VARY_ALPHA)
	    HOD.alpha = a[astart[wp.n_wp-1]+n_hod_params];
	  if(VARY_CVIR)
	    CVIR_FAC = a[astart[wp.n_wp-1]+n_hod_params+1];
	  if(VARY_ALPHA1)
	    HOD.alpha1 = a[astart[wp.n_wp-1]+n_hod_params+2];
	  wp.ncf = n_hod;

 	  switch_wp_data(i);
	  chi2+=chi2array[i]=chi2_wp_wrapper(a1);

	  printf("MMIN %d %d %e\n",i,count,HOD.M_min);

	  wp_array[i].M_min = HOD.M_min;
	  wp_array[i].sigma_logM = HOD.sigma_logM;
	  wp_array[i].M1 = HOD.M1;
	  wp_array[i].M_low = HOD.M_low;
	  wp_array[i].M_cut = HOD.M_cut;
	  wp_array[i].ngal_hod = GALAXY_DENSITY;
	}

      chi2_wp = chi2;
      printf("CHI2_WP %d %e\n",count,chi2);

      if(HOD.free[1])
	chi2 += chi2_ngal = chi2_number_density(count);
      
      switch(WP_ONLY)
	{
	case 1: break;
	case 2: 
	  chi2 += chi2_ml = chi2_ml_function(count);
	  break;
	case 3:
	  chi2 += chi2_multi = chi2_multiplicity_function(count);
	  break;
	case 0:
	  chi2 += chi2_multi = chi2_multiplicity_function(count);
	  chi2 += chi2_ml = chi2_ml_function(count);
	  break;
	}	  
      
      if(!ThisTask)
	{
	  printf("TRY_ALL %d %e %e %f %f %f %f\n",
		 count,chi2,chi2prev,OMEGA_TEMP,SIGMA_8,GAMMA,SPECTRAL_INDX);
	  fflush(stdout);
	}
      
      if(COBE_NORM)
	chi2 += chi2_pk = cobe_prior(OMEGA_TEMP);

      pcheck[pcnt]=0;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	continue;
      pcheck[pcnt]=1;

      niter++;
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
	printf("ACC_ALL %d %d %e\n",niter,count,chi2);
	for(j=0;j<n_wp;++j)
	  {
	    printf("ACC%d %d %d ",j,niter,count);
	    n_hod = n_hod_params;
	    if(j<6 && wp.n_wp==9){ n_hod = n_hod_params-1; } 
	    for(k=0,i=astart[j];i<=astart[j]+n_hod-1;++i)
	      printf("%e ",a[i]);
	    for(i=astart[wp.n_wp-1]+n_hod_params;i<=n;++i)
	      printf("%e ",a[i]);
	    if(VARY_NGAL)
	      printf("%e ",wp_array[j].ngal_hod);
	    printf("%e %e\n",chi2array[j],chi2);fflush(stdout);
	  }
	printf("CHI2_ACC %d %d %e %e %e %e %e ",niter,count,chi2,chi2_wp,chi2_ml,chi2_multi,chi2_ngal);
	if(COBE_NORM)
	  printf("%e",chi2_pk);
	printf("\n");
      }
    }

  stepfac=1.0;
  pcnt=-1;

  while(niter<imax_chain)
    {
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  stepfac = stepfac*pow(0.9,5-j);
	  if(!ThisTask)printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      stepfac = 2.4/sqrt(n);

      if(niter>NSTEP_MAX)goto SKIP_MATRIX;

      printf("HERE %d %d\n",niter,NSTEP);

      for(j=1;j<=n;++j)
	{
	  avg1[j]=0;
	  for(k=1;k<=n;++k)
	    cov1[j][k]=0;
	}
      for(i=niter-NSTEP+1;i<=niter;++i)
	{
	  for(j=1;j<=n;++j)
	    {
	      avg1[j]+=chain[i][j];
	      for(k=1;k<=n;++k)
		cov1[j][k]+=chain[i][j]*chain[i][k];
	    }
	}


      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  tmp[i][j] = cov1[i][j]/NSTEP - avg1[i]*avg1[j]/NSTEP/NSTEP;

      jacobi(tmp,n,eval,evect,&nrot);
      gaussj(evect,n,tmp1,1);

      /* Check to make sure the eigenvalues are all there.
       */
      nancheck=0;
      for(i=1;i<=n;++i)
	if(isnan(eval[i])|| eval[i]<=0)nancheck=1;


    SKIP_MATRIX:

      if(niter==NSTEP_MAX && covar_outflag)
	{
	  sprintf(fname,"covar.%d",abs(IDUM_MCMC));
	  fp=fopen(fname,"w");
	  for(i=1;i<=n;++i)
	    {
	      fprintf(fp,"%e\n",eval[i]);
	      for(j=1;j<=n;++j)
		fprintf(fp,"%e\n",evect[i][j]);
	    }
	  fclose(fp);
	  covar_outflag = 0;
	}

      if(nancheck)
	{	  
	  printf("NANCHECK %d\n",count);
	  for(i=1;i<=n;++i)
	    a[i] = gasdev(&IDUM)*sqrt(tmp[i][i])*stepfac;
	}
      else
	{
	  for(i=1;i<=n;++i)
	    atemp[i] = gasdev(&IDUM)*sqrt(eval[i])*stepfac;
	  
	  for(i=1;i<=n;++i)
	    for(a[i]=0,j=1;j<=n;++j)
	      a[i] += atemp[j]*evect[j][i];
	}

      for(i=1;i<=n;++i)
	a[i] += aprev[i];

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  i=N_HOD_PARAMS;
	  j = astart[wp.n_wp-1]+n_hod_params-1+VARY_ALPHA+VARY_CVIR+VARY_ALPHA1;
	  if(HOD.free[++i])OMEGA_TEMP = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX    = a[++j];

	  if(VARY_TF)
	    OMEGA_M = OMEGA_TEMP;
	  /*
	  if(VARY_GAMMA)
	    GAMMA = gasdev(&IDUM)*0.02 + 0.15;
	  */
	}

      chi2=0;

      /* Restrict CVIR_FAC and GAMMA
       */
      i = astart[wp.n_wp-1]+n_hod_params+VARY_ALPHA;
      if(VARY_CVIR)
	if(a[i]<0.3 || a[i]>1.2)continue;
      if(SPECTRAL_INDX>1.5 || SPECTRAL_INDX<0.5)continue;
      /*
      if(HOD.free[10])
	if(GAMMA<0.11 || GAMMA>0.19)continue;
      */

      ++count;
      for(i=n_wp-1;i>=0;--i)
	{	  
	  n_hod = n_hod_params;
	  HOD.free[5] = 1;
	  if(i<6 && wp.n_wp==9)
	    { n_hod = n_hod_params-1; HOD.free[5] = 0; HOD.sigma_logM = GLOBAL_SIGMA_LOGM; }
	  for(k=0,j=astart[i];j<=astart[i]+n_hod-1;++j)
	    {
	      a1[++k] = a[j];
	    }
	  if(VARY_ALPHA)
	    HOD.alpha = a[astart[wp.n_wp-1]+n_hod_params];
	  if(VARY_CVIR)
	    CVIR_FAC = a[astart[wp.n_wp-1]+n_hod_params+1];
	  if(VARY_ALPHA1)
	    HOD.alpha1 = a[astart[wp.n_wp-1]+n_hod_params+2];
	  wp.ncf = n_hod;

	  switch_wp_data(i);
	  chi2+=chi2array[i]=chi2_wp_wrapper(a1);

	  wp_array[i].M_min = HOD.M_min;
	  wp_array[i].sigma_logM = HOD.sigma_logM;
	  wp_array[i].M1 = HOD.M1;
	  wp_array[i].M_low = HOD.M_low;
	  wp_array[i].M_cut = HOD.M_cut;
	  wp_array[i].ngal_hod = GALAXY_DENSITY;
	}
      chi2_wp = chi2;
      printf("CHI2_WP %d %e\n",count,chi2);

      if(HOD.free[1])
	chi2 += chi2_ngal = chi2_number_density(count);

      switch(WP_ONLY)
	{
	case 1: break;
	case 2: 
	  chi2 += chi2_ml = chi2_ml_function(count);
	  break;
	case 3:
	  chi2 += chi2_multi = chi2_multiplicity_function(count);
	  break;
	case 0:
	  chi2 += chi2_multi = chi2_multiplicity_function(count);
	  chi2 += chi2_ml = chi2_ml_function(count);
	}	  

      if(!ThisTask)
	{
	  printf("TRY_ALL %d %e %e %f %f %f %f\n",
		 count,chi2,chi2prev,OMEGA_TEMP,SIGMA_8,GAMMA,SPECTRAL_INDX);
	  fflush(stdout);
	}

      if(COBE_NORM)
	chi2 += chi2_pk = cobe_prior(OMEGA_TEMP);

      pcnt++;
      pcheck[pcnt]=1;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	{
	  /*
	  for(i=1;i<=n;++i)
	    a[i] = aprev[i];
	  chi2 = chi2prev;
	  */
	  pcheck[pcnt]=0;
	  continue;
	  
	}

      /*
      pcheck[pcnt]=0;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	continue;
      pcheck[pcnt]=1;

      for(i=1;i<=n;++i)
	printf("COV %d %e %e\n",niter,avg1[i]/NSTEP,sqrt(cov1[i][i]/NSTEP));
      */

      if(NSTEP<NSTEP_MAX)NSTEP++;
      niter++;

      /*
      if(niter%NSTEP_MAX==0 && !convergence)
	{
	  convergence = 1;
	  for(i=1;i<=n;++i)
	    if(fabs(eval[i]-eval_prev[i])/eval_prev[i]>0.01)convergence = 0;
	  if(convergence)
	    printf("CONVERGENCE! %d %d\n",niter,count);
	}
      */

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
	printf("ACC_ALL %d %d %e\n",niter,count,chi2);
	for(j=0;j<wp.n_wp;++j)
	  {
	    printf("ACC%d %d %d ",j,niter,count);
	    n_hod = n_hod_params;
	    if(j<6 && wp.n_wp==9){ n_hod = n_hod_params-1; } 
	    for(k=0,i=astart[j];i<=astart[j]+n_hod-1;++i)
	      printf("%e ",a[i]);
	    for(i=astart[wp.n_wp-1]+n_hod_params;i<=n;++i)
	      printf("%e ",a[i]);
	    if(VARY_NGAL)
	      printf("%e ",wp_array[j].ngal_hod);
	    printf("%e %e\n",chi2array[j],chi2);fflush(stdout);
	  }
	printf("CHI2_ACC %d %d %e %e %e %e %e ",niter,count,chi2,chi2_wp,chi2_ml,chi2_multi,chi2_ngal);
	if(COBE_NORM)
	  printf("%e",chi2_pk);
	printf("\n");

      }
    }


  exit(0);
}

/* This routine loops over the standard wp_input function that inputs wp data,
 * and inserts each set of data into the wp_array structure.
 * Currently I have hard-wired the magnitude limits and the filename.
 */
void ml_input_data(int n_wp)
{
  int i,j,k;
  FILE *fp;
  double mstar=-19.66,phi0=0.0161,alpha=-1.21,low_magnitude=-17.0,binsize=0.5;
  float csign;

  if(n_wp==4)
    {
      binsize = 1.0;
      low_magnitude = -18.0;
    }

  /*fp=openfile("galaxy_density_array.dat");*/

  for(i=0;i<n_wp;++i)
    {
      wp_array[i].mag1 = low_magnitude-binsize*i;
      wp_array[i].mag2 = low_magnitude-binsize*(i+1);
      if(i==n_wp-1)wp_array[i].mag2 -= 1;
      wp_array[i].mean_mag = low_magnitude-binsize*(i+.5); /* TEMP */

      if(n_wp==9)
	sprintf(wp.fname_wp,"/home/tinker/cosmo/2df/data/half_mag_samples/cov.10001.%.02f.%.02f.keML.com.xi2.2.log2",
	      wp_array[i].mag2,wp_array[i].mag1);
      if(n_wp==4)
	sprintf(wp.fname_wp,"/home/tinker/cosmo/2df/data/full_mag_samples/cov.10001.%.02f.%.02f.keML.com.xi2.2.log2",
	      wp_array[i].mag2,wp_array[i].mag1);
      if(i==n_wp-1)
	sprintf(wp.fname_wp,"/home/tinker/cosmo/2df/data/threshold_samples/cov.10001.-22.00.-21.00.keML.com.xi2.2.log2");

      sprintf(wp.fname_covar,"xi_covar.%d",i);
      wp_input();

      wp_array[i].x=dvector(1,wp.np);
      wp_array[i].r=dvector(1,wp.np);
      wp_array[i].e=dvector(1,wp.np);
      wp_array[i].eigen=dvector(1,wp.np);
      wp_array[i].covar=dmatrix(1,wp.np,1,wp.np);

      for(j=1;j<=wp.np;++j)
	{
	  wp_array[i].x[j]=wp.x[j];
	  wp_array[i].r[j]=wp.r[j];
	  wp_array[i].e[j]=wp.e[j];
	  wp_array[i].eigen[j]=wp.eigen[j];
	  for(k=1;k<=wp.np;++k)
	    wp_array[i].covar[j][k] = wp.covar[j][k];
	}
      wp_array[i].np=wp.np;

      free_dvector(wp.x,1,wp.np);
      free_dvector(wp.r,1,wp.np);
      free_dvector(wp.e,1,wp.np);

      /* This is TEMP as well.
       */
      wp_array[i].ngal = qromo(func_integrate_schechter,
			      wp_array[i].mag2,wp_array[i].mag1,midpnt);
      /* wp_array[i].ngal = qromo(func_integrate_schechter, 
	 -23,wp_array[i].mag1,midpnt);*/
      /* if(i==8)wp_array[i].ngal = 1.445e-4; */
      fprintf(stdout,"%.1f %.1f %e\n",wp_array[i].mag1,wp_array[i].mag2,wp_array[i].ngal);
      //fprintf(stdout,"%e\n",qromo(func_integrate_schechter,
      //wp_array[i].mag2,wp_array[i].mag1+0.1,midpnt));
      /*fscanf(fp,"%lf",&wp_array[i].ngal);*/


      if(n_wp==4) {	
	switch(i) {
	case 0: wp_array[i].ngal = 0.0112201; wp_array[i].engal = 0.063*0.0112201; break; 
	case 1: wp_array[i].ngal = 0.0069183; wp_array[i].engal = 0.036*0.0069183; break;
	case 2: wp_array[i].ngal = 0.0018197; wp_array[i].engal = 0.044*0.0018197; break;
	case 3: wp_array[i].ngal = 0.0001445; wp_array[i].engal = 0.062*0.0001445; 
	  wp_array[i].mag2 = -23; break;
	}
      } else {
	switch(i) {
	case 0: wp_array[i].engal = 0.134*wp_array[i].ngal; break;
	case 1: wp_array[i].engal = 0.095*wp_array[i].ngal ; break;
	case 2: wp_array[i].engal = 0.064*wp_array[i].ngal ; break;
	case 3: wp_array[i].engal = 0.050*wp_array[i].ngal ; break;
	case 4: wp_array[i].engal = 0.036*wp_array[i].ngal ; break;
	case 5: wp_array[i].engal = 0.020*wp_array[i].ngal ; break;
	case 6: wp_array[i].engal = 0.044*wp_array[i].ngal ; break;
	case 7: wp_array[i].engal = 0.039*wp_array[i].ngal ; break;
	case 8: wp_array[i].engal = 0.062*wp_array[i].ngal ; break;
	}
      if(RANDOMIZE_NDENS)
	{
	  csign = drand48()*2 - 1;
	  csign = csign/fabs(csign);
	  wp_array[i].ngal += csign*wp_array[i].engal;
	}


      }
    }
  wp.x=dvector(1,wp.np);
  wp.r=dvector(1,wp.np);
  wp.e=dvector(1,wp.np);

}

/* This is just the value of the Schechter Function fit to
 * 2F at magnitude m--> for getting total number density of galaxies.
 */
double func_integrate_schechter(double m)
{
  //  double mstar=-19.66,phi0=0.0161,alpha=-1.21; Norberg 2002
  double mstar=-19.57,phi0=0.0150,alpha=-1.18; // Cole et al 2005 (z=0.1)
  return(schechter(m,mstar,phi0,alpha));
}

/* This is the luminosity-weighted Schechter function value
 * at magnitdue m--> to get the total amount of light.
 */
double func_integrate_lum_schechter(double m)
{
  //  double mstar=-19.66,phi0=0.0161,alpha=-1.21,lum; Norberg 2002
  double mstar=-19.57,phi0=0.0150,alpha=-1.18,lum; // Cole et al 2005 (z=0.1)
  lum = pow(10.0,-0.4*(m-MAGSOL));
  return(schechter(m,mstar,phi0,alpha)*lum);
}

/* Function to return the value of the Schechter function
 * given the parameters.
 */
double schechter(double mag, double mstar, double phi0,double alpha)
{
  static int n = -1;
  static double *m, *f, *y, *e;
  double aa;
  int i;
  FILE *fp;

  /* This is the Schechter function fit
   */
  //  return(0.4*log(10.0)*phi0*pow(10.0,-0.4*(mag-mstar)*(alpha+1))*
  //	 exp(-pow(10.0,-0.4*(mag-mstar))));


  /* This is the method of getting the luminosity function from the
   * stepwise-maximum likelihood luminosity function (presented in Norber 2002
   * but updated in Shaun Cole's website (link on astro-ph/0111011).
   */
  if(n<0)
    {
      fp = openfile("/home/tinker/cosmo/2df/data/2dF_lumfunc_flat.3.dat");
      n = filesize(fp);
      m = dvector(1,n);
      f = dvector(1,n);
      e = dvector(1,n);
      y = dvector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%lf %lf %lf",&m[i],&f[i],&e[i]);
	  m[i] = mabs(m[i]);
	}
      fclose(fp);
      spline(m,f,n,1.0E+30,1.0E+30,y);
    }
  splint(m,f,y,n,mabs(mag),&aa);
  return(aa);

  /* This is the Schechter function fit in Norberg 2002
   */
  return(0.4*log(10.0)*phi0*pow(10.0,-0.4*(mag-mstar)*(alpha+1))*
	 exp(-pow(10.0,-0.4*(mag-mstar))));

}


/* Routine to switch the global wp data variables with
 * the current dataset in the chi^2 loop.
 */
void switch_wp_data(int n)
{
  static float GASDEV=0;
  static long IDUM4=-555;
  int i,j;

  i=n;
  for(j=1;j<=wp_array[i].np;++j)
    {
      wp.x[j]=wp_array[i].x[j];
      wp.r[j]=wp_array[i].r[j];
      wp.e[j]=wp_array[i].e[j];
    }
  wp.np=wp_array[i].np;
  GALAXY_DENSITY2 = GALAXY_DENSITY=wp_array[i].ngal;
  if(VARY_NGAL)
    {
      GASDEV = gasdev(&IDUM4);
      GALAXY_DENSITY2 = GALAXY_DENSITY = wp_array[i].ngal_hod = GASDEV*wp_array[i].engal + wp_array[i].ngal;
      printf("SWITCH DENSITY %d %e %e %e %f\n",i,GALAXY_DENSITY,wp_array[i].ngal,wp_array[i].engal,GASDEV);
    }

  HOD.i_wp = n;
}



double ml_initialize(int nwp, int nhod, double *a, double **cov1, double *avg1, 
		     double *start_dev)
{
  int i,j=0,k,astart[10],np,n_hod;
  double x1,x2,m1[50],alpha[50],mcut[50],sig[50],*a1,chi2,chi2array[15],mmin[50];
  long IDUM = -556;
  FILE *fp;
  char fname[100];


  astart[0] = 1;
  if(wp.n_wp==9)
    {
      for(i=1;i<=6;++i)
	astart[i] = astart[i-1] + wp.ncf-1;
      for(i=7;i<=8;++i)
	astart[i] = astart[i-1] + wp.ncf;
    }
  if(wp.n_wp==4)
    {
      for(i=1;i<=wp.n_wp-1;++i)
	astart[i] = astart[i-1] + wp.ncf;
    }

  a1=dvector(1,nhod);

  fp=openfile("initial_values.dat");

  if(!HOD.free[1])
    {
      for(i=1;i<=nwp;++i)
	fscanf(fp,"%lf %lf %lf %lf",&m1[i],&alpha[i],&mcut[i],&sig[i]);
    }
  else
    {
      for(i=1;i<=nwp;++i)
	fscanf(fp,"%lf %lf %lf %lf %lf",&mmin[i],&m1[i],&alpha[i],&mcut[i],&sig[i]);
    }    
  fclose(fp);


  j=0;
  for(k=1;k<=nwp;++k)
    {
      i=0;

      if(HOD.free[++i])
	{
	  /* If M_min is free, then for now just put a place-holder value in the array.
	   */
	  /* a[++j]=log10(m1[k]/20.0); */

	  /* Okay, now we have an actual value from the input file.
	   */
	  a[++j]=log10(mmin[k]);
	  start_dev[j]=0.0001;	  
	}

      if(HOD.free[++i])
	{
	  a[++j]=log10(m1[k]);
	  start_dev[j]=0.001;
	}

      /*
       * Reserve alpha to get the same for all wp
       */
      if(HOD.free[++i])
	{
	  /*
	  a[++j]=alpha[k];
	  start_dev[j]=0.01;
	  */
	}

      if(HOD.free[++i])
	{
	  a[++j]=log10(mcut[k]);
	  start_dev[j]=0.001;
	}
      
      /* Only leave sigma_logM free for i_wp>=6
       */
      if(nwp==9) {
	if(HOD.free[++i] && k>6)
	  {
	    a[++j]=log10(sig[k]);
	    start_dev[j]=0.01;
	  }
      }
      if(nwp==4) {
	if(HOD.free[++i])
	  {
	    a[++j]=log10(sig[k]);
	    start_dev[j]=0.01;
	  }
      }
    }

  if(VARY_ALPHA)
    {
      a[++j]=HOD.alpha;
      HOD.alpha = a[j] = alpha[1];
      start_dev[j]=0.01;
    }
      
  if(VARY_CVIR)
    {
      a[++j]=CVIR_FAC;
      start_dev[j]=0.01;
    }

  if(VARY_ALPHA1)
    {
      /* HOD.alpha1 = a[++j]=alpha[1]; */
      a[++j] = HOD.alpha1;
      start_dev[j]=0.01;
    }
      
  np = j;

  /* If using Powell's method or amoeba, then stop here
   */
  for(i=0;i<nwp;++i)
    {
      printf("IV %d ",i);
      nhod = 3;
      if(i<6 && nwp==9){ nhod = nhod-1; }
      for(j=astart[i];j<=astart[i]+nhod-1;++j)
	printf("%f ",a[j]);
      printf("%f ",HOD.alpha);
      if(VARY_ALPHA1)printf("%f ",HOD.alpha1);
      printf("\n");
      fflush(stdout);
    }

  if(!MCMC)return 0;

  if(MCMC>1)
    {
      i=N_HOD_PARAMS;
      j=np;
      if(HOD.free[++i]){ a[++j]=OMEGA_TEMP; start_dev[j]=0.005; np++; }
      if(HOD.free[++i]){ a[++j]=SIGMA_8; start_dev[j]=0.01; np++; }
      if(HOD.free[++i]){ a[++j]=GAMMA; start_dev[j]=0.01; np++; }
      if(HOD.free[++i]){ a[++j]=VBIAS; start_dev[j]=0.01; np++; }
      if(HOD.free[++i]){ a[++j]=VBIAS_C; start_dev[j]=0.01; np++; }
      if(HOD.free[++i]){ a[++j]=SPECTRAL_INDX; start_dev[j]=0.01; np++; }
    }

      
  for(i=1;i<=np;++i)
    {
      avg1[i]=a[i];
      printf("BEGIN %d %f\n",i,a[i]);
      for(j=1;j<=np;++j)
	cov1[i][j]=a[i]*a[j];
    }

  chi2=0;
  n_hod = wp.ncf;
  muh(n_hod);
  for(i=nwp-1;i>=0;--i)
    {
      /* HOD.alpha = -0.05*(8-i)+1.05; */
      
      nhod = n_hod;
      HOD.free[5] = 1;
      if(i<6 && nwp==9){ nhod = n_hod-1; HOD.free[5] = 0; HOD.sigma_logM = GLOBAL_SIGMA_LOGM; }
      wp.ncf=nhod;
      for(k=0,j=astart[i];j<=astart[i]+nhod-1;++j)
	{
	  a1[++k] = a[j];
	}
      if(VARY_ALPHA)
	HOD.alpha = a[astart[nwp-1]+n_hod];
      if(VARY_CVIR)
	CVIR_FAC = a[astart[nwp-1]+n_hod+VARY_ALPHA];
      printf("ALPHA = %f CVIR = %f\n",HOD.alpha,CVIR_FAC);
      
      switch_wp_data(i);

      /* We're having problems with M_min starts away from the correct values, so start things
       * off with M_min set to give the correct number densities.
       */
      /* This makes no difference one the chain begins. Get rid of it.
       *
       GALAXY_DENSITY = wp_array[i].ngal;
       HOD.M_min = 0;
       HOD.M1 = pow(10.0,a1[2]);
       HOD.M_cut = pow(10.0,a1[3]);
       if(i>=6 || nwp==4) HOD.sigma_logM = pow(10.0,a1[4]);
       set_HOD_params();
       a1[1] = log10(HOD.M_min);
      */

      printf("GALDEN %e\n",GALAXY_DENSITY);
      chi2+=chi2array[i]=chi2_wp_wrapper(a1);
      printf("CHI %d %e\n",i,chi2array[i]);

      wp_array[i].M_min = HOD.M_min;
      wp_array[i].sigma_logM = HOD.sigma_logM;
      wp_array[i].M1 = HOD.M1;
      wp_array[i].M_low = HOD.M_low;
      wp_array[i].M_cut = HOD.M_cut;
      wp_array[i].ngal_hod = GALAXY_DENSITY;

      printf("MMIN %d %d %e\n",i,0,HOD.M_min);

      x1 = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
      
      printf("FIRST %.1f %.1f %.2e %.2e %.2e %.2f %.2f\n",wp_array[i].mag1, chi2array[i], 
	     HOD.M_min, HOD.M1, HOD.M_cut, HOD.sigma_logM, x1);
    }
  printf("CHI2_WP 0 %f\n",chi2);

  /* If the number density is free, get the chi^2 from the observed number density
   */
  if(HOD.free[1])
    chi2+=chi2_number_density(0);

  switch(WP_ONLY)
    {
    case 1: break;
    case 2: 
      chi2+=chi2_ml_function(0);
      break;
    case 3:
      chi2+=chi2_multiplicity_function(0);
	  break;
    case 0:
      chi2+=chi2_multiplicity_function(0);
      chi2+=chi2_ml_function(0);
    }	  
  if(COBE_NORM)
    chi2 += cobe_prior(OMEGA_TEMP);

  printf("INITIAL CHI2 %e\n",chi2);
  return(chi2);
}

/* This function calculates the chi2 from the number density inferred from the HOD
 * parameters and the errors on the observed number densities.
 */
double chi2_number_density(int count)
{
  int i;
  double chi2 = 0,x;

  for(i=0;i<wp.n_wp;++i)
    {
      x = (wp_array[i].ngal - wp_array[i].ngal_hod)/wp_array[i].engal;
      chi2 += x*x;
      printf("NGAL%d %d %e %e %e\n",count,i,wp_array[i].ngal,wp_array[i].ngal_hod,x*x);
    }
  printf("CHI2_NGAL%d %e\n",count,chi2);
  return(chi2);
}

/* This function takes the 2PIGG multiplicity function and returns a chi^2 value
 * based on the multplicity function implies by the current point in the chain.
 */
double chi2_multiplicity_function(int count)
{
  int i,j,k,nmass=200,nmulti=125,n_start,n_wp=9;
  double dlogm,mass,mlo,mhi,dndM,nsat,ncen,chi2a=0,chi2b=0,x,nhi,nlo;
  static double *multi,*dy,*ngal;
  FILE *fp;
  char fname[100];
  
  static FILE *fp2,*fp5;
  static int flag=0,NP;
  static float *N_2PIGG,*M1_2PIGG,*E1_2PIGG,*M2_2PIGG,*E2_2PIGG,xf;

  if(!flag)
    {
      flag = 1;
      fp = openfile("/home/tinker/cosmo/2df/2PIGG_MULTI.dat");
      NP = filesize(fp);

      N_2PIGG = vector(1,NP);
      E1_2PIGG = vector(1,NP);
      M1_2PIGG = vector(1,NP);
      E2_2PIGG = vector(1,NP);
      M2_2PIGG = vector(1,NP);

      for(i=1;i<=NP;++i)
	{
	  fscanf(fp,"%f %f %f %f %f",&N_2PIGG[i],&M2_2PIGG[i],&E2_2PIGG[i],&M1_2PIGG[i],&E1_2PIGG[i]);
	  if(M1_2PIGG[i]>-99) {
	    M1_2PIGG[i] = pow(10.0,M1_2PIGG[i]);
	    E1_2PIGG[i] = log(10.0)*M1_2PIGG[i]*E1_2PIGG[i];
	  }
	  if(M2_2PIGG[i]>-99) {
	    M2_2PIGG[i] = pow(10.0,M2_2PIGG[i]);
	    E2_2PIGG[i] = log(10.0)*M2_2PIGG[i]*E2_2PIGG[i];
	  }
	}
      fclose(fp);

      multi = dvector(1,nmulti);
      dy = dvector(1,nmulti);
      ngal = dvector(1,nmulti);
      for(i=1;i<=nmulti;++i)
	ngal[i] = i;

      sprintf(fname,"%s.multi2",Task.root_filename);
      fp2 = fopen(fname,"w");
      sprintf(fname,"%s.multi5",Task.root_filename);
      fp5 = fopen(fname,"w");
    }

  /* First do the multiplicity function for galaxies M_b<-18
   */
  n_start = 2;

  for(i=1;i<=nmulti;++i)
    multi[i] = 0;


  mlo = log(wp_array[n_start].M_low);
  mhi = log(HOD.M_max);
  dlogm = (mhi - mlo)/(nmass-1);

  for(j=1;j<=nmass;++j)
    {
      mass = exp((j-1)*dlogm + mlo);
      dndM = dndM_interp(mass);
      nsat = ncen = 0;
      for(i=n_start;i<n_wp;++i)
	{
	  HOD.i_wp = i;
	  HOD.M1 = wp_array[i].M1;
	  HOD.M_min = wp_array[i].M_min;
	  HOD.M_low = wp_array[i].M_low;
	  HOD.M_cut = wp_array[i].M_cut;
	  HOD.sigma_logM = wp_array[i].sigma_logM;
	  nsat += N_sat(mass);
	  ncen += N_cen(mass);
	}
      for(k=1;k<=nmulti;++k)
	multi[k] += poisson_prob(k-1,nsat)*dndM*dlogm*mass*ncen;
    }

  fwrite(&count,sizeof(int),1,fp2);
  for(i=2;i<=NP;++i)
    {
      if(M1_2PIGG[i]<0)continue;
      nlo = pow(10.0,N_2PIGG[i]-0.14137/2.0);
      nhi = pow(10.0,N_2PIGG[i]+0.14137/2.0);

      x=0;
      for(j=(int)nlo-1;j<=(int)nhi+1;++j)
	if(j>0 && j>=nlo && j<nhi)
	  x+=multi[j];
      x/=0.14137;

      chi2a += (x-M1_2PIGG[i])*(x-M1_2PIGG[i])/E1_2PIGG[i]/E1_2PIGG[i];
      xf = x;

      fwrite(&N_2PIGG[i],sizeof(float),1,fp2);
      fwrite(&xf,sizeof(float),1,fp2);
      fwrite(&M1_2PIGG[i],sizeof(float),1,fp2);
      fwrite(&E1_2PIGG[i],sizeof(float),1,fp2);
      printf("MULTI2 %e %e\n",N_2PIGG[i],xf);
    }

  /* First do the multiplicity function for galaxies M_b<-19.5
   */
  n_start = 5;

  for(i=1;i<=nmulti;++i)
    multi[i] = 0;

  mlo = log(wp_array[n_start].M_low);
  mhi = log(HOD.M_max);
  dlogm = (mhi - mlo)/(nmass-1);

  for(j=1;j<=nmass;++j)
    {
      mass = exp((j-1)*dlogm + mlo);
      dndM = dndM_interp(mass);
      nsat = ncen = 0;
      for(i=n_start;i<n_wp;++i)
	{
	  HOD.i_wp = i;
	  HOD.M1 = wp_array[i].M1;
	  HOD.M_min = wp_array[i].M_min;
	  HOD.M_low = wp_array[i].M_low;
	  HOD.M_cut = wp_array[i].M_cut;
	  HOD.sigma_logM = wp_array[i].sigma_logM;
	  nsat += N_sat(mass);
	  ncen += N_cen(mass);
	}
      for(k=1;k<=nmulti;++k)
	multi[k] += poisson_prob(k-1,nsat)*dndM*dlogm*mass*ncen;
    }

  fwrite(&count,sizeof(int),1,fp5);
  for(i=2;i<=NP;++i)
    {
      if(M2_2PIGG[i]<0)continue;

      nlo = pow(10.0,N_2PIGG[i]-0.14137/2.0);
      nhi = pow(10.0,N_2PIGG[i]+0.14137/2.0);

      x=0;
      for(j=(int)nlo-1;j<=(int)nhi+1;++j)
	if(j>0 && j>=nlo && j<nhi)
	  x+=multi[j];
      x/=0.14137;

      chi2b += (x-M2_2PIGG[i])*(x-M2_2PIGG[i])/E2_2PIGG[i]/E2_2PIGG[i];

      xf = x;
      fwrite(&N_2PIGG[i],sizeof(float),1,fp5);
      fwrite(&xf,sizeof(float),1,fp5);
      fwrite(&M2_2PIGG[i],sizeof(float),1,fp5);
      fwrite(&E2_2PIGG[i],sizeof(float),1,fp5);
      printf("MULTI5 %e %e\n",N_2PIGG[i],xf);
    }
  fflush(fp2);
  fflush(fp5);

  /* Sum up the chi2 values and return.
   */
  printf("CHI2_MULTI %d %e %e %e\n",count,chi2a+chi2b,chi2a,chi2b);
  muh(0);
  return(chi2a+chi2b);
}

/* This function takes in the 2PIGG (corrected) data and returns a chi^2 value
 * of the M/L-L curve 
 */
double chi2_ml_function(int count)
{
  FILE *fp;
  int i,j,k,n=100,n_wp=9;
  double *lum,*sig,mlo,mhi,dlogm,*mass,*ml_lum,*y_lum, *y_sig, *a,chi2,ml,ml_cnoc,lum_cnoc;
  double t0,t1;
  char fname[100];
  float xf;

  static FILE *fp1;
  static int flag=0,N_2PIGG,N_LENSING;
  static float *ML_2PIGG, *L_2PIGG,*E1_2PIGG,*E2_2PIGG,*ML_LENSING,*L_LENSING,*E_LENSING;

  /* return(0); */

  if(!flag)
    {
      flag = 1;
      fp = openfile("/home/tinker/cosmo/2df/2PIGG_ML_linear.dat");
      N_2PIGG = filesize(fp);

      L_2PIGG = vector(1,N_2PIGG);
      E1_2PIGG = vector(1,N_2PIGG);
      E2_2PIGG = vector(1,N_2PIGG);
      ML_2PIGG = vector(1,N_2PIGG);

      for(i=1;i<=N_2PIGG;++i)
	fscanf(fp,"%f %f %f %f",&L_2PIGG[i],&ML_2PIGG[i],&E1_2PIGG[i],&E2_2PIGG[i]);
      fclose(fp);

      fp = openfile("/home/tinker/cosmo/2df/LENSING_ML.dat");
      N_LENSING = filesize(fp);

      L_LENSING = vector(1,N_LENSING);
      E_LENSING = vector(1,N_LENSING);
      ML_LENSING = vector(1,N_LENSING);

      for(i=1;i<=N_LENSING;++i)
	fscanf(fp,"%f %f %f",&L_LENSING[i],&ML_LENSING[i],&E_LENSING[i]);
      fclose(fp);

      sprintf(fname,"%s.ml",Task.root_filename);
      fp1 = fopen(fname,"w");
    }

  n_g9=n;

  mlo=log(1.0e11);
  mhi=log(5.0e15);
  dlogm=(mhi-mlo)/(n-1);

  mass_g9=dvector(1,n);
  lum_g9=dvector(1,n);
  sig_g9=dvector(1,n);
  ml_lum=dvector(1,n);

  y_lum_g9=dvector(1,n);
  y_sig_g9=dvector(1,n);

  /* Tabulate <M/L>_M and sig_L(M) curves
   */
  for(i=1;i<=n;++i)
    {
      curlum_g9=1;
      if(i==1)curlum_g9=-1;
      mass_g9[i] = exp((i-1)*dlogm+mlo);
      lum_g9[i] = total_halo_light(mass_g9[i],n_wp,a,&sig_g9[i]);
      if(mass_g9[i]<wp_array[0].M1/2)sig_g9[i]/=2;
    }

  spline(mass_g9,lum_g9,n,1.0E+30,1.0E+30,y_lum_g9);
  spline(mass_g9,sig_g9,n,1.0E+30,1.0E+30,y_sig_g9);

  /* Calculate the (M/L)_L curve at the 2PIGG locations.
   * Here is where we use OMEGA_TEMP: scale the M/L values by
   * the new value of OMEGA.
   */
  chi2=0;
  fwrite(&count,sizeof(int),1,fp1);
  xf = chi2;
  fwrite(&xf,sizeof(float),1,fp1);

  for(i=1;i<=N_2PIGG;++i)
    {
      curlum_g9 = L_2PIGG[i];
      ml = qromo(func_ml_lum,mlo,mhi,midpnt)/qromo(func_ml_norm,mlo,mhi,midpnt);
      ml = ml*OMEGA_TEMP/OMEGA_M;
      if(ml>ML_2PIGG[i])
	chi2 += (ml-ML_2PIGG[i])*(ml-ML_2PIGG[i])/E1_2PIGG[i]/E1_2PIGG[i];
      else
	chi2 += (ml-ML_2PIGG[i])*(ml-ML_2PIGG[i])/E2_2PIGG[i]/E2_2PIGG[i];
      xf = curlum_g9;
      fwrite(&xf,sizeof(float),1,fp1);
      xf = ML_2PIGG[i];
      fwrite(&xf,sizeof(float),1,fp1);
      xf = ml;
      fwrite(&xf,sizeof(float),1,fp1);
      printf("MLC %e %e\n",curlum_g9,ml);
    }
  for(i=1;i<=N_LENSING;++i)
    {
      curlum_g9 = L_LENSING[i];
      ml = qromo(func_ml_lum,mlo,mhi,midpnt)/qromo(func_ml_norm,mlo,mhi,midpnt);
      ml = ml*OMEGA_TEMP/OMEGA_M;
      chi2 += (ml-ML_LENSING[i])*(ml-ML_LENSING[i])/E_LENSING[i]/E_LENSING[i];
      xf = curlum_g9;
      fwrite(&xf,sizeof(float),1,fp1);
      xf = ML_LENSING[i];
      fwrite(&xf,sizeof(float),1,fp1);
      xf = ml;
      fwrite(&xf,sizeof(float),1,fp1);
    }

  chi2+=chi2_CNOC_data(count,&ml_cnoc,&lum_cnoc);
  xf = lum_cnoc;
  fwrite(&xf,sizeof(float),1,fp1);
  xf = 384;
  fwrite(&xf,sizeof(float),1,fp1);
  xf = ml_cnoc;
  fwrite(&xf,sizeof(float),1,fp1);

  fflush(fp1);
  fprintf(stdout,"CHI2_ML %d %e\n",count,chi2);

  return(chi2);
}

double chi2_CNOC_data(int count, double *ml_cnoc, double *lum_cnoc)
{
  static int flag=0,n;
  static float *mass,x1,x2,lum_avg=0;
  double chi2,dx1,*ax1,lum,ml_ave=0;
  FILE *fp;
  int i;
  static FILE *fp1;

  if(!flag)
    {
      flag = 1;
      fp = openfile("/home/tinker/cosmo/2df/CNOC_ML.dat");
      n = filesize(fp);
      mass = vector(1,n);

      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f %f",&mass[i],&x1,&x2);
	  lum_avg += mass[i]/x1;
	  mass[i] /= (OMEGA_TEMP/OMEGA_M);
	}
      lum_avg=lum_avg/n/1.07*0.887;
      fclose(fp);
      
    }

  for(i=1;i<=n;++i)
    {
      lum = total_halo_light(mass[i],9,ax1,&dx1);
      ml_ave += mass[i]/lum/(OMEGA_M/OMEGA_TEMP);
    }
  ml_ave/=n;
  chi2 = (ml_ave - 382)*(ml_ave - 382)/34./34.;

  *ml_cnoc = ml_ave;
  *lum_cnoc = lum_avg;

  return(chi2);
}


/* This function converts the (M/L)_M curve to the (M/L)_L curve.
 * It assumes that the dispersion about <M/L> at fixed M is a Gaussian
 * with variance as the sum of the individual variances for each mag bin.
 *
 * NB- Each sat function is Poisson, and each cen function is Bernoulli,
 * but with nine of each we're reaching the central limit theorem. Works great
 * for M>10^14 M_sol, but more of an approximation for 10^13 M_sol.
 */

void mass2light_functions(int n_wp, double *a)
{
  int i,j,k,n=100;
  double *lum,*sig,mlo,mhi,dlogm,*mass,*ml_lum,*y_lum, *y_sig;
  double t0,t1;

  n_g9=n;

  mlo=log(5.0e10);
  mhi=log(5.0e15);
  dlogm=(mhi-mlo)/(n-1);

  mass_g9=dvector(1,n);
  lum_g9=dvector(1,n);
  sig_g9=dvector(1,n);
  ml_lum=dvector(1,n);

  y_lum_g9=dvector(1,n);
  y_sig_g9=dvector(1,n);

  t0 = clock();

  /* Tabulate <M/L>_M and sig_L(M) curves
   */
  for(i=1;i<=n;++i)
    {
      mass_g9[i] = exp((i-1)*dlogm+mlo);
      lum_g9[i] = total_halo_light(mass_g9[i],n_wp,a,&sig_g9[i]);
      if(mass_g9[i]<wp_array[0].M1/2)sig_g9[i]/=2;
      /* printf("ML %e %e %e\n",log10(mass_g9[i]),log10(lum_g9[i]),mass_g9[i]/lum_g9[i]); */
    }

  /* exit(0); */

  fprintf(stderr,"TOTAL LUM DENSITY: %e\n",qromo(func_integrate_lum_schechter,-24.0,-13.0,midpnt));

  spline(mass_g9,lum_g9,n,1.0E+30,1.0E+30,y_lum_g9);
  spline(mass_g9,sig_g9,n,1.0E+30,1.0E+30,y_sig_g9);

  /* Calculate the (M/L)_L curve
   */
  for(i=1;i<=n;++i)
    {
      curlum_g9 = lum_g9[i];
      ml_lum[i] = qromo(func_ml_lum,mlo,mhi,midpnt)/qromo(func_ml_norm,mlo,mhi,midpnt);
      ml_lum[i] *= OMEGA_TEMP/OMEGA_M;
      printf("ML %e %e %e %e %e\n",log10(mass_g9[i]),log10(lum_g9[i]),
	     mass_g9[i]/lum_g9[i],ml_lum[i],log10(sig_g9[i]));
    }
  t1 = clock();
  fprintf(stderr,"%.2f\n", difftime(t1,t0)/CLOCKS_PER_SEC);
  exit(0);
}

/* Integrand of (M/L)_L function.
 */
double func_ml_lum(double m)
{
  double s,l,x;

  m=exp(m);
  
  splint(mass_g9,lum_g9,y_lum_g9,n_g9,m,&l);
  splint(mass_g9,sig_g9,y_sig_g9,n_g9,m,&s);

  x=(curlum_g9-l)/s;
  x=exp(-x*x/2)/(RT2PI*s)*m/l*dndM_interp(m);
  return(x*m);
}

/* Integrand of (M/L)_L normalization.
 */
double func_ml_norm(double m)
{
  double s,l,x;

  m=exp(m);  
  splint(mass_g9,lum_g9,y_lum_g9,n_g9,m,&l);
  splint(mass_g9,sig_g9,y_sig_g9,n_g9,m,&s);
  x=(curlum_g9-l)/s;
  x=exp(-x*x/2)/(RT2PI*s)*dndM_interp(m);
  return(x*m);
}

/* Sums up the total amount of light in each halo, also tabulates dispersion.
 */
double total_halo_light(double mass, int n_wp, double *a, double *ss)
{
  int i,j,k;
  double ncen,nsat,mean=0,err=0,norm,corr_mean,p=1,mean_cen=0;
  static double mmin[15];

  if(OUTPUT==2)
    printf("OCC %e ",log10(mass));

  mean=err=0;
  for(j=n_wp-1;j>=0;--j)
    {
      HOD.i_wp=j;
      HOD.M1 = wp_array[j].M1;
      HOD.M_cut = wp_array[j].M_cut;
      HOD.M_min = wp_array[j].M_min;
      HOD.sigma_logM = wp_array[j].sigma_logM;
      HOD.M_low = wp_array[j].M_low;

      ncen = N_cen(mass);
      nsat = N_sat(mass);

      /* Sum up probability of no satellite galaxies in this halos.
       */
      p*=exp(-nsat);
      mean_cen+=ncen*pow(10.0,-0.4*(wp_array[j].mean_mag-MAGSOL));

      /* Sum the squares of the errors
       */
      err+=(nsat+ncen*(1-ncen))*pow(10.0,-0.8*(wp_array[j].mean_mag-MAGSOL));
      mean+=(ncen+nsat)*pow(10.0,-0.4*(wp_array[j].mean_mag-MAGSOL));

      if(OUTPUT==2)
	printf("%e %e ",ncen,nsat);
    }
  norm = nsat/func_integrate_schechter(-17.25);

  if(OUTPUT==2){
    printf("\n");
    fflush(stdout);
  }
  *ss = sqrt(err);

  /* Use the value of N_sat at -17>M_b>-17.5 to normalize the schechter function
   * for extrapolating to lower luminosities.
   */
  norm = nsat/func_integrate_schechter(-17.25);
  corr_mean = mean + qromo(func_integrate_lum_schechter,-17.0,-13.0,midpnt)*norm;
  mean = corr_mean;
  return(mean);

  /* Now remove systems with no satellite galaxies from the means.
   */
  if(p==1)p=0;
  corr_mean = (mean - mean_cen*p)/(1-p);
  mean = corr_mean; 
  fflush(stdout);

  return(mean);
}

/********************************************************************************
 * Code for getting the central occupation funciton based on the magnitude
 * bin in question, and the next highest magnitude bin. (The lower limit of the
 * brighter mag bin will define the upper limit of the fainter bin.)
 *
 * The general form of each will be pdfc==2.
 *
 * The routine assumes that the values in the global HOD structure 
 * represent the HOD parameters to be used for the ii node in the array
 * (which is relevent since M_min must be set after the other params are chosen).
 *
 * If ii==8, then we are at the largest luminosity bin, which we'll treat
 * currently as a threshold sample.
 */
double N_cen_i(double m, int ii)
{
  int i;
  double n1,n2,n,logm,n2max=0;

  logm=log10(m);
  n1=0.5*(1+erf((logm - log10(HOD.M_min))/HOD.sigma_logM));

  /* NB NB NB NB
   * set the maximum N_cen for bJ<-21 galaxies to be less than unity.
   */
  if(ii==wp.n_wp-1)return(n1*0.5);
  if(ii==wp.n_wp-1)return(n1);

  for(n2=0,i=ii+1;i<wp.n_wp;++i)
    if(m>wp_array[i].M_low)
      {
	n2=0.5*(1+erf((logm - log10(wp_array[i].M_min))/wp_array[i].sigma_logM));
	if(i+1==wp.n_wp)n2*=0.5;
	if(n2>n2max)n2max=n2;
      }
  n = n1-n2max;
  /*
  if(n2>1)n2=1;
  if(n1+n2>1) n = 1-n2;
  else n = n1;
  if(ii==6 && (logm<13.01 && logm>12.99))
    {
      printf("NCEN %e %e %e %e %e %e\n",n1,n2,HOD.M_min,
	     HOD.sigma_logM,wp_array[ii+1].M_min,wp_array[ii+1].sigma_logM);
    }
  */
  if(n<0)return(0);
  return(n);
}

/*********************************************************************************/
/*********************************************************************************
 *
 * The following functions are to calculate the multiplicity function
 * for a given magnitude threshold.
 */
/*********************************************************************************/
/*********************************************************************************/

void calc_multiplicity_function(int n_wp, int n_start)
{
  int i,j,k,nmass=2000,nmulti=200;
  double dlogm,mass,*multi,mlo,mhi,dndM,nsat,*ncumu,ncen;

  multi = dvector(1,nmulti);
  ncumu = dvector(1,nmulti);
  for(i=1;i<=nmulti;++i)
    multi[i] = 0;

  mlo = log(wp_array[n_start].M_low);
  mhi = log(HOD.M_max);
  dlogm = (mhi - mlo)/(nmass-1);

  for(j=1;j<=nmass;++j)
    {
      mass = exp((j-1)*dlogm + mlo);
      dndM = dndM_interp(mass);
      nsat = ncen = 0;
      for(i=n_start;i<n_wp;++i)
	{
	  HOD.i_wp = i;
	  HOD.M1 = wp_array[i].M1;
	  HOD.M_min = wp_array[i].M_min;
	  HOD.M_low = wp_array[i].M_low;
	  HOD.M_cut = wp_array[i].M_cut;
	  HOD.sigma_logM = wp_array[i].sigma_logM;
	  nsat += N_sat(mass);
	  ncen += N_cen(mass);
	}
      for(k=1;k<=nmulti;++k)
	{
	  multi[k] += poisson_prob(k-1,nsat)*dndM*dlogm*mass*ncen;
	}
    }
  ncumu[nmulti]=multi[nmulti];
  for(i=nmulti-1;i>0;--i)
    ncumu[i] = ncumu[i+1]+multi[i];
  for(i=1;i<=nmulti;++i)
    printf("MULTI%d %d %e %e\n",n_start,i,multi[i],ncumu[i]);
  
}


/*********************************************************************************/
/*********************************************************************************
 *
 * The following functions are for using powell/amoeba to minimze the set of
 * correlation functions.
 *
 */
/*********************************************************************************/
/*********************************************************************************/

void ml_powell()
{
  int n,niter,i,j,n_wp=9,n_hod,n_mag,k,astart[10];
  double *a,**pp,**pp2,*yy,FTOL=1.0E-3,chi2min,s1,dlogm,m,d[100],*avg1,**cov1,*start_dev,*a1,at,at2;
  FILE *fp;
  char aa[1000];

  /* OUTPUT=1; */

  fprintf(stderr,"\n\nCHI2 MINIMIZATION OF BINNED W_P(R_P) DATA..........\n");
  fprintf(stderr,    "---------------------------------------------------\n\n");

  if(POWELL)
    FTOL=1.0E-3;
  else
    FTOL=1.0E-4;

  /* In Zehavi et al, this is 40 Mpc/h,
   * but Norberg et al use 50.
   */
  wp.pi_max=50.0;

  ml_input_data(n_wp);

  Work.imodel=2;
  Work.chi2=1;
  HOD.pdfc = 9;
  MCMC=0;

  srand48(32498793);

  /* Find the number of free parameters in the minimization
   * for the real-space correlation function.
   */
  n_hod=0;
  for(n=0,i=1;i<100;++i)
    {
      if(i<=N_HOD_PARAMS)
	n+=HOD.free[i]*n_wp;
      else
	n+=HOD.free[i];
      if(i<=N_HOD_PARAMS)
	n_hod+=HOD.free[i];

      if(OUTPUT)
	printf("mcmc_min> free[%i] = %d\n",i,HOD.free[i]);
    }
  n = 21;

  wp.ncf=n_hod;
  wp.ncf_tot=n;
  a=dvector(1,n);
  start_dev=dvector(1,n);

  ml_initialize(n_wp,n_hod,a,cov1,avg1,start_dev);

  astart[0] = 1;
  for(i=1;i<=6;++i)
    astart[i] = astart[i-1] + 2;
  for(i=7;i<=8;++i)
    astart[i] = astart[i-1] + 3;

  for(i=1;i<=n;++i) 
    a[i]=pow(10.0,a[i]);

  for(i=1;i<=n;++i)
    printf("a[%d] = %e\n",i,a[i]);

  n_mag=n_wp-1;


  n_mag = 8;

  /* Loop through all magnitudes.
   */
 MAGNITUDE_LOOP:

  n = n_hod = 3;
  HOD.free[5] = 1;
  /* if(n_mag<6){ n = n_hod = 2; HOD.free[5] = 0; HOD.sigma_logM = GLOBAL_SIGMA_LOGM;} */

  /*
  HOD.free[3]=1;
  n = n_hod = 4;
  HOD.free[5] = 1;
  if(n_mag<6){ n = n_hod = 3; HOD.free[5] = 0; HOD.sigma_logM = GLOBAL_SIGMA_LOGM;}
  HOD.alpha = -0.05*(8-n_mag)+1.05;
  if(n_mag==4)POWELL=0;

  */



  n = n_hod;
  wp.ncf = n;
  a1=dvector(1,n);
  if(POWELL)
    pp=dmatrix(1,n,1,n);
  else if(n_mag==n_wp-1)
    pp=dmatrix(1,n+1,1,n);
  yy=dvector(1,n+1);

  for(k=0,j=astart[n_mag];j<=astart[n_mag]+n_hod-1;++j)
    a1[++k] = a[j];

  if(n_mag<6)a1[3] = GLOBAL_SIGMA_LOGM;

  HOD.i_wp = n_mag;
  /*
  at = a1[2];
  a1[2] = 0.9;
  for(i=3;i<n;++i)
    {
      at2 = a1[i];
      a1[i] = at;
      at = at2;
    }
  a1[n] = at;
  */

  switch_wp_data(n_mag);

  /* Make the starting stepsize 10% of the initial values.
   */
  for(i=1;i<=n;++i)
    d[i]=a1[i]*0.5;

  if(POWELL)
    {
      for(i=1;i<=n;++i)
	{
	  for(j=1;j<=n;++j)
	    {
	      pp[i][j]=0;
	      if(i==j)pp[i][j]+=d[j];
	    }
	}
    }
  else
    {
      for(j=1;j<=n;++j)
	pp[1][j]=a1[j];
      yy[1]=chi2_wp(a1);
 
      for(i=1;i<=n;++i)
	{
	  a1[i]+=d[i];
	  if(i>1)a1[i-1]-=d[i-1];
	  yy[i+1]=chi2_wp(a1);	  
	  for(j=1;j<=n;++j)
	    pp[i+1][j]=a1[j];
	}
      a1[wp.ncf_tot]-=d[wp.ncf_tot];
    }

  if(POWELL) 
    {
      if(OUTPUT)printf("wp_min> starting powell.\n");
      powell(a1,pp,n,FTOL,&niter,&chi2min,chi2_wp);
      chi2min = chi2_wp(a1);
    }
  else
    {
      if(OUTPUT)printf("wp_min> starting amoeba.\n");
      amoeba(pp,yy,n,FTOL,chi2_wp,&niter);
      for(i=1;i<=n;++i)a1[i]=pp[1][i];
      chi2min = chi2_wp(a1);
    }	
  printf("POWELL%d %e %e ",n_mag,chi2min,HOD.M_min);
  for(i=1;i<=n;++i)
    printf("%e ",a1[i]);
  printf("%e %f\n",GALAXY_BIAS,HOD.alpha);
  fflush(stdout);

  wp_array[n_mag].M_min = HOD.M_min;
  wp_array[n_mag].M_low = HOD.M_low;
  wp_array[n_mag].M_cut = HOD.M_cut;
  wp_array[n_mag].M1 = HOD.M1;
  wp_array[n_mag].sigma_logM = HOD.sigma_logM;
  n_mag--;

  if(n_mag<0)
    exit(0);

  free_dvector(a1,1,n);
  free_dvector(yy,1,n+1);
  if(POWELL)
    free_dmatrix(pp,1,n,1,n);
  else
    /* free_dmatrix(pp,1,n+1,1,n); */
  goto MAGNITUDE_LOOP;

}


/* This routine loops through all 9 magnitude bins
 * to get the total chi^2 for a set of parameters.
 */
double chi2_ml_wrapper(double *a)
{
  static double *a1;
  static int flag=1,niter=0;
  double chi2,chi2i,t0,t1;
  int i,j,nwp,nhod;
  
  t0=clock();

  nwp=9;
  nhod=wp.ncf;
  if(flag)
    {
      flag=0;
      a1=dvector(1,nhod);
    }

  chi2=0;
  ++niter;

  for(i=nwp-1;i>=0;--i)
    {
      for(j=1;j<=nhod;++j)
	a1[j] = a[i+(j-1)*nwp+1];

      switch_wp_data(i);
      chi2+=chi2i=chi2_wp(a1);

      printf("TRY%d %d ",i,niter);
      for(j=1;j<=nhod;++j)
	printf("%.4e ",a1[j]);
      printf("%e\n",chi2i);
      fflush(stdout);

      wp_array[i].M_min = HOD.M_min;
      wp_array[i].sigma_logM = HOD.sigma_logM;
    }
  t1=clock();
  printf("ITERALL %d %e %.2f\n",niter,chi2,difftime(t1,t0)/CLOCKS_PER_SEC);
  return(chi2);
}

/*********************************************************************************
 *
 * This is for Monte Carlo testing of some of the routines.
 *
 */
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/


void ml_function_test(int n_wp, double *a)
{
  FILE *fp;
  int i,j,n,imass;
  double mass;
  char aa[1000];

  /* fp=openfile("/gal/sr2/tinker/voids/LANL/halo.LANL.400"); */
  fp=openfile("/gal/sr2/tinker/voids/LANL/halo.400_z.5.dat");

  while(!feof(fp))
    {
      fscanf(fp,"%d %d",&i,&imass);
      fgets(aa,1000,fp);
      mass=imass*RHO_CRIT*OMEGA_M*pow(0.3125,3.0);
      dispersion_test(n_wp,a,mass);
      if(feof(fp))break;
    }
  exit(0);
}

/* does the same as the above function, but instead of using a mass
 * function for a simulation, it uses a mass funcation created by
 * monte carlo sampling of the analytic function.
 */
void ml_function_test_monte_carlo(int n_wp, double *a)
{
  FILE *fp;
  int i,j,n,imass;
  double mass,pmax,p;
  char aa[1000];
  long IDUMx=-442123;

  n=1E6;

  i=0;
  pmax = dndM_interp(1.0E11);
  while(i<n)
    {
      mass = (1.0e16)*ran2(&IDUMx) + 1.0e11;
      p = dndM_interp(mass)/pmax;
      if(ran2(&IDUMx)>p)continue;
      i++;
      dispersion_test(n_wp,a,mass);
    }

}

void ml_mock_catalog(int n_wp, double *a)
{
  FILE *fp,*fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9];
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[9][1000],halocnt[1000],imag;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float galarr[1000][6];
  int ngal,nsati[9],ALL_FILES=0;

  fp=openfile("/home/tinker/LANL/halo.LANL.400");
  //  fp=openfile("/home/tinker/LANL/halo.400_z.5.dat");
  for(i=j_start;i<n_wp;++i)
    {
      sprintf(aa,"%s.mock.%d",Task.root_filename,i);      
      fpa[i] = fopen(aa,"w");
      if(!ALL_FILES)continue;
      sprintf(aa,"%s.mag.%d",Task.root_filename,i);
      fpb[i] = fopen(aa,"w");
      sprintf(aa,"%s.cen.%d",Task.root_filename,i);
      fpc[i] = fopen(aa,"w");
      sprintf(aa,"%s.sat.%d",Task.root_filename,i);
      fps[i] = fopen(aa,"w");
    }


  for(i=0;i<1000;++i)
    halocnt[i]=0;
  for(j=0;j<n_wp;++j)
    for(i=0;i<1000;++i)
      galcnt[j][i]=0;

  j=j_start;
  HOD.i_wp=j;
  HOD.M_min=wp_array[j].M_min;
  HOD.M1 = wp_array[j].M1;
  HOD.M_cut = wp_array[j].M_cut;
  HOD.sigma_logM = wp_array[j].sigma_logM;
  mlo = set_low_mass();
  printf("MLO %e %e %f\n",mlo,HOD.M_min,HOD.sigma_logM);
  fflush(stdout);

  while(!feof(fp))
    {

      fscanf(fp,"%d %d %e %e %e %e %e %e %e",
	     &i,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
      mass=imass*RHO_CRIT*OMEGA_M*pow(0.3125,3.0);
      if(mass<mlo)continue;

      i1 = (int)(log10(mass)/0.1);

      for(j=j_start;j<n_wp;++j)
	{
	  /* HOD.alpha = -0.05*(8-j)+1.05; */

	  ngal=0;

	  HOD.i_wp=j;
	  HOD.M1 = wp_array[j].M1;
	  HOD.M_cut = wp_array[j].M_cut;
	  HOD.sigma_logM = wp_array[j].sigma_logM;
	  HOD.M_min=wp_array[j].M_min;
	  /* HOD.alpha = 1.0; */
	  
	  nsat = N_sat(mass);
	  if(j>j_start)
	    nc[j] = nc[j-1]+N_cen(mass);
	  else
	    nc[j] = N_cen(mass);

	  n1 = poisson_deviate(nsat);

	  nsati[j] = n1;

	  for(i=1;i<=n1;++i)
	    {
	      /* random_2df_magnitude_error(); */
	      mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
	      if(ALL_FILES)fprintf(fpb[j],"%f\n",mag);
	      r = NFW_position(mass,xg);
	      NFW_velocity(mass,vg,mag);
	      /* exponential_velocity(vgf); */
	      for(k=0;k<3;++k)
		{
		  xg[k]+=xh[k];
		  if(xg[k]<0)xg[k]+=BOX_SIZE;
		  if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
		  vg[k]+=vh[k];
		  /* vg[k] = vgf[k]; */
		  /*
		  vg[k] = gasdev(&IDUM3)*500;
		  vg[k] = non_gaussian_velocity();
		  */
		}	
	      /*
	      if(j<=3 && drand48()<0.2)
		fprintf(fpa[j],"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]);
	      if(j>3)
	      */
		fprintf(fpa[j],"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]);
	      if(ALL_FILES)
		fprintf(fps[j],"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]);
	      imag = (-mag-17)*2;
	      if(imag>=n_wp)imag = n_wp;

	      /* Bin up the galaxies by halo mass to check the HOD
	       */
	      galcnt[imag][i1]++;

	      if(ngal<1000) {
		galarr[ngal][0] = xg[0];
		galarr[ngal][1] = xg[1];
		galarr[ngal][2] = xg[2];
		galarr[ngal][3] = vg[0];
		galarr[ngal][4] = vg[1];
		galarr[ngal][5] = vg[2];
		ngal++;
	      }
	    }
	  PVD_one_halo(ngal,galarr,j,mass);
	}
      ncen=drand48();
      for(j=j_start;j<n_wp;++j)
	if(ncen<nc[j])break;
      if(j<n_wp) {
	/* Giving the central galaxies velocities
	 */
	/*
	NFW_velocity(mass,vg,mag);
	for(k=0;k<3;++k)
	  vh[k]+=vg[k];
	*/
	galarr[ngal][0] = xg[0];
	galarr[ngal][1] = xg[1];
	galarr[ngal][2] = xg[2];
	galarr[ngal][3] = vg[0];
	galarr[ngal][4] = vg[1];
	galarr[ngal][5] = vg[2];
	ngal++;

	/*
	exponential_velocity(vgf);
	*/
	/*
	for(k=0;k<3;++k)
	  vh[k] = non_gaussian_velocity();
	*/
	/*
	  vh[k] = gasdev(&IDUM3)*500;
	*/
	/*
	if(j<=3 && drand48()<0.2)
	  fprintf(fpa[j],"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
	if(j>3)
	*/
	  fprintf(fpa[j],"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
	  
	if(ALL_FILES)
	  fprintf(fpc[j],"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
	mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
	if(ALL_FILES)
	  fprintf(fpb[j],"%f\n",mag);
	imag = (-mag-17)*2;
	if(imag>=n_wp)imag = n_wp;
	galcnt[imag][i1]++;
      }
      halocnt[i1]++;	  

      if(feof(fp))break;
    }

  PVD_one_halo(-1,galarr,0,0);

  /* output the binned HOD
   */
  for(j=0;j<n_wp;++j)
    {
      sprintf(aa,"binned_HOD_err.%d",j);
      fp2=fopen(aa,"w");
      for(i=0;i<1000;++i)
	if(galcnt[j][i]>0)
	  fprintf(fp2,"%d %f %f %d %d\n",
		  i,(i+0.5)*0.1,(float)galcnt[j][i]/halocnt[i],galcnt[j][i],halocnt[i]);
      fclose(fp2);
    }
  
  exit(0);
}

void ml_mock_gadget(int n_wp, double *a)
{
  FILE *fp,*fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9];
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[9][1000],halocnt[1000],imag,igad;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float galarr[1000][6];
  int ngal,nsati[9],ALL_FILES=0;

  for(igad = 1; igad <=5 ; ++igad) {

    sprintf(aa,"/home/tinker/cosmo/gadget_runs/gamma.12/run%d/halo.02",igad);
  fp=openfile(aa);
  for(i=j_start;i<n_wp;++i)
    {
      sprintf(aa,"%s.mock%d.%d",Task.root_filename,igad,i);      
      fpa[i] = fopen(aa,"w");
      if(!ALL_FILES)continue;
      sprintf(aa,"%s.mag.%d",Task.root_filename,i);
      fpb[i] = fopen(aa,"w");
      sprintf(aa,"%s.cen.%d",Task.root_filename,i);
      fpc[i] = fopen(aa,"w");
      sprintf(aa,"%s.sat.%d",Task.root_filename,i);
      fps[i] = fopen(aa,"w");
    }


  for(i=0;i<1000;++i)
    halocnt[i]=0;
  for(j=0;j<n_wp;++j)
    for(i=0;i<1000;++i)
      galcnt[j][i]=0;

  j=j_start;
  HOD.i_wp=j;
  HOD.M_min=wp_array[j].M_min;
  HOD.M1 = wp_array[j].M1;
  HOD.M_cut = wp_array[j].M_cut;
  HOD.sigma_logM = wp_array[j].sigma_logM;
  mlo = set_low_mass();
  printf("MLO %e %e %f\n",mlo,HOD.M_min,HOD.sigma_logM);
  fflush(stdout);

  while(!feof(fp))
    {

      fscanf(fp,"%d %d %e %e %e %e %e %e %e",
	     &i,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
      mass=imass*RHO_CRIT*OMEGA_M*pow(0.702,3.0);
      if(mass<mlo)continue;

      i1 = (int)(log10(mass)/0.1);

      for(j=j_start;j<n_wp;++j)
	{
	  /* HOD.alpha = -0.05*(8-j)+1.05; */

	  ngal=0;

	  HOD.i_wp=j;
	  HOD.M1 = wp_array[j].M1;
	  HOD.M_cut = wp_array[j].M_cut;
	  HOD.sigma_logM = wp_array[j].sigma_logM;
	  HOD.M_min=wp_array[j].M_min;
	  /* HOD.alpha = 1.0; */
	  
	  nsat = N_sat(mass);
	  if(j>j_start)
	    nc[j] = nc[j-1]+N_cen(mass);
	  else
	    nc[j] = N_cen(mass);

	  n1 = poisson_deviate(nsat);

	  nsati[j] = n1;

	  for(i=1;i<=n1;++i)
	    {
	      /* random_2df_magnitude_error(); */
	      mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
	      if(ALL_FILES)fprintf(fpb[j],"%f\n",mag);
	      r = NFW_position(mass,xg);
	      NFW_velocity(mass,vg,mag);
	      /* exponential_velocity(vgf); */
	      for(k=0;k<3;++k)
		{
		  xg[k]+=xh[k];
		  if(xg[k]<0)xg[k]+=BOX_SIZE;
		  if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
		  vg[k]+=vh[k];
		  /* vg[k] = vgf[k]; */
		  /*
		  vg[k] = gasdev(&IDUM3)*500;
		  vg[k] = non_gaussian_velocity();
		  */
		}	
	      /*
	      if(j<=3 && drand48()<0.2)
		fprintf(fpa[j],"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]);
	      if(j>3)
	      */
		fprintf(fpa[j],"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]);
	      if(ALL_FILES)
		fprintf(fps[j],"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]);
	      imag = (-mag-17)*2;
	      if(imag>=n_wp)imag = n_wp;

	      /* Bin up the galaxies by halo mass to check the HOD
	       */
	      galcnt[imag][i1]++;

	      if(ngal<1000) {
		galarr[ngal][0] = xg[0];
		galarr[ngal][1] = xg[1];
		galarr[ngal][2] = xg[2];
		galarr[ngal][3] = vg[0];
		galarr[ngal][4] = vg[1];
		galarr[ngal][5] = vg[2];
		ngal++;
	      }
	    }
	  PVD_one_halo(ngal,galarr,j,mass);
	}
      ncen=drand48();
      for(j=j_start;j<n_wp;++j)
	if(ncen<nc[j])break;
      if(j<n_wp) {
	/* Giving the central galaxies velocities
	 */
	/*
	NFW_velocity(mass,vg,mag);
	for(k=0;k<3;++k)
	  vh[k]+=vg[k];
	*/
	galarr[ngal][0] = xg[0];
	galarr[ngal][1] = xg[1];
	galarr[ngal][2] = xg[2];
	galarr[ngal][3] = vg[0];
	galarr[ngal][4] = vg[1];
	galarr[ngal][5] = vg[2];
	ngal++;

	/*
	exponential_velocity(vgf);
	*/
	/*
	for(k=0;k<3;++k)
	  vh[k] = non_gaussian_velocity();
	*/
	/*
	  vh[k] = gasdev(&IDUM3)*500;
	*/
	/*
	if(j<=3 && drand48()<0.2)
	  fprintf(fpa[j],"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
	if(j>3)
	*/
	  fprintf(fpa[j],"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
	  
	if(ALL_FILES)
	  fprintf(fpc[j],"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
	mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
	if(ALL_FILES)
	  fprintf(fpb[j],"%f\n",mag);
	imag = (-mag-17)*2;
	if(imag>=n_wp)imag = n_wp;
	galcnt[imag][i1]++;
      }
      halocnt[i1]++;	  

      if(feof(fp))break;
    }
  for(i=0;i<wp.n_wp;++i)
    fclose(fpa[i]);
  }

  PVD_one_halo(-1,galarr,0,0);

  /* output the binned HOD
   */
  for(j=0;j<n_wp;++j)
    {
      sprintf(aa,"binned_HOD_err.%d",j);
      fp2=fopen(aa,"w");
      for(i=0;i<1000;++i)
	if(galcnt[j][i]>0)
	  fprintf(fp2,"%d %f %f %d %d\n",
		  i,(i+0.5)*0.1,(float)galcnt[j][i]/halocnt[i],galcnt[j][i],halocnt[i]);
      fclose(fp2);
    }
  
  exit(0);
}
void ml_mock_catalog_internal(int n_wp, double *a)
{
  FILE *fp,*fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9];
  int i,j,k,n,ii,imass,n1,j_start=0,i1,galcnt[9][1000],halocnt[1000],imag;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2;
  char aa[1000];
  float x1,xh[3],vh[3];
  long IDUM3 = -445;

  float galarr[1000][6];
  int ngal;
  static int flag = 1, nhalo;
  static float **xh_arr,**vh_arr,*mass_arr;

  if(flag)
    {
      fp=openfile("/home/tinker/LANL/halo.LANL.400");
      n=filesize(fp);
      xh_arr=matrix(1,nhalo,0,2);
      vh_arr=matrix(1,nhalo,0,2);
      mass_arr=vector(1,nhalo);

      for(nhalo=1;nhalo<=n;++nhalo)
	{
	  fscanf(fp,"%d %d %e %e %e %e %e %e %e",
		 &i,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
	  mass=imass*RHO_CRIT*OMEGA_M*pow(0.3125,3.0);
	  
	  mass_arr[nhalo]=mass;
	  for(i=0;i<3;++i)
	    {
	      xh_arr[nhalo][i]=xh[i];
	      vh_arr[nhalo][i]=vh[i];
	    }
	  if(feof(fp))break;
	}
      fclose(fp);
      flag=0;
      nhalo = n;
    }

  for(i=0;i<1000;++i)
    halocnt[i]=0;
  for(j=0;j<n_wp;++j)
    for(i=0;i<1000;++i)
      galcnt[j][i]=0;

  j=j_start;
  HOD.i_wp=j;
  HOD.M_min=wp_array[j].M_min;
  HOD.M1 = wp_array[j].M1;
  HOD.M_cut = wp_array[j].M_cut;
  HOD.sigma_logM = wp_array[j].sigma_logM;
  mlo = set_low_mass();
  printf("MLO %e %e %f\n",mlo,HOD.M_min,HOD.sigma_logM);
  fflush(stdout);


  for(ii=1;ii<=nhalo;++ii)
    {
      mass = mass_arr[nhalo]=mass;
      for(i=0;i<3;++i)
	{
	  xh[i] = xh_arr[nhalo][i];
	  vh[i] = vh_arr[nhalo][i];
	}

      if(mass<mlo)continue;

      i1 = (int)(log10(mass)/0.1);

      for(j=j_start;j<n_wp;++j)
	{
	  ngal=0;

	  HOD.i_wp=j;
	  HOD.M1 = wp_array[j].M1;
	  HOD.M_cut = wp_array[j].M_cut;
	  HOD.sigma_logM = wp_array[j].sigma_logM;
	  HOD.M_min=wp_array[j].M_min;
	  HOD.alpha = 1.0;
	  
	  nsat = N_sat(mass);
	  if(j>j_start)
	    nc[j] = nc[j-1]+N_cen(mass);
	  else
	    nc[j] = N_cen(mass);

	  n1 = poisson_deviate(nsat);

	  for(i=1;i<=n1;++i)
	    {
	      /* random_2df_magnitude_error(); */
	      mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
	      NFW_position(mass,xg);
	      NFW_velocity(mass,vg,mag);
	      for(k=0;k<3;++k)
		{
		  xg[k]+=xh[k];
		  if(xg[k]<0)xg[k]+=BOX_SIZE;
		  if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
		  vg[k]+=vh[k];
		}
	      imag = (-mag-17)*2;
	      if(imag>=n_wp)imag = n_wp;

	      /* Bin up the galaxies by halo mass to check the HOD
	       */
	      galcnt[imag][i1]++;

	      if(ngal<1000) {
		galarr[ngal][0] = xg[0];
		galarr[ngal][1] = xg[1];
		galarr[ngal][2] = xg[2];
		galarr[ngal][3] = vg[0];
		galarr[ngal][4] = vg[1];
		galarr[ngal][5] = vg[2];
		ngal++;
	      }
	    }
	  PVD_one_halo(ngal,galarr,j,mass);
	}
      ncen=drand48();
      for(j=j_start;j<n_wp;++j)
	if(ncen<nc[j])break;
      if(j<n_wp) {
	/* Giving the central galaxies velocities
	 */
	/*
	NFW_velocity(mass,vg,mag);
	for(k=0;k<3;++k)
	  vh[k]+=vg[k];
	*/
	galarr[ngal][0] = xg[0];
	galarr[ngal][1] = xg[1];
	galarr[ngal][2] = xg[2];
	galarr[ngal][3] = vg[0];
	galarr[ngal][4] = vg[1];
	galarr[ngal][5] = vg[2];
	ngal++;

	imag = (-mag-17)*2;
	if(imag>=n_wp)imag = n_wp;
	galcnt[imag][i1]++;
      }
      halocnt[i1]++;	  

      if(feof(fp))break;
    }

  PVD_one_halo(-1,galarr,0,0);
}

double random_2df_magnitude_error()
{
  double p,m;
  p=0;
  while(3*drand48()>p)
    {
      m = drand48()*2-1;
      p = 1/RT2PI/.14*exp(-m*m/2/.14/.14)*0.7 + 
	1/RT2PI/.235*exp(-pow(log(1+m),2.0)/2/.235/.235)/(1+m)*0.3;
      /*
      printf("BUH %f %f %f %f\n",m,p,
	     1/RT2PI/.14*exp(-m*m/2/.14/14),1/RT2PI/.235*exp(-pow(log(1+m),2.0)/2/.235/.235)/(1+m)*0.3 );
      */
    }
  return(m);
}


double non_gaussian_velocity()
{
  double s = 500,v,p,w,h4 = -0.1,p0;
  static long IDUM6 = -4523;

  p0 = (1+h4*3/sqrt(24.0));
  p = -1;
  while(p<drand48())
    {
      v = (1 - 2*drand48())*s*10;
      w = v/s;
      p = exp(-w*w/2)*(1+h4/sqrt(24.0)*(4*w*w*w*w - 12*w*w + 3))/p0;
    }
  return(v);
}

void exponential_velocity(float v[])
{
  int i;
  double s = 500,p;
  long IDUM = -423342;

  for(i=0;i<3;++i)
    {
      p = -1;
      while(p<drand48())
	{
	  v[i] = (1 - 2*drand48())*s*20;
	  if(v<0) p = exp(-fabs(v[i]/s)*ROOT2);
	  else p = exp(-fabs(v[i]/(s*0.9)*ROOT2));
	}
    }

  for(i=0;i<3;++i)
    v[i] = gasdev(&IDUM)*500.0;
  return;

}


void dispersion_test(int n_wp, double *a, double mass)
{
  int i,j,k,n=1,n1,n1t;
  double ml_ratio,x1,nsat,ncen,err=0,mean=0,nc[10],mag;

  for(i=1;i<=n;++i)
    {
      n1t =mean=err=x1=0;
      for(j=0;j<n_wp;++j)
	{
	  HOD.i_wp=j;
 	  HOD.M_min=wp_array[j].M_min;
 	  HOD.M_low=wp_array[j].M_low;
	  HOD.M1 = wp_array[j].M1;
	  HOD.M_cut = wp_array[j].M_cut;
	  HOD.sigma_logM = wp_array[j].sigma_logM;
	  	  
	  nsat = N_sat(mass);
	  if(j>0)
	    {
	      if(mass>HOD.M_low)
		nc[j] = nc[j-1]+N_cen(mass);
	      else
		nc[j] = nc[j-1]+0;
	    }
	  else
	    {
	      if(mass>HOD.M_low)
		nc[j] = N_cen(mass);
	      else
		nc[j] = 0;
	    }

	  n1t += n1 = poisson_deviate(nsat);
	  mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
	  x1 += (n1)*pow(10.0,-0.4*(mag-MAGSOL));

	}
      ncen=drand48();
      for(j=0;j<n_wp;++j)
	nc[j]/=nc[8];
      for(j=0;j<n_wp;++j)
	if(ncen<nc[j])break;
      if(j>=n_wp)j=n_wp-1;
      mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
      x1+=pow(10.0,-0.4*(mag-MAGSOL));
      if(x1>0) {
	HOD.i_wp = j;
	printf("LUM %f %f %f %d %f %e %d %e\n",
	       log10(x1),log10(mass),mass/x1,n1t,ncen,nc[j],j,N_cen(mass));
	fflush(stdout);
      }
    }

}

double random_luminosity(double m1, double m2)
{
  double p1,p2,pn,m;

  pn = func_integrate_schechter(m1);
  for(;;)
    {
      m = drand48()*(m2-m1)+m1;
      p1 = func_integrate_schechter(m)/pn;
      p2 = drand48();
      if(p2<p1)break;
    }
  return(m);
}


/*********************************************************************************/
/*********************************************************************************/
/*
 * For outputting the fits and the data all in one file.
 *
 */
/*********************************************************************************/
/*********************************************************************************/

void output_wp_fits(int n_wp, double *a)
{
  int i,j,k,nr=50;
  double *x[10],rlo=0.1,rhi=50,dlogr,r,x1,*a1;
  FILE *fp;
  char fname[100];

  a1 = dvector(1,4);

  nr = 50;

  dlogr = (log(rhi)-log(rlo))/(nr-1);

  for(i=0;i<n_wp;++i)
    x[i] = dvector(1,nr);

  for(i=n_wp-1;i>=0;--i)
    {
      muh(i);

      HOD.i_wp=i;
      HOD.M_low = wp_array[i].M_low;
      HOD.M_min = wp_array[i].M_min;
      a1[1] = HOD.M1 = wp_array[i].M1;
      a1[3] = HOD.M_cut = wp_array[i].M_cut;
      a1[4] = HOD.sigma_logM = wp_array[i].sigma_logM;
      GALAXY_DENSITY=wp_array[i].ngal;
      a1[2] = HOD.alpha;

      printf("BOO %d %e %e %e %e %e %e %e\n",
	     i,HOD.M_low,HOD.M_min,HOD.M1,HOD.alpha,HOD.M_cut,HOD.sigma_logM,GALAXY_DENSITY);

      //      chi2_wp(a1);

      BETA = pow(OMEGA_M,0.6)/qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)*
	GALAXY_DENSITY;
      
      RESET_KAISER++;
      RESET_FLAG_1H = RESET_FLAG_2H = 1;
 
      for(j=1;j<=nr;++j)
	{
	  r = exp(dlogr*(j-1))*rlo;
	  x[i][j] = projected_xi(r);
	  continue;

	  rlo = pow(10.0,-0.8 + (j-1)/5.);
	  rhi = pow(10.0,-0.6 + (j-1)/5.);

	  printf("%d %f %f %f %f\n",j,-0.8 + (j-1)/5.,-0.6 + (j-1)/5.,rlo,rhi);
	  x[i][j] = 0;
	  x1 = 0;
	  for(k=1;k<=100;++k)
	    {
	      r = (rhi-rlo)/100*(k-0.5) + rlo;
	      x[i][j] += projected_xi(r)*2*PI*r*(rhi-rlo)/100.0;
	      x1 += 2*PI*r*(rhi-rlo)/100.0;
	    }
	  x[i][j]/=x1;
	  printf("%f %f %f %f\n",wp_array[0].r[j],projected_xi(wp_array[0].r[j]),
		 x[i][j],projected_xi(wp_array[0].r[j])/x[i][j]);

	}
    }

  sprintf(fname,"%s.wp_fits",Task.root_filename);
  fp=fopen(fname,"w");

  for(j=1;j<=nr;++j)
    {
      if(nr!=12)
	fprintf(fp,"%e ",exp(dlogr*(j-1))*rlo);
      else
	fprintf(fp,"%e ",wp_array[0].r[j]);
      for(i=0;i<n_wp;++i)
	fprintf(fp,"%e ",x[i][j]);
      fprintf(fp,"\n");
    }
  fclose(fp);

  sprintf(fname,"%s.wp_data",Task.root_filename);
  fp=fopen(fname,"w");

  for(j=1;j<=wp_array[0].np;++j)
    {
      fprintf(fp,"%e ",wp_array[0].r[j]);
      for(i=0;i<n_wp;++i)
	fprintf(fp,"%e %e ",wp_array[i].x[j],wp_array[i].e[j]);
      fprintf(fp,"\n");
    }
  fclose(fp);

}

/***********************************************************************************
 *
 * This is to make mocks based on the Gao effect. This requires that the
 * local halo densities are read in from a file, and that M_min is
 * recalculated such that the number density is help constant.
 *
 ***********************************************************************************/

void ml_mock_catalog_gao(int n_wp, double *a)
{
  FILE *fp,*fpa[9],*fp2,*fpb[9];
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[9][1000],halocnt[1000],imag,nh,density_flag=1,ih;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,
    mmin_hi[10],mmin_lo[10],mmin_fac[10],*mmin_tmp;
  float *hmass,*hden;
  char aa[1000];
  float x1,xh[3],vh[3];
  long IDUM3 = -445;

  srand48(4323432);

  GAO_EFFECT = 1;
  DENSITY_THRESHOLD = 2.0;
  GAO.f0 = 2.0;
  GAO.s0 = 2.0;

  fprintf(stderr,"BEGINNING ML_MOCK GAO...\n");

  fp=openfile("/gal/sr2/tinker/voids/LANL/halo.LANL.400");
  fp2=openfile("/gal/sr2/tinker/voids/LANL/den.r5.LANL.400");
  for(i=j_start;i<n_wp;++i)
    {
      sprintf(aa,"%s.mock.%d",Task.root_filename,i);
      fpa[i] = fopen(aa,"w");
      sprintf(aa,"%s.mag.%d",Task.root_filename,i);
      fpb[i] = fopen(aa,"w");
    }

  j=2;
  HOD.i_wp=j;
  HOD.M_min=wp_array[j].M_min;
  HOD.M1 = pow(10.0,a[j+0*n_wp+1]);
  HOD.M_cut = pow(10.0,a[j+1*n_wp+1]);
  HOD.sigma_logM = pow(10.0,a[j+2*n_wp+1]);
  mlo = set_low_mass();
  printf("NCEN %e\n",N_avg(5.0e11));

  /* Read in the masses and the densities.
   */
  nh = filesize(fp);
  hmass = vector(1,nh);
  hden = vector(1,nh);
  for(i=1;i<=nh;++i)
    {
      fscanf(fp,"%d %d %e %e %e %e %e %e %e",
	     &i,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
      fscanf(fp2,"%e",&hden[i]);
      mass=imass*RHO_CRIT*OMEGA_M*pow(0.3125,3.0);
      hmass[i] = mass;
    }
  rewind(fp);
  fclose(fp2);

  for(i=0;i<1000;++i)
    halocnt[i]=0;
  for(j=0;j<n_wp;++j)
    for(i=0;i<1000;++i)
      galcnt[j][i]=0;


  j=j_start;
  HOD.i_wp=j;
  HOD.M_min=wp_array[j].M_min;
  HOD.M1 = pow(10.0,a[j+0*n_wp+1]);
  HOD.M_cut = pow(10.0,a[j+1*n_wp+1]);
  HOD.sigma_logM = pow(10.0,a[j+2*n_wp+1]);
  mlo = set_low_mass();
  printf("MLO %e %e %f\n",mlo,HOD.M_min,HOD.sigma_logM);
  fflush(stdout);

  for(j=n_wp-1;j>=j_start;--j)
    {
      HOD.i_wp=j;
      HOD.M_min=wp_array[j].M_min;
      HOD.M1 = pow(10.0,a[j+0*n_wp+1]);
      HOD.M_cut = pow(10.0,a[j+1*n_wp+1]);
      HOD.sigma_logM = pow(10.0,a[j+2*n_wp+1]);

      GALAXY_DENSITY = wp_array[j].ngal;
      dd_hod_functions(hmass,hden,nh);
      mmin_lo[j] = HOD.M_min_loden;
      mmin_hi[j] = HOD.M_min_loden*(1+GAO.f0*exp(-j*j/2.0/GAO.s0/GAO.s0));
      wp_array[j].M_min = HOD.M_min_loden;
      wp_array[j].M_low = HOD.M_low;

      fprintf(stdout,"MMIN %d %e %e %e\n",j,mmin_lo[j],mmin_hi[j],HOD.M_low);
    }
  fflush(stdout);

  for(ih=1;ih<=nh;++ih)
    {
      fscanf(fp,"%d %d %e %e %e %e %e %e %e",
	     &i,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
      mass=imass*RHO_CRIT*OMEGA_M*pow(0.3125,3.0);
      if(mass<mlo)continue;

      i1 = (int)(log10(mass)/0.1);
      /*
      mass = 5e11;
      hden[1] = 10;
      hden[2] = 0;
      */
      for(j=j_start;j<n_wp;++j)
	{
	  HOD.i_wp=j;
	  HOD.M_min=wp_array[j].M_min;
	  HOD.M_low=wp_array[j].M_low;
	  HOD.M1 = pow(10.0,a[j+0*n_wp+1]);
	  HOD.M_cut = pow(10.0,a[j+1*n_wp+1]);
	  HOD.sigma_logM = pow(10.0,a[j+2*n_wp+1]);
	
	  if(hden[ih]>DENSITY_THRESHOLD)
	    {
	      if(j>j_start)
		nc[j] = nc[j-1]+N_cen_hiden(mass);
	      else
		nc[j] = N_cen_hiden(mass);
	      /* printf("NCEN %d %d %e %e\n",i,ih,N_cen_hiden(mass),nc[j]);*/
	    }
	  else
	    {
	      if(j>j_start)
		nc[j] = nc[j-1]+N_cen(mass);
	      else
		nc[j] = N_cen(mass);
	      /* printf("NCEN %d %d %e %e\n",i,ih,N_cen_hiden(mass),nc[j]);*/
	    }

	  nsat = N_sat(mass);
	  n1 = poisson_deviate(nsat);

	  for(i=1;i<=n1;++i)
	    {
	      NFW_position(mass,xg);
	      NFW_velocity(mass,vg,mag);
	      for(k=0;k<3;++k)
		{
		  xg[k]+=xh[k];
		  if(xg[k]<0)xg[k]+=BOX_SIZE;
		  if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
		  vg[k]+=vh[k];
		}
	      fprintf(fpa[j],"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]); 
	      mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
	      fprintf(fpb[j],"%f\n",mag);
	      imag = (-mag-17)*2;
	      if(imag>=n_wp)imag = n_wp;

	      /* Bin up the galaxies by halo mass to check the HOD
	       */
	      galcnt[imag][i1]++;
	    }
	}
      density_flag=0;

      ncen=drand48();
      for(j=j_start;j<n_wp;++j)
	if(ncen<nc[j])break;
      if(j<n_wp) {
	/*
	if(j==4) {
	  for(k=j_start;k<n_wp;++k)
	    printf("BOO %d %d %e %e\n",i,k,nc[k],ncen);
	}
	*/
	fprintf(fpa[j],"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
	mag = random_luminosity(wp_array[j].mag1,wp_array[j].mag2);
	fprintf(fpb[j],"%f\n",mag);
	imag = (-mag-17)*2;
	if(imag>=n_wp)imag = n_wp;
	galcnt[imag][i1]++;
      }
      halocnt[i1]++;	  

      if(feof(fp))break;
    }

  /* output the binned HOD
   */
  for(j=0;j<n_wp;++j)
    {
      sprintf(aa,"binned_HOD.%d",j);
      fp2=fopen(aa,"w");
      for(i=0;i<1000;++i)
	if(galcnt[j][i]>0)
	  fprintf(fp2,"%d %f %f %d %d\n",
		  i,(i+0.5)*0.1,(float)galcnt[j][i]/halocnt[i],galcnt[j][i],halocnt[i]);
      fclose(fp2);
    }

  exit(0);
}


double N_cen_gao(double m, int ii, double *mmin)
{
  int i;
  double n1,n2,n,logm,n2max=0;

  logm=log10(m);
  n1=0.5*(1+erf((logm - log10(HOD.M_min))/HOD.sigma_logM));
  if(ii==8)return(n1);

  for(n2=0,i=ii+1;i<9;++i)
    if(m>wp_array[i].M_low)
      {
	n2=0.5*(1+erf((logm - log10(mmin[i]))/wp_array[i].sigma_logM));
	if(n2>n2max)n2max=n2;
      }
  n = n1-n2max;
  if(n<0)return(0);
  return(n);
}

double N_cen_hiden(double m)
{
  int i,ii;
  double n1,n2,n,logm,n2max=0,f;

  if(m<HOD.M_low)return(0);
  ii = HOD.i_wp;
  f=1 + GAO.f0*exp(-ii*ii/2.0/GAO.s0/GAO.s0);

  logm=log10(m);
  n1=0.5*(1+erf((logm - log10(HOD.M_min*f))/HOD.sigma_logM));
  if(ii==8)return(n1);

  for(n2=0,i=ii+1;i<9;++i)
    if(m>wp_array[i].M_low)
      {
	f=1 + GAO.f0*exp(-i*i/2.0/GAO.s0/GAO.s0);
	n2=0.5*(1+erf((logm - log10(wp_array[i].M_min*f))/wp_array[i].sigma_logM));
	if(n2>n2max)n2max=n2;
      }
  n = n1-n2max;
  if(n<0)return(0);
  return(n);
}

/*************************************************************************************88
 * This is to read in a chain and do some other calculations with it.
 */

void analyze_chain()
{
  int i,j,k,n,n_wp=9,i1,i2;
  FILE *fp,*fpa[9];
  char aa[1000];
  float x1,x2,x3,x4,x5,x6,x7;
  double xd1,xd2;

  /*
  for(i=0;i<n_wp;++i)
    {
      sprintf(aa,"grep ACC%d out.560 > acc.%d",i,i);
      fprintf(stderr,"%s\n",aa);
      system(aa);
    }
  */
  for(i=0;i<n_wp;++i)
    {
      sprintf(aa,"acc.%d",i);
      fpa[i] = openfile(aa);
    }      
  n = filesize(fpa[0]);

  for(i=1;i<=n;++i)
    {
      for(j=n_wp-1;j>=0;--j)
	{
	  if(j==n_wp-1)RESET_COSMOLOGY++;

	  if(j>=6)
	    {
	      fscanf(fpa[j],"%4s %d %d %e %e %e %e %e %e %e",aa,&i1,&i2,&x1,&x2,&x3,&x4,&x5,&x6,&x7);
	    }
	  else
	    {
	      fscanf(fpa[j],"%4s %d %d %e %e %e %e %e %e",aa,&i1,&i2,&x1,&x2,&x4,&x5,&x6,&x7);
	      x3 = log10(GLOBAL_SIGMA_LOGM);
	    }
	  printf("%d %f %f %f %f %f %f %f\n",i,x1,x2,x3,x4,x5,x6,x7);
	  
	  fgets(aa,1000,fpa[j]);
	  HOD.M1 = pow(10.0,x1);
	  HOD.M_cut = pow(10.0,x2);
	  HOD.sigma_logM = pow(10.0,x3);
	  HOD.alpha = x4;
	  CVIR_FAC = x5;
	  OMEGA_TEMP = x6;
	  SIGMA_8 = x7;

	  set_HOD_params();

	  wp_array[j].M_min = HOD.M_min;
	  wp_array[j].sigma_logM = HOD.sigma_logM;
	  wp_array[j].M1 = HOD.M1;
	  wp_array[j].M_low = HOD.M_low;
	  wp_array[j].M_cut = HOD.M_cut;
	}
      printf("SIGMA %f %f %f\n",SIGMA_8,OMEGA_TEMP,OMEGA_M);
      fflush(stdout);

      chi2_multiplicity_function(i);
      printf("CHI2_CNOC %e\n",chi2_CNOC_data(i,&xd1,&xd2));
    }

  exit(0);
 }

/*****************************************************************************
 * THis file is to reread a chain from a file so that you can start using
 * the covariance matrix from that chain.
 */

void read_chain_from_file(double **ax, int nx, int np, int *astart)
{
  int n,i,j,k,i1,i2,ii,n_hod;
  FILE *fp,*fp2;
  char aa[1000];
  float **a,alpha,ngal[9][1000],mmin[9][1000],mone[9][1000],
    mcut[9][1000],m1[9][1000],bias[9],ebias[9],bb;
  double *a1,x1,x2,x3;

  fprintf(stderr,"TOP\n");

  for(i=0;i<9;++i)
    bias[i]=ebias[i]=0;

  a1 = dvector(1,4);

  n=0;
  /* fp=openfile(RESTART_FILE); */
  fp=openfile("trunc.555");

  while(!feof(fp))
    {
      fscanf(fp,"%7s",aa);
      if(strncmp(aa,"ACC_ALL",7)!=0)
	{
	  fgets(aa,1000,fp);
	  continue;
	}
      n++;
      fgets(aa,1000,fp);
      if(feof(fp))break;
    }
  rewind(fp);

  a=matrix(1,n,1,np);

  muh(n);

  i=0;
  while(!feof(fp))
    {
      fscanf(fp,"%7s",aa);
      if(strncmp(aa,"ACC_ALL",7)!=0)
	{
	  fgets(aa,1000,fp);
	  continue;
	}

      fgets(aa,1000,fp);
      i++;
      muh(i);
      k=0;
      for(j=1;j<=6;++j)
	{
	  fscanf(fp,"%s %d %d %e %e",aa,&i1,&i2,&a[i][++k],&a[i][++k]);
	  fgets(aa,1000,fp);
	}
      for(j=1;j<=2;++j)
	{
	  fscanf(fp,"%s %d %d %e %e %e",aa,&i1,&i2,&a[i][++k],&a[i][++k],&a[i][++k]);
	  fgets(aa,1000,fp);
	}
      fscanf(fp,"%s %d %d %e %e %e",aa,&i1,&i2,&a[i][++k],&a[i][++k],&a[i][++k]);
      fscanf(fp,"%e",&alpha);

      /* &a[i][++k],&a[i][++k],&a[i][++k],&a[i][++k]);*/
      fgets(aa,1000,fp);
      /*
      sprintf(aa,"initial_values/iv.%d",i);
      fprintf(stderr,"Opening [%s]\n",aa);
      fp2=fopen(aa,"w");
      for(k=0,j=1;j<=6;++j)
	fprintf(fp2,"%e %e %e %e\n",pow(10.0,a[i][++k]),alpha,pow(10.0,a[i][++k]),GLOBAL_SIGMA_LOGM);
      for(j=1;j<=3;++j)
	fprintf(fp2,"%e %e %e %e\n",pow(10.0,a[i][++k]),alpha,pow(10.0,a[i][++k]),
		pow(10.0,a[i][++k]));
      fclose(fp2);
      */
      for(ii=8;ii>=0;--ii)
	{
	  HOD.M1 = pow(10.0,a[i][astart[ii]]);
	  HOD.M_cut = pow(10.0,a[i][astart[ii]+1]);
	  if(ii>6)HOD.sigma_logM = pow(10.0,a[i][astart[ii]+2]);
	  else HOD.sigma_logM = 0.15;
	  HOD.alpha = alpha;
	  HOD.M_min = 0;
 	  switch_wp_data(ii);
	  set_HOD_params();
	  /* chi2_wp_wrapper(a1); */
	  ngal[ii][i-1] = qromo(func_BW2,log(HOD.M_min),log(HOD.M_max),midpnt)/
	    GALAXY_DENSITY;
	  mmin[ii][i-1] = HOD.M_min;
	  mcut[ii][i-1] = HOD.M_cut;
	  m1[ii][i-1] = HOD.M1;
	  bb = qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/
	    GALAXY_DENSITY;
	  bias[ii] += bb;
	  ebias[ii] += bb*bb;
	  mone[ii][i-1] = zbrent(func_one_sat,HOD.M_min,HOD.M_max,1.0E-4);
	}

      /*
      if(i>n-nx)
	for(j=1;j<=k;++j)
	  ax[i-(n-nx)][j] = a[i][j];
      */
      if(i==1002)break;
    }
  fclose(fp);

  fp=fopen("m_one_true.dat","w");
  for(i=0;i<1000;++i)
    {
      for(j=0;j<9;++j)
	fprintf(fp,"%e %e ",mone[j][i],m1[j][i]);
      fprintf(fp,"\n");
    }
  fclose(fp);
  exit(0);

  for(j=0;j<wp.n_wp;++j)
    {
      for(i=x1=x2=x3=0;i<1000;++i)
	{
	  x1 += mmin[j][i];
	  x2 += m1[j][i];
	  x3 += mone[j][i];
	}
      printf("RATIO %f %f %e %e %f %e\n",wp_array[j].mean_mag,x2/x1,x1,x2,x3/x1,x3);
    }

  for(i=0;i<9;++i) 
    {
      bias[i]/=1000;
      ebias[i]/=1000;
      ebias[i] = sqrt(ebias[i] - bias[i]*bias[i]);
      printf("%d %f %f\n",i,bias[i],ebias[i]);
    }

  /*
  for(i=0;i<1000;++i)
    {
      k=1;
      for(j=2;j<9;++j)
	if(mcut[j][i]<0.9*mcut[j-1][i])k=0;
      if(k)printf("BOO %d\n",i);
    }
  exit(0);
  */


  fp=fopen("satfrac.dat","w");
  for(i=0;i<1000;++i)
    {
      for(j=0;j<9;++j)
	fprintf(fp,"%e %f ",ngal[j][i],log10(mmin[j][i]));
      fprintf(fp,"\n");
    }
  fclose(fp);

  fp=fopen("sat_params.dat","w");
  for(i=0;i<1000;++i)
    {
      for(j=0;j<9;++j)
	fprintf(fp,"%e %e ",m1[j][i],mcut[j][i]);
      fprintf(fp,"\n");
    }
  fclose(fp);

  exit(0);

  np=k;

  /*
  fprintf(stderr,"done with reading: %d points, %d params\n",n,k);

  for(k=0,i=1;i<=n;++i)
    { 
      ++k;
      for(j=1;j<=np;++j)
	{
	  if(j==26)printf("BOO %d %f\n",i,a[i][j]);
	}
    }
  */
  free_matrix(a,1,n,1,np);
}

void PVD_one_halo(int ngal, float ga[1000][6], int nwp, double mass)
{
  int i,j,k;
  float r,dx,dy,dz,dvz;
  static double v2[3][9], vbar[3][9],mbar[3][9];
  static int flag = 1, npair[3][9];

  if(flag)
    {
      flag=0;
      for(j=0;j<3;++j)
	for(i=0;i<9;++i)
	  mbar[j][i]=npair[j][i]=v2[j][i]=vbar[j][i]=0;
    }

  for(i=0;i<ngal;++i)
    for(j=i+1;j<ngal;++j)
      {
	dx = ga[i][0]-ga[j][0];
	dy = ga[i][1]-ga[j][1];
	dz = ga[i][2]-ga[j][2];
	r = sqrt(dx*dx+dy*dy+dz*dz);
	k=3;
	if(r>0.8788 && r<1.3573)
	  k = 2;
	if(r>0.3 && r<0.5)
	  k = 1;
	if(r<0.2)
	  k=0;
	if(k<3)
	  {
	    npair[k][nwp]++;

	    dz=fabs(ga[i][2]-ga[j][2]);
	    if(dz>0)
	      dvz=(ga[j][5]-ga[i][5])*(ga[j][2]-ga[i][2])/dz;
	    else
	      dvz=(ga[j][5]-ga[i][5]);
	    v2[k][nwp] += dvz*dvz;
	    vbar[k][nwp] += dvz;
	  }
      }
  mbar[k][nwp] += mass*ngal*(ngal-1)/2;
   
  if(ngal<0)
    {
      for(i=0;i<9;++i)
	{
	  for(j=0;j<3;++j)
	    printf("%8.2e %7.1f ",mbar[j][i]/npair[j][i],sqrt(v2[j][i]/npair[j][i]));
	  printf("\n");
	}
    }

}

void sigma_logM_compare()
{
  int i,j,k,nm=1000;
  double dlogm, mlo, mhi,m0,m1,m2,ncen,mass;

  for(j=wp.n_wp-1;j>=0;--j)
    {
      HOD.i_wp=j;
      HOD.M1 = wp_array[j].M1;
      HOD.M_cut = wp_array[j].M_cut;
      HOD.M_min = wp_array[j].M_min;
      HOD.sigma_logM = wp_array[j].sigma_logM;
      HOD.M_low = wp_array[j].M_low;

      ncen = N_cen(mass);

      mlo = log(HOD.M_low);
      mhi = log(HOD.M_max);
      dlogm = (mhi - mlo)/(nm - 1);

      m0 = m1 = m2 = 0;
      for(i=0;i<nm;++i)
	{
	  mass = exp(dlogm*i+mlo);
	  ncen = N_cen(mass);
	  m0 += ncen*dndM_interp(mass)*dlogm*mass;
	  m1 += ncen*pow(log(mass),1.0)*dndM_interp(mass)*dlogm*mass;
	  m2 += ncen*pow(log(mass),2.0)*dndM_interp(mass)*dlogm*mass;
	}
      printf("%.1f %f\n",-17-j*0.5,sqrt(m2/m0 - m1*m1/m0/m0));

    }
  exit(0);
}

double func_one_sat(double m)
{
  return(N_sat(m)-1);
}
