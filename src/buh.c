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
  curlum_g9,
  OMEGA_TEMP;
int n_g9,
  VARY_GAMMA=0,
  VARY_CVIR=0,
  VARY_ALPHA=0;


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
double chi2_CNOC_data(int count);

double chi2_ml_wrapper(double *a);

void ml_function_test(int n_wp, double *a);
void ml_function_test_monte_carlo(int n_wp, double *a);
int poisson_deviate(double nave);
double poisson_prob(int n, double nave);
void dispersion_test(int n_wp, double *a, double mass);
double random_luminosity(double m1, double m2);
void calc_multiplicity_function(int n_wp, int n_start);
double chi2_multiplicity_function(int count);
void analyze_chain(void);

void output_wp_fits(int n_wp, double *a);

/* Local Functions for simulation population.
 */
double NFW_density(double r, double rs, double ps);
double NFW_velocity(double mass, double v[]);
void NFW_position(double mass, double x[]);
void ml_mock_catalog(int n_wp, double *a);
double random_2df_magnitude_error(void);
void ml_mock_catalog_gao(int n_wp, double *a);
double N_cen_gao(double m, int ii, double *mmin);
double N_cen_hiden(double m);

/* External functions.
 */
double chi2_wp_wrapper(double *a);

void ml_minimization()
{
  int i,j,k,n_wp=9,n_hod=0,n1,n2,n3;

  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*a1,*avg1,chi2,chi2prev,*start_dev,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,chi2array[15],**chain,x,x2,x3;
  int n,nrot,niter=0,count=0,imax_chain=30000,NSTEP = 50;
  long IDUM=-555;

  int *pcheck,pcnt,ptot=10,nancheck,astart[10];



  fprintf(stderr,"\n\nCHI2 MINIMIZATION OF M/L DATA..........\n");
  fprintf(stderr,    "--------------------------------------------\n\n");

  ml_input_data(n_wp);

  Work.imodel=2;
  Work.chi2=1;
  HOD.pdfc = 9;

  /* In Zehavi et al, this is 40 Mpc/h,
   * but Norberg et al use 50.
   */
  wp.pi_max=50.0;

  srand48(32498793);

  astart[0] = 1;
  for(i=1;i<=5;++i)
    astart[i] = astart[i-1] + 2;
  for(i=6;i<=8;++i)
    astart[i] = astart[i-1] + 3;

  n = astart[8]+2;
  if(HOD.free[3]) {
    VARY_ALPHA=1;
    HOD.free[3] = 0;
    n++;
  }
  if(HOD.free[6]) {
    VARY_CVIR=1;
    HOD.free[6] = 0;
    n++;
  }
  if(HOD.free[10]==2) {
    VARY_GAMMA=1;
    HOD.free[10]=0;
  }
  
  for(i=8;i<100;++i)
    if(HOD.free[i]) { n++; MCMC++; }
  OMEGA_TEMP = 1.651693e-01;

  printf("mcmc_min> %d  free parameters\n",n);

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

  chain=dmatrix(1,imax_chain,1,n);

  chi2prev=ml_initialize(n_wp,n_hod,a,cov1,avg1,start_dev);


  if(IDUM_MCMC==-1)
    mass2light_functions(n_wp,a);
  if(IDUM_MCMC==-2)
    ml_mock_catalog(n_wp,a);
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

  niter++;
  for(i=1;i<=n;++i)
    aprev[i] = a[i];

  for(i=1;i<=n;++i)
    chain[niter][i]=a[i];

  IDUM=IDUM_MCMC;

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
      if(VARY_GAMMA)
	stepfac = 0.7;

      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  i = 7;
	  j = astart[8]+2+VARY_ALPHA+VARY_CVIR;
	  if(HOD.free[++i])OMEGA_TEMP = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])SIGV    = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX = a[++j];

	  if(VARY_GAMMA)
	    GAMMA = gasdev(&IDUM)*0.02 + 0.15;
	  printf("COSMO %d %f %f %f %d\n",count+1,OMEGA_TEMP,SIGMA_8,GAMMA,j);
	}

      chi2=0;

      /* Restrict CVIR_FAC
       */
      i = astart[8]+3+VARY_ALPHA;
      if(VARY_CVIR)
	if(a[i]<0.3 || a[i]>1.2)continue;
      if(HOD.free[10])
	if(GAMMA<0.11 && GAMMA>0.19)continue;
      
      ++count;
      for(i=n_wp-1;i>=0;--i)
	{
	  n_hod = 3;
	  HOD.free[5] = 1;
	  if(i<6){ n_hod = 2; HOD.free[5] = 0; HOD.sigma_logM = 0.2; }
	  for(k=0,j=astart[i];j<=astart[i]+n_hod-1;++j)
	    {
	      a1[++k] = a[j];
	    }
	  if(VARY_ALPHA)
	    HOD.alpha = a[astart[8]+3];
	  if(VARY_CVIR)
	    CVIR_FAC = a[astart[8]+4];
	  wp.ncf = n_hod;

	  switch_wp_data(i);
	  chi2+=chi2array[i]=chi2_wp_wrapper(a1);

	  wp_array[i].M_min = HOD.M_min;
	  wp_array[i].sigma_logM = HOD.sigma_logM;
	  wp_array[i].M1 = HOD.M1;
	  wp_array[i].M_low = HOD.M_low;
	  wp_array[i].M_cut = HOD.M_cut;

	  if(!ThisTask)
	    {
	      printf("TRY%d %d ",i,count);
	      for(j=1;j<=n_hod;++j)
		printf("%.4e ",a1[j]);
	      printf("%e\n",chi2array[i]);
	      fflush(stdout);
	    }
	}
      chi2+=chi2_ml_function(count);
      chi2+=chi2_multiplicity_function(count);

      if(!ThisTask)
	{
	  printf("TRY_ALL %d %e %e\n",count,chi2,chi2prev);
	  fflush(stdout);
	}

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
	    n_hod = 3;
	    if(j<6){ n_hod = 2; } 
	    for(k=0,i=astart[j];i<=astart[j]+n_hod-1;++i)
	      printf("%e ",a[i]);
	    for(i=astart[8]+3;i<=n;++i)
	      printf("%e ",a[i]);
	    printf("%e %e\n",chi2array[j],chi2);fflush(stdout);
	  }
      }
    }

  stepfac=1.0;
  pcnt=-1;

  while(niter<imax_chain)
    {
      pcnt++;
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  stepfac = stepfac*pow(0.9,5-j);
	  if(!ThisTask)printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      if(VARY_GAMMA)
	stepfac = 0.7;


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

      for(i=1;i<=n;++i)
	printf("%d %e\n",i,eval[i]);
	
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
	  i=7;
	  j = astart[8]+2+VARY_ALPHA+VARY_CVIR;
	  if(HOD.free[++i])OMEGA_TEMP = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])SIGV    = a[++j];

	  if(VARY_GAMMA)
	    GAMMA = gasdev(&IDUM)*0.02 + 0.15;
	  printf("COSMO %f %f %f\n",OMEGA_TEMP,SIGMA_8,GAMMA);
	}

      chi2=0;

      /* Restrict CVIR_FAC and GAMMA
       */
      i = astart[8]+3+VARY_ALPHA;
      if(VARY_CVIR)
	if(a[i]<0.3 || a[i]>1.2)continue;
      if(HOD.free[10])
	if(GAMMA<0.11 && GAMMA>0.19)continue;

      ++count;
      for(i=n_wp-1;i>=0;--i)
	{	  
	  n_hod = 3;
	  HOD.free[5] = 1;
	  if(i<6){ n_hod = 2; HOD.free[5] = 0; HOD.sigma_logM = 0.2; }
	  for(k=0,j=astart[i];j<=astart[i]+n_hod-1;++j)
	    {
	      a1[++k] = a[j];
	    }
	  if(VARY_ALPHA)
	    HOD.alpha = a[astart[8]+3];
	  if(VARY_CVIR)
	    CVIR_FAC = a[astart[8]+4];
	  wp.ncf = n_hod;

	  switch_wp_data(i);
	  chi2+=chi2array[i]=chi2_wp_wrapper(a1);

	  wp_array[i].M_min = HOD.M_min;
	  wp_array[i].sigma_logM = HOD.sigma_logM;
	  wp_array[i].M1 = HOD.M1;
	  wp_array[i].M_low = HOD.M_low;
	  wp_array[i].M_cut = HOD.M_cut;

	  if(!ThisTask)
	    {
	      printf("TRY%d %d ",i,count);
	      for(j=1;j<=n_hod;++j)
		printf("%.4e ",a1[j]);
	      printf("%e\n",chi2array[i]);
	      fflush(stdout);
	    }
	}
      chi2+=chi2_ml_function(count);
      chi2+=chi2_multiplicity_function(count);

      if(!ThisTask)
	{
	  printf("TRY_ALL %d %e %e\n",count,chi2,chi2prev);
	  fflush(stdout);
	}

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
	    n_hod = 3;
	    if(j<6){ n_hod = 2; } 
	    for(k=0,i=astart[j];i<=astart[j]+n_hod-1;++i)
	      printf("%e ",a[i]);
	    for(i=astart[8]+3;i<=n;++i)
	      printf("%e ",a[i]);
	    printf("%e %e\n",chi2array[j],chi2);fflush(stdout);
	  }
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
  double mstar=-19.66,phi0=0.0161,alpha=-1.21,low_magnitude=-17.0;

  /*fp=openfile("galaxy_density_array.dat");*/

  for(i=0;i<n_wp;++i)
    {
      wp_array[i].mag1 = low_magnitude-0.5*i;
      wp_array[i].mag2 = low_magnitude-0.5*(i+1);
      wp_array[i].mean_mag = low_magnitude-0.5*(i+.5); /* TEMP */

      /* Old data -> no PCA, plus with innermost data point
       *
      sprintf(wp.fname_wp,"/gal/sr2/tinker/ml_ratio/data/xi_june04.%.02f.%.02f.keML.txt",
	      wp_array[i].mag2,wp_array[i].mag1);
      */

      sprintf(wp.fname_wp,"/gal/sr2/tinker/voids/data/half_mag_samples/cov.10001.%.02f.%.02f.keML.com.xi2.2.log2",
	      wp_array[i].mag2,wp_array[i].mag1);
      if(i==n_wp-1)
	sprintf(wp.fname_wp,"/gal/sr2/tinker/voids/data/threshold_samples/cov.10001.-22.00.-21.00.keML.com.xi2.2.log2");


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
      if(i==8)wp_array[i].mag2 = -23;
      wp_array[i].ngal = qromo(func_integrate_schechter,
			      wp_array[i].mag2,wp_array[i].mag1,midpnt);
      if(i==8)wp_array[i].ngal = 1.25e-4;
      fprintf(stderr,"%f %e\n",wp_array[i].mean_mag,wp_array[i].ngal);
      /*fscanf(fp,"%lf",&wp_array[i].ngal);*/
    }
  /*fclose(fp);*/

  wp.x=dvector(1,wp.np);
  wp.r=dvector(1,wp.np);
  wp.e=dvector(1,wp.np);

}

/* This is just the value of the Schechter Function fit to
 * 2F at magnitude m--> for getting total number density of galaxies.
 */
double func_integrate_schechter(double m)
{
  double mstar=-19.66,phi0=0.0161,alpha=-1.21;
  return(schechter(m,mstar,phi0,alpha));
}

/* This is the luminosity-weighted Schechter function value
 * at magnitdue m--> to get the total amount of light.
 */
double func_integrate_lum_schechter(double m)
{
  double mstar=-19.66,phi0=0.0161,alpha=-1.21,lum;
  lum = pow(10.0,-0.4*(m-MAGSOL));
  return(schechter(m,mstar,phi0,alpha)*lum);
}

/* Function to return the value of the Schechter function
 * given the parameters.
 */
double schechter(double mag, double mstar, double phi0,double alpha)
{
  return(0.4*log(10.0)*phi0*pow(10.0,-0.4*(mag-mstar)*(alpha+1))*
	 exp(-pow(10.0,-0.4*(mag-mstar))));
}


/* Routine to switch the global wp data variables with
 * the current dataset in the chi^2 loop.
 */
void switch_wp_data(int n)
{
  int i,j;

  i=n;
  for(j=1;j<=wp_array[i].np;++j)
    {
      wp.x[j]=wp_array[i].x[j];
      wp.r[j]=wp_array[i].r[j];
      wp.e[j]=wp_array[i].e[j];
    }
  wp.np=wp_array[i].np;
  GALAXY_DENSITY=wp_array[i].ngal;
  
  HOD.i_wp = n;
}



double ml_initialize(int nwp, int nhod, double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0,k,astart[10],np;
  double x1,x2,m1[50],alpha[50],mcut[50],sig[50],*a1,chi2,chi2array[15];
  long IDUM = -556;
  FILE *fp;

  astart[0] = 1;
  for(i=1;i<=5;++i)
    astart[i] = astart[i-1] + 2;
  for(i=6;i<=8;++i)
    astart[i] = astart[i-1] + 3;

  a1=dvector(1,nhod);

  fp=openfile("initial_values.dat");
  for(i=1;i<=nwp;++i)
    fscanf(fp,"%lf %lf %lf %lf",&m1[i],&alpha[i],&mcut[i],&sig[i]);
  fclose(fp);


  j=0;
  for(k=1;k<=nwp;++k)
    {
      i=0;
      if(HOD.free[++i])
	{
	  /*
	   * This one will eventually be set set to one.
	   *
	   */
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
      if(HOD.free[++i] && k>5)
	{
	  a[++j]=log10(sig[k]);
	  start_dev[j]=0.01;
	}
    }

  if(VARY_ALPHA)
    {
      a[++j]=HOD.alpha;
      start_dev[j]=0.01;
    }
      
  if(VARY_CVIR)
    {
      a[++j]=CVIR_FAC;
      start_dev[j]=0.01;
    }

  np = j;

  /* If using Powell's method or amoeba, then stop here
   */
  for(i=0;i<nwp;++i)
    {
      printf("IV %d ",i);
      nhod = 3;
      if(i<6){ nhod = 2; }
      for(j=astart[i];j<=astart[i]+nhod-1;++j)
	printf("%f ",a[j]);
      printf("%f ",HOD.alpha);
      printf("\n");
      fflush(stdout);
    }

  if(!MCMC)return 0;

  if(MCMC>1)
    {
      i=7;
      j=np;
      if(HOD.free[++i]){ a[++j]=OMEGA_TEMP; start_dev[j]=0.005; np++; }
      if(HOD.free[++i]){ a[++j]=SIGMA_8; start_dev[j]=0.01; np++; }
      if(HOD.free[++i]){ a[++j]=GAMMA; start_dev[j]=0.01; np++; }
      if(HOD.free[++i]){ a[++j]=VBIAS; start_dev[j]=0.01; np++; }
      if(HOD.free[++i]){ a[++j]=VBIAS_C; start_dev[j]=0.01; np++; }
      if(HOD.free[++i]){ a[++j]=SIGV; start_dev[j]=0.01; np++; }
    }

      
  for(i=1;i<=np;++i)
    {
      avg1[i]=a[i];
      printf("BEGIN %d %f\n",i,a[i]);
      for(j=1;j<=nhod*nwp;++j)
	cov1[i][j]=a[i]*a[j];
    }

  chi2=0;
  for(i=nwp-1;i>=0;--i)
    {
      nhod = 3;
      HOD.free[5] = 1;
      if(i<6){ nhod = 2; HOD.free[5] = 0; HOD.sigma_logM = 0.2; }
      wp.ncf=nhod;
      for(k=0,j=astart[i];j<=astart[i]+nhod-1;++j)
	{
	  a1[++k] = a[j];
	}
      if(VARY_ALPHA)
	HOD.alpha = a[astart[8]+3];
      if(VARY_CVIR)
	CVIR_FAC = a[astart[8]+4];
      printf("ALPHA = %f CVIR = %f\n",HOD.alpha,CVIR_FAC);
      
      switch_wp_data(i);

      chi2+=chi2array[i]=chi2_wp_wrapper(a1);
      printf("CHI %d %e\n",i,chi2array[i]);

      wp_array[i].M_min = HOD.M_min;
      wp_array[i].sigma_logM = HOD.sigma_logM;
      wp_array[i].M1 = HOD.M1;
      wp_array[i].M_low = HOD.M_low;
      wp_array[i].M_cut = HOD.M_cut;
    }

  chi2+=chi2_ml_function(0) + chi2_multiplicity_function(0);
  printf("INITIAL CHI2 %e\n",chi2);
  return(chi2);
}

/* This function takes the 2PIGG multiplicity function and returns a chi^2 value
 * based on the multplicity function implies by the current point in the chain.
 */
double chi2_multiplicity_function(int count)
{
  int i,j,k,nmass=200,nmulti=145,n_start,n_wp=9;
  double dlogm,mass,mlo,mhi,dndM,nsat,ncen,chi2a=0,chi2b=0,x,nlo,nhi;
  static double *multi,*dy,*ngal;
  FILE *fp;

  static int flag=0,NP;
  static float *N_2PIGG,*M1_2PIGG,*E1_2PIGG,*M2_2PIGG,*E2_2PIGG;

  if(!flag)
    {
      flag = 1;
      fp = openfile("/gal/sr2/tinker/ml_ratio/2PIGG_MULTI.dat");
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
      printf("MULTI%d %d %e %e %e %e\n",n_start,count,N_2PIGG[i],x,M1_2PIGG[i],E1_2PIGG[i]);
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
      printf("MULTI%d %d %e %e %e %e\n",n_start,count,N_2PIGG[i],x,M2_2PIGG[i],E2_2PIGG[i]);
    }

  /* Sum up the chi2 values and return.
   */
  printf("CHI2_MULTI %d %e %e %e\n",count,chi2a+chi2b,chi2a,chi2b);
  return(chi2a+chi2b);
}

/* This function takes in the 2PIGG (corrected) data and returns a chi^2 value
 * of the M/L-L curve 
 */
double chi2_ml_function(int count)
{
  FILE *fp;
  int i,j,k,n=100,n_wp=9;
  double *lum,*sig,mlo,mhi,dlogm,*mass,*ml_lum,*y_lum, *y_sig, *a,chi2,ml;
  double t0,t1;

  static int flag=0,N_2PIGG,N_LENSING;
  static float *ML_2PIGG, *L_2PIGG,*E1_2PIGG,*E2_2PIGG,*ML_LENSING,*L_LENSING,*E_LENSING;

  /* return(0); */

  if(!flag)
    {
      flag = 1;
      fp = openfile("/gal/sr2/tinker/ml_ratio/2PIGG_ML_linear.dat");
      N_2PIGG = filesize(fp);
      muh(N_2PIGG);

      L_2PIGG = vector(1,N_2PIGG);
      E1_2PIGG = vector(1,N_2PIGG);
      E2_2PIGG = vector(1,N_2PIGG);
      ML_2PIGG = vector(1,N_2PIGG);

      for(i=1;i<=N_2PIGG;++i)
	fscanf(fp,"%f %f %f %f",&L_2PIGG[i],&ML_2PIGG[i],&E1_2PIGG[i],&E2_2PIGG[i]);
      fclose(fp);

      fp = openfile("/gal/sr2/tinker/ml_ratio/LENSING_ML.dat");
      N_LENSING = filesize(fp);
      muh(N_LENSING);

      L_LENSING = vector(1,N_LENSING);
      E_LENSING = vector(1,N_LENSING);
      ML_LENSING = vector(1,N_LENSING);

      for(i=1;i<=N_LENSING;++i)
	fscanf(fp,"%f %f %f",&L_LENSING[i],&ML_LENSING[i],&E_LENSING[i]);
      fclose(fp);

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
  for(i=1;i<=N_2PIGG;++i)
    {
      curlum_g9 = L_2PIGG[i];
      ml = qromo(func_ml_lum,mlo,mhi,midpnt)/qromo(func_ml_norm,mlo,mhi,midpnt);
      ml = ml*OMEGA_TEMP/OMEGA_M;
      if(ml>ML_2PIGG[i])
	chi2 += (ml-ML_2PIGG[i])*(ml-ML_2PIGG[i])/E1_2PIGG[i]/E1_2PIGG[i];
      else
	chi2 += (ml-ML_2PIGG[i])*(ml-ML_2PIGG[i])/E2_2PIGG[i]/E2_2PIGG[i];
      printf("MLC %d %e %e %e\n",count,curlum_g9,ML_2PIGG[i],ml);
    }
  for(i=1;i<=N_LENSING;++i)
    {
      curlum_g9 = L_LENSING[i];
      ml = qromo(func_ml_lum,mlo,mhi,midpnt)/qromo(func_ml_norm,mlo,mhi,midpnt);
      ml = ml*OMEGA_TEMP/OMEGA_M;
      chi2 += (ml-ML_LENSING[i])*(ml-ML_LENSING[i])/E_LENSING[i]/E_LENSING[i];
      printf("MLC %d %e %e %e\n",count,curlum_g9,ML_LENSING[i],ml);
    }

  printf("ML_CHI2 %d %e\n",count,chi2);
  chi2+=chi2_CNOC_data(count);
  return(chi2);
}

double chi2_CNOC_data(int count)
{
  static int flag=0,n;
  static float *mass,x1,x2,lum_avg=0;
  double chi2,dx1,*ax1,lum,ml_ave=0;
  FILE *fp;
  int i;

  if(!flag)
    {
      flag = 1;
      fp = openfile("/gal/sr2/tinker/ml_ratio/CNOC_ML.dat");
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
  printf("MLC %d %e %e %e %f\n",count,lum_avg,ml_ave,382.0,OMEGA_TEMP/OMEGA_M);
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

  mlo=log(1.0e11);
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
  if(ii==8)return(n1);

  for(n2=0,i=ii+1;i<9;++i)
    if(m>wp_array[i].M_low)
      {
	n2=0.5*(1+erf((logm - log10(wp_array[i].M_min))/wp_array[i].sigma_logM));
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
  double *a,**pp,*yy,FTOL=1.0E-3,chi2min,s1,dlogm,m,d[100],*avg1,**cov1,*start_dev,*a1;
  FILE *fp;
  char aa[1000];

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
      if(i<=7)
	n+=HOD.free[i]*n_wp;
      else
	n+=HOD.free[i];
      if(i<=7)
	n_hod+=HOD.free[i];

      if(OUTPUT)
	printf("mcmc_min> free[%i] = %d\n",i,HOD.free[i]);
    }
  wp.ncf=n_hod;
  wp.ncf_tot=n;
  a=dvector(1,n);
  start_dev=dvector(1,n);

  ml_initialize(n_wp,n_hod,a,cov1,avg1,start_dev);

  astart[0] = 1;
  for(i=1;i<=5;++i)
    astart[i] = astart[i-1] + 2;
  for(i=6;i<=8;++i)
    astart[i] = astart[i-1] + 3;

  for(i=1;i<=n;++i) 
    a[i]=pow(10.0,a[i]);

  for(i=1;i<=n;++i)
    printf("a[%d] = %e\n",i,a[i]);

  n_mag=n_wp-1;



  /* Loop through all magnitudes.
   */
 MAGNITUDE_LOOP:

  n = n_hod = 3;
  HOD.free[5] = 1;
  if(n_mag<6){ n = n_hod = 2; HOD.free[5] = 0; HOD.sigma_logM = 0.2;}

  n = n_hod;
  a1=dvector(1,n);
  if(POWELL)
    pp=dmatrix(1,n,1,n);
  else
    pp=dmatrix(1,n+1,1,n);
  yy=dvector(1,n+1);

  for(k=0,j=astart[n_mag];j<=astart[n_mag]+n_hod-1;++j)
    a1[++k] = a[j];

  switch_wp_data(n_mag);
  
  /* Make the starting stepsize 10% of the initial values.
   */
  for(i=1;i<=n;++i)
    d[i]=a1[i]*0.1;

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
      yy[1]=chi2_ml_wrapper(a1);
    
      for(i=1;i<=n;++i)
	{
	  a1[i]+=d[i];
	  if(i>1)a1[i-1]-=d[i-1];
	    yy[i+1]=chi2_ml_wrapper(a1);	  
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
  printf("%e\n",GALAXY_BIAS);

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
    free_dmatrix(pp,1,n+1,1,n);

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
  FILE *fp,*fpa[9],*fp2,*fpb[9];
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[9][1000],halocnt[1000],imag;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2;
  char aa[1000];
  float x1,xh[3],vh[3];
  long IDUM3 = -445;

  fp=openfile("/gal/sr2/tinker/voids/LANL/halo.LANL.400");
  for(i=j_start;i<n_wp;++i)
    {
      sprintf(aa,"%s.mock.%d",Task.root_filename,i);
      fpa[i] = fopen(aa,"w");
      sprintf(aa,"%s.mag.%d",Task.root_filename,i);
      fpb[i] = fopen(aa,"w");
    }

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

  while(!feof(fp))
    {
      fscanf(fp,"%d %d %e %e %e %e %e %e %e",
	     &i,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
      mass=imass*RHO_CRIT*OMEGA_M*pow(0.3125,3.0);
      if(mass<mlo)continue;

      i1 = (int)(log10(mass)/0.1);

      for(j=j_start;j<n_wp;++j)
	{
	  HOD.i_wp=j;
	  HOD.M_min=wp_array[j].M_min;
	  HOD.M1 = pow(10.0,a[j+0*n_wp+1]);
	  HOD.M_cut = pow(10.0,a[j+1*n_wp+1]);
	  HOD.sigma_logM = pow(10.0,a[j+2*n_wp+1]);
	  HOD.alpha = 1.0;
	  
	  nsat = N_sat(mass);
	  if(j>j_start)
	    nc[j] = nc[j-1]+N_cen(mass);
	  else
	    nc[j] = N_cen(mass);

	  n1 = poisson_deviate(nsat);
	  for(i=1;i<=n1;++i)
	    {
	      NFW_position(mass,xg);
	      NFW_velocity(mass,vg);
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
	      /* random_2df_magnitude_error(); */
	      imag = (-mag-17)*2;
	      if(imag>=n_wp)imag = n_wp;

	      /* Bin up the galaxies by halo mass to check the HOD
	       */
	      galcnt[imag][i1]++;
	    }
	}
      ncen=drand48();
      for(j=j_start;j<n_wp;++j)
	if(ncen<nc[j])break;
      if(j<n_wp) {
	fprintf(fpa[j],"%e %e %e %e %e %e %d\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2],j);
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

void NFW_position(double mass, double x[])
{
  double r,pr,max_p,costheta,sintheta,phi1,signs,rvir,rs,cvir;
  
  cvir=halo_concentration(mass)*CVIR_FAC;
  rvir=pow(3*mass/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rs=rvir/cvir;
  max_p=NFW_density(rs,rs,1.0)*rs*rs*4.0*PI;

  for(;;) {
    r=drand48()*rvir;
    pr=NFW_density(r,rs,1.0)*r*r*4.0*PI/max_p;
    
    if(drand48()<=pr)
      {
	costheta=2.*(drand48()-.5);
	sintheta=sqrt(1.-costheta*costheta);
	signs=2.*(drand48()-.5);
	costheta=signs*costheta/fabs(signs);
	phi1=2.0*PI*drand48();
	
	x[0]=r*sintheta*cos(phi1);
	x[1]=r*sintheta*sin(phi1);
	x[2]=r*costheta;
	return;
      }
  }
}

/* This is the NFW density profile
 */
double NFW_density(double r, double rs, double ps)
{
  return(ps*rs/(r*(1+r/rs)*(1+r/rs)));
}

/* This sets the velocity to be isotropic Gaussian.
 */
double NFW_velocity(double mass, double v[])
{
  static long IDUM2=-455;
  static double fac = -1;
  double sigv;
  int i;

  if(fac<0)
      fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS;
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

int poisson_deviate(double nave)
{
  static int flag=0;
  double p,pp;
  int n;

  p=0;
  pp=1;
  if(!flag){
    srand48(42398);
    flag++;
  }

  while(p<pp)
    {
      if(nave<1)
	n=(int)(drand48()*20);
      else
	n=(int)(drand48()*30*nave);
      p=poisson_prob(n,nave);
      pp=drand48();
    }
  return(n);
}

double poisson_prob(int n, double nave)
{
  int i;
  double fac=1;

  if(n>0)
    for(i=1;i<=n;++i)
      fac*=nave/i;

  return((float)(fac*exp(-nave)));
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
  int i,j,nr=50;
  double *x[10],rlo=0.1,rhi=50,dlogr,r;
  FILE *fp;
  char fname[100];

  dlogr = (log(rhi)-log(rlo))/(nr-1);

  for(i=0;i<n_wp;++i)
    x[i] = dvector(1,nr);

  for(i=0;i<n_wp;++i)
    {
      HOD.i_wp=i;
      HOD.M_low = wp_array[i].M_low;
      HOD.M_min = wp_array[i].M_min;
      HOD.M1 = wp_array[i].M1;
      HOD.M_cut = wp_array[i].M_cut;
      HOD.sigma_logM = wp_array[i].sigma_logM;
      GALAXY_DENSITY=wp_array[i].ngal;
      /*
      printf("%d %e %e %e %e %e %e\n",
	     i,HOD.M_low,HOD.M_min,HOD.M1,HOD.alpha,HOD.M_cut,HOD.sigma_logM);
      */
      RESET_FLAG_1H = RESET_FLAG_2H = 1;
 
      for(j=1;j<=nr;++j)
	{
	  r = exp(dlogr*(j-1))*rlo;
	  x[i][j] = projected_xi(r);
	}
    }

  sprintf(fname,"%s.wp_fits",Task.root_filename);
  fp=fopen(fname,"w");

  for(j=1;j<=nr;++j)
    {
      fprintf(fp,"%e ",exp(dlogr*(j-1))*rlo);
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
	      NFW_velocity(mass,vg);
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
	      x3 = log10(0.2);
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
      printf("CHI2_CNOC %e\n",chi2_CNOC_data(i));
    }

  exit(0);
}
