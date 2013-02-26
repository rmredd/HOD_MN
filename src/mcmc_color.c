#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"


/* External functions from wp_minimization.c
 */
void wp_input(void);
void wp_color_input(void);
double chi2_wp_color(double *a);
double chi2_wp_color_wrapper(double *a);

/* inteernal functions
 */
double mcmc_color_initialize(double *a, double **cov1, double *avg1, double *start_dev);


void mcmc_color_minimization()
{
  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*avg1,chi2,chi2prev,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,**chain,*start_dev,*eval_prev;
  int n,i,j,k,nrot,niter=0,count=0,imax_chain=100000,NSTEP=50,NSTEP_MAX=10000,convergence=0;
  long IDUM=-555;

  int *pcheck,pcnt,ptot=20,firstflag=1,*iweight,total_weight;
  double t0,tprev,temp,chi2a,chi2b;

  opar=dvector(1,100);

  MCMC=Task.MCMC;

  pcheck=calloc(ptot,sizeof(int));

  Work.imodel=2;
  Work.chi2=1;

  
  OUTPUT=1;
  

  fprintf(stderr,"\n\nMCMC OF W_P(R_P) COLOR  DATA..........\n");
  fprintf(stderr,    "--------------------------------------------\n\n");

  HOD.blue_fraction = 0.5857; /* <- Millenium fraction (sloan);; SDSS fraction -> 0.565; */
  HOD.blue_fraction = 0.6555; /* <- Millenium fraction (B-V>0.8) */
  HOD.blue_fraction = 0.565; /* SDSS -19,-20 */
  HOD.blue_fraction = 0.492; /* SDSS -20,-21 */
  HOD.blue_fraction = 0.379; /* SDSS -21 */

  wp_color.ON = 1;

  wp_color_input();


  srand48(32498793);

  /* Find the number of free parameters in the minimization
   * for the real-space correlation function.
   */
  for(n=0,i=1;i<12;++i)
    {
      n+=HOD.free[i];
      /* if(i>N_HOD_PARAMS && HOD.free[i])MCMC=3;*/
      if(OUTPUT)
	printf("mcmc_min> free[%i] = %d\n",i,HOD.free[i]);
    }

  // add 3 parameters for the color stuff
  n+=3;
  wp.ncf=n;

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

  //chi2prev=mcmc_initialize(a,cov1,avg1,start_dev);
  chi2prev=mcmc_color_initialize(a,cov1,avg1,start_dev);
  //initial_color_values(a,pp,yy);

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
      /* stepfac=0.7; */
      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];

      
      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
	  i=N_HOD_PARAMS;
	  if(HOD.free[++i])OMEGA_M = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];
	  /* if(HOD.free[++i])SIGV    = a[++j]; */
	}
      if(VBIAS_C<0)continue;

      /* Hard-wire CVIR variation
       */
      if(HOD.free[6])
	CVIR_FAC = a[3];

      chi2=chi2_wp_color_wrapper(a);	  

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
	printf("%e\n",chi2);fflush(stdout);
	printf("HSTATS %d %e %e %e %e\n",niter,HOD.M_min,number_weighted_halo_mass(),
	       number_weighted_central_mass(),
	       qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
	printf("FSAT %d %e %e %e %e\n",niter,HOD.M_min,wp.fsat_red,wp.fsat_blue,wp.fsat_all);

      }

    }
 RESTART_POINT:

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
	  //stepfac=1.6/sqrt(n);
	  //stepfac=0.6/sqrt(n);
	  //	  stepfac = stepfac*pow(0.9,6-j);
	  stepfac = 0.25;
	  stepfac = 0.5;
	  stepfac=1.6/sqrt(n);
	  if(!ThisTask)printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      stepfac=1.6/sqrt(n);
      //stepfac = 0;

      if(convergence)goto SKIP_MATRIX;
      // if(niter>NSTEP_MAX && niter%NSTEP_MAX!=0)goto SKIP_MATRIX;

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
      if(RESTART==4)convergence = 1;

      if(RESTART && count==0)stepfac=0;
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

      // CHeck that the chi2 for the last point in the restart chain
      // is recovered.
      /*
      if(RESTART && !count)
	for(i=1;i<=n;++i)
	  a[i] = aprev[i];
      */
      if(RESTART==5)
	{
	  j = count+1;
	  if(j>4000)exit(0);
	  for(i=1;i<=n;++i)
	    a[i] = chain[j][i];
	}

      /*
      if(!ThisTask)
	for(i=1;i<=n;++i)
	  {
	    printf("COV %d %d %e ",count,i,sqrt(eval[i]));
	    for(j=1;j<=n;++j)
	      printf("%e ",evect[j][i]);
	    printf("\n");
	  }
      */

      /* Using only variances
       */
      //for(i=1;i<=n;++i) 
      //	a[i] = aprev[i] + gasdev(&IDUM)*sqrt(tmp[i][i])*stepfac;

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=N_HOD_PARAMS;++i)if(HOD.free[i])j++;
	  i=N_HOD_PARAMS;
	  if(HOD.free[++i])OMEGA_M = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	  if(HOD.free[++i])GAMMA   = a[++j];
	  if(HOD.free[++i])SPECTRAL_INDX   = a[++j];
	  /* if(HOD.free[++i])SIGV    = a[++j]; */
	}
      if(VBIAS_C<0)continue;
      
      /* Hard-wire CVIR variation
       */
      if(HOD.free[6])
	CVIR_FAC = a[3];

      chi2=chi2_wp_color_wrapper(a);	  

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
	printf("%e\n",chi2);fflush(stdout);
	
	printf("FSAT %d %e %e %e %e\n",niter,HOD.M_min,wp.fsat_red,wp.fsat_blue,wp.fsat_all);
	if(MCMC==1)
	  {
	    printf("HSTATS %d %e %e %e %e\n",niter,HOD.M_min,number_weighted_halo_mass(),
		   number_weighted_central_mass(),
		   qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
	    
	    fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	  }
      }

    }
}

double chi2_wp_color_wrapper(double *a)
{
  static int flag=1;
  static double *b;
  int i,j;

  if(flag)
    {
      b=dvector(1,100);
      flag=0;
    }

  for(j=0,i=1;i<=N_HOD_PARAMS;++i) {
    if(HOD.free[i] && i!=5) { 
      if(a[++j]<=0) { printf("NEG %d %d %e\n",i,j,a[j]); return(1.0E7); } }
    if(HOD.free[i] && i==5) {
      ++j; }
  }

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

  // now the 3 color params
  ++j;
  b[j] = a[j];
  ++j;
  b[j] = a[j];
  ++j;
  b[j] = a[j];

  muh(1);

  return(chi2_wp_color(b));
}

double mcmc_color_initialize(double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0;
  double x1,x2,omega_m;
  long IDUM = -556;

  i=0;j=0;
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_min); start_dev[j]=0.001; }
  if(HOD.free[++i]){ a[++j]=log10(HOD.M1); start_dev[j] = 0.001; }
  if(HOD.free[++i]){ a[++j]=HOD.alpha; start_dev[j] = 0.03; }
  if(HOD.free[++i]){ a[++j]=log10(HOD.M_cut); start_dev[j] = 0.01; }
  if(HOD.free[++i]){ a[++j]=log10(HOD.sigma_logM); start_dev[j] = 0.01; }
  if(HOD.free[++i]){ a[++j]=CVIR_FAC; start_dev[j] = 0.02; }
  if(HOD.pdfc>=7){
    if(HOD.free[++i]){ a[++j]=log10(HOD.M_cen_max); start_dev[j] = 0.001; }}
  else {
    if(HOD.free[++i]) { a[++j]=log10(HOD.MaxCen); start_dev[j] = 0.001; }}
  
  a[++j] = HOD.fblue0_sat;
  start_dev[j] = 0.01;
  a[++j] = HOD.sigma_fblue_sat;
  start_dev[j] = 0.01;
  a[++j] = HOD.sigma_fblue_cen;
  start_dev[j] = 0.01;

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


  x1=chi2_wp_color_wrapper(a);	  
  x2 = 0;

  if(!ThisTask) {
    printf("TRY 0 ");
    for(i=1;i<=wp.ncf;++i)
      printf("%.4e ",a[i]);
    printf("%e\n",x1+x2);fflush(stdout);
    printf("INITIAL CHI2: %e %e\n",x1,x2);
    fflush(stdout);
  }
  return(x1+x2);
}

