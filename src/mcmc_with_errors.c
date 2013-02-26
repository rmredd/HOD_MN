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
double mcmc_initialize(double *a, double **cov1, double *avg1);


/* Internal functions.
 */
void choose_bias_fit(void);
double chi2_wp_wrapper(double *);
void choose_dndM_fit(void);

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
 * [7] ->  OMEGA_M
 * [8] ->  SIGMA_8
 * [9] ->  VBIAS
 * 
 */
void mcmc_minimization()
{
  int EJENK=1,EBIAS=0;
  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*avg1,chi2,chi2prev,
    **evect,*eval,*aprev,*atemp,**tmp1;
  int n,i,j,nrot,niter=0,count=0;
  long IDUM=-555;

  float original_jenkins_a,
    original_jenkins_b,
    original_jenkins_c;

  int *pcheck,pcnt,ptot=10;

  original_jenkins_a=JENKINS_A;
  original_jenkins_b=JENKINS_B;
  original_jenkins_c=JENKINS_C;

  pcheck=calloc(ptot,sizeof(int));

  if(MCMC>1)
    wp.esys=0.08;
    
  wp_input();

  Work.imodel=2;
  Work.chi2=1;
  MCMC=Task.MCMC;

  OUTPUT=0;

  srand48(32498793);

  /* Find the number of free parameters in the minimization
   * for the real-space correlation function.
   */
  for(n=0,i=1;i<=10;++i)
    {
      n+=HOD.free[i];
      if(OUTPUT)
	printf("mcmc_min> free[%i] = %d\n",i,HOD.free[i]);
    }
  wp.ncf=n;

  if(OUTPUT)
    printf("mcmc_min> %d  free parameters\n",n);

  a=dvector(1,n);
  aprev=dvector(1,n);
  atemp=dvector(1,n);
  cov1=dmatrix(1,n,1,n);
  avg1=dvector(1,n);

  tmp=dmatrix(1,n,1,n);
  tmp1=dmatrix(1,n,1,1);
  evect=dmatrix(1,n,1,n);
  eval=dvector(1,n);

  chi2prev=mcmc_initialize(a,cov1,avg1);
  niter++;
  for(i=1;i<=n;++i)
    aprev[i] = a[i];

  IDUM=IDUM_MCMC;

  pcnt=0;
  pcheck[pcnt]=1;

  stepfac=1;
  while(niter<10)
    {
      pcnt++;
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  stepfac = stepfac*pow(0.9,5-j);
	  printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*0.004*stepfac)*aprev[i];

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=6;++i)if(HOD.free[i])j++;
	  i=6;
	  if(HOD.free[++i])OMEGA_M = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	}

      /* Take the parameters of the Jenkins
       * mass function from Gaussian distributions.
       */
      /*
      RESET_COSMOLOGY++;
      if(EJENK) 
	choose_dndM_fit();
      if(EBIAS)
	choose_bias_fit();
      */

      chi2=chi2_wp_wrapper(a);
      if(MCMC>1)chi2+=chi2_zspace(a);

      printf("TRY %d ",++count);
      for(i=1;i<=n;++i)
	printf("%.4e ",a[i]);
      printf("%e\n",chi2);fflush(stdout);

      pcheck[pcnt]=0;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	continue;
      pcheck[pcnt]=1;

      niter++;
      for(i=1;i<=n;++i)
	avg1[i] += a[i];
      for(i=1;i<=n;++i)
	aprev[i] = a[i];
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cov1[i][j] += a[i]*a[j];
      chi2prev=chi2;

      printf("ACCEPT %d %d ",niter,count);
      for(i=1;i<=n;++i)
	printf("%e ",a[i]);
      printf("%e\n",chi2);fflush(stdout);

    }
  stepfac=1.0;
  pcnt=-1;
  while(error>tolerance)
    {
      pcnt++;
      if(pcnt==ptot)
	{
	  for(j=i=0;i<ptot;++i)j+=pcheck[i];
	  stepfac = stepfac*pow(0.9,5-j);
	  printf("STEPFAC %f %d %d\n",stepfac,j,count);
	  pcnt=0;
	}
      stepfac=1;

      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  tmp[i][j] = cov1[i][j]/niter - avg1[i]*avg1[j]/niter/niter;

      jacobi(tmp,n,eval,evect,&nrot);
      gaussj(evect,n,tmp1,1);

      for(i=1;i<=n;++i)
	atemp[i] = gasdev(&IDUM)*sqrt(eval[i])*stepfac;

      for(i=1;i<=n;++i)
	for(a[i]=0,j=1;j<=n;++j)
	  a[i] += atemp[j]*evect[j][i];

      for(i=1;i<=n;++i)
	a[i] += aprev[i];

      printf("BOO %d %f %f\n",niter,a[1],a[2]);
      fflush(stdout);

      if(MCMC>1)
	{
	  RESET_COSMOLOGY++;
	  j=0;
	  for(i=1;i<=6;++i)if(HOD.free[i])j++;
	  i=6;
	  if(HOD.free[++i])OMEGA_M = a[++j];
	  if(HOD.free[++i])SIGMA_8 = a[++j];
	  if(HOD.free[++i])VBIAS   = a[++j];
	  if(HOD.free[++i])VBIAS_C = a[++j];
	}

      /* Take the parameters of the Jenkins
       * mass function from Gaussian distributions.
       */
      
      RESET_COSMOLOGY++;
      if(EJENK) 
	{
	  /*
	  JENKINS_A = original_jenkins_a*(1+gasdev(&IDUM)*sqrt(3.919662e-07));
	  JENKINS_B = original_jenkins_b*(1+gasdev(&IDUM)*sqrt(9.265636e-06));
	  JENKINS_C = original_jenkins_c*(1+gasdev(&IDUM)*sqrt(2.365370e-03));
	  choose_dndM_fit();
	  */
	}
      if(EBIAS)
	choose_bias_fit();

      chi2=chi2_wp_wrapper(a);
      if(MCMC>1)chi2+=chi2_zspace(a);

      printf("TRY %d ",++count);
      for(i=1;i<=n;++i)
	printf("%.4e ",a[i]);
      printf("%e\n",chi2);fflush(stdout);

      pcheck[pcnt]=0;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	continue;
      pcheck[pcnt]=1;

      niter++;
      for(i=1;i<=n;++i)
	avg1[i] += a[i];
      for(i=1;i<=n;++i)
	aprev[i] = a[i];
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cov1[i][j] += a[i]*a[j];
      chi2prev=chi2;

      printf("ACCEPT %d %d ",niter,count);
      for(i=1;i<=n;++i)
	printf("%e ",a[i]);
      printf("%e\n",chi2);fflush(stdout); 
    }
}

double chi2_wp_wrapper(double *a)
{
  static int flag=1;
  static double *b;
  int i,j;

  if(flag)
    {
      b=dvector(1,100);
      flag=0;
    }

  for(j=0,i=1;i<=6;++i)
    if(HOD.free[i])
      if(a[++j]<=0)return(1.0E7);

  i=0;j=0;
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);}
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);}
  if(HOD.free[++i]){j++;b[j]=a[j];}
  if(HOD.free[++i]){j++;b[j]=pow(10.0,a[j]);}
  if(HOD.free[++i]){j++;b[j]=a[j];}
  if(HOD.free[++i]){j++;b[j]=a[j];}

  return(chi2_wp(b));
}

double mcmc_initialize(double *a, double **cov1, double *avg1)
{
  int i,j=0;
  double x1,x2;
  long IDUM = -556;

  i=0;j=0;
  if(HOD.free[++i])a[++j]=log10(HOD.M_min);
  if(HOD.free[++i])a[++j]=log10(HOD.M1);
  if(HOD.free[++i])a[++j]=HOD.alpha;
  if(HOD.free[++i])a[++j]=log10(HOD.M_cut);
  if(HOD.free[++i])a[++j]=HOD.sigma_logM;
  if(HOD.free[++i])a[++j]=CVIR_FAC;
  if(MCMC>1)
    {
      if(HOD.free[++i])a[++j]=OMEGA_M;
      if(HOD.free[++i])a[++j]=SIGMA_8;
      if(HOD.free[++i])a[++j]=VBIAS;
      if(HOD.free[++i])a[++j]=VBIAS_C;
    }
  printf("INITIAL VALUES: ");
  for(i=1;i<=wp.ncf;++i)printf("%e ",a[i]);
  printf("\n");
  
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
      for(i=1;i<=6;++i)if(HOD.free[i])j++;
      i=6;
      if(HOD.free[++i])OMEGA_M = a[++j];
      if(HOD.free[++i])SIGMA_8 = a[++j];
      if(HOD.free[++i])VBIAS   = a[++j];
      if(HOD.free[++i])VBIAS_C = a[++j];
    }

  x1=chi2_wp_wrapper(a);
  
  if(MCMC>1)
    x2=chi2_zspace(a);
  else
    x2=0;
  
  printf("TRY 0 ");
  for(i=1;i<=wp.ncf;++i)
    printf("%.4e ",a[i]);
  printf("%e\n",x1+x2);fflush(stdout);

  printf("INITIAL CHI2: %e %e\n",x1,x2);
  fflush(stdout);
  return(x1+x2);
}

void choose_bias_fit()
{
  static long IDUM1=-444;
  static int flag=1,n;
  static float *a,*b,*c;

  FILE *fp;
  int i;
  char string[1000];

  if(flag)
    {
      flag=0;
      fp=openfile("/home/tinker/TABLES/bias_errors.dat");
      n=filesize(fp);
      a=vector(1,n);
      b=vector(1,n);
      c=vector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f %f",&a[i],&b[i],&c[i]);
	  fgets(string,1000,fp);
	}
    }

  i=(int)(ran1(&IDUM1)*n) + 1;
  BIAS_A = a[i];
  BIAS_B = b[i];
  BIAS_C = c[i];
  /*
  printf("BIAS %f %f %f\n",a[i],b[i],c[i]);
  */
}

void choose_dndM_fit()
{
  static long IDUM1=-444;
  static int flag=1,n;
  static float *a,*b,*c,*d,*e;

  FILE *fp;
  int i;
  char string[1000];

  if(flag)
    {
      flag=0;
      fp=openfile("/home/tinker/TABLES/dndM_errors.dat");
      n=filesize(fp);
      a=vector(1,n);
      b=vector(1,n);
      c=vector(1,n);
      d=vector(1,n);
      e=vector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f %f %f %f",&a[i],&b[i],&c[i],&d[i],&e[i]);
	  fgets(string,1000,fp);
	}
    }

  i=(int)(ran1(&IDUM1)*n) + 1;
  DNDM_PARAMS[1] = a[i];
  DNDM_PARAMS[2] = b[i];
  DNDM_PARAMS[3] = c[i];
  DNDM_PARAMS[4] = d[i];
  DNDM_PARAMS[5] = e[i];
}
