#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* external functions
 */
double func_findmshift(double x);
double func_drg1(double m);
double func_drg1a(double m);
double chi2wtheta_covar(float *d, float *x);

/* Internal functions.
 */
double chi2_wtheta_wrapper(double *a);
double mcmc_wtheta_initialize(double *a, double **cov1, double *avg1, double *start_dev);
double chi2_wtheta(void);

int USE_IWEIGHT1 = 0;

/******************************************************************
 *
 * HOD.free[] also controls which variables will be held constant/vary
 * during MCMC minimization. Since this routine will also so z-space
 * minimization if requested, indices>6 are cosmological.
 *
 *  i     variable
 * ---    --------
 * [1] ->  Mshift
 * [2] ->  fsat
 * [3] ->  drg_number_density
 * 
 */
void mcmc_wtheta()
{
  double stepfac=1;
  double error=1,tolerance=0,**cov1,**tmp,*a,*avg1,chi2,chi2prev,
    **evect,*eval,*aprev,*atemp,**tmp1,*opar,x1,fsat,**chain,*start_dev,*eval_prev;
  int n,i,j,k,nrot,niter=0,count=0,imax_chain=100000,NSTEP=50,NSTEP_MAX=10000,convergence=0;
  long IDUM=-555;

  int *pcheck,pcnt,ptot=20,firstflag=1,*iweight,total_weight;
  double t0,tprev,temp,chi2a,chi2b;

  // initialize various things
  NGAL_DRG = 6.5e-4;

  // get mean mass of centrals (ALL)
  MALL = qromo(func_drg1,log(HOD.M_low),log(HOD.M_max),midpnt)/
    qromo(func_drg1a,log(HOD.M_low),log(HOD.M_max),midpnt);

  //reset the central HOD
  HOD.pdfc = 10;
  //HOD.shift_alpha = 1.0;

  opar=dvector(1,100);

  MCMC=Task.MCMC;

  pcheck=calloc(ptot,sizeof(int));

  srand48(32498793);

  if(ARGC>3)
    IDUM = atoi(ARGV[3]);

  n = 3;
  wp.ncf=n;

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

  //IDUM=IDUM_MCMC;


  chi2prev=mcmc_wtheta_initialize(a,cov1,avg1,start_dev);
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
      stepfac=1;
      for(i=1;i<=n;++i)
	a[i] = (1+gasdev(&IDUM)*start_dev[i]*stepfac)*aprev[i];

      chi2=chi2_wtheta_wrapper(a);
      if(chi2>9e6)continue;

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
	  if(USE_IWEIGHT1)
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
	printf("%e %e ",HOD.mass_shift,GALAXY_DENSITY);
	printf("%e\n",chi2);fflush(stdout);
      }

    }

  stepfac=1.6/sqrt(n);
  pcnt=-1;
  t0 = second();

  NSTEP = niter;

  while(niter<imax_chain)
    {
      stepfac=2.5/sqrt(n);

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

      for(i=1;i<=n;++i)
	atemp[i] = gasdev(&IDUM)*sqrt(eval[i])*stepfac;

      for(i=1;i<=n;++i)
	for(a[i]=0,j=1;j<=n;++j)
	  a[i] += atemp[j]*evect[j][i];

      for(i=1;i<=n;++i) 
	a[i] += aprev[i];

      chi2=chi2_wtheta_wrapper(a);
      if(chi2>9e6)continue;


      tprev = t0;
      t0 = second();
      ++count;
      printf("TRY %d ",count);
      for(i=1;i<=n;++i)
	printf("%.4e ",a[i]);
      printf("%e %.2f\n",chi2,
	     timediff(tprev,t0));fflush(stdout);
      
      pcheck[pcnt]=0;
      if(!(chi2<chi2prev || drand48() <= exp(-(chi2-chi2prev)/2)))
	{
	  if(USE_IWEIGHT1)
	    iweight[niter+1]++;
	  continue;
	}
      pcheck[pcnt]=1;

      niter++;
      if(!convergence)NSTEP = niter;
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

      if(!ThisTask) {
	printf("ACCEPT %d %d ",niter,count);
	for(i=1;i<=n;++i)
	  printf("%e ",a[i]);
	printf("%e %e ",HOD.mass_shift,GALAXY_DENSITY);
	printf("%e\n",chi2);fflush(stdout);
      }

    }
}

double chi2_wtheta_wrapper(double *a)
{
  static int flag=1;
  static double *b;
  int i,j;

  if(a[1]<0)return 1E7;
  if(a[2]<0)return 1E7;
  if(a[3]<0)return 1E7;
  if(a[2]>1 || a[3]>1)return 1E7;

  MSHIFT = pow(10.0,a[1]);
  HOD.freds = a[2];
  HOD.fredc = a[3];
  HOD.mass_shift = zbrent(func_findmshift,0.0,10.0,1.0E-4);

  GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
  fmuh(GALAXY_DENSITY);


  // do a 2-sigma clipping on the number density
  //if(fabs(GALAXY_DENSITY - 0.00065)/(0.13*.00065)>2)return 2E7;

  RESET_FLAG_1H++;
  RESET_FLAG_2H++;

  return(chi2_wtheta());
}

double chi2_wtheta()
{
  int k,i,n=40,j,i1,i2,ngal;
  double dlogr,r,mass,xk,pnorm,psp,rm,sig,t0,t1,pnorm1,mp,mlo,mhi,dr,
    slo,shi,dsdM,rlo,rhi,dlogm,f1,f2,fac,cvir,rvir,rs,r1,r2,mtot=0,m,
    chi2lf,chi2w,x0,mlow;
  FILE *fp,*fp2;
  float x1,x2,x3,x4;
  char aa[1000],fname[1000],**argv;

  float *xx,*yy,*zz,*vx,*vy,*vz;
  int magcnt[100];

  double tmin, tmax, dtheta, theta, delta_mag, maglo, mred, mblue;


  static int nwdata, nlf;
  static float *wdata,*rdata,*edata,*lfdata,*lferr,*lfmag,*wmodel;
  static int niter = -1;

  BOX_SIZE = 160.;

  if(niter==-1)
    {
      fp = openfile("wtheta_corrected.data");
      nwdata = filesize(fp);
      wdata = vector(1,nwdata);
      rdata = vector(1,nwdata);
      edata = vector(1,nwdata);
      wmodel = vector(1,nwdata);
      for(i=1;i<=nwdata;++i)
	fscanf(fp,"%f %f %f",&rdata[i],&wdata[i],&edata[i]);
      fclose(fp);
      
      fp = openfile("LF_DRG.data");
      nlf = filesize(fp);
      lfmag = vector(1,nlf);
      lfdata = vector(1,nlf);
      lferr = vector(1,nlf);
      for(i=1;i<=nlf;++i)
	fscanf(fp,"%f %f %f",&lfmag[i],&lfdata[i],&lferr[i]);
      fclose(fp);
      delta_mag = 0.45;
      maglo = -23.86-delta_mag/2;
    }
  niter ++;


  // calculate the chi^2 for the wtheta values. ( append values on the end for plotting purposes)
  chi2w = 0;
  printf("WTH %d %e %e %e %e\n",niter,1.0,0.0,wtheta(1.0),1.0);
  for(i=1;i<=nwdata;++i)
    {
      x0 = wtheta(rdata[i]);
      wmodel[i] = x0;
      printf("WTH %d %e %e %e %e\n",niter,rdata[i],wdata[i],x0,edata[i]);
      chi2w += (wdata[i]-x0)*(wdata[i]-x0)/(edata[i]*edata[i]);
    }

  //use covariance matrix
  chi2w = chi2wtheta_covar(wdata,wmodel);

  // what's mass of blue-fraction
  HOD.freds = 1-HOD.freds;
  HOD.pdfc = 11;
  mblue = qromo(func_drg1,log(HOD.M_low),log(HOD.M_max),midpnt)/
    qromo(func_drg1a,log(HOD.M_low),log(HOD.M_max),midpnt);
  HOD.freds = 1-HOD.freds;
  HOD.pdfc = 10;
  mred = qromo(func_drg1,log(HOD.M_low),log(HOD.M_max),midpnt)/
    qromo(func_drg1a,log(HOD.M_low),log(HOD.M_max),midpnt);
  

  printf("WTH %d %e %e %e %e\n",niter,1200.0,0.0,wtheta(1200.),1.0);
	  
  // use the number density as the other constraint
  if(ARGC>4)
    {
      chi2lf = (GALAXY_DENSITY-NGAL_DRG)*(GALAXY_DENSITY-NGAL_DRG)/(0.13*0.13*NGAL_DRG*NGAL_DRG);
      
      printf("CHI %d %f %f %f %e %e %f %e %f\n",niter,HOD.freds,HOD.fredc,HOD.mass_shift,chi2lf,chi2w,MSHIFT,GALAXY_DENSITY,mred/mblue);
      fflush(stdout);
      
      return chi2w+chi2lf;
    }


  // make the sham-galaxy distribution
  sprintf(fname,"sham ../../../SHAM/halosub_0.284 hod_mshift %f %f %f %f > galtemp",HOD.freds,HOD.fredc,HOD.mass_shift,GALAXY_DENSITY);
  fprintf(stdout,"[%s]\n",fname);
  system(fname);
  
  // calculate the luminosity function
  sprintf(fname,"galtemp");
  fp = openfile(fname);
  n = filesize(fp);
  for(i=1;i<=nlf;++i)
    magcnt[i] = 0;
  for(i=1;i<=n;++i)
    {
      for(k=1;k<=10;++k) fscanf(fp,"%f",&x1);// printf("%e\n",x1); }
      k = (x1-maglo)/delta_mag + 1;
      if(k>=1 && k<=nlf)magcnt[k]+=1;
      fscanf(fp,"%f",&x1);
    }
  fclose(fp);

  sprintf(fname,"rm -f galtemp");
  system(fname);

  // calculate the chi^2 for the luminosity function
  chi2lf = 0;
  for(i=1;i<=nlf;++i)
    {
      if(i==nlf)
	x0 = log10(magcnt[i]/pow(BOX_SIZE/0.7,3.0)/delta_mag);
      else
	x0 = log10(magcnt[i]/pow(BOX_SIZE/0.7,3.0)/delta_mag);
      printf("LF %d %f %f %f\n",niter,x0,lfdata[i],lferr[i]);
      chi2lf += (x0 - lfdata[i])*(x0 - lfdata[i])/(lferr[i]*lferr[i]);
    }
  

  printf("CHI %d %f %f %f %e %e %f %e %f\n",niter,HOD.freds,HOD.fredc,HOD.mass_shift,chi2lf,chi2w,MSHIFT,GALAXY_DENSITY,mred/mblue);
  fflush(stdout);
  
  return(chi2lf+chi2w);

}

double mcmc_wtheta_initialize(double *a, double **cov1, double *avg1, double *start_dev)
{
  int i,j=0;
  double x1,x2,omega_m;
  long IDUM = -556;

  a[1] = 0.1; //for WTHETA+LF
  a[2] = 0.8;
  a[3] = 0.55;

  if(ARGC>4)
    {
      a[1] = 0.2; //for WTHETA only
      a[2] = 0.6;
      a[3] = 0.9;
    }

  start_dev[1] = 0.2;
  start_dev[2] = 0.2;
  start_dev[3] = 0.2;

  printf("INITIAL VALUES: ");
  for(i=1;i<=wp.ncf;++i)printf("%e ",a[i]);
  printf("\n");

  for(i=1;i<=wp.ncf;++i)
    {
      avg1[i]=a[i];
      for(j=1;j<=wp.ncf;++j)
	cov1[i][j]=a[i]*a[j];
    }
  x1=chi2_wtheta_wrapper(a);
 
  printf("TRY 0 ");
  for(i=1;i<=wp.ncf;++i)
    printf("%.4e ",a[i]);
  printf("%e\n",x1);fflush(stdout);
  printf("INITIAL CHI2: %e\n",x1);
  fflush(stdout);
  return(x1);
}
