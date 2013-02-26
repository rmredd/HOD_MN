#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

double NGAL;
int USE_COVAR = 1;

double wtheta(double theta);
void nden_check(void);

void output_drg_lf(int iout);
double lf_all(double mag);
double lf_number_density(double mag1, double mag2);
double mag_at_fixed_ngal(double ngal);
double sham_func1(double mag);
double chi2wtheta_covar(float *d, float *x);

double func_drg1(double m);
double func_drg1a(double m);
double func_findmshift(double x);
double func_findmshift(double x)
{
  double m;
  HOD.mass_shift = x;
  m = qromo(func_drg1,log(HOD.M_low),log(HOD.M_max),midpnt)/
    qromo(func_drg1a,log(HOD.M_low),log(HOD.M_max),midpnt)/MALL;
  //printf("B %f %f\n",m,MSHIFT);
  return m - MSHIFT;
}
double func_findmshift2(double x);
double func_findmshift2(double x)
{
  double m;
  HOD.mass_shift = x;
  printf("C %f %e %e\n",x,qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt),NGAL_DRG);
  return qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt) - NGAL_DRG;
}

double func_findfredc(double f);
double func_findfredc(double f)
{
  HOD.fredc = f;
  printf("%e %e\n",qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt),NGAL_DRG);
  return qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt) - NGAL_DRG;
}


void drg_model()
{
  int k,i,n=40,j,i1,i2,ngal;
  double dlogr,r,mass,xk,pnorm,psp,rm,sig,t0,t1,pnorm1,mp,mlo,mhi,dr,
    slo,shi,dsdM,rlo,rhi,dlogm,f1,f2,fac,cvir,rvir,rs,r1,r2,mtot=0,m,
    chi2lf,chi2w,x0,mlow;
  FILE *fp,*fp2;
  float x1,x2,x3,x4;
  char aa[1000],fname[1000],**argv;

  float *xx,*yy,*zz,*vx,*vy,*vz;
  float *wdata,*rdata,*edata,*lfdata,*lferr,*lfmag,*wmodel;
  int magcnt[100];
  int nwdata, nlf, j1, ngrid = 10, ingal;

  double tmin, tmax, dtheta, theta, delta_mag, maglo;

  /* FITTING RIK'S DATA
   */
  rlo = log(1.0);
  rhi = log(1000);
  dtheta = (rhi-rlo)/(19);
  for(i=0;i<20;++i)
    {
      theta = exp(rlo+i*dtheta);
      printf("ALL %e %e\n",theta,wtheta(theta));
      fflush(stdout);
    }

  sprintf(Task.root_filename,"all");
  tasks(ARGC,ARGV);

  RESET_FLAG_1H++;
  RESET_FLAG_2H++;

  GALAXY_DENSITY = 4.9e-4;
  NGAL_DRG = 4.9e-4;
  HOD.pdfc = 10; // set up centrals for RED
  HOD.freds = 0.3; // for RED
  HOD.fredc = 1;
  HOD.mass_shift = zbrent(func_findmshift2,0.0,2.434682,1.0E-4); //spaced in fcen0
  fprintf(stdout,"fsat = %.2f fcen = %.2f mu = %.2f ngal= %e %e\n",HOD.freds,HOD.fredc, HOD.mass_shift,GALAXY_DENSITY,qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt));

  for(i=0;i<20;++i)
    {
      theta = exp(rlo+i*dtheta);
      printf("RED %e %e\n",theta,wtheta(theta));
      fflush(stdout);
    }

  sprintf(Task.root_filename,"red");
  tasks(ARGC,ARGV);

  RESET_FLAG_1H++;
  RESET_FLAG_2H++;

  GALAXY_DENSITY = 1.9e-3;
  HOD.pdfc = 11; // set up centrals for BLUE
  HOD.freds = 1 - HOD.freds; // for BLUE
  printf("ngal = %e\n",qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt));

  for(i=0;i<20;++i)
    {
      theta = exp(rlo+i*dtheta);
      printf("BLUE %e %e\n",theta,wtheta(theta));
      fflush(stdout);
    }

  sprintf(Task.root_filename,"blue");
  tasks(ARGC,ARGV);

  exit(0);



  BOX_SIZE = 160.;

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

  MSTAR = mstar();

  tmin = log(1.0);
  tmax = log(1000);
  dtheta = (tmax-tmin)/(19);

  //HOD.M_min *= 1.5;
  //HOD.M_low *= 1.5;

  // get mean mass of centrals (ALL)
  MALL = qromo(func_drg1,log(HOD.M_low),log(HOD.M_max),midpnt)/
    qromo(func_drg1a,log(HOD.M_low),log(HOD.M_max),midpnt);

  //new model
  HOD.mass_shift = 0.5;
  HOD.pdfc = 10;
  //HOD.shift_alpha = 1;

  mlow = 1;
  mhi = 2;
  dlogm = log(mhi/mlow)/ngrid;

  //analytic_sham();

  /**********************************
   */
  // let's do the prediction for the blue galaxies.
  // best-fit model (alpha=1) fit10 is 12 28
  // best-fit model (alpha=1) fit9 is 30 19

  j1 = 12;
  j = 19;

  // set up HOD as DRGs
  // for stepping in mshift
  ERROR_FLAG = 0;
  MSHIFT = exp((j1-0.5)*dlogm)*mlow;
  HOD.fredc = 1;
  HOD.mass_shift = zbrent(func_findmshift,0.0,10.0,1.0E-4);
  HOD.freds = (j-0.5)/ngrid;
  HOD.fredc = zbrent(func_findfredc,0.0,1.0,1.0E-4);

  // for stepping in fsat0
  HOD.fredc = (j1-0.5)/ngrid;
  HOD.freds = (j-0.5)/ngrid;
  HOD.mass_shift = zbrent(func_findmshift2,0.0,1.434682,1.0E-4); //spaced in fcen0

  // 9panel c2 5.9
  //CHI 993 0.978682 0.660204 0.181224 2.119014e-01 8.793237e+00 1.237535 6.664499e-04
  HOD.mass_shift = 0.181;
  HOD.fredc = 0.66;
  HOD.freds = 0.3;//0.979;


  HOD.freds = 1 - HOD.freds; // for NON-DRGs
  HOD.pdfc = 11; // set up centrals as 1-f(DRG)

  //HOD.pdfc = 7; // square-well blues at low-mass end
  //HOD.M_min = 7e11;
  //HOD.M_cen_max = 1.5e12;
  //HOD.M_low = 7e11;
 
  GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
 fprintf(stdout,"fsat = %.2f fcen = %.2f ngal= %e\n",HOD.freds,HOD.fredc, GALAXY_DENSITY);
  RESET_FLAG_1H++;
  RESET_FLAG_2H++;

  for(i=0;i<20;++i)
    {
      theta = exp(tmin+i*dtheta);
      fprintf(stdout,"WTH %e %e\n",theta,wtheta(theta));
    }
  sprintf(Task.root_filename,"BX");
  //tasks(2,argv);
  exit(0);
  /*
  ****************************************/

  if(ARGC>3)
    ingal = atoi(ARGV[3]);
  else
    ingal = 1;
      
  NGAL_DRG = 6.5e-4;//*(1+(ingal-5)/2.5*0.13);


  // do an outer loop of the mass shifts
  for(j1=1;j1<=ngrid;++j1)
    {
      // go in steps of mean mass shift from 1 to 2.
      MSHIFT = exp((j1-0.5)*dlogm)*mlow;
      HOD.fredc = 1;
      HOD.mass_shift = zbrent(func_findmshift,0.0,10.0,1.0E-4);

      //HOD.mass_shift = 1.116703; //figure 3

      // doing equally spaced in fcen0
      //HOD.fredc = (j1-0.5)/ngrid;

      for(j=1;j<=ngrid;++j)
	{
	  ERROR_FLAG = 0;
	  HOD.freds = (j-0.5)/ngrid;
	  HOD.fredc = zbrent(func_findfredc,0.0,1.0,1.0E-4);
	  if(ERROR_FLAG)
	    HOD.fredc = 1.0;
	  ERROR_FLAG = 0;

	  /*
	  HOD.mass_shift = zbrent(func_findmshift2,0.0,1.434682,1.0E-4); //spaced in fcen0
	  if(ERROR_FLAG)
	    {
	      chi2lf = chi2w = 1.0e6;
	      printf("CHI %d %d %f %f %f %e %e %f\n",j1,j,HOD.freds,HOD.fredc,HOD.mass_shift,chi2lf,chi2w,MSHIFT);
	      fflush(stdout);
	      continue;
	    }
	  */


	  //best fit model from chains (wth+n only)
	  //1.648589e-01 7.732068e-01 9.796977e-01
	  MSHIFT = pow(10.0,0.164);
	  HOD.freds = 0.7732;
	  HOD.fredc = 0.977;
	  HOD.mass_shift = zbrent(func_findmshift,0.0,10.0,1.0E-4);
	  
	  // third column 7.10
	  //CHI 7 10 0.950000 1.000000 0.661607 9.897945e+00 1.109673e+01 1.569168 6.500000e-04 5.905374e-04
	  HOD.mass_shift = 0.6616;
	  HOD.fredc = 1.00;
	  HOD.freds = 0.95;
	  j1 = 7;
	  j = 10;

	  // 9panel c2 5.9
	  //CHI 993 0.978682 0.660204 0.181224 2.119014e-01 8.793237e+00 1.237535 6.664499e-04
	  j1 = 5;
	  j = 9;
	  HOD.mass_shift = 0.181;
	  HOD.fredc = 0.66;
	  HOD.freds = 0.979;

	  //for 9panel c1 1.1
	  j1 = 1;
	  j = 1;
	  HOD.mass_shift = 0;
	  HOD.fredc = NGAL_DRG/GALAXY_DENSITY;
	  HOD.freds = NGAL_DRG/GALAXY_DENSITY;

	  //best-fit model without LF
	  HOD.mass_shift = 0.48;
	  HOD.fredc = 0.99;
	  HOD.freds = 0.69;

	  //best-fit model with LF (bottom row of 6panel)
	  //8.234815e-02 9.011035e-01 6.467542e-01 
	  j1 = 7;
	  j = 10;
	  HOD.mass_shift = 1.505979e-01 ;
	  HOD.fredc = 6.467542e-01 ;
	  HOD.freds = 9.011035e-01 ;

	  GALAXY_DENSITY = qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
	  fprintf(stdout,"fsat = %.1f fcen = %f ngal= %e\n",HOD.freds,HOD.fredc, GALAXY_DENSITY);
	  RESET_FLAG_1H++;
	  RESET_FLAG_2H++;

	  /*
	  for(i=0;i<20;++i)
	    {
	      theta = exp(tmin+i*dtheta);
	      fprintf(stdout,"WTH %e %e\n",theta,wtheta(theta));
	    }
	  exit(0);
	  */

	  //populate_simulation();
	  //exit(0);

	  // get qudari model points
	  //fp = openfile("q8m.dat");
	  

	  // calculate the chi^2 for the wtheta values.
	  chi2w = 0;
	  for(i=1;i<=nwdata;++i)
	    {
	      x0 = wtheta(rdata[i]);
	      wmodel[i] = x0;
	      //q8 model
	      //fscanf(fp,"%f %f",&x1,&wmodel[i]);
	      //x0 = wmodel[i];
	      printf("XX %e %e %e %e\n",rdata[i],wdata[i],x0,edata[i]);
	      chi2w += (wdata[i]-x0)*(wdata[i]-x0)/(edata[i]*edata[i]);
	    }
	  //fclose(fp);

	  fmuh(chi2w);
	  if(USE_COVAR)
	    chi2w = chi2wtheta_covar(wdata,wmodel);
	  fmuh(chi2w);
	  //exit(0);
	  
	  sprintf(fname,"wth_mshift_%d.%d",j1,j);
	  fp = fopen(fname,"w");
	  for(i=0;i<20;++i)
	    {
	      theta = exp(tmin+i*dtheta);
	      fprintf(fp,"%e %e\n",theta,wtheta(theta));
	    }
	  fclose(fp);
	  //continue;
	  
	  // do the real-space clustering and HOD
	  sprintf(Task.root_filename,"mshift_%d.%d",j1,j);
	  tasks(2,argv);
	  
	  // now make the luminosity function for this model
	  //output_drg_lf(j);      
	  sprintf(fname,"sham ../../SHAM/halosub_0.284 hod_mshift_%d.%d %f %f %f %f > gal_mshift_%d.%d",j1,j,HOD.freds,HOD.fredc,HOD.mass_shift,GALAXY_DENSITY,j1,j);
	  //sprintf(fname,"sham ../SHAM_120/halosub_0.3323.dat hod_mshift_%.1f.%d %f %f %f > gal_mshift_%.1f.%d",HOD.mass_shift,j,HOD.freds,HOD.fredc,HOD.mass_shift,HOD.mass_shift,j);
	  fprintf(stderr,"[%s]\n",fname);
	  system(fname);
	  
	  // calculate the clustering of this drg sample
	  sprintf(fname,"covar3 0.1 15 12 160 0 160 1 gal_mshift_%d.%d a 0 1 auto > xi.mshift_%d.%d",j1,j,j1,j);
	  //fprintf(stderr,"[%s]\n",fname);
	  //system(fname);
	  

	  // calculate the luminosity function
	  sprintf(fname,"gal_mshift_%d.%d",j1,j);
	  fp = openfile(fname);
	  n = filesize(fp);
	  for(i=1;i<=nlf;++i)
	    magcnt[i] = 0;
	  for(i=1;i<=n;++i)
	    {
	      for(k=1;k<=10;++k) fscanf(fp,"%f",&x1);// printf("%e\n",x1); }
	      k = (x1-maglo)/delta_mag + 1;
	      if(k>=1 && k<=nlf)magcnt[k]+=1;
	      //printf("%d %d %f %f %f\n",i,k,x1,maglo,delta_mag);
	      fscanf(fp,"%f",&x1);
	    }
	  fclose(fp);
	  
	  // calculate the chi^2 for the luminosity function
	  chi2lf = 0;
	  for(i=1;i<=nlf;++i)
	    {
	      if(i==nlf)
		x0 = log10(magcnt[i]/pow(BOX_SIZE/0.7,3.0)/delta_mag);
	      else
		x0 = log10(magcnt[i]/pow(BOX_SIZE/0.7,3.0)/delta_mag);
	      printf("LF %d %d %f %f %f\n",j1,j,x0,lfdata[i],lferr[i]);
	      chi2lf += (x0 - lfdata[i])*(x0 - lfdata[i])/(lferr[i]*lferr[i]);
	    }
	  
	  printf("CHI %d %d %f %f %f %e %e %f %e %e\n",j1,j,HOD.freds,HOD.fredc,HOD.mass_shift,chi2lf,chi2w,MSHIFT,NGAL_DRG,GALAXY_DENSITY);
	  fflush(stdout);
	  exit(0);
	}
    }

  exit(0);
}

void analytic_sham()
{
  int i,j,k,nlogm=400,imag,nmag,jmin;
  double dlogm,halocnt[401],magcnt[200],subcnt[401],m,msub,n1,mlow,dmag,magmin,nsub;

  double *hmass, *mag, *yy;

  dmag = 0.2;
  magmin = -18.0;
  nmag = 5/dmag;
  
  for(i=0;i<200;++i)
    magcnt[i] = 0;
  hmass = dvector(1,nlogm);
  mag = dvector(1,nlogm);
  yy = dvector(1,nlogm);
  
  mlow = HOD.M_low;

  for(i=1;i<=nlogm;++i)
    halocnt[i] = subcnt[i] = 0;
  
  // go through the mass function to find TOTAL number of halos
  dlogm = log(HOD.M_max/mlow)/nlogm;
  
  // get all HOD
  HOD.mass_shift = 0;
  HOD.fredc = 0.44;

  HOD.mass_shift = 0.6616;
  HOD.fredc = 1.00;

  HOD.mass_shift = 0.2;
  HOD.fredc = 0.66;

  for(i=1;i<=nlogm;++i)
    {
      m = exp(dlogm*(i-0.5))*mlow;
      halocnt[i] += dndM_interp(m)*m*dlogm;
      
      n1 = 0;
      for(j=1;j<=nlogm;++j)
	{
	  msub = exp(dlogm*(j-0.5))*mlow;
	  if(msub>m/3)continue; // no subhalo greater than half total mass
	  //n1 = pow(m,0.5)/2.6*pow(msub,-1.5)*exp(-pow(msub/m/0.7,6))*dlogm*m/2; //original +2
	  //n1 = pow(m,0.7)/2.6*pow(msub,-1.7)*exp(-pow(msub/m/0.7,6))*dlogm*m/3.16;
	  n1 = pow(msub,-2)*m*m*dlogm*15.7*HOD.M1;
	  halocnt[0] += n1;
	  subcnt[0] += N_cen(msub)*n1;	  
	}
      //printf("%e %e %e\n",m,N_sat(m),n1);
    }
  printf("%f\n",subcnt[0]/halocnt[0]);
  exit(0);

}

void output_drg_lf(int iout)
{
  FILE *fp;
  char fname[100];

  int i,j,k,nlogm=400,imag,nmag,jmin;
  double dlogm,halocnt[401],magcnt[200],m,msub,n1,mlow,dmag,magmin,nsub;

  double *hmass, *mag, *yy;

  dmag = 0.2;
  magmin = -18.0;
  nmag = 5/dmag;
  
  for(i=0;i<200;++i)
    magcnt[i] = 0;
  hmass = dvector(1,nlogm);
  mag = dvector(1,nlogm);
  yy = dvector(1,nlogm);
  
  mlow = HOD.M_low/5.0;

  for(i=1;i<=nlogm;++i)
    halocnt[i] = 0;
  
  // go through the mass function to find TOTAL number of halos
  dlogm = log(HOD.M_max/mlow)/nlogm;
  
  for(i=1;i<=nlogm;++i)
    {
      m = exp(dlogm*(i-0.5))*mlow;
      halocnt[i] += dndM_interp(m)*m*dlogm;
      
      for(j=1;j<=nlogm;++j)
	{
	  msub = exp(dlogm*(j-0.5))*mlow;
	  if(msub>m/2)continue; // no subhalo greater than half total mass
	  
	  // charlie's fit (with extra factor of log(10)?)
	  halocnt[j] += 0.02*log(10)*pow(msub/m,-0.9)*exp(-msub/(2*m))*dlogm*dndM_interp(m)*m*dlogm;
	  
	  // gao et al 2004-- gives ~4% at M=1e12 z=0.
	  //halocnt[i] += 6.3e-4*pow(msub,-1.9)*m*msub*dlogm*dndM_interp(m)*m*dlogm;
	}	  
    }
  
  // make it cumulative
  // and associate a magnitude at every mass.
  mag[nlogm] = mag_at_fixed_ngal(halocnt[nlogm]);
  hmass[nlogm] = exp(dlogm*(nlogm-0.5)*mlow);

  for(i=nlogm-1;i>=1;--i)
    {
      halocnt[i] += halocnt[i+1];
      mag[i] = mag_at_fixed_ngal(halocnt[i]);
      hmass[i] = exp(dlogm*(i-0.5))*mlow;
      //printf("CNT %e %e %e\n",hmass[i],halocnt[i],lf_number_density(-27.9,-2.0));
    }
  
  // prepare for spline interpolation
  spline(hmass,mag,nlogm,1.0E+30,1.0E+30,yy);
  
  // now go through HOD and create the DRG luminosity function
  for(i=1;i<=nlogm-1;++i)
    {
      m = hmass[i];
      if(m<HOD.M_low)continue;
      imag = (magmin-mag[i])/dmag+1;
      //printf("%e\n",magcnt[i]);
      magcnt[imag] += N_cen(m)*dndM_interp(m)*dlogm*m;

      //printf("%e %f %d %e %e %e %e %e\n",m,mag[i],imag,magcnt[imag],N_cen(m),m,dlogm,dndM_interp(m));
      //continue;
      
      // find the minimum subhalo mass
      nsub = 0;
      for(j=nlogm;j>=1;--j)
	{
	  msub = hmass[j];
	  if(msub>m/2)continue;
	  nsub += 0.02*log(10)*pow(msub/m,-0.9)*exp(msub/(2*m))*dlogm;
	  if(nsub>N_sat(m)/HOD.freds) break;
	}
      if(j==nlogm)continue;
      jmin = j;
      if(jmin>=0)
	{
	  printf("ERR %d %e %e %e %e %f\n",j,m,msub,nsub,N_sat(m)/HOD.freds,mag[jmin]);
	}

      // go through the satellites
      for(j=jmin;j<=nlogm;++j)
	{
	  msub = hmass[j];
	  if(msub>m/2)break;
	  imag = (magmin-mag[j])/dmag+1;
	  magcnt[imag] += HOD.freds*0.02*log(10)*pow(msub/m,-0.9)*exp(msub/(2*m))*dlogm*dndM_interp(m)*m*dlogm;
	}
    }

  // print out the results
  sprintf(fname,"drg_lf.%d",iout);
  fp = fopen(fname,"w");

  n1 = 0;
  for(i=0;i<200;++i)
    {
      if(-(i-0.5)*dmag+magmin >-21.6)continue;
      magcnt[i]/=dmag;
      if(magcnt[i]>0)
	fprintf(fp,"%f %e\n",-(i-0.5)*dmag + magmin,magcnt[i]);
      n1 += magcnt[i]*dmag;
    }
  printf("TOTAL ngal = %e\n",n1);
  fclose(fp);
}

double mag_at_fixed_ngal(double ngal)
{
  NGAL = ngal;
  return zbrent(sham_func1,-27.9,-2.0,1.0E-5);
}
double sham_func1(double mag)
{
  //printf("%e %e %f\n",NGAL,lf_number_density(-28.0,mag),mag);
  return NGAL  - lf_number_density(-28.00,mag);
}
double lf_number_density(double mag1, double mag2)
{
  return qromo(lf_all,mag1,mag2,midpnt);
}
double lf_all(double mag)
{
  return 0.4*2.30258*exp(-pow(10.0,-0.4*(mag+22.67)))*(10.65e-4*pow(10.0,-0.4*(mag+22.67)*(-1.01+1)));
}

// this is to calcualte the mean number of DRG central galaxies
double func_drg1a(double m)
{
  m = exp(m);
  return N_cen(m)*m*dndM_interp(m);
}
// this is to calculate the numer-weighted mass of DRG galaxies
double func_drg1(double m)
{
  m = exp(m);
  return N_cen(m)*m*dndM_interp(m)*m;
}

double chi2wtheta_covar(float *d, float *x)
{
  static int flag = 1, n=12;
  static float **cij;
  double **tmp, **tmp2, chi2=0, x1;
  int i,j;
  FILE *fp;

  if(flag)
    {
      flag = 0;
      fp = openfile("covar.dat");
      cij = matrix(1,n,1,n);
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  fscanf(fp,"%e",&cij[i][j]);
      fclose(fp);

      tmp=dmatrix(1,n,1,1);
      tmp2=dmatrix(1,n,1,n);
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  tmp2[i][j]=cij[i][j];
      gaussj(tmp2,n,tmp,1);
      for(i=1;i<=n;++i)
	for(j=1;j<=n;++j)
	  cij[i][j]=tmp2[i][j];
      free_dmatrix(tmp,1,n,1,1);
      free_dmatrix(tmp2,1,n,1,n);
    }

  for(i=1;i<=n;++i)
    for(j=1;j<=n;++j)
      {
	x1=(x[i]-d[i])*(x[j]-d[j])*cij[i][j];
	chi2+=x1;
      }
  return chi2;
}
