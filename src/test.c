#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

void interp_hod_params(double *a, double *b);
void grid_check(void);
void grid_check2(void);
double tfunc1(double);
void xi_void_model(int argc, char **argv);
void mass2number(void);
int mcmc_chain_post_processing(int ARGC, char **ARGV);

double delta_pdf(double delta, double r);

double mass_gt,delta_gt,g20_ncen;

double func_fsat(double m)
{
  m=exp(m);
  return(dndM_interp(m)*N_cen(m)*m);
}
double func_BW(double m)
{
  double n;
  m=exp(m);
  n=N_sat(m)-1;
  if(n<0)n=0;
  return(dndM_interp(m)*n*m);
}
double func_ngh(double m)
{
  double n;
  m=exp(m);
  n=N_sat(m)+N_cen(m);
  if(n>1)n=1;
  return(dndM_interp(m)*n*m);
}
double func_BW1(double m)
{
  double n,func_BW2();
  HOD.M1=exp(m);
  n=qromo(func_BW2,log(HOD.M_min),log(HOD.M_max),midpnt);
  return(n-GALAXY_DENSITY);
}
double func_BW2(double m)
{
  double n;
  m=exp(m);
  n=N_sat(m);
  return(dndM_interp(m)*n*m);
}
double func_BW3(double m)
{
  double n;
  m=exp(m);
  n=N_cen(m);
  return(dndM_interp(m)*n*m);
}
double func_maxcen(double x)
{
  double n;
  HOD.MaxCen=x;
  n=qromo(func_BW3,log(HOD.M_min*1.0E-3),log(HOD.M_min*1.0E+3),midpnt);
  return(n-g20_ncen);
}
double fixfunc1(double m);

double func_total_bias(double m)
{
  m=exp(m);
  return(halo_mass_function(m)*bias(m)*m*m/(RHO_CRIT*OMEGA_M));
}
double func_total_mass(double m)
{
  m=exp(m);
  return(halo_mass_function(m)*m*m/(RHO_CRIT*OMEGA_M));
}
double func_find_m1(double m)
{
  double x;
  m=exp(m);
  HOD.M1 = m;
  x = qromo(func_galaxy_density,log(HOD.M_min),log(HOD.M_max),midpnt);
  /* printf("> %e %e %e\n",m,x,GALAXY_DENSITY); */
  return(x-GALAXY_DENSITY);
}

double wtheta(double theta);
void nden_check(void);
void drg_model(void);

double func_findfredc(double f);

void test(int argc, char **argv)
{
  int k,i,n=40,j,i1,i2,ngal;
  double dlogr,r,mass,xk,pnorm,psp,rm,sig,t0,t1,pnorm1,mp,mlo,mhi,dr,
    slo,shi,dsdM,rlo,rhi,dlogm,f1,f2,fac,cvir,rvir,rs,r1,r2,mtot=0;
  FILE *fp,*fp2;
  float x1,x2,x3,x4;
  char aa[1000],fname[1000];

  float *xx,*yy,*zz,*vx,*vy,*vz;

  double ml_minimization(),ml_powell();

  double func_fsat(),func_BW(),func_BW1(),func_BW2(),func_maxcen(),func_BW3(),
    func_ngh(),func_total_bias(),func_total_mass(),func_find_m1(),func_halo_density();

  double tmin, tmax, dtheta, theta;

  mcmc_chain_post_processing(argc, argv);
  exit(0);

  m2n_mcmc();
  exit(0);
  drg_model();

  tmin = log(1.0);
  tmax = log(1000);
  dtheta = (tmax-tmin)/(19);
  for(i=0;i<20;++i)
    {
      theta = exp(tmin+i*dtheta);
      printf("ALL %e %e\n",theta,wtheta(theta));
      fflush(stdout);
    }

  RESET_FLAG_1H++;
  RESET_FLAG_2H++;

  HOD.freds = 0.5; // for NON-DRGs
  HOD.pdfc = 10; // set up centrals as 1-f(DRG)
  HOD.mass_shift = 0.3;

  for(i=0;i<20;++i)
    {
      theta = exp(tmin+i*dtheta);
      printf("RED %e %e\n",theta,wtheta(theta));
      fflush(stdout);
    }


  exit(0);
  
  mcmc_wtheta();


  mass = 1.41*1000*2.775e11*pow(500./2140,3.0)*0.25;
  mass = 2.86e11;
  mass = 2.10e14;
  fmuh(mass);
  fmuh(sigmac_interp(mass));
  fmuh(bias(mass));
  fmuh(bias_interp(mass,-1));
  for(i=0;i<=180;++i)
    {
      r = pow(10.0,i/100.0);
      printf("%e %e\n",r,bias_interp(mass,r)*bias_interp(mass,r)*xi_interp(r));
    }
  exit(0);


  /* bias echeck
   */
  for(i=0;i<=12;++i)
    {
      mass = pow(2.0,i+0.5)*2.775e11*pow(500./2140,3.0)*0.25*82;
      sig = sigmac_interp(mass);
      printf("%e %f %f %f %f\n",mass,bias(mass),sig,log10(bias(mass)),log10(1/sig));
    }
  exit(0);



  drg_model();



  //printf("%f %f\n",log10(1/sigmac_interp(HOD.M_min)),log10(bias_interp(HOD.M_min,-1)));
  //printf("%f %f\n",log10(1/sigmac_interp(4.5e12)),log10(bias_interp(4.5e12,-1)));
  //exit(0);


  //nden_check();


  REDSHIFT = 0.1;
  SIGMA_8 = 0.8*growthfactor(0.1);
  RESET_COSMOLOGY++;
  fmuh(halo_mass_function(1.0e14));

  REDSHIFT = 0.3;
  SIGMA_8 = 0.8*growthfactor(0.3);
  RESET_COSMOLOGY++;
  fmuh(halo_mass_function(1.0e14));
  exit(0);

  output_halo_mass_function();
  exit(0);

  for(i=0;i<=13;++i)
    {
      mass = pow(RESOLUTION,3.0)*100*pow(2.0,i+0.5)*RHO_CRIT*OMEGA_M;
      sprintf(Task.root_filename,"mr%d",i);
      output_halo_correlation_function(mass);
    }
  exit(0);

  output_halo_correlation_function(5.7e11);
  exit(0);    

  /* this is for fitting the DR7 data:
   * iterate through several values of sigma8
   */
  for(i=1;i<=9;++i)
    {
      RESET_COSMOLOGY++;
      SIGMA_8 = 0.6 + (i-1.)*0.4/8.0;
      printf("SIGMA8= %f\n",SIGMA_8);
      sprintf(Task.root_filename,"s8_%.2f",SIGMA_8);
      wp_minimization(argv[1]);
      mass2number();

      // output the real-space correlation function
      sprintf(fname,"%s.r_space",Task.root_filename);
      fp=fopen(fname,"w");
      dr=(log(70.0)-log(0.01))/49.0;
      for(j=0;j<50;++j)
	{
	  r=exp(j*dr+log(0.01));
	  x1=one_halo_real_space(r);
	  x2=two_halo_real_space(r);
	  x3=projected_xi(r);
	  fprintf(fp,"%f %e %e %e %e\n",r,x1,x2,x1+x2,x3);
	  fflush(fp);
	}
      fclose(fp);

    }
  exit(0);


  for(i=1;i<=2000;++i)
    {
      r = (double)i/100;
      printf("BIAS %f %f\n",r,bias_interp(HOD.M_min,r)/bias_interp(HOD.M_min,-1));
    }
  exit(0);


  /* Make a grid of gal files.
   */
  slo = 0.6;
  shi = 1.4;
  
  //slo = 0;
  //shi = 0.5;

  set_HOD_params();
  ngal = 2*GALAXY_DENSITY*BOX_SIZE*BOX_SIZE*BOX_SIZE;
  fprintf(stderr,"%d %f %f\n",ngal,GALAXY_DENSITY,BOX_SIZE);

  xx = malloc(ngal*sizeof(float));
  yy = malloc(ngal*sizeof(float));
  zz = malloc(ngal*sizeof(float));
  vx = malloc(ngal*sizeof(float));
  vy = malloc(ngal*sizeof(float));
  vz = malloc(ngal*sizeof(float));

  for(i=1;i<=1;++i)
    {
      //VBIAS = (shi-slo)/20.0*(i-0.5)+slo;
      //VBIAS_C = (shi-slo)/20.0*(i-0.5)+slo;
      for(j=1;j<=10;++j)
	{
	  IDUM_MCMC = 555 + j;
	  sprintf(Task.root_filename,"test_%d.%d",i,j);
	  printf("GRID %s.mock %d %d %f %f\n",Task.root_filename,i,j,VBIAS,VBIAS_C);
	  internal_populate_simulation(xx,yy,zz,1.0,1,vx,vy,vz);
	}
    }
  exit(0);

  for(i=90;i<=160;++i)
    {
      mass = pow(10.0,i/10.0);
      mtot += mass*mass*dndM_interp(mass)/RHO_CRIT/OMEGA_M*0.1*log(10);
      printf("%e %e %e\n",mass, mass*dndM_interp(mass)/RHO_CRIT/OMEGA_M,mtot);
    }
  exit(0);

  /*
  for(i=-400;i<=400;++i)
    {
      xk = pow(10.0,i/100.0);
      x1 = transfnc(xk);
      printf("%e %e %e %e %e\n",xk,x1,x1,x1,x1);
    }
  exit(0);

  output_halo_mass_function();
  exit(0);
  
  for(i=-10;i<=14;++i)
    printf("%f %f\n",pow(10.0,i/10.0),nbody_xi(pow(10.0,i/10.0)));
  exit(0);

  covar_test();

  XCORR=1;
  HOD2.M_min = 1.0e12;
  HOD2.M_cen_max = 4.25e12;
  HOD2.M_max = 1.0E16;
  HOD2.M1 = pow(10.0,13.38);
  HOD2.alpha = 1.16;
  HOD2.pdfs = 1;
  HOD2.pdfc = 7;

  set_HOD2_params();
  one_halo_real_space(1.0);
  exit(0);
  */

  /*
  fit_dispersion_model();
  exit(0);
  pairwise_velocity_dispersion();
  exit(0);

  printf("%e %e\n",N_sat(1.0e13),one_halo(0.2,0.2));
  exit(0);

  for(i=100;i<=155;++i)
    {
      mass = pow(10.0,i/10.0);
      cvir = halo_concentration(mass);
      r1 = pow(3*mass/(4*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
      r2 = pow((3*halo_mass_conversion2(mass,cvir,180.0,334.0)/
		(4*PI*334.0*RHO_CRIT*OMEGA_M)),1.0/3.0);
      printf("%e %e\n",r1,r2);
    }
  exit(0);
  */

  

  /* Test cobenorm
   */
  //  cobe_prior(argv[3]);

  /* TESTING TESTING TESTING
   */  
  if(argc>3)
    IDUM_MCMC=atoi(argv[3]);
  if(IDUM_MCMC==101)goto OUTPUT_PSPEC;
  if(IDUM_MCMC>0)
    ml_powell();
  ml_minimization();
  exit(0);

 OUTPUT_PSPEC:
  for(k=-400;k<=400;++k)
    {
      xk=pow(10.0,k/100.0);
      printf("%e %e %e 0.0 0.0\n",xk,transfnc(xk),transfnc(xk));
      /* printf("%f %e %e\n",xk,linear_power_spectrum(xk),nonlinear_power_spectrum(xk)); */
    }
  exit(0);



  /*
  for(i=-30;i<=20;++i)
    {
      xk = pow(10.0,i/10.0);
      printf("%f %e\n",log10(xk),log10(linear_power_spectrum(xk)/(xk*xk*xk)*TWOPI*TWOPI));
    }
  exit(0);
  */

  /*
  vpf();

  printf("%e\n",qromo(func_halo_density,log(3.9e11),log(HOD.M_max),midpnt));

  for(i=1;i<=1;++i)
    {
      HOD.M_cut = HOD.M_min = 4.1e11*pow(1.02,i-1);
      GALAXY_DENSITY = 0.01;
      HOD.alpha = 0.75;
      mass = zbrent(func_find_m1,log(HOD.M_min),log(HOD.M_max*100),1.0E-5);
      printf("%e %e\n",HOD.M_min,exp(mass));
    }
  exit(0);
  */

  /*
  vpf();

  for(i=150;i>=100;--i)
    {
      mass = pow(10.0,i/10.0);
      x2 = qromo(func_total_mass,log(mass),log(1.0e17),midpnt);
      printf("%f %e\n",i/10.0,x2);
    }      
  exit(0);

  fprintf(stderr,"%e %e\n",halo_mass_function(MSTAR), bias(MSTAR));

  mass = 1.0e8;
  x1 = qromo(func_total_bias,log(mass),log(1.0e17),midpnt);
  x2 = qromo(func_total_mass,log(mass),log(1.0e17),midpnt);

  fprintf(stderr,"%f\n",(1-x1)/(1-x2));

  printf("%e\n",qromo(func_total_bias,log(mass),log(1.0e17),midpnt));
  printf("%e\n",qromo(func_total_mass,log(mass),log(1.0e17),midpnt));
  exit(0);

  mass = 1.0e14;
  cvir = halo_concentration(mass);
  rvir = pow(4*mass/(3*PI*OMEGA_M*RHO_CRIT*DELTA_HALO),1.0/3.0);
  rs = rvir/cvir;
  fac = sqrt(4.499E-48/2.0)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sig = fac*pow(mass,0.33333333333);
  */


  for(i=50;i<=150;++i)
    {
      CVIR_FAC = i/100.0;
      printf("%f %f %f\n",sig,rvir,velocity_dispersion_profile(mass,cvir*CVIR_FAC,rvir,rvir*0.5));
    }
  exit(0);
  fp=openfile(argv[3]);
  n=filesize(fp);
  f1 = 50.0/(n-200);
  j=0;
  for(i=1;i<=n;++i)
    {
      if(!(drand48()<f1 && i>200))continue;
      j++;
      sprintf(aa,"bat_tail.%d",j);
      fp2=fopen(aa,"w");
      fscanf(fp,"%6s %d %d %e %e %e %e",aa,&i1,&i2,&x1,&x2,&x3,&x4);
      HOD.M1 = pow(10.0,x1);
      HOD.M_cut = pow(10.0,x2);
      HOD.sigma_logM = pow(10.0,x3);
      set_HOD_params();
      fprintf(fp2,"M_min %e\n",HOD.M_min);
      fprintf(fp2,"M1 %e\n",HOD.M1);
      fprintf(fp2,"M_cut %e\n",HOD.M_cut);
      fprintf(fp2,"sigma_logM %e\n",HOD.sigma_logM);
      fclose(fp2);
      printf("%d %d %e\n",i,j,x4);
    }
  exit(0);

  f1=qromo(func_ngh,log(HOD.M_low),log(HOD.M_max),midpnt);
  printf("%e %f\n",f1,f1/GALAXY_DENSITY);
  exit(0);

  set_HOD_params();
  f1=qromo(func_fsat,log(HOD.M_low),log(HOD.M_max),midpnt);
  printf("%e %f %f\n",HOD.M_min,f1/GALAXY_DENSITY,1-f1/GALAXY_DENSITY);
  exit(0);  

  /* Output some velocity dispersions for halos of different
   * masses.
   */
  fac=sqrt(4.499E-48/2.0)*
    pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  mass = pow(10.0,12.7);
  sig=fac*pow(mass,0.33333333333);
  printf("%e\n",sig);

  mass = pow(10.0,13.15);
  sig=fac*pow(mass,0.33333333333);
  printf("%e\n",sig);
  exit(0);


  /* THis was to determine M1 for the mag-bin models
   * with varying cen width (but same MaxCen).
   */
  HOD.sigma_logM=0.22;
  HOD.M_min=5.0E+12;
  HOD.M1=HOD.M_min*20;
  HOD.alpha=1;
  HOD.pdfs=1;
  HOD.pdfc=6;
  HOD.MaxCen=1.0;
  OUTPUT=0;
  set_HOD_params();
  g20_ncen=qromo(func_BW3,log(HOD.M_min*1.0E-3),log(HOD.M_min*1.0E+3),midpnt);
  printf("%e %e\n",g20_ncen,g20_ncen/GALAXY_DENSITY);
  HOD.sigma_logM=0.09;
  for(i=2;i<=25;++i)
    {
      HOD.M1=-1;
      HOD.sigma_logM+=0.01;
      set_HOD_params();
      f1=qromo(func_BW3,log(HOD.M_min*1.0E-3),log(HOD.M_min*1.0E+3),midpnt);
      printf("%e %f %e %f\n",HOD.M_min,HOD.sigma_logM,HOD.M1,f1/GALAXY_DENSITY);
    }
  exit(0);


  /* THis was to determine MaxCen for the mag-bin models
   * with varying cen width. (Same N_sat)
   */
  HOD.sigma_logM=0.2;
  HOD.M_min=5.0E+12;
  HOD.M1=HOD.M_min*23;
  HOD.alpha=1;
  HOD.pdfs=1;
  HOD.pdfc=6;
  HOD.MaxCen=1.0;
  OUTPUT=1;
  set_HOD_params();
  g20_ncen=qromo(func_BW3,log(HOD.M_min*1.0E-3),log(HOD.M_min*1.0E+3),midpnt);
  printf("%e %e\n",g20_ncen,g20_ncen/GALAXY_DENSITY);
  for(i=2;i<=5;++i)
    {
      HOD.sigma_logM=i/10.0;
      HOD.MaxCen=zbrent(func_maxcen,1.0E-5,1.0,1.0E-5);
      printf("%e %f %f\n",HOD.M_min,i/10.0,HOD.MaxCen);
    }
  exit(0);

  set_HOD_params();
  f1=qromo(func_fsat,log(HOD.M_low),log(HOD.M_max),midpnt);
  printf("%e %f\n",HOD.M_min,1-f1/GALAXY_DENSITY);
  exit(0);  
  
  HOD.M_min=2.8e11;
  for(i=0;i<=15;++i)
    {
      HOD.M_min+=0.1e11;
      HOD.M1=0;
      set_HOD_params();
      f1=qromo(func_fsat,log(HOD.M_low),log(HOD.M_max),midpnt);
      printf("%e %f\n",HOD.M_min,1-f1/GALAXY_DENSITY);
    }
  exit(0);


  HOD.pdfs=1;
  HOD.pdfc=1;
  f1=GALAXY_DENSITY=0.01;
  HOD.M_low=HOD.M_min = 0.7e11;
  HOD.alpha=0.5;
  for(i=0;i<=3;++i)
    {
      HOD.M_low=HOD.M_min=2*HOD.M_min;
      HOD.M1 = exp(zbrent(func_BW1,log(1.0E9),log(1.0e15),1.0E-5));
      f1 = qromo(func_BW,log(HOD.M_min),log(HOD.M_max),midpnt);
      printf("%e %e %e %e\n",f1,f1/GALAXY_DENSITY,HOD.M_min,HOD.M1);
    }
  exit(0);

  xi_void_model(argc,argv);
  /* grid_check2(); */

  /* This is for the void paper. Do the f(>d) calculation analyically
   */
  n=10;
  mlo = log(3.0E+11);
  mhi = log(1.0E+13);
  dlogm = (mhi-mlo)/(n-1);
  delta_gt = 0.5;
  for(i=1;i<=n;++i)
    {
      mass_gt = exp((i-1)*dlogm+mlo);
      f1 = qromo(tfunc1,log(1.0e-5),log(1.0e+5),midpnt);
      f2 = qromo(tfunc1,log(delta_gt),log(1.0E+5),midpnt);
      printf("%e %f\n",mass_gt,1-f2/f1);
    }
  exit(0);

  fp=openfile(argv[3]);
  n=filesize(fp);
  BEST_FIT=1;
  RESET_COSMOLOGY++;
  pnorm = SIGMA_8/sigmac(8.0);
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%f %f %f %f",&x1,&x2,&x3,&x4);
      mass=x1;
      r=pow(3.0*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
      sig = pnorm*sigmac(r);

      mlo=0.99*mass;
      mhi=1.01*mass;
      rlo=pow(3.0*mlo/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
      rhi=pow(3.0*mhi/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
      slo=pnorm*sigmac(rlo);
      shi=pnorm*sigmac(rhi);
      dsdM=-(shi-slo)/(mhi-mlo)*OMEGA_M*RHO_CRIT/mass/sig;
      printf("%e %e %e %e\n",1/sig,x2/dsdM,x3/dsdM,dndM_interp(mass)/dsdM);
    }
  exit(0);

  if(argc<8)
    {
      fprintf(stderr,"not enough parameters\n");
      exit(0);
    }
  fprintf(stderr,"%s\n",argv[3]);
  fp=fopen(argv[3],"r");
  n=filesize(fp);

  RESET_COSMOLOGY++;
  SIGMA_8 = atof(argv[4]);
  GAMMA = atof(argv[5]);
  ITRANS = atoi(argv[6]);
  sprintf(Files.TF_file,"CMBFAST_trans.dat");
  j = atoi(argv[7]);

  pnorm=SIGMA_8/sigmac(8.0);
  switch(j){
  case 400:
    mp=2.78e11*.3*pow(400/1280.0,3.0); break;
  case 768:
    mp=2.78e11*.3*pow(768/1024.0,3.0); break;
  case 384:
    mp=2.78e11*.3*pow(384/1024.0,3.0); break;
  case 253:
    mp=2.78e11*.3*pow(253/360.,3.0); break;
  }
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%d %e %e %e",&k,&x1,&x2,&x3);
      x1 = x1 - 0.75*sqrt(x1);
      x1*=mp;
      rm = pow(3*x1/(4*PI*OMEGA_M*RHO_CRIT),THIRD);
      printf("%e %f %f\n",x1,1/(pnorm*sigmac(rm)),rm);
    }      
  exit(0);




  OUTPUT=1;

  printf("%e\n",xi2d_interp(0.1, 0.1, 1.1, 1.1));

  printf("%e\n",integrated_bin(0.1,0.1,1.0,1.0,2));
  printf("%e\n",integrated_bin(0.1,0.1,1.0,1.0,4));
  printf("%e\n",integrated_bin(0.1,0.1,1.0,1.0,10));
  printf("%e\n",integrated_bin(0.1,0.1,1.0,1.0,16));
  printf("%e\n",integrated_bin(0.1,0.1,1.0,1.0,32));
  exit(0);

  fprintf(stderr,"%e\n",sigmac(1.0)/sigmac(100.0));
  exit(0);

  mp = RHO_CRIT*OMEGA_M*pow(400/1280.,3.)*30;
  pnorm = SIGMA_8/sigmac(8.0);
  for(i=0;i<=15;++i)
    {
      mass = pow(2.0,i+0.5)*mp;
      rm = pow(3*mass/(4*PI*OMEGA_M*RHO_CRIT),THIRD);
      printf("%e %f %f\n",mass,pnorm*sigmac(rm),rm);
    }
  exit(0);

  two_halo(1.,1.);
  t0=second();
  for(i=1;i<=100;++i)
    printf("%f %e\n",i/10.0,two_halo(2.,i/10.0));
  t1=second();
  fprintf(stdout,"%.2f\n",timediff(t0,t1));
  exit(0);

  pnorm=SIGMA_8/sigmac(8.0);
  pnorm1=SIGMA_8/nonlinear_sigmac(8.0);
  pnorm1 = pnorm*sigmac(80.0)/nonlinear_sigmac(80.0);
  dlogr = (log(40)-log(0.1))/(n-1);
  for(k=1;k<=n;++k)
    {
      r=exp(dlogr*(k-1))*0.1;
      printf("boo %f %f %f %f\n",r,pnorm*sigmac(r),pnorm1*nonlinear_sigmac(r),pnorm/pnorm1);
    }
  exit(0);


  for(k=-200;k<=200;++k)
    {
      xk=pow(10.0,k/100.0);
      printf("%f %e %e\n",xk,linear_power_spectrum(xk),nonlinear_power_spectrum(xk));
    }
  exit(0);

  dlogr = (log(40)-log(0.1))/(n-1);
  for(k=1;k<=n;++k)
    {
      r=exp(dlogr*(k-1))*0.1;
      printf("boo %f %e %e\n",r,xi_interp(r),xi_linear_interp(r));
    }
  exit(0);

}


double tfunc1(double x)
{
  double p,b;
  x=exp(x)-1;
  b=bias_interp(mass_gt,-1.0);
  p=delta_pdf(x,10.0);
  return((1+b*x)*p);
}

void xi_void_model(int argc, char **argv)
{
  int i,j,n,n1=50;
  double *nv,*vpf,*rv,r1,r2,r,dlogr,dr,y1,fac,x1,x2,x3,fv,fac2;
  FILE *fp;

  fprintf(stderr,"here %s\n",argv[3]);

  fp=openfile(argv[3]);
  n=filesize(fp);

  vpf=dvector(1,n);
  nv =dvector(1,n);
  rv =dvector(1,n);


  for(i=1;i<=n;++i)
    fscanf(fp,"%lf %lf",&rv[i],&nv[i]);
  dr = rv[2]-rv[1];

  r1=log(0.1);
  r2=log(30.0);
  dlogr=(r2-r1)/(n1-1);
  for(i=0;i<n1;++i)
    {
      r=exp(dlogr*i+r1);
      x1 = one_halo_real_space(r);
      x2 = two_halo_real_space(r);
      fac=0;
      fac2=0;
      for(j=1;j<=n;++j)
	{
	  fv = 4./3.*rv[j]*rv[j]*rv[j]*nv[j]*100;
	  if(r<=2*rv[j])
	    y1 = exp(fv/2*(pow(r/2/rv[j],3.0) - 3*r/2/rv[j]+2));
	  else
	    y1 = 1;
	  fac+=y1*dr*nv[j];
	  fac2+=dr*nv[j];
	}	  
      fac/=fac2;
      /* fac = qromo(func_void,1.0,r,midpnt); */
      x3 = fac*(1+x2)-1;
      printf("%e %e %e %e %e\n",r,x1+x2,x1+x3,x3,fac);
    }
  exit(0);
}


