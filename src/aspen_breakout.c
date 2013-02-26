#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

void output_hod(char fname[]);
void output_wp(char fname[]);
double func_find_mone(double m);
double func_find_mone2(double m);
double func_find_alpha(double a);
double fsat_g1;

void aspen_breakout()
{
  int i,j,k,i1;
  double x1,x2,x3,x4,x5,x6,r,dr,fsat,m1,mcut,mmin,dlogm;
  char fname[100];
  FILE *fp;

  //goto loop_m1;
  //  goto loop_cvir;
  //goto loop_alpha;
  //goto loop_mcut2;

  /* Vary M_min, find M1 such that fsat is fixed.
   */
  OUTPUT=3;
  HOD.alpha = 1.0;
  for(i=1;i<=5;++i)
    {
      muh(0);
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      RESET_KAISER++;

      sprintf(Task.root_filename,"fsat.%d",i);

      HOD.M_min = pow(10.0,i+9.0);
      HOD.M_cut = 4*HOD.M_min;
      if(i==1) {
	HOD.M1 = 30*HOD.M_min;
	set_HOD_params();
	fsat = fsat_g1 = qromo(func_satfrac,log(HOD.M_min),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
	fprintf(stderr,"%e %e %e %f %f\n",HOD.M_min,HOD.M1,GALAXY_DENSITY,fsat,HOD.M1/HOD.M_min);
	muh(2);
      }
      set_HOD_params();
      HOD.M1 = exp(zbrent(func_find_mone2,log(HOD.M_low/10),log(HOD.M_max),1.0E-4));
      fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
      
      fprintf(stdout,"MUH %e %e %e %f %f\n",HOD.M_min,HOD.M1,GALAXY_DENSITY,fsat,HOD.M1/HOD.M_min);

      sprintf(fname,"fsat_1h.%d",i);
      output_wp(fname);
      sprintf(fname,"fsat_wp.%d",i);
      output_wp(fname);
      sprintf(fname,"fsat_hod.%d",i);
      output_hod(fname);

    }
  exit(0);

  /* Vary sigma8 at fixed HOD/ngal (so change mmin)
   */
  for(i=0;i<=3;++i)
    {
      SIGMA_8 = 0.9*growthfactor((double)i);
      fprintf(stderr,"z=%d, sigma8= %f\n",i,SIGMA_8);
      RESET_COSMOLOGY++;
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      RESET_KAISER++;
      //HOD.M_min = 0;
      set_HOD_params();
      
      sprintf(fname,"SIGMA8a.%d",i);
      output_wp(fname);
    }
  exit(0);
  SIGMA_8=0.9;
  RESET_COSMOLOGY++;

 loop_m1:

  /* Vary fsat by varying M1
   */

  m1 = HOD.M1;
  i1 = 0;
  for(i=-2;i<=2;++i)
    {
      HOD.M1 = m1*pow(10.0,i/10.0);
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      RESET_KAISER++;
      HOD.M_min = HOD.M_low = 0;
      set_HOD_params();
      
      fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
      fprintf(stderr,"M1= %e fsat=%f\n",HOD.M1,fsat);
      
      sprintf(fname,"M1_wp.%d",++i1);
      output_wp(fname);
      sprintf(fname,"M1_hod.%d",i1);
      output_hod(fname);
    }
  HOD.M1 = m1;
  exit(0);

 loop_mcut1:

  /* Vary M_cut, fix M1 and fsat
   */
  mcut = HOD.M_cut;
  HOD.M_min = 0;
  set_HOD_params();

  dlogm = (log(HOD.M1) - log(HOD.M_min))/4;
  mmin = log(HOD.M_min);

  fsat_g1 = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  i1 = 0;
  for(i=0;i<5;++i)
    {
      HOD.M_cut = exp(dlogm*i+mmin);
      HOD.M1 = exp(zbrent(func_find_mone,log(HOD.M_min),log(HOD.M_max),1.0E-4));
      
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      RESET_KAISER++;
      HOD.M_min = HOD.M_low = 0;
      set_HOD_params();

      fprintf(stderr,"M_cut= %e M1= %e %f\n",HOD.M_cut,HOD.M1,HOD.M1/HOD.M_cut);
      
      sprintf(fname,"Mcut1_wp.%d",++i1);
      output_wp(fname);
      sprintf(fname,"Mcut1_hod.%d",i1);
      output_hod(fname);
    }

 loop_mcut2:

  /* Vary M_cut/fsat, keep M_cut/M1 = 1
   */
  mcut = 3.0e13;
  m1 = 3.0e13;

  i1 = 0;
  dlogm = (log(3.0e13) - log(10.0e12))/4;
  for(i=0;i<5;++i)
    {
      HOD.M_cut = HOD.M1 = exp(dlogm*i)*10.0e12;

      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      RESET_KAISER++;
      HOD.M_min = HOD.M_low = 0;
      set_HOD_params();

      fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
      fprintf(stderr,"M_min= %e M1= %e fsat=%f\n",HOD.M_min,HOD.M1,fsat);
      

      sprintf(fname,"Mcut2_wp.%d",++i1);
      output_wp(fname);
      sprintf(fname,"Mcut2_hod.%d",i1);
      output_hod(fname);
    }

 loop_alpha:
  
  /* Vary Mcut as above, but fix f_sat by varying alpha
   */
  mcut = 3.0e13;
  m1 = 3.0e13;

  i1 = 0;
  mmin = log(6.0e12);
  dlogm = (log(3.0e13) - mmin)/4;
  HOD.M_cut = HOD.M1 = exp(mmin);
  
  HOD.alpha = 0.1;
  HOD.M_min = HOD.M_low = 0;
  set_HOD_params();

  fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  fprintf(stderr,"fsat = %f %e\n",fsat,HOD.M_min);
  for(i=0;i<5;++i)
    {
      HOD.M_cut = HOD.M1 = exp(dlogm*i + mmin);

      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      RESET_KAISER++;

      HOD.alpha = zbrent(func_find_alpha,0.05,2.0,1.0E-4);
      fprintf(stderr,"M_cut= %e alpha = %f\n",HOD.M_cut,HOD.alpha);
      

      sprintf(fname,"alpha_wp.%d",++i1);
      output_wp(fname);
      sprintf(fname,"alpha_hod.%d",i1);
      output_hod(fname);
    }
  
 loop_cvir:

  /* Vary cvir with central one above
   */
  HOD.M_cut = HOD.M1 = exp(dlogm*2 + mmin);
  
  HOD.alpha = zbrent(func_find_alpha,0.05,2.0,1.0E-4);
  fprintf(stderr,"M_cut= %e alpha = %f\n",HOD.M_cut,HOD.alpha);

  i1 = 0;
  for(i=0;i<5;++i)
    {
      CVIR_FAC = pow(3.0,i-2);
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      RESET_KAISER++;

      sprintf(fname,"cvir_wp.%d",++i1);
      output_wp(fname);
      sprintf(fname,"cvir_hod.%d",i1);
      output_hod(fname);
    }

  exit(0);
}

double func_find_alpha(double a)
{
  HOD.alpha = a;
  return qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt) - GALAXY_DENSITY;
}

double func_find_mone(double m)
{
  double fsat;
  HOD.M1=exp(m);
  HOD.M_min = 0;
  set_HOD_params();
  fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  return fsat-fsat_g1;
}
double func_find_mone2(double m)
{
  double fsat;
  HOD.M1=exp(m);
  set_HOD_params();
  fsat = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY;
  return fsat-fsat_g1;
}

void output_hod(char fname[])
{
  FILE *fp;
  double dlogm,sig,m;
  int i;

  fp=fopen(fname,"w");
  dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
  for(i=1;i<=100;++i)
    {
      m=exp((i-1)*dlogm)*HOD.M_low;
      sig = N_sat(m)*N_sat(m) + N_cen(m)*(1-N_cen(m));
      fprintf(fp,"%e %e %e %e %e %e\n",
	      m,N_cen(m),N_sat(m),N_avg(m),sig,sig/(N_avg(m)*N_avg(m)));
    }
  fclose(fp);  
}

void output_wp(char fname[])
{
  FILE *fp;
  double r,dr,x1,x2,x3,x4,x5;
  int j;

  fp=fopen(fname,"w");
  dr=(log(70.0)-log(0.05))/49.0;
  for(j=0;j<50;++j)
    {
      r=exp(j*dr+log(0.05));
      x1=one_halo_real_space(r);
      x2=two_halo_real_space(r);
      x3=projected_xi(r);
      //x4 = projected_xi1h(r);
      //x5 = projected_xi2h(r);
      fprintf(fp,"%f %e %e %e %e\n",r,x1,x2,x1+x2,x3);
      fflush(fp);
    }
  fclose(fp);
}
