#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

#define MAGSOL_R 4.76

float sloan_densities(int n);
float SAM_densities(int n);
void switch_params(int n, float hod[10][6]);
double func_m_one(double m);

void ml_ratios()
{
  int i,j,k,n,nh,nprev,i1,i2,ibin,total=0,imag,jmag,imax,nm,nwp_start,ALPHA=1;
  FILE *fp,*fp1;
  char aa[1000],fname[100];
  float x,*h,dx,xmin,xmax,x1,x2,x3,x4,x5,hmax,hod[10][6],chimin,
    mmin,mmax,dlogm,mass,luminosity,magnitude,n1,n2,n1cen,n2cen,
    luminosity_cen;
  double hod_mean[10][6],hod_err[10][6];

  dx = 0.05;
  xmin = 0.3;
  xmax = 2.0;

  dx = 0.2;
  xmin = -1.5;
  xmax = 1.5;

  nh = (xmax - xmin)/dx+1;
  h = vector(1,nh);
  for(i=1;i<=nh;++i)
    h[i] = 0;
      fprintf(stderr,"here\n");

  for(jmag=3, imag=190; imag <= 220; imag+=5, jmag++)
    {
      //if(imag==220) ALPHA=0;

      for(i=0;i<10;++i)
	for(j=0;j<5;++j)
	  hod_mean[i][j] = hod_err[i][j] = 0;
      

      sprintf(fname,"acc.%d",imag);
      fp = openfile(fname);
      n = filesize(fp);
      i2 = 0;
      for(i=1;i<=n;++i)
	{
	  nprev = i2;
	  if(ALPHA)
	    //fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x,&x4);
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4);
	  else
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4);
	    
	    //if(i<200)continue;
	  ibin = (x - xmin)/dx + 1;
	  if(ibin>nh)continue;
	  h[ibin] ++;
	  total ++;
	  //      h[ibin] += n1 - nprev;
	  //      total += n1 - nprev;
	}
      hmax = 0;
      for(i=1;i<=nh;++i)
	{
	  x = (i-0.5)*dx + xmin;
	  h[i] /= (total*dx);
	  if(h[i]>hmax){ imax = i; hmax = h[i]; }
	  // if(h[i]>hmax && (i-0.5)*dx + xmin < log10(0.7)){ imax = i; hmax = h[i]; }
	}

      rewind(fp);

      chimin = 1e6;
      for(j=0,i=1;i<=n;++i)
	{
	  nprev = i1;
	  if(ALPHA)
	    //fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x4,&x);
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4);
	  else
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x3,&x4,&x);
	  if(i<0)continue;
	  // NO MCUT::
	  x3 = 6;
	  
	  ibin = (x4 - xmin)/dx + 1;
	  //if((x4)>log10(0.5))continue;
	  
	  j++;
	  hod_mean[jmag][0] += pow(10.0,x1);
	  hod_err[jmag][0] += pow(10.0,2*x1);
	  if(ALPHA) {
	    hod_mean[jmag][1] += x2;
	    hod_err[jmag][1] +=  x2*x2;
	  } else {
	    hod_mean[jmag][1] += (pow(10.0,x1)+pow(10.0,x3));
	    hod_err[jmag][1] +=  (pow(10.0,x1)+pow(10.0,x3))*(pow(10.0,x1)+pow(10.0,x3));
	  }
	  hod_mean[jmag][2] += pow(10.0,x3);
	  hod_err[jmag][2] += pow(10.0,2*x3);
	  hod_mean[jmag][3] += pow(10.0,x4);
	  hod_err[jmag][3] += pow(10.0,2*x4);

	  HOD.M1 = pow(10.0,x1);
	  if(ALPHA) HOD.alpha = x2;
	  else HOD.alpha = 1;
	  HOD.M_cut = pow(10.0,x3);
	  HOD.sigma_logM = pow(10.0,x4);
	  GALAXY_DENSITY = sloan_densities(jmag);
	  set_HOD_params();

	  x5 = exp(zbrent(func_m_one,log(HOD.M_low),log(HOD.M_max),1.0E-4));
	  hod_mean[jmag][4] += x5;
	  hod_err[jmag][4] += x5*x5;

	  //	  if(ibin!=imax)continue;
	  if(x<chimin) {
	    hod[jmag][0] = x1;
	    hod[jmag][1] = 1;
	    if(ALPHA) hod[jmag][1] = x2;
	    hod[jmag][2] = x3;
	    hod[jmag][3] = x4;
	    hod[jmag][4] = x5;
	    chimin = x;
	  }
	}
      
      /*
      if(jmag==9)
	{
	  hod[9][0] = 1.585574e+01;
	  hod[9][1] = 9.044981e-01;
	  hod[9][2] =  1.337405e+01;
	  hod[9][3] = -2.145606e-01;
	}
      */
      fclose(fp);
      
      for(i=0;i<=4;++i)
	hod_err[jmag][i] = sqrt(hod_err[jmag][i]/j - hod_mean[jmag][i]*hod_mean[jmag][i]/j/j);
      x1 = sqrt(hod_err[jmag][0]*hod_err[jmag][0] + hod_err[jmag][2]*hod_err[jmag][2]);

      fprintf(stderr,"\nMin Model %d (chi^2 = %f):\n",imag,chimin);
      fprintf(stderr,"> %d elements\n",j);
      fprintf(stderr,"> M1 = %e +/- %e\n",pow(10.0,hod[jmag][0]),hod_err[jmag][0]);
      fprintf(stderr,"> M_one = %e +/- %e\n",hod[jmag][4],hod_err[jmag][4]);
      fprintf(stderr,"> alpha = %e +/- %e\n",hod[jmag][1],hod_err[jmag][1]);
      fprintf(stderr,"> M_cut = %e +/- %e\n",pow(10.0,hod[jmag][2]),hod_err[jmag][2]);
      fprintf(stderr,"> sigma_logM = %e +/- %e\n",pow(10.0,hod[jmag][3]),hod_err[jmag][3]);
      
      printf("BOO1 %e %e %e %e %e\n",hod[jmag][0],hod[jmag][1],hod[jmag][2],hod[jmag][3],hod[jmag][4]);
      printf("BOO2 %e %e %e %e %e\n",hod_err[jmag][0],hod_err[jmag][1],hod_err[jmag][2],hod_err[jmag][3],hod_err[jmag][4]);
    }

  nwp_start = 3;

  switch_params(nwp_start,hod);
  nm = 100;
  mmin = log(HOD.M_low);
  mmax = log(3.0e15);
  dlogm = (mmax - mmin)/(nm-1);

  luminosity = 0;
  for(j=nwp_start;j<=9;++j)
    {
      switch_params(j,hod);
      n1 = sloan_densities(j);
      n2 = sloan_densities(j+1);
      magnitude = -17.5 - j*0.5;
      luminosity += (n1-n2)*pow(10.0,-0.4*(magnitude - 0.25 - MAGSOL_R));
      fprintf(stderr,"mag = %.1f M_min= %e\n",magnitude,HOD.M_min);
    }
  fprintf(stderr,"%e %f\n",luminosity,2.78e11*.3/(luminosity));

  /* HODs
   */
  fp1 = fopen("Q1.dat","w");
  for(i=1;i<=nm;++i)
    {
      mass = exp(dlogm*(i-1) + mmin);
      luminosity = 0;

      fprintf(fp1,"%e ",mass);

      // every other bin
      for(j=nwp_start;j<=9;j+=1)
	{
	  switch_params(j,hod);
	  fprintf(fp1,"%e %e ",N_sat(mass),N_cen(mass));
	}
      fprintf(fp1,"\n");
    }
  fclose(fp1);



  /* Mass to light ratios.
   */
  fp1 = fopen("Q5.dat","w");
  for(i=1;i<=nm;++i)
    {
      mass = exp(dlogm*(i-1) + mmin);
      luminosity = 0;
      luminosity_cen = 0;
      for(j=5;j<=9;++j) //start wtih -20
	{
	  switch_params(j,hod);
	  magnitude = -17.5 - j*0.5;
	  n1 = N_sat(mass);
	  n1cen = N_cen(mass);
	  if(j<9)switch_params(j+1,hod);
	  if(j==9)n2 = 0;
	  else n2 = N_sat(mass);
	  if(j==9)n2cen = 0;
	  else n2cen = N_cen(mass);
	  if(n1-n2>=0)
	    luminosity += pow(10.0,-0.4*(magnitude - MAGSOL_R))*(n1 - n2);
	  if(n1cen-n2cen>=0)
	    luminosity += pow(10.0,-0.4*(magnitude - MAGSOL_R))*(n1cen - n2cen);
	  if(n1cen-n2cen>=0)
	    luminosity_cen += pow(10.0,-0.4*(magnitude - MAGSOL_R))*(n1cen - n2cen);
	  //	  else
	  //fprintf(stdout,"%e %d %f %f\n",mass,j,n1,n2);
	}
      fprintf(fp1,"%e %e %e\n",mass,luminosity,luminosity_cen/luminosity);
      fflush(fp1);
    }
  fclose(fp1);

  /* CLFs
   */
  fp1 = fopen("Q4.dat","w");
  for(j=nwp_start;j<=9;++j)
    {
      magnitude = -17.5 - j*0.5;
      fprintf(fp1,"%e ",magnitude);
      mass = 1.0e13;
      switch_params(j,hod);
      n1 = N_avg(mass);
      if(j<9)switch_params(j+1,hod);
      if(j==9)n2 = 0;
      else n2 = N_avg(mass);
      fprintf(fp1,"%e ",(n1-n2)*0.5);

      mass = 1.0e14;
      switch_params(j,hod);
      n1 = N_avg(mass);
      if(j<9)switch_params(j+1,hod);
      if(j==9)n2 = 0;
      else n2 = N_avg(mass);
      fprintf(fp1,"%e\n",(n1-n2)*0.5);
    }
  fclose(fp1);

  exit(0);

}


void switch_params(int n, float hod[10][6])
{
  HOD.M1 = pow(10.0,hod[n][0]);
  HOD.alpha = hod[n][1];
  HOD.M_cut = pow(10.0,hod[n][2]);
  HOD.sigma_logM = pow(10.0,hod[n][3]);
  GALAXY_DENSITY = sloan_densities(n);
  HOD.M_min = 0;
  set_HOD_params();
}

/*
float sloan_densities(int n)
{
  switch(n) {
  case 1: return 0.02692;
  case 2: return 0.02060;
  case 3: return 0.01507;
  case 4: return 0.01015;
  case 5: return 0.00574;
  case 6: return 0.00308;
  case 7: return 0.00117;
  case 8: return 0.00031;
  case 9: return 0.00006;
  }
  return 0;
}
*/

// these are actually the SAM densities.
float sloan_densities(int n)
{
  switch(n) {
  case 1: return 0.02692;
  case 2: return 0.02060;
  case 3: return 0.01873;
  case 4: return 0.01234;
  case 5: return 0.00773;
  case 6: return 0.00413;
  case 7: return 0.00173;
  case 8: return 0.0005457;
  case 9: return 0.000120;
  }
  return 0;
}

double func_m_one(double m)
{
  m = exp(m);
  return(1 - N_sat(m));
}
