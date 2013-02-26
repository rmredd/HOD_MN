#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

//#define MAGSOL_R -20.44//4.76
#define MAGSOL_R 4.76

float sloan_densities(int n);
void switch_params(int n, float hod[10][6]);
double func_m_one(double m);

void populate_simulation_clf()
{
  int argc;
  char **argv;

  int i,j,k,n,nh,nprev,i1,i2,ibin,total=0,imag,jmag,imax,nm,nwp_start,ALPHA=1;
  FILE *fp,*fp1,*fp2;
  char aa[1000],fname[100],command[100];
  float x,*h,dx,xmin,xmax,x1,x2,x3,x4,x5,hmax,hod[10][6],chimin,
    mmin,mmax,dlogm,luminosity,magnitude,ngal,
    luminosity_cen;
  double hod_mean[10][6],hod_err[10][6];

  int nchi_arr[10],nchi_tot=0;
  float *hod_arr[10][6];
  int ii[10],nrand,ncombos=0,imag_start=3,bad_cnt[10];

  float prior_hi[10][6],prior_lo[10][6];
  double *mbar_cen,*sigma_logm,*mag_arr,*yy1,*yy2;
  double m1,m2,dlogL_dlogM,mass,n1cen,n2cen,n1,n2;
  double lum,norm,lum_bar,lum_hi,lum_lo,dlum,lum_2bar,ntot;

  // stuff for errors part
  float ***ncen_arr,***nsat_arr;
  float ncen_max,nsat_max,nsat_bar,ncen_bar,nsat_min,ncen_min,min,max;

  float **m1arr,**lum1arr,**s1arr,
    *m1bar,*s1bar,*lum1bar,
    *m1err,*s1err,*lum1err;
  int ichimin[10];
  double xchimin[10];

  // stuff for populate sim
  int SO_FILE = 0;
  float nsat,mag,r,b2,ncen,nc[10];
  int n_wp=9, imass, haloid, insat;
  float xh[3],vh[3];
  double xg[3],vg[3];

  mbar_cen = dvector(1,7);
  sigma_logm = dvector(1,7);
  mag_arr = dvector(1,7);
  yy1 = dvector(1,7);
  yy2 = dvector(1,7);

  for(j=0;j<10;++j)
    {
      for(i=0;i<6;++i)
	{
	  prior_hi[j][i] = 1e6;
	  prior_lo[j][i] = -2;
	}
    }

  for(j=0;j<10;++j)
    bad_cnt[j]=0;

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
  

  // UBER COMPARISON
  imag_start = 3;
  // ERIN ML
  //imag_start = 2;
  // IRO ML
  //imag_start = 1;
  // POPULATING SIMULATION FOR MG2 COMPARISON
  imag_start = 1;

  for(jmag=imag_start, imag=180+5*(imag_start-1); imag <= 220; imag+=5, jmag++)
    {
      //if(imag==220) ALPHA=0; // USE FOR 4params HOD

      chimin = 1.0e6;

      // find chi2 min

      sprintf(fname,"acc.%d",imag);
      fp = openfile(fname);
      n = filesize(fp);
      i2 = 0;
      for(i=1;i<=n;++i)
	{
	  nprev = i2;
	  if(ALPHA)
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4); // 3 free (UBER)
	    //fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x,&x4); //4param free (IRO)
	  else
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4);

	  // PRIORS
	  //if(imag == 210)
	  //if(x<-0.4)continue;
	  //  if(x1-x3<0.5)continue;

	  if(x4<chimin) { chimin=x4; }
	}
      fprintf(stderr,"%d %f\n",imag,chimin);

      // find the total number with Delta\chi^2<1

      rewind(fp);
      nchi_tot=0;
      for(i=1;i<=n;++i)
	{
	  if(ALPHA)
	    //fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x,&x4); // 4param free (IRO)
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4); // 3 param free (UBER)
	  else
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4);
	  // PRIORS
	  //if(imag == 195)
	  // if(x1-x3<0.5)continue;
	  //if(imag == 210)
	  //  if(x<-0.4)continue;

	  if(x4-1<chimin)nchi_tot++;
	}

      nchi_arr[jmag] = nchi_tot;
      hod_arr[jmag][1] = vector(1,nchi_tot);
      hod_arr[jmag][2] = vector(1,nchi_tot);
      hod_arr[jmag][3] = vector(1,nchi_tot);
      hod_arr[jmag][4] = vector(1,nchi_tot);
      hod_arr[jmag][5] = vector(1,nchi_tot);

      fprintf(stderr,"found > %d %d models\n",imag,nchi_tot);

      // read in all acceptable models
      rewind(fp);

      nchi_tot=0;
      for(i=1;i<=n;++i)
	{
	  if(ALPHA)
	    //fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x,&x4);
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4);
	  else
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x3,&x,&x4);
	  if(!ALPHA)x2=1.0;
	  //x3 = 6; //UBER!!
	  
	  // PRIORS
	  // if(imag == 195)
	  //if(x1-x3<0.5)continue;
	  //if(imag == 210)
	  //if(x<-0.4)continue;

	  if(x4-1<chimin) {
	    if(x4-chimin<1.0e-3) { ichimin[jmag] = nchi_tot; xchimin[jmag] = chimin; }
	    nchi_tot++;
	    hod_arr[jmag][1][nchi_tot] = x1;
	    hod_arr[jmag][2][nchi_tot] = x2;
	    hod_arr[jmag][3][nchi_tot] = x3;
	    hod_arr[jmag][4][nchi_tot] = x;
	    hod_arr[jmag][5][nchi_tot] = x4;
	    //printf("MAG%d %e %e %e %e %e\n",imag,x1,x2,x3,pow(10.0,x),x4);
	  }
	    
	}
    }

  fflush(stdout);
  muh(0);
  /*
  for(i=nwp_start;i<=9;++i)
    {
      x1 = x2 = x3 = x4 = 0;
      for(j=1;j<=nchi_tot[i];++j)
	{
	  x1 += hod_arr
  */

  nm = 50;
  mmin = log(HOD.M_low);
  mmax = log(3.0e14);
  dlogm = (mmax - mmin)/(nm-1);

  nwp_start = imag_start;
  goto STATS;


  //goto STATS;

  // go through ALL combinations until a good combo has been found

  /*
  nrand = 1000;
  for(i1=1;i1<=nrand;++i1)
    {
      for(j=nwp_start;j<=9;++j)
	ii[j] = drand48()*nchi_arr[j] + 1;
  
      //printf("%d %d %d\n",ii[3],ii[4],ii[5]);
      */
  //goto ERRORS;

  // new seed
  srand48(23498023);

  for(i1=1;;i1++)
    {
    RESTART_LOOP:

      for(j=nwp_start;j<=9;++j)
	ii[j] = drand48()*nchi_arr[j] + 1;
  
      for(j=nwp_start;j<=9;++j)
	for(k=0;k<4;++k)
	  hod[j][k] = hod_arr[j][k+1][ii[j]];      

      //if(i1==12)
      //goto STATS;

      for(i=1;i<=nm;++i)
	{
	  mass = exp(dlogm*(i-1) + mmin);
	  for(j=nwp_start;j<=9;++j)
	    {
	      switch_params(j,hod);
	      //printf("%d> %e %e %e %e\n",j,hod[j][0],hod[j][1],hod[j][2],hod[j][3]);
	      n1 = N_sat(mass);
	      n1cen = N_cen(mass);
	      if(j<9)switch_params(j+1,hod);
	      if(j==9)n2 = 0;
	      else n2 = N_sat(mass);
	      if(j==9)n2cen = 0;
	      else n2cen = N_cen(mass);
	      //if(n2>n1)break;
	      if(n2>n1 && n1>0.1)break;
	      if(n2cen>n1cen && n1cen>0.1)break;
	    }
	  //muh(j);
	  if(j<10) { bad_cnt[j]++; break; }
	}
      if(i==nm+1){ ncombos++; fprintf(stdout,"%d %d %d> %d\n",ii[3],ii[4],ii[5],ncombos); 
	fprintf(stdout,"BUH %d> %e\n",ncombos, HOD.M1);
	fprintf(stdout,"BUH %d> %d %e\n",ncombos,ii[3],hod[3][1]);
	switch_params(3,hod);
	fprintf(stdout,"BUH %d> %e\n",ncombos,HOD.M1);
	fflush(stdout); goto STATS; }
    }
  printf("%d %d\n",ncombos,nrand);
  for(j=nwp_start;j<=9;++j)
    printf("BAD: %d %d\n",j,bad_cnt[j]);
  exit(0);


  for(j=nwp_start;j<=9;++j)
    for(k=0;k<4;++k)
      hod[j][k] = hod_arr[j][k+1][1];
  

 STATS:

  srand48(23423);
  drand48();
  
  // i think this is for using the minimum chi^2 fits
  for(j=nwp_start;j<=9;++j)
    {
      for(k=0;k<4;++k)
	hod[j][k] = hod_arr[j][k+1][ichimin[j]];      
      muh(ichimin[j]);
      fmuh(xchimin[j]);
    }
  //exit(0);

  // read in the simulation
  fp=openfile(Files.HaloFile);
  sprintf(aa,"%s.mock",Task.root_filename);      
  fp2 = fopen(aa,"w");

  while(!feof(fp))
    {
      if(SO_FILE)
	{
	  fscanf(fp,"%d %lf %f %f %f %f %f %f %f %f",
		 &i,&mass,&x1,&x1,&xh[0],&xh[1],&xh[2],&vh[0],&vh[1],&vh[2]);
	  haloid = i;
	}
      else
	{
	  fscanf(fp,"%d %d %e %e %e %e %e %e %e",
		 &i,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
	  mass=imass*RHO_CRIT*OMEGA_M*pow(RESOLUTION,3.0);
	}
      haloid = i;
      if(mass>HOD.M_max)continue;


      // get the distribution of centrals for this halo mass.
      for(j=nwp_start;j<=9;++j)
	{
	  switch_params(j,hod);
	  n1 = N_cen(mass);
	  if(j<9) 
	    {
	      switch_params(j+1,hod);
	      n2 = N_cen(mass);
	    }
	  else
	    n2 = 0;

	  if(j>nwp_start)
	    nc[j] = nc[j-1]+(n1-n2);
	  else
	    nc[j] = (n1-n2);
	  //printf("%d %e\n",j,nc[j]);
	}
      
      ncen=drand48();
      for(j=nwp_start;j<=n_wp;++j)
	if(ncen<nc[j])break;
      switch_params(j,hod);
      if(j<=n_wp) {
	magnitude = -17.5 - j*0.5;
	fprintf(fp2,"%e %e %e %e %e %e %.1f %e %d\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2],magnitude,mass,1);
      }
      //exit(0);
      //continue;

      /* CLFs
       */
      for(j=nwp_start;j<=9;++j)
	{
	  magnitude = -17.5 - j*0.5;
	  switch_params(j,hod);
	  n1 = N_sat(mass);
	  if(j<9)switch_params(j+1,hod);
	  if(j==9)n2 = 0;
	  else n2 = N_sat(mass);
	  
	  nsat = n1 - n2;
	  insat = poisson_deviate(nsat);
	
	  for(i=0;i<insat;++i)
	    {
	      r = NFW_position(mass,xg);
	      NFW_velocity(mass,vg,mag);
	      for(k=0;k<3;++k)
		{
		  xg[k]+=xh[k];
		  if(xg[k]<0)xg[k]+=BOX_SIZE;
		  if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
		  vg[k]+=vh[k];
		}	
	      fprintf(fp2,"%e %e %e %e %e %e %.1f %e %d\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2],magnitude,mass,0);
	    }	

	}
    }

}
