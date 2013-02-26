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

void uber_comparison()
{
  int argc;
  char **argv;

  int i,j,k,n,nh,nprev,i1,i2,ibin,total=0,imag,jmag,imax,nm,nwp_start,ALPHA=1,SLOGM=1;
  FILE *fp,*fp1;
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

  for(jmag=imag_start, imag=180+5*(imag_start-1); imag <= 220; imag+=5, jmag++)
    {
      //if(imag==220) ALPHA=0; // USE FOR 4params HOD

      chimin = 1.0e6;

      SLOGM=1;
      if(imag<205)SLOGM=0;

      // find chi2 min

      sprintf(fname,"acc.%d",imag);
      fp = openfile(fname);
      n = filesize(fp);
      i2 = 0;
      for(i=1;i<=n;++i)
	{
	  nprev = i2;
	  if(ALPHA)
	    {
	      if(!SLOGM)
		{
		  fscanf(fp,"%s %d %d %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x4); //4param free (slogm)
		}
	      else
		{
		  //fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4); // 3 free (UBER)
		  fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x,&x4); //4param free UBERSAM
		}
	    }
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
	    {
	      if(!SLOGM)
		{
		  fscanf(fp,"%s %d %d %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x4); //4param free (slogm)
		}
	      else
		{
		  //fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4); // 3 free (UBER)
		  fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x,&x4); //4param free UBERSAM
		}
	    }
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
	    {
	      if(!SLOGM)
		{
		  fscanf(fp,"%s %d %d %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x4); //4param free (slogm)
		}
	      else
		{
		  //fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4); // 3 free (UBER)
		  fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x3,&x,&x4); //4param free UBERSAM
		}
	    }
	  else
	    fscanf(fp,"%s %d %d %f %f %f %f %f",aa,&i1,&i2,&x1,&x2,&x,&x4);

	  if(!ALPHA)x2=1.0;
	  if(!SLOGM)x = log10(0.15);
	  //x3 = 6; //UBER!!
	  
	  // PRIORS
	  // if(imag == 195)
	  //if(x1-x3<0.5)continue;
	  //if(imag == 210)
	  //if(x<-0.4)continue;

	  if(x4-1<chimin) {
	    if(x4-chimin<1.0e-3) { ichimin[jmag] = nchi_tot; xchimin[jmag] = chimin; }
	    nchi_tot++;
	    hod_arr[jmag][1][nchi_tot] = x1; //m1
	    hod_arr[jmag][2][nchi_tot] = x2; //alpha
	    hod_arr[jmag][3][nchi_tot] = x3; //mcut
	    hod_arr[jmag][4][nchi_tot] = x;  //slogm
	    hod_arr[jmag][5][nchi_tot] = x4; //chi2
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
  //goto STATS;


  //goto STATS;
  //goto ERRORS;

  // go through ALL combinations until a good combo has been found

  /*
  nrand = 1000;
  for(i1=1;i1<=nrand;++i1)
    {
      for(j=nwp_start;j<=9;++j)
	ii[j] = drand48()*nchi_arr[j] + 1;
  
      //printf("%d %d %d\n",ii[3],ii[4],ii[5]);
      */

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
	  //printf("%d %d %e %e %e %e %e\n",i,j,mass,n2,n1,n2cen,n1cen);
	  if(j<10) { bad_cnt[j]++; break; }
	}
      if(i==nm+1){ ncombos++; fprintf(stdout,"%d %d %d> %d\n",ii[3],ii[4],ii[5],ncombos); 
	fprintf(stdout,"BUH %d> %e\n",ncombos, HOD.M1);
	fprintf(stdout,"BUH %d> %d %e\n",ncombos,ii[3],hod[3][1]);
	switch_params(3,hod);
	fprintf(stdout,"BUH %d> %e\n",ncombos,HOD.M1);
	fprintf(stdout,"BOO %d %d %d %d %d %d %d\n",ii[3],ii[4],ii[5],ii[6],ii[7],ii[8],ii[9]);
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

  
  // make sure to remove this when doing the errors
  /*
  for(j=nwp_start;j<=9;++j)
    {
      for(k=0;k<4;++k)
	hod[j][k] = hod_arr[j][k+1][ichimin[j]];      
      muh(ichimin[j]);
      fmuh(xchimin[j]);
    }
  */
  //exit(0);

  /* Output the correlation functions.
   */
  
  for(j=nwp_start;j<=9;++j)
    {
      switch_params(j,hod);

      printf("%e %e %e %e %e\n",HOD.M_min,HOD.M1,HOD.alpha,HOD.sigma_logM,GALAXY_DENSITY);
      RESET_FLAG_1H++;
      RESET_FLAG_2H++;
      sprintf(Task.root_filename,"M%d",j);
      //if(j==7)wp_minimization(fname);
      //tasks(argc,argv);
    }
  //exit(0);
  

  switch_params(nwp_start,hod);
  nm = 100;
  mmin = log(HOD.M_low);
  mmax = log(1.0e15);
  dlogm = (mmax - mmin)/(nm-1);

  /* Average satellite luminosity.
   */  
  fprintf(stderr,"Starting Q8 (Mean sat magnitude)...\n");

  fp1 = fopen("Q8.dat","w");
  for(i=1;i<=nm;++i)
    {
      ntot = 0;
      mass = exp(dlogm*(i-1) + mmin);
      luminosity = 0;
      luminosity_cen = 0;
      //for(j=5;j<=9;++j) //start wtih -20
      for(j=nwp_start;j<=9;++j) //start wtih -19
	{
	  switch_params(j,hod);
	  magnitude = -17.5 - j*0.5;
	  n1 = N_sat(mass);
	  n1cen = N_cen(mass);
	  if(j==4)ngal = n1+n1cen;
	  if(j<9)switch_params(j+1,hod);
	  if(j==9)n2 = 0;
	  else n2 = N_sat(mass);
	  if(j==9)n2cen = 0;
	  else n2cen = N_cen(mass);
	  if(n1-n2>=0)
	    luminosity += (magnitude - 0.2)*(n1 - n2);
	  if(n1-n2>=0)
	    ntot += (n1-n2);
	}
      fprintf(fp1,"%e %e %e %e\n",mass,luminosity/ntot);
      fflush(fp1);
    }
  fclose(fp1);
  //exit(0);

  for(i=nwp_start;i<=9;++i)
    fmuh(hod[i][1]);

  luminosity = 0;
  for(j=nwp_start;j<=9;++j)
    {
      switch_params(j,hod);
      n1 = sloan_densities(j);
      n2 = sloan_densities(j+1);
      magnitude = -17.5 - j*0.5;
      if(n1-n2>0)
	luminosity += (n1-n2)*pow(10.0,-0.4*(magnitude - 0.25 - MAGSOL_R));
      fprintf(stderr,"mag = %.1f M_min= %e\n",magnitude,HOD.M_min);
      fprintf(stdout,"XX %d %f %f %f\n",ncombos, (10.0,-0.4*(magnitude - 0.25 - MAGSOL_R)),log10(HOD.M_min),
	      HOD.sigma_logM);
      fflush(stdout);
    }
  fprintf(stderr,"%e %f\n",luminosity,2.78e11*.3/(luminosity));

  /* Satellite fraction
   */
  fp1 = fopen("Q6.dat","w");
  for(j=nwp_start;j<=9;j+=1)
    {
      switch_params(j,hod);
      n1 = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt);
      x1 = GALAXY_DENSITY;
      if(j<9) {
	switch_params(j+1,hod);
	n2 = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt);
	x2 = GALAXY_DENSITY;
      } else {
	n2 = 0;
	x2 = 0;
      }
      switch_params(j,hod);
      fprintf(fp1,"%f %f %f\n",-17.5-0.5*j,
	      (n1-n2)/(x1-x2),
	      qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
    }
  fclose(fp1);

  /* Satellite fraction - thresholds
   */
  fp1 = fopen("Q7.dat","w");
  for(j=nwp_start;j<=9;j+=1)
    {
      switch_params(j,hod);
      n1 = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt);
      x1 = GALAXY_DENSITY;
      if(j<-9) {
	switch_params(j+1,hod);
	n2 = qromo(func_satfrac,log(HOD.M_low),log(HOD.M_max),midpnt);
	x2 = GALAXY_DENSITY;
      } else {
	n2 = 0;
	x2 = 0;
      }
      switch_params(j,hod);
      fprintf(fp1,"%f %f %f\n",-17.5-0.5*j,
	      (n1-n2)/(x1-x2),
	      qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt)/GALAXY_DENSITY);
    }
  fclose(fp1);

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

  fprintf(stderr,"Starting Q5 (M/L)...\n");

  /* Mass to light ratios.
   */
  fp1 = fopen("Q5.dat","w");
  for(i=1;i<=nm;++i)
    {
      mass = exp(dlogm*(i-1) + mmin);
      luminosity = 0;
      luminosity_cen = 0;
      //for(j=5;j<=9;++j) //start wtih -20
      for(j=2;j<=9;++j) //start wtih -18.5
	{
	  switch_params(j,hod);
	  magnitude = -17.5 - j*0.5;
	  n1 = N_sat(mass);
	  n1cen = N_cen(mass);
	  if(j==4)ngal = n1+n1cen;
	  if(j<9)switch_params(j+1,hod);
	  if(j==9)n2 = 0;
	  else n2 = N_sat(mass);
	  if(j==9)n2cen = 0;
	  else n2cen = N_cen(mass);
	  if(n1-n2>=0)
	    luminosity += pow(10.0,-0.4*(magnitude - 0.2 - MAGSOL_R))*(n1 - n2);
	  if(n1cen-n2cen>=0)
	    luminosity += pow(10.0,-0.4*(magnitude - 0.2 -  MAGSOL_R))*(n1cen - n2cen);
	  if(n1cen-n2cen>=0)
	    luminosity_cen += pow(10.0,-0.4*(magnitude - 0.2 - MAGSOL_R))*(n1cen - n2cen);
	  //	  else
	  //fprintf(stdout,"%e %d %f %f\n",mass,j,n1,n2);
	}
      //fprintf(fp1,"%e %e %e\n",mass,luminosity,luminosity_cen/luminosity);
      fprintf(fp1,"%e %e %e %e\n",mass,ngal,mass/luminosity,mass/(luminosity*pow(10.0,0.4*0.4)));
      fflush(fp1);
    }
  fclose(fp1);


  fprintf(stderr,"Starting Q4 (CLFs)...\n");


  /* CLFs
   */
  fp1 = fopen("Q4.dat","w");
  for(j=nwp_start;j<=9;++j)
    {
      magnitude = -17.5 - j*0.5;
      fprintf(fp1,"%e ",magnitude);
      mass = 1.0e13;
      switch_params(j,hod);
      n1 = N_sat(mass);
      if(j<9)switch_params(j+1,hod);
      if(j==9)n2 = 0;
      else n2 = N_sat(mass);
      fprintf(fp1,"%e ",(n1-n2)/0.5);

      mass = 1.0e14;
      switch_params(j,hod);
      n1 = N_sat(mass);
      if(j<9)switch_params(j+1,hod);
      if(j==9)n2 = 0;
      else n2 = N_sat(mass);
      fprintf(fp1,"%e\n",(n1-n2)/0.5);
    }
  fclose(fp1);

  /* L_cen, sigma_logL
   */
  fp1 = fopen("Q3.dat","w");
  for(j=nwp_start;j<=9;++j) 
    {
      magnitude = -17.5 - j*0.5;
      switch_params(j,hod);
      mbar_cen[j-2] = log10(HOD.M_min);
      sigma_logm[j-2] = HOD.sigma_logM;
      mag_arr[j-2] = -0.4*(magnitude - MAGSOL_R);
    }
      
  spline(mbar_cen,mag_arr,7,1.0E30,1.0E30,yy1);
  
  for(i=1,j=nwp_start;j<=9;++j,++i)
    {
      n1cen = mbar_cen[i] - 0.02;
      n2cen = mbar_cen[i] + 0.02;
      splint(mbar_cen,mag_arr,yy1,7,n1cen,&n1);
      splint(mbar_cen,mag_arr,yy1,7,n2cen,&n2);
      dlogL_dlogM = (n2-n1)/(n2cen-n1cen);
      sigma_logm[i] = sigma_logm[i]/ROOT2*dlogL_dlogM;
      printf("%e %e %e\n",mag_arr[j-2],mbar_cen[j-2],sigma_logm[j-2]);
      fprintf(fp1,"%e %e %e\n",mag_arr[j-2],mbar_cen[j-2],sigma_logm[j-2]);
    }

  fclose(fp1);
  //exit(0); // exit here if only doing for one model.

  for(i=1;i<=7;++i)
    {
      sprintf(command,"mv Q%d.dat tempfiles2/Q%d.%d",i,i,ncombos);
      fprintf(stderr,"command: %s\n",command);
      system(command);
    }
  if(ncombos<100)
    goto RESTART_LOOP;

  //  exit(0);

 ERRORS:

  fp = fopen("tempfiles2/Q1.1","r");
  fscanf(fp,"%e",&mmin);
  fclose(fp);

  nm = 100;
  //mmin = log(1.499343e+11);
  mmin = log(mmin);
  mmax = log(1.0e15);
  dlogm = (mmax - mmin)/(nm-1);
  ncombos = 100;
 
  /* ERRORS ERRORS ERRORS ERRORS ERRORS ERRORS
   */

  /* Mass to Light ratios
   */
  m1arr = matrix(1,nm,1,ncombos);
  m1bar = vector(1,nm);
  m1err = vector(1,nm);
  s1arr = matrix(1,nm,1,ncombos);
  s1bar = vector(1,nm);
  s1err = vector(1,nm);
  
  for(i=1;i<=ncombos;++i)
    {
      sprintf(fname,"tempfiles2/Q5.%d",i);
      fp = fopen(fname,"r");
      for(j=1;j<=nm;++j)
	fscanf(fp,"%f %f %f",&x1,&m1arr[j][i],&s1arr[j][i]);
      fclose(fp);
    }
  for(i=1;i<=nm;++i)
    m1bar[i] = m1err[i] = 0;
  for(i=1;i<=nm;++i)
    s1bar[i] = s1err[i] = 0;

  for(i=1;i<=nm;++i)
    {
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  m1bar[i] += m1arr[i][j]/ncombos;
	  if(m1arr[i][j]>max)max = m1arr[i][j];
	  if(m1arr[i][j]<min)min = m1arr[i][j];
	}
      m1err[i] = 0.5*(max-min);
  
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  s1bar[i] += s1arr[i][j]/ncombos;
	  if(s1arr[i][j]>max)max = s1arr[i][j];
	  if(s1arr[i][j]<min)min = s1arr[i][j];
	}
      s1err[i] = 0.5*(max-min);
    }  

  fp1 = fopen("Q5.dat","w");
  for(i=1;i<=nm;++i)
    fprintf(fp1,"%e %e %e %e %e\n",exp(dlogm*(i-1)+mmin),m1bar[i],s1bar[i],m1err[i],s1err[i]);
  fclose(fp1);
  
  free_matrix(m1arr,1,nm,1,ncombos);
  free_matrix(s1arr,1,nm,1,ncombos);
  free_vector(s1bar,1,nm);
  free_vector(s1err,1,nm);
  free_vector(m1bar,1,nm);
  free_vector(m1err,1,nm);


  /* satellite fraction
   */
  m1arr = matrix(1,7,1,ncombos);
  m1bar = vector(1,7);
  m1err = vector(1,7);
  lum1arr = matrix(1,7,1,ncombos);
  lum1bar = vector(1,7);
  lum1err = vector(1,7);
  s1arr = matrix(1,7,1,ncombos);
  s1bar = vector(1,7);
  s1err = vector(1,7);
  
  for(i=1;i<=ncombos;++i)
    {
      sprintf(fname,"tempfiles2/Q6.%d",i);
      fp = fopen(fname,"r");
      for(j=1;j<=7;++j)
	fscanf(fp,"%f %f %f",&x1,&m1arr[j][i],&s1arr[j][i]);
      fclose(fp);
    }
  for(i=1;i<=7;++i)
    m1bar[i] = m1err[i] = 0;
  for(i=1;i<=7;++i)
    s1bar[i] = s1err[i] = 0;

  for(i=1;i<=7;++i)
    {
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  m1bar[i] += m1arr[i][j]/ncombos;
	  if(m1arr[i][j]>max)max = m1arr[i][j];
	  if(m1arr[i][j]<min)min = m1arr[i][j];
	}
      m1err[i] = 0.5*(max-min);
  
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  s1bar[i] += s1arr[i][j]/ncombos;
	  if(s1arr[i][j]>max)max = s1arr[i][j];
	  if(s1arr[i][j]<min)min = s1arr[i][j];
	}
      s1err[i] = 0.5*(max-min);
    }  

  fp1 = fopen("Q6.dat","w");
  for(i=1;i<=7;++i)
    fprintf(fp1,"%e %e %e %e %e\n",-18.5-i*0.5,m1bar[i],s1bar[i],m1err[i],s1err[i]);
  fclose(fp1);
  
  /* threshold fsats
   */
  for(i=1;i<=ncombos;++i)
    {
      sprintf(fname,"tempfiles2/Q7.%d",i);
      fp = fopen(fname,"r");
      for(j=1;j<=7;++j)
	fscanf(fp,"%f %f %f",&x1,&m1arr[j][i],&s1arr[j][i]);
      fclose(fp);
    }
  for(i=1;i<=7;++i)
    m1bar[i] = m1err[i] = 0;
  for(i=1;i<=7;++i)
    s1bar[i] = s1err[i] = 0;

  for(i=1;i<=7;++i)
    {
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  m1bar[i] += m1arr[i][j]/ncombos;
	  if(m1arr[i][j]>max)max = m1arr[i][j];
	  if(m1arr[i][j]<min)min = m1arr[i][j];
	}
      m1err[i] = 0.5*(max-min);
  
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  s1bar[i] += s1arr[i][j]/ncombos;
	  if(s1arr[i][j]>max)max = s1arr[i][j];
	  if(s1arr[i][j]<min)min = s1arr[i][j];
	}
      s1err[i] = 0.5*(max-min);
    }  

  fp1 = fopen("Q7.dat","w");
  for(i=1;i<=7;++i)
    fprintf(fp1,"%e %e %e %e %e\n",-18.5-i*0.5,m1bar[i],s1bar[i],m1err[i],s1err[i]);
  fclose(fp1);
  


  /* CLFs
   */
  for(i=1;i<=ncombos;++i)
    {
      sprintf(fname,"tempfiles2/Q4.%d",i);
      fp = fopen(fname,"r");
      for(j=1;j<=7;++j)
	fscanf(fp,"%f %f %f",&x1,&m1arr[j][i],&s1arr[j][i]);
      fclose(fp);
    }
  for(i=1;i<=7;++i)
    m1bar[i] = m1err[i] = 0;
  for(i=1;i<=7;++i)
    s1bar[i] = s1err[i] = 0;

  for(i=1;i<=7;++i)
    {
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  m1bar[i] += m1arr[i][j]/ncombos;
	  if(m1arr[i][j]>max)max = m1arr[i][j];
	  if(m1arr[i][j]<min)min = m1arr[i][j];
	}
      m1err[i] = 0.5*(max-min);
  
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  s1bar[i] += s1arr[i][j]/ncombos;
	  if(s1arr[i][j]>max)max = s1arr[i][j];
	  if(s1arr[i][j]<min)min = s1arr[i][j];
	}
      s1err[i] = 0.5*(max-min);
    }  

  fp1 = fopen("Q4.dat","w");
  for(i=1;i<=7;++i)
    fprintf(fp1,"%f %e %e %e %e\n",-18.5-i*0.5,m1bar[i],s1bar[i],m1err[i],s1err[i]);
  fclose(fp1);


  /* Lcen, sigmaLogL
   */
  for(i=1;i<=ncombos;++i)
    {
      sprintf(fname,"tempfiles2/Q3.%d",i);
      fp = fopen(fname,"r");
      for(j=1;j<=7;++j)
	fscanf(fp,"%f %f %f",&lum1arr[j][i],&m1arr[j][i],&s1arr[j][i]);
      fclose(fp);
    }
  for(i=1;i<=7;++i)
    m1bar[i] = m1err[i] = 0;
  for(i=1;i<=7;++i)
    s1bar[i] = s1err[i] = 0;
  for(i=1;i<=7;++i)
    lum1bar[i] = lum1err[i] = 0;
  
  for(i=1;i<=7;++i)
    {
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  m1bar[i] += m1arr[i][j]/ncombos;
	  if(m1arr[i][j]>max)max = m1arr[i][j];
	  if(m1arr[i][j]<min)min = m1arr[i][j];
	}
      m1err[i] = 0.5*(max-min);
  
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  s1bar[i] += s1arr[i][j]/ncombos;
	  if(s1arr[i][j]>max)max = s1arr[i][j];
	  if(s1arr[i][j]<min)min = s1arr[i][j];
	}
      s1err[i] = 0.5*(max-min);
  
      min = 1E6;
      max = 0;
      for(j=1;j<=ncombos;++j)
	{
	  //if(i==1)printf("%d %f\n",j,lum1arr[i][j]);
	  lum1bar[i] += lum1arr[i][j]/ncombos;
	  if(lum1arr[i][j]>max)max = lum1arr[i][j];
	  if(lum1arr[i][j]<min)min = lum1arr[i][j];
	}
      lum1err[i] = 0;
    }
  
  fp1 = fopen("Q3.dat","w");
  for(i=1;i<=7;++i)
    fprintf(fp1,"%f %f %f %f %f %f\n",m1bar[i],lum1bar[i],s1bar[i],m1err[i],lum1err[i],s1err[i]);
  fclose(fp1);

  /* HOD
   */  
  ncen_arr = f3tensor(1,ncombos,1,4,1,nm);
  nsat_arr = f3tensor(1,ncombos,1,4,1,nm);

  for(i=1;i<=ncombos;++i)
    {
      sprintf(fname,"tempfiles2/Q1.%d",i);
      fp = fopen(fname,"r");
      for(j=1;j<=nm;++j) {
	fscanf(fp,"%f %f %f %f %f %f %f %f %f",&x1,
	       &nsat_arr[i][1][j],&ncen_arr[i][1][j],
	       &nsat_arr[i][2][j],&ncen_arr[i][2][j],
	       &nsat_arr[i][3][j],&ncen_arr[i][3][j],
	       &nsat_arr[i][4][j],&ncen_arr[i][4][j]);
      }
      fclose(fp);
    }
  fp1 = fopen("Q1.dat","w");
  for(i=1;i<=nm;++i)
    {
      mass = (exp(dlogm*(i-1)+mmin));
      fprintf(fp1,"%e ",mass);
      fflush(fp1);
      for(k=1;k<=4;++k)
	fprintf(fp1,"%e %e ",nsat_arr[1][k][i],ncen_arr[1][k][i]);
      for(k=1;k<=4;++k)
	{
	  ncen_min = 1e6;
	  ncen_max = 0;
	  nsat_min = 1E6;
	  nsat_max = 0;
	  ncen_bar = 0;
	  nsat_bar = 0;
	  for(j=1;j<=ncombos;++j)
	    {
	      ncen_bar += ncen_arr[j][k][i]/ncombos;
	      if(ncen_arr[j][k][i]>ncen_max)ncen_max = ncen_arr[j][k][i];
	      if(ncen_arr[j][k][i]<ncen_min)ncen_min = ncen_arr[j][k][i];
	      nsat_bar += nsat_arr[j][k][i]/ncombos;
	      if(nsat_arr[j][k][i]>nsat_max)nsat_max = nsat_arr[j][k][i];
	      if(nsat_arr[j][k][i]<nsat_min)nsat_min = nsat_arr[j][k][i];
	    }
	  fprintf(fp1,"%e %e %e %e ",nsat_min,ncen_min,nsat_max,ncen_max);
	}
      fprintf(fp1,"\n");
    }
  fclose(fp1);

  exit(0);

}
