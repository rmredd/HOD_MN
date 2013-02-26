#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

#define OUTPUTZ 1
#define ZFORMAT 2

/* Local functions.
 */
void initialize_chi2_zspace(void);
void input_datafiles(void);
double esys_mono(double r);
double esys_quad(double r);
double esys_half(double r);

double chi2_zspace(double *a)
{
  static double **xiz;
  int i,j,na=0,nb=0,nc=0,i1,i2,j1,j2;
  double chi2d=0,chi2a=0,chi2b=0,chi2c=0,chi2,e,rhalf[100],**tmp,**tmp2;
  double x,rslo,rplo,drs,drp,dx1,rshi,rphi,x1,x2,psend=0,precv;
  static int iter=0,flag=1;


  if(flag)
    {
      flag=0;
      COVARZ=0;
      Work.chi2=1;
      Work.SysErrFlag=0;
      initialize_chi2_zspace();
      input_datafiles();
      if(Work.izspace)
	xiz = dmatrix(1,Work.n_z,1,Work.n_z);

      if(Work.ihalf_covar)
	{
	  tmp=dmatrix(1,Work.n_half,1,1);
	  tmp2=dmatrix(1,Work.n_half,1,Work.n_half);
	  for(i=0;i<Work.n_half;++i)
	    for(j=0;j<Work.n_half;++j)
	      {
		tmp2[i+1][j+1]=Work.covar_h[i][j];
	      }
	  gaussj(tmp2,Work.n_half,tmp,1);
	  for(i=1;i<=Work.n_half;++i)
	    for(j=1;j<=Work.n_half;++j)
	      {
		Work.covar_h[i-1][j-1]=tmp2[i][j];
		if(i==j && !ThisTask)
		  printf("DIAG HALF %d %d %e\n",i,j,Work.covar_q[i-1][j-1]);
	      }
	  free_dmatrix(tmp,1,Work.n_half,1,1);
	  free_dmatrix(tmp2,1,Work.n_half,1,Work.n_half);
	}
      if(Work.iquad_covar)
	{
	  tmp=dmatrix(1,Work.n_quad,1,1);
	  tmp2=dmatrix(1,Work.n_quad,1,Work.n_quad);
	  for(i=0;i<Work.n_quad;++i)
	    for(j=0;j<Work.n_quad;++j)
	      {
		tmp2[i+1][j+1]=Work.covar_q[i][j];
	      }
	  gaussj(tmp2,Work.n_quad,tmp,1);
	  for(i=1;i<=Work.n_quad;++i)
	    for(j=1;j<=Work.n_quad;++j)
	      {
		Work.covar_q[i-1][j-1]=tmp2[i][j];
		if(i==j && !ThisTask)
		  printf("DIAG QUAD %d %d %e\n",i,j,Work.covar_q[i-1][j-1]);
	      }
	  free_dmatrix(tmp,1,Work.n_quad,1,1);
	  free_dmatrix(tmp2,1,Work.n_quad,1,Work.n_quad);
	}
      if(Work.imono_covar)
	{
	  tmp=dmatrix(1,Work.n_mono,1,1);
	  tmp2=dmatrix(1,Work.n_mono,1,Work.n_mono);
	  for(i=0;i<Work.n_mono;++i)
	    for(j=0;j<Work.n_mono;++j)
	      {
		tmp2[i+1][j+1]=Work.covar_m[i][j];
	      }
	  gaussj(tmp2,Work.n_mono,tmp,1);
	  for(i=1;i<=Work.n_mono;++i)
	    for(j=1;j<=Work.n_mono;++j)
	      {
		Work.covar_m[i-1][j-1]=tmp2[i][j];
		if(i==j && !ThisTask)
		  printf("DIAG MONO %d %d %e\n",i,j,Work.covar_m[i-1][j-1]);
	      }
	  free_dmatrix(tmp,1,Work.n_mono,1,1);
	  free_dmatrix(tmp2,1,Work.n_mono,1,Work.n_mono);
	}
    }

  if(!LINEAR_PSP)
    nonlinear_sigmac(8.0);

  RESET_PVZ=1;

  Work.em_abs=Work.eq_abs=0;
  Work.percentage_error=0;

  if(Work.imono || Work.iquad)
    {
      xi_multipoles();
       
      if(Work.imono && !Work.imono_covar)
	{
	  for(i=0;i<Work.n_mono;++i)
	    {
	      if(Work.r_mono[i]>=Work.rmlo && Work.r_mono[i]<=Work.rmhi)
		{
		  e=0;//esys_mono(Work.r_mono[i])*Work.xi_mono[i];
		  chi2a+=(Work.xi_mono[i]-Work.data_m[i])*(Work.xi_mono[i]-Work.data_m[i])/
		    (Work.err_m[i]*Work.err_m[i] + e*e);
		  if(!ThisTask && OUTPUTZ)
		    printf("CHIMONO%d %d %f %e %e %e %e\n",Work.imodel,iter,Work.r_mono[i],
			   Work.xi_mono[i],Work.data_m[i],Work.err_m[i],e);
		}
	    }
	}
	  
      if(Work.imono && Work.imono_covar)
	{
	  for(i=0;i<Work.n_mono;++i)
	    for(j=0;j<Work.n_mono;++j)
	      {
		chi2a+=(Work.xi_mono[i]-Work.data_m[i])*(Work.xi_mono[j]-Work.data_m[j])*
		  Work.covar_m[i][j];
		if(!ThisTask && OUTPUTZ && i==j)
		  printf("CHIMONO %d %f %e %e %e\n",iter,Work.r_mono[i],
			 Work.xi_mono[i],Work.data_m[i],Work.err_m[i]);
	      }
	}

      if(Work.iquad && !Work.iquad_covar)
	{
	  for(i=0;i<Work.n_quad;++i)
	    {
	      if(Work.r_quad[i]>=Work.rqlo && Work.r_quad[i]<=Work.rqhi)
		{
		  e=esys_quad(Work.r_quad[i]);
		  if(Work.r_quad[i]<3)e*=Work.xi_quad[i];

		  /* Try percentage error on (1+abs(Q)) -- originally 0.02
		   */
		  e = 0.013*(1+mabs(Work.xi_quad[i]));
		  if(Work.xi_quad[i]>0.2)
		    e = 0.013*(Work.xi_quad[i]);

		  x1=(Work.xi_quad[i]-Work.data_q[i])*
		    (Work.xi_quad[i]-Work.data_q[i])/
		    (Work.err_q[i]*Work.err_q[i] + e*e);
		  chi2b+=x1;
		  if(!ThisTask && OUTPUTZ) {
		    printf("CHIQUAD%d %d %f %e %e %e %e %e\n",
			   Work.imodel,iter,Work.r_quad[i],
			   Work.xi_quad[i],Work.data_q[i],Work.err_q[i],e,x1);
		    fflush(stdout); }
		}
	    }
	}      
    }
  if(Work.iquad && Work.iquad_covar)
    {

      for(i=0;i<Work.n_quad;++i)
	for(j=0;j<Work.n_quad;++j)
	  {
	    chi2b+= (Work.xi_quad[i]-Work.data_q[i])*(Work.xi_quad[j]-Work.data_q[j])*
	      (Work.covar_q[i][j]);
	    if(!ThisTask && OUTPUTZ && i==j)
	      printf("CHIQUAD %d %f %e %e %e\n",
		     iter,Work.r_quad[j],
		     Work.xi_quad[j],Work.data_q[j],Work.err_q[i]);
	  }
    }

  if(Work.ihalf && !Work.ihalf_covar)
    {
      if(!ThisTask)printf("HERE %d\n",Work.n_half);
      calc_rhalf(Work.r_half,rhalf,Work.n_half);
      for(i=0;i<Work.n_half;++i)
	{
	  // if(Work.r_half[i]<Work.rhlo || Work.r_half[i]>Work.rhhi)continue;
	  e=esys_half(Work.r_half[i])*rhalf[i];
	  e=0;
	  /*rhalf=small_scale_measure(Work.r_half[i]);*/
	  chi2c+=(rhalf[i]-Work.data_h[i])*(rhalf[i]-Work.data_h[i])/
	    (Work.err_h[i]*Work.err_h[i] + e*e);
	  if(!ThisTask && OUTPUTZ)
	    printf("CHIHALF%d %d %f %e %e %e %e\n",Work.imodel,iter,Work.r_half[i],
		   rhalf[i],Work.data_h[i],Work.err_h[i],e);
	}
    }
  if(Work.ihalf && Work.ihalf_covar)
    {
      calc_rhalf(Work.r_half,rhalf,Work.n_half);
      for(i=0;i<Work.n_half;++i)
	for(j=0;j<Work.n_half;++j)
	  {
	    chi2c+=(rhalf[i]-Work.data_h[i])*(rhalf[j]-Work.data_h[j])*
	      (Work.covar_h[i][j]);
	    if(!ThisTask && OUTPUTZ && i==j)
	      printf("CHIHALF %d %f %e %e %e\n",
		     iter,Work.r_half[j],
		     rhalf[j],Work.data_h[j],Work.err_h[i]);
	    /*
	    if(OUTPUTZ && i==j)
	      printf("CPUHALF%02d %d %f %e %e %e\n",ThisTask,
		     iter,Work.r_half[j],
		     rhalf[j],Work.data_h[j],Work.err_h[i]);
	    */
	  }
    }


  if(Work.izspace)
    {
      for(i=0;i<Work.n_z;++i)
	for(j=0;j<Work.n_z;++j)
	  {
	    if(!i)continue;
	    if(i>12 || j>12)continue;

	    rplo=pow(10.0,-1.12+0.2*j);
	    if(!j)rplo=0;
	    rphi=pow(10.0,-1.12+0.2*(j+1));
	    rslo=pow(10.0,-1.12+0.2*i);
	    if(!i)rslo=0;
	    rshi=pow(10.0,-1.12+0.2*(i+1));
	    /*
	    x2 = integrated_bin(rslo,rplo,rshi-rslo,rphi-rplo,10);
	    x=one_halo(Work.rsigma[i][j],Work.rpi[i][j]) + 
	      two_halo(Work.rsigma[i][j],Work.rpi[i][j]);
	    */
	    x = xi2d_interp(rslo,rplo,rshi,rphi);
	    e = 0.03*x;

	      chi2d+=(x-Work.data_z[i][j])*(x-Work.data_z[i][j])/
		(Work.err_z[i][j]*Work.err_z[i][j] + e*e);
	    if(!ThisTask && OUTPUTZ){
	      printf("CHIZSPACE %d %f %f %e %e %e %e %e\n",iter,Work.rsigma[i][j],
		     Work.rpi[i][j],Work.data_z[i][j],Work.err_z[i][j],x,x2,x1);
	      fflush(stdout); }
	  }
    }

  chi2=chi2a+chi2b+chi2c+chi2d;
  ++iter;
  if(!ThisTask)
    {
      printf("ITERZ %d %e %e %e %e %e\n",iter,chi2,chi2a,chi2b,chi2c,chi2d);
      fflush(stdout);
    }
  //  printf("CPUCHIZ%02d %d %e %e %e %e %e\n",ThisTask,iter,chi2,chi2a,chi2b,chi2c,chi2d);
  //fflush(stdout);
  return(chi2);
}

void input_datafiles()
{
  char a[1000];
  FILE *fp;
  int i,n,j,i1,i2,nn,*indx;
  float x1,x2,x3,x4,x5,x6;
  double **tmp,**tmp2,dd,*col;

  if(Work.imono)
    {
      if(!(fp=fopen(Work.monofile,"r")))
	{
	  fprintf(stdout,"ERROR: Could not open [%s]\n",Work.monofile);
	  endrun("Exiting with file read error.");
	}
      i=-1;
      while(!feof(fp))
	{
	  i++;
	  fgets(a,1000,fp);
	}
      Work.n_mono=i;
      Work.nrad=i;
      rewind(fp);
      for(i=0;i<Work.n_mono;++i)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  Work.rad[i]=Work.r_mono[i]=x1;
	  Work.data_m[i]=x2;
	  Work.err_m[i]=x3;
	}
      fclose(fp);

      if(Work.imono_covar)
	{
	  fp = openfile(Work.monocovarfile);
	  n = filesize(fp);
	  if(n!=Work.n_mono*Work.n_mono)
	    {
	      fprintf(stdout,"ERROR: filesize mismatch monofiles %d vs %d\n",n,Work.n_mono);
	      exit(0);
	    }
	  for(i=0;i<Work.n_mono;++i)
	    for(j=0;j<Work.n_mono;++j)
	      fscanf(fp,"%lf",&Work.covar_m[i][j]);
	  fclose(fp);
	}
    }

  if(Work.iquad)
    {
      if(!(fp=fopen(Work.quadfile,"r")))
	{
	  fprintf(stdout,"ERROR: Could not open [%s]\n",Work.quadfile);
	  endrun("Exiting with file read error.");
	}
      n = filesize(fp);
      /*
      Work.n_quad = 0;
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f",&x1);
	  fgets(a,1000,fp);
	  if(x1>=Work.rqlo && x1<=Work.rqhi)
	    Work.n_quad++;
	}
      */
      Work.n_quad = n;
      Work.nrad=Work.n_quad;
      rewind(fp);
      for(i=1;i<=n-Work.n_quad;++i)
	fgets(a,1000,fp);

      for(i=0;i<Work.n_quad;++i)
	{
	  if(ZFORMAT==1)
	    fscanf(fp,"%f %f %f",&x1,&x2,&x3); /* n-body format */
	  if(ZFORMAT==2)
	    fscanf(fp,"%f %f %f",&x1,&x2,&x3); /* n-body format */
	    /* fscanf(fp,"%f %f %f %f %f",&x1,&x4,&x5,&x2,&x3);*/ /* SDSS format */ 
	  Work.rad[i]=Work.r_quad[i]=x1;
	  Work.data_q[i]=x2;
	  Work.err_q[i]=x3;
	}
      fclose(fp);
    }

  if(Work.iquad_covar)
    {
      if(!(fp=fopen(Work.quadcovarfile,"r")))
	{
	  fprintf(stdout,"ERROR: Could not open [%s]\n",Work.quadcovarfile);
	  endrun("Exiting with file read error.");
	}
      i = filesize(fp);
      if(Work.n_quad!=(int)sqrt(i))
	{
	  fprintf(stderr,"ERROR filesize mismatch between quad/covar\n");
	  exit(0);
	}
      for(i=0;i<Work.n_quad;++i)
	for(j=0;j<Work.n_quad;++j)
	{
	  fscanf(fp,"%f",&x1);
	  Work.covar_q[i][j] = x1;
	}
      fclose(fp);
    }

  if(Work.ihalf)
    {
      if(!(fp=fopen(Work.halffile,"r")))
	{
	  fprintf(stdout,"ERROR: Could not open [%s]\n",Work.halffile);
	  endrun("Exiting with file read error.");
	}
      Work.n_half = filesize(fp);
      for(i=0;i<Work.n_half;++i)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  Work.r_half[i]=x1;
	  Work.data_h[i]=x2;
	  Work.err_h[i]=x3;
	}
      fclose(fp);
    }

  if(Work.ihalf_covar)
    {
      fp = openfile(Work.halfcovarfile);
      n = filesize(fp);
      if((int)sqrt(n) != Work.n_half)
	{
	  fprintf(stderr,"ERROR: file size mismatch rhallf+covar files %d %d\n",n,Work.n_half);
	  exit(0);
	}
      for(i=0;i<Work.n_half;++i)
	for(j=0;j<Work.n_half;++j)
	{
	  fscanf(fp,"%f",&x1);
	  Work.covar_h[i][j] = x1;
	}
      fclose(fp);
    }      

  if(!ThisTask)
    printf("Done reading data. %d %d %d\n",Work.n_mono,Work.n_quad,Work.n_half);

  if(Work.use_asymptotic_values)
    {
      for(x1=j=0,i=Work.n_quad-1;i>=0;--i)
	{
	  if(Work.r_quad[i]<=Work.rqhi)
	    {
	      j++;
	      x1+=Work.data_q[i];
	    }
	  if(j==4)break;
	}
      Work.eq_abs=Work.percentage_error*(x1/j);

      for(x1=j=0,i=Work.n_mono-1;i>=0;--i)
	{
	  if(Work.r_mono[i]<=Work.rmhi)
	    {
	      j++;
	      x1+=Work.data_m[i];
	    }
	  if(j==4)break;
	}
      Work.em_abs=Work.percentage_error*(x1/j);

      printf("ABSOLUTE ERRORS: eq= %f em= %f\n",Work.eq_abs,Work.em_abs);
    }
  fflush(stdout);


  if(Work.izspace)
    {
      if(!(fp=fopen(Work.zfile,"r")))
	{
	  fprintf(stdout,"ERROR: Could not open [%s]\n",Work.zfile);
	  endrun("Exiting with file read error.");
	}
      i=-1;
      while(!feof(fp))
	{
	  i++;
	  fgets(a,1000,fp);
	}
      Work.n_z=(int)sqrt(i);
      rewind(fp);
      for(i=0;i<Work.n_z;++i)
	for(j=0;j<Work.n_z;++j)
	  {	
	    /* This is idit's standard format
	     */
	    fscanf(fp,"%d %d %f %f %f %f %f",&i1,&i2,&x1,&x2,&x3,&x4,&x5);
	    Work.rsigma[i][j]=x2;
	    Work.rpi[i][j]=x1;
	    Work.data_z[i][j]=x3;
	    Work.err_z[i][j]=x4;
	    continue;

	    /* THis is my format for truncating the files.
	     */
	    fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	    Work.rsigma[i][j]=x1;
	    Work.rpi[i][j]=x2;
	    Work.data_z[i][j]=x3;
	    continue;

	    /*
	    fscanf(fp,"%f %f",&x1,&x2);
	    Work.data_z[i][j]=x1;
	    Work.err_z[i][j]=x2;
	    */
	  }
      fclose(fp);

    }
}


/* Until I put all this stuff into the batch file, I'll simply
 * hard-wire it.
 */
void initialize_chi2_zspace()
{
  /* This top section is now the control panel for the pm fake data.
   */
  Work.rmlo = 0;
  Work.rmhi = 100;

  Work.ihalf = 1;
  Work.iquad = 1;
  Work.imono = 1;
  Work.i_quad2mono = 1;
  Work.i_monopole = 1;

  Work.imono_covar = 1;
  Work.iquad_covar = 1;
  Work.ihalf_covar = 1;

  /* For the vbias fits.
   */

  Work.ihalf = 0;
  Work.iquad = 0;
  Work.imono = 1;
  Work.i_quad2mono = 0;
  Work.i_monopole = 1;

  Work.imono_covar = 0;
  Work.iquad_covar = 0;
  Work.ihalf_covar = 0;


  /* These are the Millennium Rud data (local)
   */
  sprintf(Work.monofile,
	  "xi_mono.data");
  sprintf(Work.quadfile,
	  "xi_20.data");
  sprintf(Work.halffile,
	  "rhalf.data");

  sprintf(Work.monocovarfile,
	  "xi_mono.covar");
  sprintf(Work.quadcovarfile,
	  "xi_20.covar");
  sprintf(Work.halfcovarfile,
	  "rhalf.covar");
  return;

  /* THese are for the actual SDSS data,
   * The "+_lin" files are using absolute errors, rather than propto errors (from PM mocks)
   */
  sprintf(Work.monofile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/xi_mono.dat");
  sprintf(Work.quadfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/xi_20_lin.dat");
  sprintf(Work.halffile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/rhalf.dat");

  sprintf(Work.monocovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/xi_mono.covar");
  sprintf(Work.quadcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/xi_20_lin.covar");
  sprintf(Work.halfcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/rhalf.covar");
  return;

  /* THese are the actual SDSS data where I've cut out any r<1 data.
   * also _20 is linear.
   */
  sprintf(Work.monofile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/xi_mono_LS.dat");
  sprintf(Work.quadfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/xi_20_LS.dat");
  sprintf(Work.halffile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/rhalf_LS.dat");

  sprintf(Work.monocovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/xi_mono_LS.covar");
  sprintf(Work.quadcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/xi_20_LS.covar");
  sprintf(Work.halfcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data/rhalf_LS.covar");
  return;


  
  /* THese are fiber-collision WMAP3 mock data
   */
  sprintf(Work.monofile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data2/xi_mono.dat");
  sprintf(Work.quadfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data2/xi_20.dat");
  sprintf(Work.halffile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data2/rhalf.dat");

  sprintf(Work.monocovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data2/xi_mono.covar");
  sprintf(Work.quadcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data2/xi_20.covar");
  sprintf(Work.halfcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/dry_runs/data2/rhalf.covar");
  return;
  
  /* These are for the 384 box, new WMAP3 cosmology
   */
  sprintf(Work.monofile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_WMAP/xi_mono.data");
  sprintf(Work.quadfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_WMAP/xi_20.data");
  sprintf(Work.halffile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_WMAP/rhalf.data");

  sprintf(Work.monocovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_WMAP/xi_mono.covar");
  sprintf(Work.quadcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_WMAP/xi_20.covar");
  sprintf(Work.halfcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_WMAP/rhalf.covar");

  return;


  /* THese are for the 400 box, concordance cosmology
   */
  sprintf(Work.monofile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_results/xi_mono.data");
  sprintf(Work.quadfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_results/xi_20.data");
  sprintf(Work.halffile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_results/rhalf.data");

  sprintf(Work.monocovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_results/xi_mono.covar");
  sprintf(Work.quadcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_results/xi_20.covar");
  sprintf(Work.halfcovarfile,
	  "/home/tinker/cosmo/SDSS/zspace_analysis/pm_mock_data/normalized_results/rhalf.covar");
  return;
  

  /* Here is the old stuff.
   */
  Work.ihalf=1;
  Work.imono=1;
  Work.iquad=1;
  Work.izspace=0;
  Work.iout=1;

  Work.imodel=3;

  sprintf(Work.monofile,"xi_mono.0%d.data",Work.imodel);
  sprintf(Work.quadfile,"xi_quad.0%d.data",Work.imodel);
  sprintf(Work.halffile,"xi_half.0%d.data",Work.imodel);
  sprintf(Work.zfile,"xi_z.data");
  sprintf(Work.covarz_file,"xi_z.0%d.covar",Work.imodel);
  sprintf(Work.esysfile,"esys.dat");

  Work.rmlo=0.1;
  Work.rmhi=30;
  Work.rqlo=1.0;
  Work.rqhi=40;
  Work.rhlo=0.0;
  Work.rhhi=30.0;
  Work.zlo=0.1;

  Work.use_asymptotic_values=0;
  Work.em_per=0.080;
  Work.eq_abs=0.112;
  Work.eh_per=0.111;

  Work.em_per=0.036;
  Work.eq_abs=0.080;
  Work.eh_per=0.090;

  /* If using n-body data, stop here
   */
  if(ZFORMAT==1)return;

  /* This is the setp for SDSS work
   */
  Work.SDSS_bins = 1;
  Work.i_quad2mono = 0;

  Work.ihalf = 1;
  Work.ihalf_covar = 1;
  Work.iquad = 0;
  Work.iquad_covar = 0;
  Work.imono = 1;
  Work.imono_covar = 1;
  Work.izspace = 0;
  
  sprintf(Work.halffile,"/home/tinker/cosmo/SDSS/zData/jack_21.0/rhalf.dat");
  sprintf(Work.halfcovarfile,"/home/tinker/cosmo/SDSS/zData/jack_21.0/rhalf.covar");

  sprintf(Work.quadfile,"/home/tinker/cosmo/SDSS/zData/Q_21.0.dat");
  /* sprintf(Work.quadfile,"quad.dat");  */
  sprintf(Work.quadfile,"/home/tinker/cosmo/SDSS/zData/quad.21.data");
  sprintf(Work.quadcovarfile,"/home/tinker/cosmo/SDSS/zData/quad.21.covar");

  /* Mock data.
   */
  sprintf(Work.halffile,"/home/tinker/cosmo/SDSS/zspace_analysis/mock_data/rhalf.dat");
  sprintf(Work.halfcovarfile,"/home/tinker/cosmo/SDSS/zspace_analysis/mock_data/rhalf.covar");

  sprintf(Work.quadfile,"/home/tinker/cosmo/SDSS/zspace_analysis/mock_data/quad.dat");
  sprintf(Work.quadcovarfile,"/home/tinker/cosmo/SDSS/zspace_analysis/mock_data/quad.covar");

  if(Work.i_quad2mono)
    {
      sprintf(Work.quadfile,"/home/tinker/cosmo/SDSS/zspace_analysis/mock_data/xi_20.dat");
      sprintf(Work.quadcovarfile,"/home/tinker/cosmo/SDSS/zspace_analysis/mock_data/xi_20.covar");
    }
  sprintf(Work.monofile,"/home/tinker/cosmo/SDSS/zspace_analysis/mock_data/mono.dat");
  sprintf(Work.monocovarfile,"/home/tinker/cosmo/SDSS/zspace_analysis/mock_data/mono.covar");

}


double esys_mono(double r)
{
  static int flag=1,n;
  static double *x,*y,*y2;
  FILE *fp;
  int i;
  double e,m,b;
  float x1,x2,x3,x4,x5;
  
  if(Work.SysErrFlag==0)return(0);

  return 0.03;

  if(flag)
    {
      flag=0;
      if(!(fp=fopen(Work.esysfile,"r")))
	{
	  fprintf(stdout,"ERROR opening [%s]\n",Work.esysfile);
	  endrun("error opening esys file");
	}
      n=filesize(fp);
      x=dvector(1,n);
      y=dvector(1,n);
      y2=dvector(1,n);

      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f %f %f %f",&x1,&x2,&x3,&x4,&x5);
	  x[i]=x1;
	  y[i]=x2;
	  /*
	  y[i]*=5./3.;
	  if(y[i]<0.07)y[i]=0.07;
	  */
	}
      spline(x,y,n,1.0E+30,1.0E+30,y2);
      fclose(fp);
    }
  /*
  for(i=1;i<=n;++i)
    if(r>x[i])break;
  if(i==1)return(y[1]);
  if(i>=n)return(y[n]);
  m=(y[i-1]-y[i])/(x[i-1]-x[i]);
  b=y[i]-m*x[i];
  return(m*r+b);
  */
  splint(x,y,y2,n,r,&e);   
  return(e);
}

double esys_quad(double r)
{
  static int flag=1,n;
  static double *x,*y,*y2;
  FILE *fp;
  int i;
  double e,m,b;
  float x1,x2,x3,x4,x5;

  if(Work.SysErrFlag==0)return(0);

  if(r<3)return(0.03);
  return 0.005;

  if(r<3)return(0.04);
  return 0.006;

  if(flag)
    {
      flag=0;
      if(!(fp=fopen(Work.esysfile,"r")))
	{
	  fprintf(stdout,"ERROR opening [%s]\n",Work.esysfile);
	  endrun("error opening esys file");
	}
      n=filesize(fp);
      x=dvector(1,n);
      y=dvector(1,n);
      y2=dvector(1,n);

      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f %f %f %f",&x1,&x2,&x3,&x4,&x5);
	  x[i]=x1;
	  y[i]=x3;
	}
      spline(x,y,n,1.0E+30,1.0E+30,y2);
      fclose(fp);
    }
  /*
  for(i=1;i<=n;++i)
    if(r>x[i])break;
  if(i==1)return(y[1]);
  if(i>=n)return(y[n]);
  m=(y[i-1]-y[i])/(x[i-1]-x[i]);
  b=y[i]-m*x[i];
  return(m*r+b);
  */

  splint(x,y,y2,n,r,&e);   
  return(e);
}

double esys_half(double r)
{
  static int flag=1,n;
  static double *x,*y,*y2;
  FILE *fp;
  int i;
  double e,m,b;
  float x1,x2,x3,x4,x5;

  if(Work.SysErrFlag==0)return(0);

  return 0.03;
  return 0.04;

  if(flag)
    {
      flag=0;
      if(!(fp=fopen(Work.esysfile,"r")))
	{
	  fprintf(stdout,"ERROR opening [%s]\n",Work.esysfile);
	  endrun("error opening esys file");
	}
      n=filesize(fp);
      x=dvector(1,n);
      y=dvector(1,n);
      y2=dvector(1,n);

      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f %f %f %f",&x1,&x2,&x3,&x4,&x5);
	  x[i]=x4;
	  y[i]=x5;
	  /*
	  y[i]*=5./3.;
	  if(y[i]<0.07)y[i]=0.07;
	  */
	}
      spline(x,y,n,1.0E+30,1.0E+30,y2);
      fclose(fp);
    }
  /*
  for(i=1;i<=n;++i)
    if(r>x[i])break;
  if(i==1)return(y[1]);
  if(i>=n)return(y[n]);
  m=(y[i-1]-y[i])/(x[i-1]-x[i]);
  b=y[i]-m*x[i];
  return(m*r+b);
  */
  splint(x,y,y2,n,r,&e);  
  return(e);
}
