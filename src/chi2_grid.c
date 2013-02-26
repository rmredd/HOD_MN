#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

void wp_input(void);

void chi2_grid(int argc, char**argv)
{
  int i1,i2,i3,i4,n,nhod = 3,i,j;
  double alo[4],ahi[5],*a;
  float avector[100][3],chi2,bias,mmin;
  FILE *fp;
  char fname[100];

  /* 
   * OMEGA_M
   * SIGMA_8
   * VBIAS
   * VBIAS_C
   */

  n = 10;

  alo[0] = 0.18;
  ahi[0] = 0.26;
  alo[1] = 0.84;
  ahi[1] = 0.94;
  alo[2] = 0.5;
  ahi[2] = 1.0;
  alo[3] = 0.0;
  ahi[3] = 0.5;

  a = dvector(1,nhod);

  wp_input();

  if(!ThisTask)
    {
      sprintf(fname,"%s.grid",Task.root_filename);
      fp = fopen(fname,"w");
    }
  MCMC=1;

  for(i3=1;i3<=n;++i3)
    {
      VBIAS = (i3-1)*(ahi[2] - alo[2])/(n-1) + alo[2];
      for(i4=1;i4<=n;++i4)
	{
	  VBIAS_C = (i4-1)*(ahi[3] - alo[3])/(n-1) + alo[3];
	  RESET_COSMOLOGY++;
	  one_halo_real_space(1);
	  two_halo_real_space(1);

	  a[1] = VBIAS;
	  a[2] = VBIAS_C;

	  chi2 = chi2_zspace(a);
	  if(!ThisTask)
	    {
	      printf("GRID %d %d %f %f %e\n",i3,i4,VBIAS,VBIAS_C,chi2);
	      fprintf(fp,"%d %d %f %f %e\n",i3,i4,VBIAS,VBIAS_C,chi2);
	    }
	}
    }
  exit(0);
  
  for(i=66;i<=98;++i)
    {
      sprintf(fname,"xi_%d.fit",i);
      fp = openfile(fname);
      fscanf(fp,"%e %e %e %e %e %e",&chi2,&mmin,&avector[i][0],&avector[i][1],&avector[i][2],&bias);
      fclose(fp);
    }
  wp.ncf = 3;

  for(i1=1;i1<=n;++i1)
    {
      OMEGA_M = (i1-1)*(ahi[0] - alo[0])/(n-1) + alo[0];
      for(i2=1;i2<=n;++i2)
	{
	  SIGMA_8= (i2-1)*(ahi[1] - alo[1])/(n-1) + alo[1];
	  for(i3=1;i3<=n;++i3)
	    {
	      VBIAS = (i3-1)*(ahi[2] - alo[2])/(n-1) + alo[2];
	      for(i4=1;i4<=n;++i4)
		{
		  VBIAS_C = (i4-1)*(ahi[3] - alo[3])/(n-1) + alo[3];
		  RESET_COSMOLOGY++;
		  
		  j = 100.001*SIGMA_8;
		  a[1] = avector[j][0]*OMEGA_M/0.25;
		  a[2] = avector[j][1];
		  a[3] = avector[j][2]*OMEGA_M/0.25;

		  if(!ThisTask)
		    {
		      printf("GOO %f %d %e %e %e\n",SIGMA_8,j,a[1],a[2],a[3]);
		      fflush(stdout);
		    }

		  chi2 = chi2_wp(a);
		  chi2 += chi2_zspace(a);
		  if(!ThisTask)
		    {
		      printf("GRID %d %d %d %d %f %f %f %f %e\n",i1,i2,i3,i4,OMEGA_M,SIGMA_8,VBIAS,VBIAS_C,chi2);
		      fprintf(fp,"%d %d %d %d %f %f %f %f %e\n",i1,i2,i3,i4,OMEGA_M,SIGMA_8,VBIAS,VBIAS_C,chi2);
		    }
		}
	    }
	}
    }
  if(!ThisTask)
    fclose(fp);
  exit(0);
      

}
