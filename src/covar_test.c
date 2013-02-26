#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

void covar_test()
{
  int i,j,k,n,nd=2,nrot;
  float x[5],x1,x2,x3,x4,x0;
  double **cov,**tmp,*avg,**tmp1,*eval,**evect,*a,*atemp;
  FILE *fp;

  fp = openfile("acc.out");
  n = filesize(fp);

  a = dvector(1,nd);
  atemp = dvector(1,nd);

  cov = dmatrix(1,nd,1,nd);
  evect = dmatrix(1,nd,1,nd);
  avg = dvector(1,nd);
  eval = dvector(1,nd);
  tmp = dmatrix(1,nd,1,nd);
  tmp1 = dmatrix(1,nd,1,1);

  for(j=1;j<=nd;++j)
    {
      avg[j] = 0;
      tmp1[j][1] = 0;
      for(k=1;k<=nd;++k)
	cov[j][k] = 0;
    }

  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%f %d %f %f %f %f",&x0,&j,&x[1],&x[2],&x[3],&x[4]);
      x[1] = log(x[1]);
      x[3] = log(x[3]);
      x[4] = log(x[4]);
      for(j=1;j<=nd;++j)
	{
	  avg[j] += x[j];
	  for(k=1;k<=nd;++k)
	    cov[j][k] += x[j]*x[k];
	}
    }

  for(j=1;j<=nd;++j)
    for(k=1;k<=nd;++k)
      cov[j][k] = cov[j][k]/n - avg[j]*avg[k]/(n*n);

  for(j=1;j<=nd;++j)
    {
      printf("cov %d> ",j);
      for(k=1;k<=nd;++k)
	printf("%f ",cov[j][k]);
      printf("\n");
    }
      printf("\n");

  n = nd;

  jacobi(cov,n,eval,evect,&nrot);

  for(j=1;j<=nd;++j)
    {
      printf("jac %d> ",j);
      for(k=1;k<=nd;++k)
	printf("%f ",evect[j][k]);
      printf("\n");
    }
  printf("\n");

  for(k=1;k<=nd;++k)
    printf("%f ",eval[k]);
  printf(" %d \n",nrot);
  printf("\n");
  

  gaussj(evect,n,tmp1,1);

  for(j=1;j<=nd;++j)
    {
      printf("gj %d> ",j);
      for(k=1;k<=nd;++k)
	printf("%f ",evect[j][k]);
      printf("\n");
    }
      printf("\n");
  

      for(i=1;i<=n;++i)
	atemp[i] = sqrt(eval[i]);
      
      for(i=1;i<=n;++i)
	for(a[i]=0,j=1;j<=n;++j)
	  a[i] += atemp[j]*evect[j][i];
      

  for(k=1;k<=nd;++k)
    printf("%f ",a[k]);
  printf("\n\n");


  exit(0);
}

