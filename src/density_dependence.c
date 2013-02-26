#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header.h"

#define NBINLOOKUP2 10000000
double *xm1,*xf1,*dx1,*mfunc1,*mfy1,*xm2;
int inx;

double N_cen_hiden(double m);
double func_ng6(double m);

void dd_hod_functions(float *hmass, float *hdensity, int nhalo)
{
  static int flag=0;

  int i,j,k,nh,i1,i2,ibin,*nbin,*fbin;
  float *den,*mbin,mass,density,ng1,ng2;
  FILE *fp,*fp2;
  char a[1000];

  double func_ng3(),func_ng4();

  /* Variables for the log bins.
   */
  float dlogm,dmax=5.0e16,dmin=1.0e9,*dupp,binfac;
  int idmax=1000,kbin,*binlookup;

  if(flag)goto SKIP1;
  flag=1;

  /*********************** 
   *initializing the logarithmic bins 
   */  
  idmax = 100;

  dlogm=log(dmax/dmin)/(float)(idmax-1) ;
  dupp=malloc(idmax*sizeof(float));
  binlookup=(int *)calloc(NBINLOOKUP2+2,sizeof(int)) ;
  ibin=0 ;
  for (i=0;i<=NBINLOOKUP2;i++)  {
    mass=dmax*i/NBINLOOKUP2 ;
    if (mass>0)  {
      kbin=(int)floor(log(mass/dmin)/dlogm+1.0) ;
    }
    else {
      kbin=0 ;
    }
    if (kbin<0) kbin=0 ;
    if (kbin>ibin)  {
      dupp[ibin]=mass ;
      ibin=kbin ;
    }
    binlookup[i]=kbin ;
  }
  binlookup[NBINLOOKUP2+1]=idmax ;
  binfac=NBINLOOKUP2/dmax ;
  dupp[idmax-1]=dmax;

  nbin=calloc(idmax,sizeof(int));
  mbin=calloc(idmax,sizeof(double));
  fbin=calloc(idmax,sizeof(int));

  for(i=1;i<=nhalo;++i)
    {
      mass=hmass[i];
      j=binlookup[(int)(mass*binfac)];
      mbin[j]+=mass;
      nbin[j]++;
      if(GAO_EFFECT) {
	if(hdensity[i]>DENSITY_THRESHOLD)
	  fbin[j]++;
      } else {
	if(hdensity[i]<DENSITY_THRESHOLD)
	  fbin[j]++;
      }
    }

  for(j=i=0;i<idmax;++i)
    if(nbin[i])j++;

  xm1 = dvector(1,j);
  xm2 = dvector(1,j);
  xf1 = dvector(1,j);
  dx1 = dvector(1,j);
  mfunc1=dvector(1,j);
  mfy1=dvector(1,j);
  inx = j;

  for(j=i=0;i<idmax;++i)
    if(nbin[i])
      {
	xm1[++j] = mbin[i]/nbin[i];
	xf1[j]   = (double)fbin[i]/nbin[i];
	mfunc1[j] = log(nbin[i]/pow(BOX_SIZE,3.0)/(dlogm*xm1[j]));
	xm2[j] = log(xm1[j]);
	printf("DENFRAC %e %e %e %e %d %e\n",
	       xm1[j],xf1[j],mfunc1[j],dndM_interp(xm1[j]),nbin[i],dlogm);
	fflush(stdout);
      }
  spline(xm1,xf1,inx,1.0E+30,1.0E+30,dx1);
  spline(xm2,mfunc1,inx,1.0E+30,1.0E+30,mfy1);

 SKIP1:

  if(wp_color.ON)
    return;

  fprintf(stdout,"Old M_min: %e\n",HOD.M_min);
  fflush(stdout);
  if(GAO_EFFECT)
    {
      HOD.M_min_loden = exp(zbrent(func_ng4,log(HOD.M_min/10.0),log(HOD.M_min*10),1.0E-5));
      HOD.M_min = HOD.M_min_loden;
    }
  else
    {
      HOD.M_min_hiden = exp(zbrent(func_ng4,log(HOD.M_min/10.0),log(HOD.M_min*10),1.0E-5));
      HOD.M_min_loden = HOD.M_min_hiden*HOD.M_min_fac;
      HOD.M_min = HOD.M_min_hiden;
    }
  fprintf(stdout,"New M_min: %e %e %e\n",HOD.M_min,HOD.M_min_loden,HOD.M_min_fac);
  fflush(stdout);
}

double func_ng3(double m)
{
  double x,dn;
  splint(xm2,mfunc1,mfy1,inx,m,&dn);
  dn=exp(dn);
  m=exp(m);
  //  dn = dndM_interp(m);
  splint(xm1,xf1,dx1,inx,m,&x);
  if(x>1)x=1;
  if(x<0)x=0;
  //  x=x*dn*(N_sat(m) + N_cen_hiden(m))*m; <-- this is for GAO EFFECT
  x=x*dn*N_avg(m)*m;
  return(x);
}

double func_ng5(double m)
{
  double x,dn;
  splint(xm2,mfunc1,mfy1,inx,m,&dn);
  dn=exp(dn);
  m=exp(m);
  // dn = dndM_interp(m);
  splint(xm1,xf1,dx1,inx,m,&x);
  if(x>1)x=1;
  if(x<0)x=0;
  x=(1-x)*dn*N_avg(m)*m;
  return(x);
}

double func_ng4(double m)
{
  double mlo,x1,x2,func_ng5(),func_mlow();

  HOD.M_min=exp(m);
  if(SOFT_CENTRAL_CUTOFF)
    mlo = (zbrent(func_mlow,log(HOD.M_min*1.0E-4),log(HOD.M_min),1.0E-5));
  else
    mlo = m;
  HOD.M_low=exp(mlo);
  x1 = qromo(func_ng5,mlo,log(HOD.M_max),midpnt);
  HOD.M_min=HOD.M_min*HOD.M_min_fac;
  x2 = qromo(func_ng3,mlo,log(HOD.M_max),midpnt);
  //  printf("%e %e %e\n",x1,x2,GALAXY_DENSITY);
  return((x1+x2)-GALAXY_DENSITY);
}


double density_fraction(double m)
{
  double x;
  splint(xm1,xf1,dx1,inx,m,&x);
  return(x);
}

double func_ng_hiden(double m)
{
  return(0);
}
double func_ng_loden(double m)
{
  return(0);
}

/* Color-dependent models.
 */
double dd_func_red_fraction(double x)
{
  double n;
  HOD.fblue0_cen=pow(10.0,x);
  //HOD.fblue0_cen = 0.67;
  n = qtrap(func_ng6,log(HOD.M_low),log(HOD.M_max),1.0E-5);
  printf("%f %e %e %e %e %e\n",HOD.fblue0_cen,n,wp_color.ngal_red,N_avg(8.0e11),HOD.M_min,HOD.M_low); 
  //  exit(0);
  return n - wp_color.ngal_red;
}

double func_ng6(double m)
{
  double x,x1,dn,xt;
  splint(xm2,mfunc1,mfy1,inx,m,&dn);
  dn=exp(dn);
  m=exp(m);
  //dn = dndM_interp(m);
  splint(xm1,xf1,dx1,inx,m,&x);
  if(x>1)x=1;
  if(x<0)x=0;
  //  printf("%e %e\n",m,x);
  x1=(1-x)*dn*N_avg(m)*m; // 1-x is fraction above critical density.
  xt = HOD.fblue0_cen;
  HOD.fblue0_cen = (1. - HOD.M_min_fac*(1.0 - HOD.fblue0_cen));
  x1+=x*dn*N_avg(m)*m;
  HOD.fblue0_cen = xt;
  return(x1);
}
