#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "header.h"

#define NBINLOOKUP 10000

/* Loocal functions.
 */
void nbody_covar(double *rad, double *xi, int nr);

/* External functions.
 */
void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart,int ip);
void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);

double nbody_xi(double r)
{
  static int flag=1,n=30;
  static double *x,*y,*y2;
  double a;
  int i;

  if(flag || RESET_FLAG_1H)
    {
      if(flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=0;
      if(OUTPUT)
	printf("Calculating nbody xi(r)...\n");
      nbody_covar(x,y,n);
      if(OUTPUT)
	printf("finished\n");
      for(i=1;i<=n;++i)
	y[i]=log(y[i]);
      spline(x,y,n,1.0E+30,1.0E+30,y2);
      RESET_FLAG_1H=0;
    }
  r=log(r);
  splint(x,y,y2,n,r,&a);
  return(exp(a));
}


void nbody_covar(double *rad, double *xi, int nr)
{
  float rmax,rmin,lrstep,binfac,rcube,weight0,fac,rlow,density,weightrandom,vol,r;
  int ibin,kbin,nbin,i,j,k,*binlookup,*npair,ngal,ngal_temp;

  float *rupp,*rsqr,reduction_factor,*vxt,*vyt,*vzt;
  int *meshparts, ***meshstart,nmesh,meshfac,nbrmax,*indx;
  double *rbar,galden;

  static float *xg,*yg,*zg,*xt,*yt,*zt;

  static int flag = 1;

  rmin=0.1;
  rmax=30.0;
  nbin=nr;
  rcube=BOX_SIZE;

  if(flag) {
    galden = GALAXY_DENSITY;
    if(HOD.color == 1) galden = wp_color.ngal_blue;
    if(HOD.color == 2) galden = wp_color.ngal_red;
    ngal = galden*rcube*rcube*rcube*2;
    reduction_factor = 1;
    if(galden > 0.005)reduction_factor = 0.2;
    if(wp_color.ON)reduction_factor = 0.5;
    xg = malloc(reduction_factor*ngal*sizeof(float));
    yg = malloc(reduction_factor*ngal*sizeof(float));
    zg = malloc(reduction_factor*ngal*sizeof(float));

    xt = malloc(ngal*sizeof(float));
    yt = malloc(ngal*sizeof(float));
    zt = malloc(ngal*sizeof(float));
    flag = 1;
  }
  
  printf("Starting population...\n");
  ngal = internal_populate_simulation(xt,yt,zt,1.0,0,vxt,vyt,vzt);
  printf("Done. [%d] galaxies [%.0f]\n",ngal,galden*rcube*rcube*rcube);
  fflush(stdout);

  ngal_temp=0;
  for(i=0;i<ngal;++i)
    {
      if(drand48()>reduction_factor)continue;
      xg[ngal_temp] = xt[i];
      yg[ngal_temp] = yt[i];
      zg[ngal_temp] = zt[i];
      ngal_temp++;
    }
  ngal = ngal_temp;

  indx=malloc(ngal*sizeof(int));
  rsqr=malloc(ngal*sizeof(float));

  /*********************** 
   * initializing the logarithmic bins 
   */  
  rupp = (float *) calloc(nbin+1,sizeof(float)) ;
  lrstep=log(rmax/rmin)/(float)(nbin-1) ;
  binlookup=(int *)calloc(NBINLOOKUP+2,sizeof(int)) ;
  ibin=0 ;
  for (i=0;i<=NBINLOOKUP;i++)  {
    r=rmax*i/NBINLOOKUP ;
    if (r>0)  {
      kbin=(int)floor(log(r/rmin)/lrstep+1.0) ;
    }
    else {
      kbin=0 ;
    }
    if (kbin<0) kbin=0 ;
    if (kbin>ibin)  {
      rupp[ibin]=r ;
      ibin=kbin ;
    }
    binlookup[i]=kbin ;
  }
  binlookup[NBINLOOKUP+1]=nbin ;
  rupp[nbin-1]=rmax ;
  rupp[nbin]=rmax ;
  binfac=NBINLOOKUP/rmax ;

  rbar=calloc(nbin,sizeof(double));
  npair=calloc(nbin,sizeof(int));

  nmesh=0;
  meshlink2(ngal,&nmesh,0.0,rcube,rmax,xg,yg,zg,&meshparts,&meshstart,meshfac);

  for(i=0;i<ngal;++i)
    {

      nbrmax=ngal;
      nbrsfind2(0.0,rcube,rmax,nmesh,xg[i],yg[i],zg[i],&nbrmax,indx,rsqr,xg,yg,zg,
		meshparts,meshstart,i);
      for(j=0;j<nbrmax;++j)
	{
	  r=sqrt(rsqr[j]) ;
	  kbin=binlookup[(int)(binfac*r)] ;
	  if(kbin>=nbin)continue;
	  npair[kbin]++;
	  rbar[kbin]+=r;
	}
    }

  density=ngal/(rcube*rcube*rcube) ;

  rlow=0;
  for (kbin=0;kbin<nbin;kbin++)  {
    weight0=(float) npair[kbin];
    
    if (weight0>0.0)  {
      fac=1./weight0 ;
      rbar[kbin] *= fac ;
    }
    else  {					/* avoid errors in empty bins */
      rbar[kbin]=(rupp[kbin]+rlow)*0.5;
    }
    
    /* compute xi, dividing summed weight by that expected for a random set */
    vol=4.*PI/3.*(rupp[kbin]*rupp[kbin]*rupp[kbin]-rlow*rlow*rlow) ;
    weightrandom=ngal*density*vol ;
    
    rad[kbin+1]=log(rbar[kbin]);
    xi[kbin+1]=weight0/weightrandom-1 ;
    rlow=rupp[kbin] ;    

    if(OUTPUT>1)
      {
	fprintf(stdout,"nbody_xi> %f %f %.0f\n",rbar[kbin],xi[kbin+1],weight0);
	fflush(stdout);
      }
  }


  free(indx);
  free(rsqr);
  free(rupp);
  free(binlookup);
  free(npair);
  free(rbar);

  free_i3tensor(meshstart,0,nmesh-1,0,nmesh-1,0,nmesh-1);
} 

int internal_populate_simulation(float *x, float *y, float *z, float reduction_factor, int ivel, 
				 float *vx, float *vy, float *vz)
{
  FILE *fp,*fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9];
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[1000],halocnt[1000],imag;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,high_density_f0,low_density_f0;
  char aa[1000],fname[100];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM = -445, IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1,ii,*indx;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax,temp1,temp2,fac,sigma;
  float **haloarr,**temp_pos,*temp_den,**temp_vel;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0;
  

  static int 
    flag = 1,
    nhalo,
    SO_FILE = 0;

  FILE *fp3;

  static float
    **halo_pos,
    **halo_vel,
    *hdensity,
    *halo_mass;


  if(flag){
    
    fp = openfile(Files.HaloFile);
    nhalo = filesize(fp);

    if(DENSITY_DEPENDENCE)
      {
	fp2 = openfile(Files.HaloDensityFile);
	if(nhalo != filesize(fp2))
	  {
	    fprintf(stderr,"ERROR: filesize mismatch with [%s] and [%s]\n",
		    Files.HaloFile, Files.HaloDensityFile);
	    exit(0);
	  }
      }

    halo_pos = matrix(1,nhalo,0,2);
    halo_mass = vector(1,nhalo);
    temp_pos = matrix(1,nhalo,0,2);
    halo_vel = matrix(1,nhalo,0,2);
    temp_vel = matrix(1,nhalo,0,2);
    indx = ivector(1,nhalo);

    fprintf(stderr,"HERE %d halos\n",nhalo);

    if(DENSITY_DEPENDENCE)
      {
	temp_den = vector(1,nhalo);
	hdensity = vector(1,nhalo);
      }

    for(i=1;i<=nhalo;++i)
      {
	if(DENSITY_DEPENDENCE)
	  fscanf(fp2,"%f",&temp_den[i]);
	if(SO_FILE)
	  {
	    fscanf(fp,"%d %lf %f %f %f %f %f %f %f %f",
		   &j,&mass,&x1,&x1,&xh[0],&xh[1],&xh[2],&vh[0],&vh[1],&vh[2]);
	    halo_mass[i] = -mass;
	  }
	else
	  {
	    fscanf(fp,"%d %d %e %e %e %e %e %e %e",
		   &j,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
	    mass=imass*RHO_CRIT*OMEGA_M*pow(RESOLUTION,3.0);
	    halo_mass[i]= - imass*RHO_CRIT*OMEGA_M*pow(RESOLUTION,3.0);
	  }
	for(j=0;j<3;++j) {
	  temp_pos[i][j] = xh[j];
	  temp_vel[i][j] = vh[j]; }
	indx[i] = i;
      }
    sort2(nhalo,halo_mass,indx);

    for(i=1;i<=nhalo;++i)
      {
	halo_mass[i] *= -1;
	for(j=0;j<3;++j) {
	  halo_pos[i][j] = temp_pos[indx[i]][j];
	  halo_vel[i][j] = temp_vel[indx[i]][j]; }
	if(DENSITY_DEPENDENCE)
	  hdensity[i] = temp_den[indx[i]];
      }
    free_matrix(temp_pos,1,nhalo,0,2);
    free_matrix(temp_vel,1,nhalo,0,2);
    free_ivector(indx,1,nhalo);
    flag = 0;

    fclose(fp);
    if(DENSITY_DEPENDENCE)
      fclose(fp2);
  }

  if(wp_color.ON)
    {
      temp1 = HOD.M_min;
      temp2 = HOD.M_low;
    }
  HOD.M_min = 0;
  set_HOD_params();
  if(DENSITY_DEPENDENCE)
    dd_hod_functions(halo_mass,hdensity,nhalo);
  if(wp_color.ON)
    {
      HOD.M_min = temp1;
      HOD.M_low = temp2;
      fprintf(stderr,"true M_min = %e\n",HOD.M_min);
    }
  HOD.M_low = set_low_mass();
  mlo = HOD.M_low;
  if(HOD.M_min_fac<0)
    mlo *= HOD.M_min_fac;

  high_density_f0 = HOD.fblue0_cen;
  low_density_f0 = (1 - HOD.M_min_fac*(1-HOD.fblue0_cen));

  fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;

  fprintf(stderr,"%e %e %e\n",HOD.M_min,HOD.M_low,halo_mass[1]);

  sprintf(fname,"%s.mock_halo",Task.root_filename);
  fp2 = fopen(fname,"w");

  for(ii=1;ii<=nhalo;++ii)
    {
      //      printf("%d\n",ii);
      //fflush(stdout);
      mass=halo_mass[ii];
      if(mass<mlo)break;
      if(ran2(&IDUM)>reduction_factor)continue;

      for(i=0;i<3;++i)
	{
	  xh[i] = halo_pos[ii][i];
	  if(xh[i]<0)xh[i]+=BOX_SIZE;
	  if(xh[i]>BOX_SIZE)xh[i]-=BOX_SIZE;
	  if(ivel)
	    vh[i] = halo_vel[ii][i];
	}

      nsat = N_sat(mass);
      if(nsat>250)
	n1 = gasdev(&IDUM3)*sqrt(nsat) + nsat;
      else
	n1 = poisson_deviate(nsat);      
      
      for(i=1;i<=n1;++i)
	{
	  r = NFW_position(mass,xg);
	  for(k=0;k<3;++k)
	    {
	      xg[k]+=xh[k];
	      if(xg[k]<0)xg[k]+=BOX_SIZE;
	      if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
	    }
	  if(ivel) {
	    NFW_velocity(mass,vg,mag);
	    vx[id2] = vg[0]+vh[0];
	    vy[id2] = vg[1]+vh[1];
	    vz[id2] = vg[2]+vh[2];
	  }
	  x[id2] = xg[0];
	  y[id2] = xg[1];
	  z[id2] = xg[2];
	  id2++;
	  fprintf(fp2,"%d\n",ii);

	}
      if(DENSITY_DEPENDENCE && !wp_color.ON)
	{
	  //	  printf("%e %e\n",mass,hdensity[ii]);
	  HOD.M_min=HOD.M_min_hiden;
	  if(hdensity[ii]<DENSITY_THRESHOLD)
	    HOD.M_min=HOD.M_min*HOD.M_min_fac;
	}
      if(DENSITY_DEPENDENCE && wp_color.ON && HOD.color == 2)
	{
	  if(hdensity[ii]<DENSITY_THRESHOLD)
	    HOD.fblue0_cen = low_density_f0;
	  else
	    HOD.fblue0_cen = high_density_f0;
	}
      if(DENSITY_DEPENDENCE && wp_color.ON && HOD.color == 1)
	{
	  if(hdensity[ii]<DENSITY_THRESHOLD)
	    HOD.fblue0_cen = low_density_f0;
	  else
	    HOD.fblue0_cen = high_density_f0;
	}

      ncen=N_cen(mass);
      if(ran2(&IDUM)>ncen)continue;
      sigma=fac*pow(mass,1.0/3.0)/sqrt(2.0);
      x[id2] = xh[0];
      y[id2] = xh[1];
      z[id2] = xh[2];
      if(ivel) {
	vx[id2] = vh[0] + gasdev(&IDUM)*VBIAS_C*sigma;
	vy[id2] = vh[1] + gasdev(&IDUM)*VBIAS_C*sigma;
	vz[id2] = vh[2] + gasdev(&IDUM)*VBIAS_C*sigma;
      }
      id2++;
      fprintf(fp2,"%d\n",ii);
    }
  fclose(fp2);

  muh(999);
  //  ivel = 1;
  if(ivel) {
    sprintf(aa,"%s.mock",Task.root_filename);      
    fp = fopen(aa,"w");
    for(i=0;i<id2;++i)
      fprintf(fp,"%e %e %e %e %e %e\n",x[i],y[i],z[i],vx[i],vy[i],vz[i]);
      //  fprintf(fp,"%e %e %e\n",x[i],y[i],z[i]);
    fclose(fp);
  }

  return id2;

}
