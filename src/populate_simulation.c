#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "header.h"

#define NBINLOOKUP 10000

/* This function reads in a halo file and creates a mock galaxy distribution 
 * by populating the dark matter halos with galaxies with the currently specified
 * HOD functions.
 *
 * The format of the halo file needs to be in ASCII:
 * col 1 - halo ID (not used)
 * col 2 - number of particles in halo
 * cols 3-5 - x,y,z position of halo, in Mpc/h
 * col 6 - space velocity of halo (not used)
 * cols 7-9 - vx,vy,vz velocity of halo, in km/s
 *
 * The values of RESOLUTION, BOX_SIZE, OMEGA_M from the batch file are used
 * to convert particle number into mass, and to box-wrap galaxies.
 * If you want another file format, by all means edit.
 *
 * Output: mock galaxies are printed to [root].mock
 * in format: x,y,z [Mpc/h] vz,vy,vz [km/s]
 *
 * Output: HOD calculated from populated halos (for cross-checking with true HOD)
 * [root].binned_HOD in format: bin id, log10(mass), <N_gal>, N_halo
 * 
 * Satellite positions are chosen from a random sampling of the NFW profile
 * for the mass of the halo being populated. If CVIR_FAC is specified, then that value
 * will be used to adjust the NFW profile. Satellite velocities are chosen
 * from a Gaussin distribution with width = to virial dispersion of halo (plus
 * halo center-of-mass).
 *
 * NB-- If you are using the code to calculate M_min from a desired space density,
 * be sure that the linear matter power spectrum is the same as that used in the
 * simulation, or else the space densities won't match. [Mass function is used for this
 * purpose].
 *
 */
double NFW_central_velocity(double mass, double v[], double mag);
void calc_nbody_two_halo(float **gal, int *id, int ngal);

double *g21_rad, *g21_xi;
int g21_nrad;

/* External functions.
 */
void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart,int ip);
void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);


void populate_simulation()
{
  FILE *fp,*fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9],*fpt,*fpxtr;
  int j,k,n,n1,j_start=0,i1,galcnt[1000],halocnt[1000],imag;
  int i, imass;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax;
  float **haloarr;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0,haloid;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;

  int SO_FILE = 0,
    JEANS_DISPERSION = 0;
    
  //int SO_FILE = 1,
  //  JEANS_DISPERSION = 0;
    
  //  if(RESOLUTION > 1.6)
  // WARREN_MASS_CORRECTION = 1;

  TRACK_GALAXIES=1;

  //  srand48(IDUM_MCMC);

  if(DENSITY_DEPENDENCE)
    {
      ngal = GALAXY_DENSITY*BOX_SIZE*BOX_SIZE*BOX_SIZE*1.4;
      xt = vector(0,ngal);
      yt = vector(0,ngal);
      zt = vector(0,ngal);
      vxt = vector(0,ngal);
      vyt = vector(0,ngal);
      vzt = vector(0,ngal);
      internal_populate_simulation(xt,yt,zt,1.0,1,vxt,vyt,vzt);
      return;
    }
  galarr = matrix(1,1000000,0,5);
  haloarr = matrix(1,1000000,0,6);
  galid = ivector(1,1000000);

  fp=openfile(Files.HaloFile);
  sprintf(aa,"%s.mock",Task.root_filename);      
  fp2 = fopen(aa,"w");
  sprintf(aa,"%s.mock_halo",Task.root_filename);      
  fpt = fopen(aa,"w");
  sprintf(aa,"%s.mock_extra",Task.root_filename);
  fpxtr = fopen(aa,"w");

  for(i=0;i<1000;++i)
    halocnt[i]=0;
  for(i=0;i<1000;++i)
    galcnt[i]=0;

  set_HOD_params();
  mlo = HOD.M_low;
  printf("MLO %e %e %f\n",mlo,HOD.M_min,HOD.sigma_logM);
  printf("BOX_SIZE %f\n",BOX_SIZE);
  fflush(stdout);

int nlines = filesize(fp);
fprintf(stderr,"Number of lines: %d\n",nlines);                      

    float extra;
  
  //  int ii=0;
  //  for(ii<nlines; ii++;)
    while(!feof(fp))
    {
      if(SO_FILE)
	{
	  fscanf(fp,"%d %lf %f %f %f %f %f %f %f %f",
		 &i,&mass,&x1,&x1,&xh[0],&xh[1],&xh[2],&vh[0],&vh[1],&vh[2]);
	  fgets(aa,1000,fp);
	  /*
	  fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
	  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
	  printf("%e %e %f\n",mass,sigv,VBIAS);
	  if(i==2)
	    one_halo(1,1);	    
	  */
	  haloid = i;
	}
      else
	{
	  fscanf(fp,"%d %d %e %e %e %e %e %e %lf %e",
		 &i,&imass,&xh[0],&xh[1],&xh[2],&vh[0],&vh[1],&vh[2],&mass,&extra);
	  //mass=imass*RHO_CRIT*OMEGA_M*pow(RESOLUTION,3.0);
      //fprintf(stderr, "%e %e %e\n", xh[0], vh[2], mass);
	}
      haloid = i;
      if(mass>HOD.M_max)continue;

      if(WARREN_MASS_CORRECTION)
	mass = imass*(1-pow(imass,-0.6))*RHO_CRIT*OMEGA_M*pow(RESOLUTION,3.0);


      if(mass<mlo)continue;

      for(i=0;i<3;++i)
	{
	  if(xh[i]<0)xh[i]+=BOX_SIZE;
	  if(xh[i]>BOX_SIZE)xh[i]-=BOX_SIZE;
	}

      i1 = (int)(log10(mass)/0.1);
      halocnt[i1]++;	  
      ncen=N_cen(mass);
      if(drand48()>ncen)goto SATELLITES;
      if(VBIAS_C>0)
	{
	  NFW_central_velocity(mass,vg,mag);
	  for(i=0;i<3;++i)
	    vh[i]+=vg[i];
	}
      fprintf(fp2,"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
      fprintf(fpt,"%d\n",haloid);
      fprintf(fpxtr,"%e -1\n",mass);
        
      // THis mocks up the ~22 columns of the MR mocks
      //      fprintf(fpt,"%d %d %d %e %e %e 0.0 0.0 0.0 %e %e %e 0.0 0.0 0.0 0.0 0.0 0.0 %e 0.0 0.0\n",
      //	      0,0,0,xh[0],xh[1],xh[2],vh[0],vh[1],vh[2],log10(mass));

      if(VBIAS_C>0)
	{
	  for(i=0;i<3;++i)
	    vh[i]-=vg[i];
	}
      
      galcnt[i1]++;

      if(TRACK_GALAXIES)
	{
	  id2++;
	  galid[id2] = haloid;
	  galarr[id2][0] = xh[0];
	  galarr[id2][1] = xh[1];
	  galarr[id2][2] = xh[2];
	  galarr[id2][3] = vh[0];
	  galarr[id2][4] = vh[1];
	  galarr[id2][5] = vh[2];
	}

      if(TRACK_GALAXIES)
	{
	  id1++;
	  haloarr[id1][0] = xh[0];
	  haloarr[id1][1] = xh[1];
	  haloarr[id1][2] = xh[2];
	  haloarr[id1][3] = vh[0];
	  haloarr[id1][4] = vh[1];
	  haloarr[id1][5] = vh[2];
	  haloarr[id1][6] = mass;
	}

    SATELLITES:
      nsat = N_sat(mass);
      if(nsat>250)
	n1 = gasdev(&IDUM3)*sqrt(nsat) + nsat;
      else
	n1 = poisson_deviate(nsat);      
      /*
      if(nsat>1)
	fprintf(stdout,"BUH %d %f %e %d %e %f\n",haloid,nsat,mass,n1,HOD.M1,pow(mass/HOD.M1,HOD.alpha));
      */
      
      for(i=1;i<=n1;++i)
	{
	  r = NFW_position(mass,xg);
	  if(JEANS_DISPERSION)
	    jeans_dispersion(mass,r,vg);
	  else
	    NFW_velocity(mass,vg,mag);
	  for(k=0;k<3;++k)
	    {
	      xg[k]+=xh[k];
	      if(xg[k]<0)xg[k]+=BOX_SIZE;
	      if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
	      vg[k]+=vh[k];
	      /* vg[k] = vgf[k]; */
	      /*
		vg[k] = gasdev(&IDUM3)*500;
		vg[k] = non_gaussian_velocity();
	      */
	    }	
	  fprintf(fp2,"%e %e %e %e %e %e\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2]);
      fprintf(fpt,"%d\n",haloid);
      fprintf(fpxtr,"%e 1\n",mass);

	  //	  fprintf(fpt,"%d %d %d %e %e %e 0.0 0.0 0.0 %e %e %e 0.0 0.0 0.0 0.0 0.0 0.0 %e 0.0 0.0\n",1,1,1,xg[0],xg[1],xg[2],vg[0],vg[1],vg[2],log10(mass));

	  if(TRACK_GALAXIES)
	    {
	      id2++;
	      galid[id2] = haloid;
	      galarr[id2][0] = xg[0];
	      galarr[id2][1] = xg[1];
	      galarr[id2][2] = xg[2];
	      galarr[id2][3] = vg[0];
	      galarr[id2][4] = vg[1];
	      galarr[id2][5] = vg[2];
	    }

	  /* Bin up the galaxies by halo mass to check the HOD
	   */
	  galcnt[i1]++;	  
	}
      
      if(feof(fp))break;
    }
  fclose(fp2);
  fclose(fpt);
  fclose(fpxtr);

  /* output the binned HOD
   */
  sprintf(aa,"%s.binned_HOD",Task.root_filename);      
  fp2=fopen(aa,"w");
  for(i=0;i<1000;++i)
    if(galcnt[i]>0)
      fprintf(fp2,"%d %f %f %d %d\n",
	      i,(i+0.5)*0.1,(float)galcnt[i]/halocnt[i],galcnt[i],halocnt[i]);
  fclose(fp2);

  /* Calculate the two-halo term
   */
  if(TRACK_GALAXIES)
    {
      fprintf(stderr,"Calling nbody_2halo...\n");
      //calc_nbody_two_halo(galarr,galid,id2);
    }
  return ;
  
  /* Track the close pairs - for Non-Tinkers, don't worry about this.
   */
  rmin = 1.5;
  rmax = 2.5;
  for(i=1;i<=id2;++i)
    for(j=i+1;j<=id2;++j)
      {
	if(galid[i]==galid[j])continue;
	dx = galarr[i][0] - galarr[j][0];
	if(dx>BOX_SIZE/2)dx = BOX_SIZE-dx;
	if(dx>rmax)continue;
	dy = galarr[i][1] - galarr[j][1];
	if(dy>BOX_SIZE/2)dy = BOX_SIZE-dy;
	if(dy>rmax)continue;
	dz = galarr[i][2] - galarr[j][2];
	if(dz>BOX_SIZE/2)dz = BOX_SIZE-dz;
	if(dz>rmax)continue;
	dr = sqrt(dx*dx+dy*dy+dz*dz);
	if(dr<rmin || dr>rmax)continue;

	i1 = galid[i];
	j1 = galid[j];
	dx = haloarr[i1][0] - haloarr[j1][0];
	if(dx>BOX_SIZE/2)dx = BOX_SIZE-dx;
	dy = haloarr[i1][1] - haloarr[j1][1];
	if(dy>BOX_SIZE/2)dy = BOX_SIZE-dy;
	dz = haloarr[i1][2] - haloarr[j1][2];
	if(dz>BOX_SIZE/2)dz = BOX_SIZE-dz;
	drh = sqrt(dx*dx+dy*dy+dz*dz);
	
	rv1 = pow(3*haloarr[i1][6]/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
	rv2 = pow(3*haloarr[j1][6]/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
	if(dr < rv1+rv2)
	  {
	    dx = galarr[i][3] - galarr[j][3];
	    dy = galarr[i][4] - galarr[j][4];
	    dz = galarr[i][5] - galarr[j][5];
	    printf("%e %e %e %e %e %e %e\n",dr,rv1+rv2,haloarr[i1][6],haloarr[j1][6],dx,dy,dz);
	    fflush(stdout);
	  }
      }

  exit(0);
}

/* Generate a random integer based on a Poisson distribution 
 * with mean given as input.
 */
int poisson_deviate(double nave)
{
  static int flag=0;
  double p,pp;
  int n;

  p=0;
  pp=1;

  while(p<pp)
    {
      if(nave<1)
	n=(int)(drand48()*20);
      else
	n=(int)(drand48()*30*nave);
      p=poisson_prob(n,nave);
      pp=drand48();
    }
  return(n);
}

/* Poisson probability of n given n_average
 */
double poisson_prob(int n, double nave)
{
  int i;
  double fac=1;

  if(n>0)
    for(i=1;i<=n;++i)
      fac*=nave/i;

  return((float)(fac*exp(-nave)));
}

/* Randomy generates a position away from the origin with 
 * a probability given by the NFW profile for a halo of the input
 * mass (and including the CVIR_FAC)
 */
double NFW_position(double mass, double x[])
{
  double r,pr,max_p,costheta,sintheta,phi1,signs,rvir,rs,cvir;
  
  cvir=halo_concentration(mass)*CVIR_FAC;
  rvir=pow(3*mass/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rs=rvir/cvir;
  max_p=NFW_density(rs,rs,1.0)*rs*rs*4.0*PI;

  for(;;) {
    r=drand48()*rvir;
    pr=NFW_density(r,rs,1.0)*r*r*4.0*PI/max_p;
    
    if(drand48()<=pr)
      {
	costheta=2.*(drand48()-.5);
	sintheta=sqrt(1.-costheta*costheta);
	signs=2.*(drand48()-.5);
	costheta=signs*costheta/fabs(signs);
	phi1=2.0*PI*drand48();
	
	x[0]=r*sintheta*cos(phi1);
	x[1]=r*sintheta*sin(phi1);
	x[2]=r*costheta;
	return r;
      }
  }
}

/* This is the NFW density profile
 */
double NFW_density(double r, double rs, double ps)
{
  return(ps*rs/(r*(1+r/rs)*(1+r/rs)));
}

/* This sets the velocity to be isotropic Gaussian.
 */
double NFW_velocity(double mass, double v[], double mag)
{
  static long IDUM2=-455;
  static double fac = -1;
  double sigv,vbias=1;
  int i;

  if(fac<0)
      fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS;
  return(0);
}

/* This sets the velocity to be isotropic Gaussian.
 */
double NFW_central_velocity(double mass, double v[], double mag)
{
  static long IDUM2=-455;
  static double fac = -1;
  double sigv,vbias=1;
  int i;

  if(fac<0)
      fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS_C;
  return(0);
}


/* This is the routine where we tabulate the 2-halo term and
 * interpolate as needed.
 */
double nbody_two_halo(double r)
{
  static int flag=1,n=30;
  static double *x,*y,*y2;
  double a;
  int i;

  if(flag || RESET_FLAG_1H)
    {
      populate_simulation();
      if(flag)
	{
	  n = g21_nrad;
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=0;
      for(i=1;i<=n;++i)
	{
	  x[i] = g21_rad[i];
	  y[i] = g21_xi[i];
	}
      /*
      if(OUTPUT)
	printf("Calculating nbody xi(r)...\n");
      calc_nbody_two_halo(x,y,n);
      if(OUTPUT)
	printf("finished\n");
      for(i=1;i<=n;++i)
	y[i]=log(y[i]);
      */
      spline(x,y,n,1.0E+30,1.0E+30,y2);
      RESET_FLAG_1H=0;
    }
  //  r=log(r);
  if(r<0.2)return(-1);
  splint(x,y,y2,n,r,&a);
  return(a);
  //  return(exp(a));
}

/* This is the two-halo term.
 */
void calc_nbody_two_halo(float **gal, int *id, int ngal)
{
  float rmax,rmin,lrstep,binfac,rcube,weight0,fac,rlow,density,weightrandom,vol,r;
  int ibin,kbin,nbin,i,j,k,*binlookup,*npair,ngal_temp;

  float *rupp,*rsqr,reduction_factor,*vxt,*vyt,*vzt;
  int *meshparts, ***meshstart,nmesh,meshfac,nbrmax,*indx;
  double *rbar,galden;

  float *xg,*yg,*zg,*xt,*yt,*zt;

  static float *rad, *xi;
  static int nr, flag = 0;

  char fname[100];
  FILE *fp;

  sprintf(fname,"%s.nbody_2halo",Task.root_filename);
  fp = fopen(fname,"w");

  rmax = 30;
  rmin = 0.1;
  nbin = nr = 20;

  rad = vector(1,nr);
  xi = vector(1,nr);
  g21_rad = dvector(1,nr);
  g21_xi = dvector(1,nr);
  g21_nrad = nr;

  xg = malloc(ngal*sizeof(float));
  yg = malloc(ngal*sizeof(float));
  zg = malloc(ngal*sizeof(float));

  ngal_temp=0;
  for(i=1;i<=ngal;++i)
    {
      xg[ngal_temp] = gal[i][0];
      yg[ngal_temp] = gal[i][1];
      zg[ngal_temp] = gal[i][2];
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

  rcube = BOX_SIZE;

  nmesh=0;
  meshlink2(ngal,&nmesh,0.0,rcube,rmax,xg,yg,zg,&meshparts,&meshstart,meshfac);


  for(i=0;i<ngal;++i)
    {

      nbrmax=ngal;
      nbrsfind2(0.0,rcube,rmax,nmesh,xg[i],yg[i],zg[i],&nbrmax,indx,rsqr,xg,yg,zg,
		meshparts,meshstart,i);
      for(j=0;j<nbrmax;++j)
	{
	  k = indx[j];

	  // if gals in same halo, skip
	  if(id[i+1]==id[k+1])continue;

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
    
    g21_rad[kbin+1] = rad[kbin+1]=log(rbar[kbin]);
    g21_xi[kbin+1] = xi[kbin+1]=weight0/weightrandom-1 ;
    rlow=rupp[kbin] ;    

    fprintf(fp,"%e %e %.0f\n",rbar[kbin],xi[kbin+1],weight0);

    if(OUTPUT>1)
      {
	fprintf(stdout,"nbody_xi> %f %f %.0f\n",rbar[kbin],xi[kbin+1],weight0);
	fflush(stdout);
      }
  }
  fclose(fp);

  free(indx);
  free(rsqr);
  free(rupp);
  free(binlookup);
  free(npair);
  free(rbar);

  free_i3tensor(meshstart,0,nmesh-1,0,nmesh-1,0,nmesh-1);
} 


/**************************
 */
void populate_sampled_simulation()
{
  FILE *fp,*fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9],*fpt;
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[1000],halocnt[1000],imag;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax;
  float **haloarr;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0,haloid;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;

  int SO_FILE = 0,
    JEANS_DISPERSION = 0;

  int ioct, prev_haloid, npart, ibuf, ihalo, overall_halo_count;
  float xbuf[6], x0, y0, z0,rhalf, BUFFER = 4;
  int NSAMPLE_MAX=100000;
  float halopart[NSAMPLE_MAX][7];
  int OUTPUT_HALOID=1;
  char fname[1000];
  FILE *fpp;

  rhalf = BOX_SIZE/2;

  TRACK_GALAXIES=0;

  galarr = matrix(1,1000000,0,5);
  haloarr = matrix(1,1000000,0,6);
  galid = ivector(1,1000000);

  fp=openfile(Files.HaloFile);
  sprintf(aa,"%s.mock",Task.root_filename);      
  fp2 = fopen(aa,"w");

  if(OUTPUT_HALOID)
    {
      sprintf(aa,"%s.mock_halo",Task.root_filename);      
      fpt = fopen(aa,"w");
    }

  for(i=0;i<1000;++i)
    halocnt[i]=0;
  for(i=0;i<1000;++i)
    galcnt[i]=0;

  set_HOD_params();
  mlo = HOD.M_low;
  printf("MLO %e %e %f\n",mlo,HOD.M_min,HOD.sigma_logM);
  printf("BOX_SIZE %f\n",BOX_SIZE);
  fflush(stdout);

  /* Open the particle file.
   */
  ioct = 1;
  sprintf(fname,"part.%d",ioct);
  fprintf(stderr,"OPENING [%s]\n",fname);
  fpp = openfile(fname);
  fscanf(fpp,"%d %f %f %f %f %f %f",&ibuf,&xbuf[0],&xbuf[1],
	 &xbuf[2],&xbuf[3],&xbuf[4],&xbuf[5]);
  //  ftread(&ibuf,sizeof(int),1,fpp);
  //ftread(&xbuf,sizeof(float),6,fpp);

  x0 = -BUFFER;
  y0 = -BUFFER;
  z0 = -BUFFER;
  
  if(ioct%2==0)x0+=rhalf;
  if(ioct==3 || ioct==4 || ioct==7 || ioct==8)y0+=rhalf;
  if(ioct>4)z0+=rhalf;
  
  srand48(IDUM_MCMC);

  prev_haloid = 0;
  overall_halo_count=0;
  while(!feof(fp))
    {
      if(SO_FILE)
	{
	  fscanf(fp,"%d %lf %f %f %f %f %f %f %f %f",
		 &i,&mass,&x1,&x1,&xh[0],&xh[1],&xh[2],&vh[0],&vh[1],&vh[2]);
	}
      else
	{
	  fscanf(fp,"%d %d %e %e %e %e %e %e %e",
		 &i,&imass,&xh[0],&xh[1],&xh[2],&x1,&vh[0],&vh[1],&vh[2]);
	  mass=imass*RHO_CRIT*OMEGA_M*pow(RESOLUTION,3.0);
	}
      if(feof(fp))break;

      overall_halo_count++;
      haloid = i;

      /* Check to see if we've moved on to another octant
       */
      if(haloid<prev_haloid)
	{
	  ioct++;
	  if(ioct>8)
	    {
	      fprintf(stderr,"ERROR: we're on octant 9 at halo %d %d\n",
		      haloid,overall_halo_count);
	      exit(0);
	    }

	  x0 = -BUFFER;
	  y0 = -BUFFER;
	  z0 = -BUFFER;
	  
	  if(ioct%2==0)x0+=rhalf;
	  if(ioct==3 || ioct==4 || ioct==7 || ioct==8)y0+=rhalf;
	  if(ioct>4)z0+=rhalf;
	      
	  fclose(fpp);
	  sprintf(fname,"part.%d",ioct);
	  fprintf(stderr,"OPENING [%s] at halo %d %d\n",fname,haloid,overall_halo_count);
	  fpp = openfile(fname);
	  fscanf(fpp,"%d %f %f %f %f %f %f",&ibuf,&xbuf[0],&xbuf[1],
		 &xbuf[2],&xbuf[3],&xbuf[4],&xbuf[5]);
	  //ftread(&ibuf,sizeof(int),1,fpp);
	  //ftread(&xbuf,sizeof(float),6,fpp);
	  ihalo=ibuf;
	}
      prev_haloid = haloid;

      /* Read in all the particles for this halo.
       */
    READ_LOOP:
      ihalo = ibuf;
      npart = 0;
      while(ihalo==ibuf)
	{
	  if(npart<NSAMPLE_MAX)
	    {
	      halopart[npart][0] = xbuf[0]+x0;
	      halopart[npart][1] = xbuf[1]+y0;
	      halopart[npart][2] = xbuf[2]+z0;

	      for(i=3;i<6;++i)
		halopart[npart][i] = xbuf[i];
	      halopart[npart][6] = 1;
	    }
	  npart++;
	  //if(npart==NSAMPLE_MAX)
	  // endrun("too many sampled particle for a halo");
	  fscanf(fpp,"%d %f %f %f %f %f %f",&ibuf,&xbuf[0],&xbuf[1],
		 &xbuf[2],&xbuf[3],&xbuf[4],&xbuf[5]);
	  //ftread(&ibuf,sizeof(int),1,fpp);
	  //ftread(&xbuf,sizeof(float),6,fpp);
	  //fprintf(stderr,"%d %d\n",ihalo,ibuf);
	  if(feof(fpp))break; 
	}
      if(npart>NSAMPLE_MAX)
	{
	  fprintf(stderr,"Warning, %d particles were read for halo %d\n",npart,ihalo);
	  npart = NSAMPLE_MAX;
	}

      if(ihalo<haloid)goto READ_LOOP;
      if(ihalo>haloid)
	{
	  fprintf(stderr,"ERROR: particles at halo %d, file at halo %d\n",ihalo,haloid);
	  exit(0);
	}
      
      if(mass>HOD.M_max)continue;

      for(i=0;i<3;++i)
	{
	  if(xh[i]<0)xh[i]+=BOX_SIZE;
	  if(xh[i]>BOX_SIZE)xh[i]-=BOX_SIZE;
	}

      halocnt[i1]++;	  
      ncen=N_cen(mass);
      if(drand48()>ncen)continue;

      fprintf(fp2,"%e %e %e %e %e %e\n",xh[0],xh[1],xh[2],vh[0],vh[1],vh[2]);
      if(OUTPUT_HALOID)
	fprintf(fpt,"%d\n",haloid);
 
      i1 = (int)(log10(mass)/0.1);
      galcnt[i1]++;

      if(TRACK_GALAXIES)
	{
	  id2++;
	  galid[id2] = haloid;
	  galarr[id2][0] = xh[0];
	  galarr[id2][1] = xh[1];
	  galarr[id2][2] = xh[2];
	  galarr[id2][3] = vh[0];
	  galarr[id2][4] = vh[1];
	  galarr[id2][5] = vh[2];
	}

      if(TRACK_GALAXIES)
	{
	  id1++;
	  haloarr[id1][0] = xh[0];
	  haloarr[id1][1] = xh[1];
	  haloarr[id1][2] = xh[2];
	  haloarr[id1][3] = vh[0];
	  haloarr[id1][4] = vh[1];
	  haloarr[id1][5] = vh[2];
	  haloarr[id1][6] = mass;
	}

      nsat = N_sat(mass);
      if(nsat>250)
	n1 = gasdev(&IDUM3)*sqrt(nsat) + nsat;
      else
	n1 = poisson_deviate(nsat);      

      /* Check to see if we have enough particles.
       */
      if(npart<n1)
	{
	  printf("ERROR: halo %d %e wants %d sats, only %d in file\n",
		 haloid,mass,n1,npart);
	  exit(0);
	}
      
      i=0;
      while(i<n1)
	{
	  j = drand48()*npart;
	  if(halopart[j][6]<0)continue;
	  
	  if(halopart[j][0]>BOX_SIZE)halopart[j][0] -= BOX_SIZE;
	  if(halopart[j][1]>BOX_SIZE)halopart[j][1] -= BOX_SIZE;
	  if(halopart[j][2]>BOX_SIZE)halopart[j][2] -= BOX_SIZE;
	  if(halopart[j][0]<=0)halopart[j][0] += BOX_SIZE;
	  if(halopart[j][1]<=0)halopart[j][1] += BOX_SIZE;
	  if(halopart[j][2]<=0)halopart[j][2] += BOX_SIZE;

	  if(VBIAS!=1)
	    {
	      NFW_velocity(mass,vg,mag);
	      for(k=0;k<3;++k)
		halopart[j][k+3] = vh[k] + vg[k];
	    }
	  fprintf(fp2,"%e %e %e %e %e %e\n",halopart[j][0],halopart[j][1],halopart[j][2],
		  halopart[j][3],halopart[j][4],halopart[j][5]);
	  halopart[j][6] = -1;

	  if(OUTPUT_HALOID)
	    fprintf(fpt,"%d\n",haloid);
	  i++;

	  if(TRACK_GALAXIES)
	    {
	      id2++;
	      galid[id2] = haloid;
	      galarr[id2][0] = xg[0];
	      galarr[id2][1] = xg[1];
	      galarr[id2][2] = xg[2];
	      galarr[id2][3] = vg[0];
	      galarr[id2][4] = vg[1];
	      galarr[id2][5] = vg[2];
	    }

	  /* Bin up the galaxies by halo mass to check the HOD
	   */
	  galcnt[i1]++;	  
	} 
      if(feof(fp))
	break;
    }      
  fclose(fp2);
  fclose(fpt);
  fclose(fpp);

  /* output the binned HOD
   */
  sprintf(aa,"%s.binned_HOD",Task.root_filename);      
  fp2=fopen(aa,"w");
  for(i=0;i<1000;++i)
    if(galcnt[i]>0)
      fprintf(fp2,"%d %f %f %d %d\n",
	      i,(i+0.5)*0.1,(float)galcnt[i]/halocnt[i],galcnt[i],halocnt[i]);
  fclose(fp2);

  /* Calculate the two-halo term
   */
  if(TRACK_GALAXIES)
    {
      fprintf(stderr,"Calling nbody_2halo...\n");
      calc_nbody_two_halo(galarr,galid,id2);
    }
  return ;
  
}
