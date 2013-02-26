
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* This is a simple integrator to calculate the correlation function multipoles.
 * The xi(sigma,pi) will be calculated on a polar grid in log r, phi, and then
 * simple midpoint method to integrate over angle at each r.
 * There will be 16 angular bins and maybe 20 radial bins.
 *
 * I suppose I could interpolate over angle from the 16 calculations, but this is the
 * most direct way of comparing to the n-body simulations, which are calculated exactly
 * this way.
 */

int nrad_g6,nphi_g6,NR,nrad_g6_temp=50;
double *phi_g6,*rr_g6,*xi_mono,**xi_g6,*xi_mono_temp,*rr_g6_temp;
double spherically_averaged_xi(double s);
double spherically_averaged_xi_temp(double s);
double monopole(double mu);
double quadrupole(double mu);

void xi_multipoles()
{
  static int flag=1;
  static double *xi_quad;

  int i,j,k,n,istep=1,istart=1,jstep=1,jstart=1,sdss_ilimit=3;
  double rlo,rhi,rs,rp,delta_phi,phi1,t1h=0,t2h=0,xi,
    x1,x2,x3,x1h,x2h,a,b,t0,t1,*send,*recv,dlogr,xi_bar,drbin,
    philo,phihi,rp_min=0;
  FILE *fp;
  char fname[100],aa[100];
  static int flag1 = 1;

#ifdef PARALLEL
  istep = NTask;
  istart = ThisTask + 1;
  jstart = 1;
  jstep = 1;
#endif

  nrad_g6=50;
  if(Work.chi2)nrad_g6=Work.nrad;

  nphi_g6=36;
  delta_phi=PI/2./(nphi_g6);

  rlo=log(0.1);
  rhi=log(40.0);
  dlogr=(rhi-rlo)/(nrad_g6-1);

  if(Work.SDSS_bins)
    {
      nrad_g6 = Work.nrad;
      rlo = 1.0;
      rhi = 40.0;
      dlogr = 1.0;

      // By-hand stuff
      nrad_g6 = 11;
      nphi_g6 = 18;
    }

  if(flag) {
    flag=0;
    rr_g6=dvector(1,nrad_g6);
    phi_g6=dvector(1,nphi_g6);
    xi_mono=dvector(1,nrad_g6);
    xi_quad=dvector(1,nrad_g6);
    xi_g6=dmatrix(1,nrad_g6,1,nphi_g6);
  }

  recv=dvector(1,nrad_g6);

  two_halo(1,1);

  for(i=1;i<=nrad_g6;++i)
    xi_mono[i]=xi_quad[i]=0;

  x1=one_halo_real_space(1)+two_halo_real_space(1);

  /* Check to see if we need to do internal bin sums
   */
  if(MCMC)
    {
      sdss_ilimit=0;
      for(i=0;i<nrad_g6;++i)
	if(Work.rad[i]<1.0)sdss_ilimit++;
      if(!ThisTask && flag1) {
	printf("SDSS_ILIMIT = %d (%f)\n",sdss_ilimit,Work.rad[sdss_ilimit-1]);
	flag1 = 0;
      }
      rp_min = 0.12;
    }

  /* Initialize the xi2d_interp grid
   */
  if(MCMC && sdss_ilimit)
    xi2d_interp_polar(0.2,0.23,PI/4,PI/4*1.1);

  t0=second();

  /* Calculate the redshift space correlation function on the polar grid.
   */
  for(i=istart;i<=nrad_g6;i+=istep)
    {
      NR=i;
      xi_mono[i]=0;
      xi_quad[i]=0;

      rr_g6[i]=exp((i-1)*dlogr+rlo);
      if(Work.chi2)rr_g6[i]=Work.rad[i-1];
      if(Work.SDSS_bins)rr_g6[i]=Work.rad[i-1];

      for(j=nphi_g6;j>=1;--j)
	{
	  phi_g6[j]=(j-0.5)/(nphi_g6)*PI/2.0;
	  rs=rr_g6[i]*cos(phi_g6[j]);
	  rp=rr_g6[i]*sin(phi_g6[j]);

	  /* If MCMC-ing z-space data,
	   * use polar bins for small-r multipoles.
	   * 
	   * upper i-limit is r=~1 Mpc/h
	   */	  
	  if(MCMC>1 && i<=sdss_ilimit)
	    {

	      /* New log bins -- assumes innermost bin is 0.31 (i==2)
	       */
	      rlo = pow(10.0,((i+2)*.2 - 1.09));
	      rhi = pow(10.0,((i+3)*.2 - 1.09));
	      philo = (PI/2 - asin(rp_min/rlo))/nphi_g6*(j-1);
	      phihi = (PI/2 - asin(rp_min/rlo))/nphi_g6*(j);
	      phi_g6[j] = 0.5*(philo + phihi);
	      delta_phi = (PI/2 - asin(rp_min/rlo))/nphi_g6;
	      rr_g6[i] = sqrt(rlo*rhi);
	      
	      xi = xi2d_interp_polar(rlo,rhi,philo,phihi);
	    }
	  else
	    {	      

	      t0=second();
	      x1h=one_halo(rs,rp);
	      t1=second();
	      t1h+=timediff(t0,t1);
	      
	      t0=second();
	      x2h=two_halo(rs,rp);
	      t1=second();
	      t2h+=timediff(t0,t1);
	      
	      xi=x1h+x2h;
	    }

	  xi_g6[i][j] = xi;

	  t1=second();

	  /* Calculate the monopole and quadrupole for this bin.
	   */
	  phi1=PI/2 - phi_g6[j];
	  
	  xi_mono[i]+=xi*sin(phi1)*delta_phi;
	  xi_quad[i]+=5*(3*cos(phi1)*cos(phi1)-1)/2.0*sin(phi1)*delta_phi*xi; 

	}

      /* Use spline interpolation to integrate the multipoles
       * unless in the small-r regime.
       */
      for(j=1;j<=nphi_g6;++j)
	phi_g6[j]=cos(PI/2 - phi_g6[j]);

      if(i<=sdss_ilimit && MCMC)
	x1=xi_mono[i];
      else
	x1=qtrap(monopole,0.0,cos(rp_min/rr_g6[i]),1.0E-4);

      if(i<=sdss_ilimit && MCMC)
	x2=xi_quad[i];
      else
	x2=qtrap(quadrupole,0.0,cos(rp_min/rr_g6[i]),1.0E-4);
      
      xi_mono[i]=x1;
      xi_quad[i]=x2;
      
      for(xi_bar=0,j=1;j<=i;++j)
	xi_bar+=xi_mono[j]*rr_g6[j]*rr_g6[j];
      xi_bar*=(3/pow(rr_g6[i]+0.5,3.0));

      if(OUTPUT && !ThisTask)
	{
	  fprintf(stdout,"MULTIPOLES%d %f %e %e %e %.2f %.2f\n",
		  ThisTask,rr_g6[i],xi_mono[i],xi_quad[i],-xi_quad[i]/(xi_bar - xi_mono[i]),t1h,t2h);
	  fflush(stdout);
	}
    }

#ifdef PARALLEL
  MPI_Allreduce(&xi_mono[1],&recv[1],nrad_g6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(i=1;i<=nrad_g6;++i)
    xi_mono[i]=recv[i];
  MPI_Allreduce(&xi_quad[1],&recv[1],nrad_g6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(i=1;i<=nrad_g6;++i)
    xi_quad[i]=recv[i];
#endif

  /* Now output all the information to a file (if we're not
   * in the middle of a chi2 fitting.
   */  
  sprintf(fname,"%s.multipoles",Task.root_filename);
  if(!ThisTask && !Work.chi2)
    fp=fopen(fname,"w");

  for(i=1;i<=nrad_g6;++i)
    {
      rr_g6[i]=exp((i-1)*dlogr+rlo);
      if(Work.chi2)rr_g6[i]=Work.rad[i-1];
    }

  dlogr = (log(rr_g6[2]) - log(rr_g6[1]));

  for(i=1;i<=nrad_g6;++i)
    {
      x1=one_halo_real_space(rr_g6[i])+two_halo_real_space(rr_g6[i]);

      for(j=1;j<=nrad_g6;++j)
	{
	  xi_mono[j] = log(xi_mono[j]);
	  rr_g6[j] = log(rr_g6[j]);
	}

      if(Work.chi2==-1)
	x2=qtrap(spherically_averaged_xi_temp,0.05,rr_g6[i],1.0E-3)*3/rr_g6[i]/rr_g6[i]/rr_g6[i];
      else
	x2=qtrap(spherically_averaged_xi,log(0.1),rr_g6[i],1.0E-3);


      x2*=3/(pow(exp(rr_g6[i]),3.0)-0.1*0.1*0.1);

      if(Work.i_monopole)
	Work.xi_mono[i-1]=exp(xi_mono[i]);
	//Work.xi_mono[i-1] = qtrap(spherically_averaged_xi,rr_g6[i]-0.5*dlogr,rr_g6[i]+0.5*dlogr,1.0E-4)*3/(pow(exp(rr_g6[i]+0.5*dlogr),3.0) - pow(exp(rr_g6[i]-0.5*dlogr),3.0));
      else
	Work.xi_mono[i-1]=xi_mono[i]/x1;

      for(j=1;j<=nrad_g6;++j)
	{
	  xi_mono[j] = exp(xi_mono[j]);
	  rr_g6[j] = exp(rr_g6[j]);
	}

      if(Work.i_quad2mono)
	Work.xi_quad[i-1]=xi_quad[i]/(xi_mono[i]);
      else
	Work.xi_quad[i-1]=xi_quad[i]/(xi_mono[i]-x2);

      /* Work.xi_rad[i-1]=rr_g6[i-1];*/

      if(!ThisTask && !Work.chi2)
	{
	  fprintf(fp,"%f %e %e %e %e %e\n",rr_g6[i],xi_mono[i]/x1,
		  xi_quad[i]/(xi_mono[i]-x2),xi_mono[i],xi_quad[i],x2);
	  fflush(fp);
	}
      if(!ThisTask && !Work.chi2)
	{
	  fprintf(stdout,"MULTI %f %e %e %e %e %e\n",rr_g6[i],xi_mono[i]/x1,
		  xi_quad[i]/(xi_mono[i]-x2),xi_mono[i],xi_quad[i],x2);
	  fflush(stdout);
	}
    }
  if(!ThisTask && !Work.chi2)
    fclose(fp);

}

/* Use the tabluated values of xi(s,p) and calculate the monpole.
 */
double monopole(double mu)
{
  static double *y2;
  static int flag=0, prev_nr=0;
  double a;

  if(prev_nr!=NR)
    {
      y2=dvector(1,nphi_g6);
      spline(phi_g6,xi_g6[NR],nphi_g6,1.0E+30,1.0E+30,y2);
    }
  prev_nr=NR;
  splint(phi_g6,xi_g6[NR],y2,nphi_g6,mu,&a);
  return(a);
}

/* Use the tabluated values of xi(s,p) and calculate the quadrupole.
 */
 
double quadrupole(double mu)
{
  static double *y2;
  static int flag=0, prev_nr=0;
  double a;
  if(prev_nr!=NR)
    {
      y2=dvector(1,nphi_g6);
      spline(phi_g6,xi_g6[NR],nphi_g6,1.0E+30,1.0E+30,y2);
    }
  prev_nr=NR;
  splint(phi_g6,xi_g6[NR],y2,nphi_g6,mu,&a);
  return(2.5*(3*mu*mu-1)*a);
}


/* For the quadrupole to monopole ratio, calculate the spherically averaged
 * monopole as in Hawkins etal. (2003)
 * For Kaiser model, FLAG=0. For HOD model, FLAG=1
 */
double spherically_averaged_xi(double s)
{
  static double *y2,*k2;
  static int flag=0,prev=0;
  double a;

  if(!flag || prev!=RESET_COSMOLOGY)
    {
      if(!flag) {
	y2=dvector(1,nrad_g6);
	k2=dvector(1,nrad_g6); 
      }
      flag++;
      spline(rr_g6,xi_mono,nrad_g6,1.0E+30,1.0E+30,y2);
      prev=RESET_COSMOLOGY;
    }
  splint(rr_g6,xi_mono,y2,nrad_g6,s,&a);    
  s=exp(s);
  return(s*s*s*exp(a));
}

double spherically_averaged_xi_temp(double s)
{
  static double *y2,*k2;
  static int flag=0,prev=0;
  double a;

  if(!flag || prev!=RESET_COSMOLOGY)
    {
      if(!flag) {
	y2=dvector(1,nrad_g6_temp);
	k2=dvector(1,nrad_g6_temp); 
      }
      flag++;
      spline(rr_g6_temp,xi_mono_temp,nrad_g6_temp,1.0E+30,1.0E+30,y2);
      prev=RESET_COSMOLOGY;
    }
  splint(rr_g6_temp,xi_mono_temp,y2,nrad_g6_temp,s,&a);    
  return(s*s*a);
}









