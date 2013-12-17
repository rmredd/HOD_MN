#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "header.h"

/***********************************************************************
 * This is the main function for controlling the activities of the code.
 * If you ask for any sort of minimization or MCMC analysis, the code is
 * redirected to another routine, but once finished will complete the rest
 * of the tasks checked off (so if you want to minimize and then print out
 * related statistics, the code will do that automatically if you ask it).
 *
 *
 * Each task generates an output file using the "root_filename" in the hod.bat
 * file and tacking on an extension at the end:
 *
 *  [filename].r_space  --> real space correlation function
 *  [filename].r_half   --> r_xi/2 statistic (z-space)
 *  [filename].qp_ratio --> multipole diagnostics of z-space xi
 *  [filename].linlin   --> 2-d z-space xi with linear bins
 *  [filename].loglog   --> 2-d z-space xi with logarithmic bins
 *  [filename].HOD      --> prints out the mean N(M)
 *  [filename].M2N		--> prints out M/N (against richness or mass)
 */


// A couple functions to make the M/N with richness easier to calculate
double p_rich(double, double);
double func_msum(double);
double func_nsatsum(double);
double func_ntotsum(double);

void tasks(int argc, char **argv)
{
  int i,j,nsize,ix,iy,nr,i1,i2,n;
  float x1,x2,x3,x4,x5,x6,err,npairs;
  double r,rs,rp,dx1,dx2,dx3,dx4,dx5,**xi2d,*xx2d,*yy2d,**xi2d_data,**xi2d_kaiser,
    xi0_m,xi2_m,xi0_k,xi2_k,xi0_d,xi2_d,xi_r,rlo,rhi,rr[50],rhalf[50],dr,
    rphi,rshi,rslo,rplo,dlogm,m,sig;
  FILE *fp,*fp2,*fp1;
  char fname[100];
  


  /* This is for chi^2 minimization of data for the projected correlation
   * function.
   */
  if(Task.wp_minimize && !HOD.color)
    wp_minimization(argv[1]);
  if(Task.wp_minimize && HOD.color)
    fit_color_samples();

  /* This is for chi^2 minimization of data for both the projected correlation
   * function and the redshift-space diagnostics.
   */
  if(Task.zspace_minimize)
    zspace_minimization(argv[1]);

  /* This is for Monte-Carlo Markov Chain exploration of the posterior
   * distribution of the parameters, either real-space or redshift-space,
   * depending on what MCMC is set to.
   */
muh(Task.MCMC);
  if(Task.MCMC==1) 
    mcmc_minimization();
  if(Task.MCMC>1)
    m2n_mcmc();

  /* This is to output the shape of the mean occupation function and the 
   * scatter about the mean. 
   * File columns are:
   *  1 - halo mass (M_sol/h)
   *  2 - N_cen(M)
   *  3 - N_sat(M)
   *  4 - N_tot(M)
   *  5 - <N(N-1)> 
   *  6 - ratio of scatter to Poisson scatter (often called "alpha")
   *  7 - mass fn dn(m)/dm
   */
  if(Task.HOD || Task.All)
    {
      sprintf(fname,"%s.HOD",Task.root_filename);
      fp=fopen(fname,"w");
      dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
      for(i=1;i<=100;++i)
	{
	  m=exp((i-1)*dlogm)*HOD.M_low;
	  sig = N_sat(m)*N_sat(m) + N_cen(m)*(1-N_cen(m));
/*	  fprintf(fp,"%e %e %e %e %e %e\n",
		  m,N_cen(m),N_sat(m),N_avg(m),sig,sig/(N_avg(m)*N_avg(m)));*/
        fprintf(fp,"%e %e %e %e %e %e %e\n",
            m,N_cen(m),N_sat(m),N_avg(m),sig,sig/(N_avg(m)*N_avg(m)),m*dndM_interp(m));
	}
      fclose(fp);
    }

  /* This sets in motion the calculation and tabulation of the real-space
   * correlation function and outputs it to a file.
   * File columns are:
   *  1 - radius (Mpc/h)
   *  2 - one-halo term (real-space)
   *  3 - two-halo term (real-space)
   *  4 - full xi(r)
   *  5 - projected correlation function (1st column is now r_p)
   */
  if(Task.real_space_xi || Task.All)
    {
      fprintf(stderr,"\n\nCALCULATING REAL-SPACE CORRELATION FUNCTION.\n");
      fprintf(stderr,    "--------------------------------------------\n\n");

      sprintf(fname,"%s.r_space",Task.root_filename);
      fp=fopen(fname,"w");
      dr=(log(70.0)-log(0.01))/49.0;
      for(i=0;i<50;++i)
	{
	  r=exp(i*dr+log(0.01));
	  x1=one_halo_real_space(r);
	  x2=two_halo_real_space(r);
	  x3=projected_xi(r);
	  x4 = projected_xi_rspace(r);
	  /*
	  rlo=pow(10.0,-1.+0.176*i);
	  rhi=pow(10.0,-1.+0.176*(i+1));
	  if(i<15)
	    x4 = qtrap(integrated_wp_bin,rlo,rhi,1.0E-3)/(0.5*(rhi*rhi-rlo*rlo));
	  else
	    x4 = 0;
	  */
	  fprintf(fp,"%f %e %e %e %e %e\n",r,x1,x2,x1+x2,x3,x4);
	  fflush(fp);
	}
      fclose(fp);
      
    }

  /* This calculates and outputs the angular correlation function,
   * and outputs it to a file.
   * File columns are:
   * 1 - theta (degrees)
   * 2 - w(theta)
   *
   * Note that this also reads in the n(z) file if that hasn't been done
   */
  if(Task.angular_xi || Task.All)
    {
      fprintf(stderr,"\n\nCALCULATING ANGULAR CORRELATION FUNCTION.\n");
      fprintf(stderr,    "--------------------------------------------\n\n");
      
      //First, check to see if the n(z) has already been read in, and, if
      //not, read it in
      if(wp.np_nz == 0) {
	if(!(fp=fopen(wp.fname_nz,"r"))) {
	  fprintf(stdout, "Error opening [%s]\n",wp.fname_nz);
	  endrun("error in Task.angular_xi");
	}
	wp.np_nz=filesize(fp);
	wp.z = dvector(1,wp.np_nz);
	wp.nz = dvector(1,wp.np_nz);

	for(i=1; i<=wp.np_nz; ++i); {
	  fscanf(fp,"%e %e",&x1,&x2);
	  wp.z[i] = x1;
	  wp.nz[i] = x2;
	}
	fclose(fp);
	fprintf(stderr,"Done reading %d lines from [%s]\n",wp.np_nz,wp.fname_nz);
      }

      sprintf(fname,"%s.wtheta",Task.root_filename);
      fp = fopen(fname,"w");
      dr = (log(100.)-log(0.01))/49.0; //This is actually dlntheta
      for(i=0; i<50; ++i) {
	r = exp(i*dr + log(0.01));
	x1 = wtheta_fit(r);

	fprintf(fp,"%f %e",r,x1);
	fflush(fp);
      }
      fclose(fp);

    }
	
  /* This runs the calculation of M/N from the HOD and outputs it to a file
   * File columns are:
   * 1 - halo mass or richness
   * 2 - M/N*ndens/rho_C
   * 3 - M/N
   */
   
  if(Task.M2N || Task.All)
    {
	  fprintf(stderr,"\n\nCALCULATING M/N.\n");
      fprintf(stderr,    "--------------------------------------------\n\n");

      sprintf(fname,"%s.m2n",Task.root_filename);
      fp=fopen(fname,"w");
	  fprintf(stderr,"   Type: %d \n",Task.M2N_type);
	  if(Task.M2N_type==0 || Task.M2N_type==2) {
		//Calculate for mass
	  
	    dlogm=(log(1E16)-log(1E13))/49.0;
	    for(i=0;i<50;++i)
	      {
		m=exp(i*dlogm+log(1E13));
		x1 = N_sat(m);  //nsat
		if (Task.M2N_type==2) x1 += N_cen(m);  //ntot
		x2 = m/x1;
		
		fprintf(fp,"%e %e %e %e\n",log10(m),x2*GALAXY_DENSITY/RHO_CRIT,x2,x1);
		fflush(fp);
		
	      }
	  } else {
	    //For richness instead
	    double x1,x2;
	    double binsize=1;
	    for(i=5;i<150/binsize;++i) {
	      x1=0;x2=0;
	      M2N_mass.curr_n200=i*binsize;
	      x1 = qromo(func_msum,log(1e11),log(1e16),midpnt);
	      if (Task.M2N_type==1) x2 = qromo(func_nsatsum,log(1e11),log(1e16),midpnt);
	      if (Task.M2N_type==3) x2 = qromo(func_ntotsum,log(1e11),log(1e16),midpnt);
	      fprintf(fp,"%e %e %e\n",i*binsize,x1/x2*GALAXY_DENSITY/RHO_CRIT,x1/x2);
	      fflush(fp);
	    }
	  }
	
      fclose(fp);
    }

  /* This takes a halofile from a simulation and populates the halos
   * with galaxies according the HOD specified in the batch file.
   */ 
  if(Task.populate_sim==1)
    populate_simulation();
  if(Task.populate_sim==2)
    populate_sampled_simulation();

  if(!ThisTask)
    OUTPUT=1;

  /* Using the analytic velocity model, calculate the pairwise velocity 
   * statistics of the galaxy/HOD.
   *
   * The format of the output file [root].PVD is:
   *   1 - r [Mpc/h]
   *   2 - sigma_v(radial) [km/s]
   *   3 - one-halo sigma_v (radial)
   *   4 - two-halo sigma_v (radial)
   *   3 - sigma_v(tangentail) [km/s]
   */
  if(Task.PVD)
    pairwise_velocity_dispersion();

  if(Task.massfunc)
    output_halo_mass_function();

  /* Output the linear and non-linear dark matter power spectrum.
   * Non-linear P(k) calculated using Smith et al.
   *
   * Format of file [root].matter_pk
   *   1 - k [h/Mpc]
   *   2 - linear P(k) [Mpc/h]^3
   *   3 - non-linear P(k) [Mpc/h]^3
   */
  if(Task.matter_pk)
    output_matter_power_spectrum();

  /* Output the linear and non-linear dark matter power spectrum.
   * Non-linear xi(r) is Fourier transform of Smith et al (above)
   *
   * Format of file [root].matter_pk
   *   1 - r [Mpc/h]
   *   2 - linear xi(r)
   *   3 - non-linear xi(r)
   */
  if(Task.matter_xi)
    output_matter_correlation_function();

  /* Output the matter variance as a function of scale.
   *
   * Format of file [root].sigma_r 
   *  1 - r [Mpc/h]
   *  2 - sigma(r) [linear]
   *  3 - sigma(r) [non-linear, using Smith]
   *  4 - mass [M_sol/h] mass = (4/3)*PI*r^3*rho_bar
   */
  if(Task.sigma_r)
    output_matter_variance();

  /* Output the halo concentrations using Bullock et al (2001) model
   * 
   * Format of file [root].civr
   *  1 - mass [Msol/h] --> "virial mass" using Bryan & Norman function for DELTA_VIR
   *  2 - mass [Msol/h] --> mass specified by DELTA_HALO (input file).
   *  3 - halo concentration.
   */
  if(Task.cvir)
    output_halo_concentrations();

  /* This redirects the program to calculate the multipoles of the redshift-space
   * correlation function. The function xi_multipoles will output a file with the
   * extension ".qp_ratio". File columns are:
   *  1 - radius (Mpc/h)
   *  2 - monopole to real-space ratio 
   *  3 - quadrupole moment (see definition in Hawkins et al or Tinker et al).
   *  4 - monopole
   *  5 - quadrupole
   */
  if(Task.multipoles)
    xi_multipoles();

  /* This calls functions to calculate and tabulate the r_xi/2 statistic.
   * The format of the output file is:
   *  1 - r_\sigma
   *  2 - r_xi/2
   */
  if(Task.r_half || Task.All)
    {
      rlo=-1;
      rhi=1.6;
      nr=50;
      for(i=0;i<nr;++i)
	rr[i]=pow(10.0,(i*1.0)/(nr-1.)*(rhi-rlo)+rlo);

      calc_rhalf(rr,rhalf,nr);
      if(!ThisTask)
	{
	  sprintf(fname,"%s.r_half",Task.root_filename);
	  fp2=fopen(fname,"w");
	  for(i=0;i<nr;++i)
	    {
	      fprintf(fp2,"%f %f\n",rr[i],rhalf[i]);
	      fflush(fp2);
	    }
	  fclose(fp2);
	}
    }


  /* This calls functions to tabulate and output the full two-dimensional
   * correlation function in redshift space, xi(r_sigma, r_pi). It is outputted
   * in two files with two different bin spacings.
   *
   * [filename].linlin -> linearly spaced bins of 2 Mpc/h per side. This format is
   * good for plotting two-dimensional contour plots. The size of the bin can be 
   * altered in the file linlin_bins.c. Because these bins are large, the 
   * correlation function must be calculated through a volume-weighted integral over 
   * the volume of the bin. Therefore columns 4 and 5 are only estimates of the true
   * one-halo and two-halo terms for that bin. Column 3 holds the most accurate value.
   *
   * The format of the file is:
   *  1 - r_sigma
   *  2 - r_pi
   *  3 - xi(rs,rp)
   *  4 - one-halo term (at center of bin)
   *  5 - two-halo term (at center of bin)
   *
   * [filename].loglog -> Logarithmically-spaced bins where the edges are defined by:
   *    r=10.0^(-1.12+0.2*j) for j>0 and r=0 for j=0. 
   *
   * Columns of output file are the same. Correlation function is also volume-weighted
   * in the 3rd column as well.
   */
  if(Task.z_space_xi==1 || Task.All)
    {
      if(!ThisTask) {
	fprintf(stderr,"\n\nCALCULATING REDSHIFT-SPACE CORRELATION FUNCTION.\n");
	fprintf(stderr,    "------------------------------------------------\n\n");


      linlin_bins();

      sprintf(fname,"%s.loglog",Task.root_filename);
      fp2=fopen(fname,"w");

      for(i=0;i<13;++i)
	for(j=0;j<13;++j)
	  {
	    rplo=pow(10.0,-1.12+0.2*j);
	    if(!j)rplo=0;
	    rphi=pow(10.0,-1.12+0.2*(j+1));
	    rslo=pow(10.0,-1.12+0.2*i);
	    if(!i)rslo=0;
	    rshi=pow(10.0,-1.12+0.2*(i+1));

	    rs=0.5*(rslo+rshi);
	    rp=0.5*(rplo+rphi);

	    dx1=one_halo(rs,rp);
	    dx2=two_halo(rs,rp);
	    dx5=xi2d_interp(rslo,rplo,rshi,rphi);
	    fprintf(fp2,"%e %e %e %e %e\n",rs,rp,dx5,dx1,dx2);
	    fprintf(stdout,"%e %e %e %e %e\n",rs,rp,dx5,dx1,dx2);
	    fflush(stdout);
	    fflush(fp2);
	  }
      }
    }

  if(Task.z_space_xi==2 || Task.All)
    {
      if(!ThisTask) {
	fprintf(stderr,"\n\nCALCULATING REDSHIFT-SPACE CORRELATION FUNCTION.\n");
	fprintf(stderr,    "------------------------------------------------\n\n");

	fp1=fopen("xi2d.02.g212.avg","r");
	n=12;

      sprintf(fname,"%s.linlin",Task.root_filename);
      fp2=fopen(fname,"w");

      rslo = 0;
      for(i=0;i<n;++i)
	{
	  rplo=0;
	  for(j=0;j<n;++j)
	    {
	      fscanf(fp1,"%f %f %f %f %f %f",&x1,&x2,&x3,&x4,&x5,&x6);
	      rphi=x6;
	      rshi=x5;
		
	      rs=0.5*(rslo+rshi);
	      rp=0.5*(rplo+rphi);
	      
	      dx5=xi2d_interp(rslo,rplo,rshi,rphi);
	      fprintf(fp2,"%e %e %e %e %e %e %e\n",dx5,x1,x4,rplo,rphi,rslo,rshi);
	      fprintf(stdout,"%e %e %e %e %e %e %e\n",dx5,x1,x4,rplo,rphi,rslo,rshi);
	      fflush(stdout);
	      fflush(fp2);
	      rplo = rphi;
	    }
	  rslo = rshi;
	}
      }
    }

  //endrun("finished with tasks");
}
