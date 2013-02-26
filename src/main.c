#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* test file routines.
 */
void test(int argc, char **argv);
void chi2_grid(int argc, char **argv);
void fit_scale_bias(int argc, char **argv);
void aspen_breakout(void);
void populate_simulation_clf(void);

int main(int argc, char **argv)
{
  double s1;
  int i;

#ifdef PARALLEL
  printf("STARTING>>>\n");
  fflush(stdout);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  printf("TASK %d reporting for duty.\n",ThisTask);
  fflush(stdout);
#endif

  ARGC = argc;
  ARGV = argv;

  OUTPUT=0;
  HOD.fredc = HOD.freds = 1.0;

  for(i=1;i<=99;++i)
    HOD.free[i]=0;
  wp.esys=0;

  Work.chi2=0;
  Work.imodel=1;

  USE_ERRORS = 0;
  ITRANS=4;
  HUBBLE=0.7;
  BEST_FIT = 0;
  HOD.M_sat_break = 1.0e14;
  HOD.alpha1 = 1.0;

  if(argc==1)
    endrun("./HOD.x hod.bat_file [MCMC output]");

  read_parameter_file(argv[1]);
  
  if(argc>2) {
		sprintf(Task.mcmcfilename,"%s",argv[2]);
	}
	else {
		sprintf(Task.mcmcfilename,"%s.mcmc",Task.root_filename);
	}
  fprintf(stderr,"Output MCMC file: %s \n",Task.mcmcfilename);

  /* If there's no cross-correlation function,
   * set the second number density equal to the first
   */
  if(!XCORR)
    GALAXY_DENSITY2 = GALAXY_DENSITY;

  /* Initialize the non-linear power spectrum.
   */
  nonlinear_sigmac(8.0);
  sigmac_interp(1.0E13);
  sigmac_radius_interp(1.0);

  if(REDSHIFT>0)
    {
      SIGMA_8 = SIGMA_8*growthfactor(REDSHIFT);
      fprintf(stderr,"SIGMA8(z)= %f\n",SIGMA_8);
      RESET_COSMOLOGY++;
    }

  if(argc>2)
    IDUM_MCMC=atoi(argv[2]);
  //if(MCMC)m2n_mcmc();

  /* Get the galaxy bias factor
   */
  s1=qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt);
  GALAXY_BIAS=s1/GALAXY_DENSITY;
  if(OUTPUT)
    fprintf(stdout,"Galaxy Bias bg= %f\n",GALAXY_BIAS);
  fflush(stdout);

  /* Get the galaxy satellite fraction
   */
  s1=qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/
    GALAXY_DENSITY;
  if(OUTPUT)
    fprintf(stdout,"fsat %e\n",s1);
  fflush(stdout);

  /* Mean halo mass.
   */
  if(OUTPUT)
    fprintf(stdout,"M_eff %e\n",number_weighted_halo_mass());
  fflush(stdout);

  /* Set up BETA for wp integration.
   */
  BETA = pow(OMEGA_M,0.6)/GALAXY_BIAS;
  if(OUTPUT)
    printf("BETA = %f\n",BETA);

  //two_halo_real_space(1.0);
  //two_halo_real_space(35.0);

  /* neg_binomial();*/

  /* Check for extra commands:
   * arg==999 goes to the test program, superceding tasks.
   * arg<0 supercedes the MCMC random number in the batfile.
   */
  if(argc>2)
    {
      if(atoi(argv[2])==999)
	test(argc,argv);
      //if(atoi(argv[2])==991)
      //chi2_grid(argc,argv);
      if(atoi(argv[2])==990)
	fit_scale_bias(argc,argv);
      if(atoi(argv[2])==998)
	output_velocity_moments(atoi(argv[3]));
       if(atoi(argv[2])<0)
	IDUM_MCMC=atoi(argv[2]);
      if(atoi(argv[2])==997)
	{
	  VBIAS_MASS_THRESHOLD=pow(10.0,atof(argv[3])/10.0);
	  if(argc>4)
	    VBIAS=atof(argv[4]);
	  if(argc>5)
	    sprintf(Task.root_filename,"%s",argv[5]);
	  fprintf(stdout,"VBIAS_MASS_THRESHOLD= %e\n",VBIAS_MASS_THRESHOLD);
	  fprintf(stdout,"VBIAS= %e\n",VBIAS);
	  fprintf(stdout,"root filename [%s]\n",Task.root_filename);
	}	  
      if(atoi(argv[2])==998)
	{
	  VBIAS_SLOPE=atof(argv[3]);
	  if(argc>4)
	    sprintf(Task.root_filename,"%s",argv[4]);
	  fprintf(stdout,"VBIAS= %e\n",VBIAS);
	  fprintf(stdout,"VBIAS_SLOPE= %e\n",VBIAS_SLOPE);
	  fprintf(stdout,"root filename [%s]\n",Task.root_filename);
	}
      if(atoi(argv[2])==991)
	//populate_simulation_clf();
	uber_comparison();
      if(atoi(argv[2])==9910)
	ml_ratios();
      if(atoi(argv[2])==988)
	aspen_breakout();
	 
    }

  fprintf(stderr,"M_low = %e \n",HOD.M_low);

  tasks(argc,argv);

  fprintf(stderr,"M_low = %e \n",HOD.M_low);

}

