#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "header.h"

/* This function is taken from Volker Springel's GADGET code and
 * modified for the parameters used for the HOD.x code.
 */

/*
 *  This function parses the parameterfile in a simple way.
 *  Each paramater is defined by a keyword (`tag'), and can be
 *  either of type douple, int, or character string.
 *  The routine makes sure that each parameter appears 
 *  exactly once in the parameterfile.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define CHAR 4
#define LONG 4
#define MAXTAGS 300

  FILE *fd,*fdout;

  char buf[200],buf1[200],buf2[200],buf3[200],tempchar;
  int  i,j,nt,ii,nn,ctemp;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][200];
  int  errorFlag=0;
  int IDUM_MCMC_TEMP=-555;

  nt=0;
 
  fprintf(stdout,"Tags are %i of %d\n",nt,MAXTAGS);
    
  strcpy(tag[nt],"RESTART");
  addr[nt]=&RESTART;
  id[nt++]=INT;

  strcpy(tag[nt],"RESTART_FILE");
  addr[nt]=&RESTART_FILE;
  id[nt++]=STRING;

  strcpy(tag[nt],"GAMMA");
  addr[nt]=&GAMMA;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"OMEGA_M");
  addr[nt]=&OMEGA_M;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"OMEGA_TEMP");
  addr[nt]=&OMEGA_TEMP;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"OMEGA_B");
  addr[nt]=&OMEGA_B;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"SIGMA_8");
  addr[nt]=&SIGMA_8;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"RHO_CRIT");
  addr[nt]=&RHO_CRIT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"HUBBLE");
  addr[nt]=&HUBBLE;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"REDSHIFT");
  addr[nt]=&REDSHIFT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"ITRANS");
  addr[nt]=&ITRANS;
  id[nt++]=INT;
  
  strcpy(tag[nt],"LINEAR_PSP");
  addr[nt]=&LINEAR_PSP;
  id[nt++]=INT;
  
  strcpy(tag[nt],"KAISER");
  addr[nt]=&KAISER;
  id[nt++]=INT;
  
  strcpy(tag[nt],"BETA");
  addr[nt]=&BETA;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"SIGV");
  addr[nt]=&SIGV;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"IDUM_MCMC");
  addr[nt]=&IDUM_MCMC_TEMP;
  id[nt++]=INT;
  
  strcpy(tag[nt],"SPECTRAL_INDX");
  addr[nt]=&SPECTRAL_INDX;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"DELTA_CRIT");
  addr[nt]=&DELTA_CRIT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"DELTA_HALO");
  addr[nt]=&DELTA_HALO;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"BOX_SIZE");
  addr[nt]=&BOX_SIZE;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"RESOLUTION");
  addr[nt]=&RESOLUTION;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"VBIAS");
  addr[nt]=&VBIAS;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"VBIAS_C");
  addr[nt]=&VBIAS_C;
  id[nt++]=DOUBLE;

  //Added scale-dep bias parameters
  strcpy(tag[nt],"HBIAS_C1");
  addr[nt]=&HBIAS_C1;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"HBIAS_C2");
  addr[nt]=&HBIAS_C2;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"HBIAS_D1");
  addr[nt]=&HBIAS_D1;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"HBIAS_D2");
  addr[nt]=&HBIAS_D2;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"TF_file"); 
  addr[nt]=Files.TF_file;
  id[nt++]=STRING;

  strcpy(tag[nt],"M1");
  addr[nt]=&HOD.M1;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"M_min");
  addr[nt]=&HOD.M_min;
  id[nt++]=DOUBLE;
        
  strcpy(tag[nt],"M_cen_max");
  addr[nt]=&HOD.M_cen_max;
  id[nt++]=DOUBLE;
        
  strcpy(tag[nt],"M_cut");
  addr[nt]=&HOD.M_cut;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"M_max");
  addr[nt]=&HOD.M_max;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"sigma_logM");
  addr[nt]=&HOD.sigma_logM;
  id[nt++]=DOUBLE;

  HOD.MaxCen=1;
  strcpy(tag[nt],"MaxCen");
  addr[nt]=&HOD.MaxCen;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"alpha");
  addr[nt]=&HOD.alpha;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"alpha1");
  addr[nt]=&HOD.alpha1;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"M_sat_break");
  addr[nt]=&HOD.M_sat_break;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"M_cen_lin");
  addr[nt]=&HOD.M_cen_lin;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"pdfc");
  addr[nt]=&HOD.pdfc;
  id[nt++]=INT;

  strcpy(tag[nt],"pdf");
  addr[nt]=&tempchar;
  id[nt++]=CHAR;
    
  strcpy(tag[nt],"pdfs");
  addr[nt]=&HOD.pdfs;
  id[nt++]=INT;

  strcpy(tag[nt],"M_min_fac");
  addr[nt]=&HOD.M_min_fac;
  id[nt++]=DOUBLE;

  /* Paramaters for the second HOD function
   * (for x-corr)
   */
  strcpy(tag[nt],"XCORR");
  addr[nt]=&XCORR;
  id[nt++]=INT;
  XCORR=0;

  strcpy(tag[nt],"GALAXY_DENSITY2");
  addr[nt]=&GALAXY_DENSITY2;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"GALDENS_ERR2");
  addr[nt]=&GALDENS_ERR2;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"HOD2.M1");
  addr[nt]=&HOD2.M1;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"HOD2.M_min");
  addr[nt]=&HOD2.M_min;
  id[nt++]=DOUBLE;
        
  strcpy(tag[nt],"HOD2.M_cen_max");
  addr[nt]=&HOD2.M_cen_max;
  id[nt++]=DOUBLE;
        
  strcpy(tag[nt],"HOD2.M_cut");
  addr[nt]=&HOD2.M_cut;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"HOD2.M_max");
  addr[nt]=&HOD2.M_max;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"HOD2.sigma_logM");
  addr[nt]=&HOD2.sigma_logM;
  id[nt++]=DOUBLE;

  HOD2.MaxCen=1;
  strcpy(tag[nt],"HOD2.MaxCen");
  addr[nt]=&HOD2.MaxCen;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"HOD2.alpha");
  addr[nt]=&HOD2.alpha;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"HOD2.alpha1");
  addr[nt]=&HOD2.alpha1;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"HOD2.M_sat_break");
  addr[nt]=&HOD2.M_sat_break;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"HOD2.pdfc");
  addr[nt]=&HOD2.pdfc;
  id[nt++]=INT;
  HOD2.pdfc=-1;

  strcpy(tag[nt],"HOD2.pdfs");
  addr[nt]=&HOD2.pdfs;
  id[nt++]=INT;
  HOD2.pdfs=-1;
    
  /* Finished with HOD2 params
   */

  strcpy(tag[nt],"color");
  addr[nt]=&HOD.color;
  id[nt++]=INT;
  HOD.color=0;

  strcpy(tag[nt],"fblue0_cen");
  addr[nt]=&HOD.fblue0_cen;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"sigma_fblue_cen");
  addr[nt]=&HOD.sigma_fblue_cen;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"fblue0_sat");
  addr[nt]=&HOD.fblue0_sat;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"sigma_fblue_sat");
  addr[nt]=&HOD.sigma_fblue_sat;
  id[nt++]=DOUBLE;
    
    
  strcpy(tag[nt],"GALAXY_DENSITY");
  addr[nt]=&GALAXY_DENSITY;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"GALDENS_ERR");
  addr[nt]=&GALDENS_ERR;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"EXCLUSION");
  addr[nt]=&EXCLUSION;
  id[nt++]=INT;

  strcpy(tag[nt],"FIX_PARAM");
  addr[nt]=&FIX_PARAM;
  id[nt++]=INT;

  strcpy(tag[nt],"POWELL");
  addr[nt]=&POWELL;
  id[nt++]=INT;

  strcpy(tag[nt],"OUTPUT");
  addr[nt]=&OUTPUT;
  id[nt++]=INT;

  for(i=0;i<=21;++i)
    {
      sprintf(tag[nt],"free[%d]",i);
      addr[nt]=&HOD.free[i];
      id[nt++]=INT;
    }

  strcpy(tag[nt],"All");
  addr[nt]=&Task.All;
  id[nt++]=INT;
    
  strcpy(tag[nt],"real_space_xi");
  addr[nt]=&Task.real_space_xi;
  id[nt++]=INT;
    
  strcpy(tag[nt],"z_space_xi");
  addr[nt]=&Task.z_space_xi;
  id[nt++]=INT;
    
  strcpy(tag[nt],"kaiser_xi");
  addr[nt]=&Task.kaiser_xi;
  id[nt++]=INT;
    
  strcpy(tag[nt],"multipoles");
  addr[nt]=&Task.multipoles;
  id[nt++]=INT;
    
  strcpy(tag[nt],"r_half");
  addr[nt]=&Task.r_half;
  id[nt++]=INT;
    
  strcpy(tag[nt],"PVD");
  addr[nt]=&Task.PVD;
  id[nt++]=INT;
  Task.PVD = 0;

  strcpy(tag[nt],"cvir");
  addr[nt]=&Task.cvir;
  id[nt++]=INT;
  Task.cvir = 0;
    
  strcpy(tag[nt],"massfunc");
  addr[nt]=&Task.massfunc;
  id[nt++]=INT;
  Task.cvir = 0;
    
  strcpy(tag[nt],"matter_xi");
  addr[nt]=&Task.matter_xi;
  id[nt++]=INT;
  Task.matter_xi = 0;
    
  strcpy(tag[nt],"matter_pk");
  addr[nt]=&Task.matter_pk;
  id[nt++]=INT;
  Task.matter_pk = 0;
    
  strcpy(tag[nt],"sigma_r");
  addr[nt]=&Task.sigma_r;
  id[nt++]=INT;
  Task.sigma_r = 0;
    
  strcpy(tag[nt],"wp_minimize");
  addr[nt]=&Task.wp_minimize;
  id[nt++]=INT;
  Task.wp_minimize=0;
  
  strcpy(tag[nt],"m2n_minimize");
  addr[nt]=&Task.m2n_minimize;
  id[nt++]=INT;
  Task.m2n_minimize=0;

  strcpy(tag[nt],"zspace_minimize");
  addr[nt]=&Task.zspace_minimize;
  id[nt++]=INT;
  Task.zspace_minimize=0;

  strcpy(tag[nt],"MCMC");
  addr[nt]=&Task.MCMC;
  id[nt++]=INT;
  Task.MCMC=0;

  strcpy(tag[nt],"populate_sim");
  addr[nt]=&Task.populate_sim;
  id[nt++]=INT;
  Task.populate_sim=0;

  strcpy(tag[nt],"HaloFile");
  addr[nt]=&Files.HaloFile;
  id[nt++]=STRING;

  strcpy(tag[nt],"HaloDensityFile");
  addr[nt]=&Files.HaloDensityFile;
  id[nt++]=STRING;
  
  //Added to supply input covariance matrix for MCMC
  strcpy(tag[nt],"CovarMCMCFile");
  addr[nt]=&Files.CovarMCMCFile;
  id[nt++]=STRING;
  Files.i_CovarMCMC=1;

  strcpy(tag[nt],"DENSITY_DEPENDENCE");
  addr[nt]=&DENSITY_DEPENDENCE;
  id[nt++]=INT;

  strcpy(tag[nt],"DENSITY_THRESHOLD");
  addr[nt]=&DENSITY_THRESHOLD;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"WP_ONLY");
  addr[nt]=&WP_ONLY;
  id[nt++]=INT;  

  strcpy(tag[nt],"M2N");
  addr[nt]=&Task.M2N;
  id[nt++]=INT;
  Task.M2N=0;  
  
  strcpy(tag[nt],"M2N_type");
  addr[nt]=&Task.M2N_type;
  id[nt++]=INT;
  Task.M2N_type=0;  

  strcpy(tag[nt],"HOD");
  addr[nt]=&Task.HOD;
  id[nt++]=INT;
  Task.HOD=0;

  strcpy(tag[nt],"COVAR");
  addr[nt]=&COVAR;
  id[nt++]=INT;
  COVAR=1;

  strcpy(tag[nt],"MN_COVAR");
  addr[nt]=&MN_COVAR;
  id[nt++]=INT;
  MN_COVAR=1;
    
  strcpy(tag[nt],"PCA");
  addr[nt]=&PCA;
  id[nt++]=INT;
  PCA=0;
    
  strcpy(tag[nt],"wp_npca");
  addr[nt]=&wp.npca;
  id[nt++]=INT;
  wp.npca=0;
    
  strcpy(tag[nt],"DEPROJECTED");
  addr[nt]=&DEPROJECTED;
  id[nt++]=INT;
    
  strcpy(tag[nt],"fname_covar");
  addr[nt]=&wp.fname_covar;
  id[nt++]=STRING;
    
  strcpy(tag[nt],"mn_fname_covar");
  addr[nt]=&M2N_mass.fname_covar;
  id[nt++]=STRING;

  wp.pi_max=40;
  strcpy(tag[nt],"pi_max");
  addr[nt]=&wp.pi_max;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"esys");
  addr[nt]=&wp.esys;
  id[nt++]=DOUBLE;
 
  strcpy(tag[nt],"fname_wp");
  addr[nt]=&wp.fname_wp;
  id[nt++]=STRING;

  wp.format=1;
  strcpy(tag[nt],"wp_format");
  addr[nt]=&wp.format;
  id[nt++]=INT;

  wp.n_wp=9;
  strcpy(tag[nt],"n_wp");
  addr[nt]=&wp.n_wp;
  id[nt++]=INT;

  strcpy(tag[nt],"root_filename");
  addr[nt]=&Task.root_filename;
  id[nt++]=STRING;
   
  strcpy(tag[nt],"CVIR_FAC");
  addr[nt]=&CVIR_FAC;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"JENKINS_A");
  addr[nt]=&JENKINS_A;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"JENKINS_C");
  addr[nt]=&JENKINS_C;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"JENKINS_B");
  addr[nt]=&JENKINS_B;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"BEST_FIT");
  addr[nt]=&BEST_FIT;
  id[nt++]=INT;


  /* M2N filename and params
   */
  strcpy(tag[nt],"m2n_filename");
  addr[nt]=&M2N_mass.filename;
  id[nt++]=STRING;

  /* Input filename and parameters for w(theta) calculations
   */
  strcpy(tag[nt],"FIT_WTHETA");
  addr[nt]=&Task.FIT_WTHETA;
  id[nt++]=INT;
  Task.FIT_WTHETA=0;

  strcpy(tag[nt],"fname_nz");
  addr[nt]=&wp.fname_nz;
  id[nt++]=STRING;
  wp.np_nz=0;

  strcpy(tag[nt],"angular_xi");
  addr[nt]=&Task.angular_xi;
  id[nt++]=INT;
  Task.angular_xi = 0;

  /* Input parameters for correctly handling redshift range for number density
   */
  strcpy(tag[nt],"REDSHIFT_MIN");
  addr[nt]=&REDSHIFT_MIN;
  id[nt++]=DOUBLE;
  REDSHIFT_MIN=-1;

  strcpy(tag[nt],"REDSHIFT_MAX");
  addr[nt]=&REDSHIFT_MAX;
  id[nt++]=DOUBLE;
  REDSHIFT_MAX=-1;

  if((fd=fopen(fname,"r")))
    {
      nn=filesize(fd);
      sprintf(buf,"%s","hod-usedvalues");
      if(!(fdout=fopen(buf,"w")))
	{
	  fprintf(stdout,"error opening file '%s' \n",buf);
	  errorFlag=1; 
	}
      else
	{
	  /*while(!feof(fd))*/
	  for(ii=1;ii<=nn;++ii)
	    {
	      fgets(buf,200,fd);
	      if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
		continue;
	      
	      if(buf1[0]=='%')
		continue;
	      
	      for(i=0,j=-1;i<nt;i++)
		if(strcmp(buf1,tag[i])==0)
		  {
		    j=i;
		    tag[i][0]=0;
		    break;
		  }
	      
	      if(j>=0)
		{
		  switch(id[j])
		    {
		    case DOUBLE:
		      *((double*)addr[j])=atof(buf2); 
		      fprintf(fdout,"%-35s%g\n",buf1,*((double*)addr[j]));
		      break;
		    case STRING:
		      strcpy(addr[j],buf2);
		      fprintf(fdout,"%-35s%s\n",buf1,buf2);
		      break;
		    case INT:
		      *((int*)addr[j])=atoi(buf2);
		      fprintf(fdout,"%-35s%d\n",buf1,*((int*)addr[j]));
		      break;
		    case CHAR:
		      *((char*)addr[j])=buf2[0];
		      fprintf(fdout,"%-35s%c\n",buf1,*((int*)addr[j]));
		      break;
		    }
		}
	      else
		{
		  fprintf(stderr,"Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			  fname,buf1);
		  errorFlag=1;
		}
	    }
	}
      fclose(fd);
      fclose(fdout);

    }
  else
    {
      fprintf(stderr,"Parameter file %s not found.\n", fname);
      exit(1);
    }


  /* Check the params for some basic interpretation
   */
  MASS_PER_PARTICLE=pow(RESOLUTION,3.0)*RHO_CRIT*OMEGA_M;
  if(HOD.M_min<1.0E4 && HOD.M_min>0)
    {
      HOD.M_min=pow(RESOLUTION,3.0)*RHO_CRIT*OMEGA_M*HOD.M_min;
      fprintf(stderr,"HOD.M_min= %e\n",HOD.M_min);
    }
  if(HOD.M1<1.0E5)
    {
      HOD.M1=HOD.M1*MASS_PER_PARTICLE;
      fprintf(stderr,"HOD.M1= %e\n",HOD.M1);
    }


  /* SCale the M_max value by OMEGA_M=0.3
   */
  /*
  HOD.M_max*=(OMEGA_M/0.3);
  */

  for(i=0;i<nt;i++)
    {      
      if(!strcmp(tag[i],"DENSITY_DEPENDENCE"))continue;
      if(!strcmp(tag[i],"DENSITY_THRESHOLD") && !DENSITY_DEPENDENCE)continue;
      if(!strcmp(tag[i],"M_min_fac") && !DENSITY_DEPENDENCE)continue;
      if(!strcmp(tag[i],"HaloDensityFile") && !DENSITY_DEPENDENCE)continue;
      
      if(!strcmp(tag[i],"HOD2.M_min") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.M_max") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.MaxCen") && !(HOD2.pdfc>=4 && HOD.pdfc<=6))continue;
      if(!strcmp(tag[i],"HOD2.M_cut") && HOD2.pdfs<2)continue;
      if(!strcmp(tag[i],"HOD2.M1") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.M_cen_max") && HOD2.pdfc!=7)continue;
      if(!strcmp(tag[i],"HOD2.sigma_logM") && (!XCORR || (HOD2.pdfc==1 || HOD2.pdfc==7)))continue;
      if(!strcmp(tag[i],"HOD2.pdfs") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.pdfc") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.alpha") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.alpha1") && HOD2.pdfs!=4 && HOD2.pdfs!=5)continue;
      if(!strcmp(tag[i],"HOD2.M_sat_break") && HOD2.pdfs!=4 && HOD2.pdfs!=5)continue;
      if(!strcmp(tag[i],"GALDENS_ERR")) continue;
      if(!strcmp(tag[i],"GALAXY_DENSITY2") && !XCORR)continue;
      if(!strcmp(tag[i],"GALDENS_ERR2") && !XCORR)continue;
      if(!strcmp(tag[i],"XCORR"))continue;

      if(!strcmp(tag[i],"alpha1"))continue;
      if(!strcmp(tag[i],"M_sat_break"))continue;
      if(!strcmp(tag[i],"M_cen_lin"))continue;
      
      if(!strcmp(tag[i],"OMEGA_TEMP"))
	{
	  OMEGA_TEMP = OMEGA_M;
	  continue;
	}
      if(!strcmp(tag[i],"pdf"))continue;
      if(!strcmp(tag[i],"REDSHIFT"))continue;
      if(!strcmp(tag[i],"WP_ONLY"))continue;
      if(!strcmp(tag[i],"RESTART"))continue;
      if(!strcmp(tag[i],"RESTART_FILE"))continue;
      if(!strcmp(tag[i],"DEPROJECTED"))continue;
      if(!strcmp(tag[i],"POWELL"))continue;
      if(!strcmp(tag[i],"LINEAR_PSP"))continue;
      if(!strcmp(tag[i],"BOX_SIZE"))continue;
      if(!strcmp(tag[i],"RESOLUTION"))continue;
      if(!strcmp(tag[i],"DELTA_HALO"))continue;
      if(!strcmp(tag[i],"VBIAS"))continue;
      if(!strcmp(tag[i],"COVAR"))continue;
      if(!strcmp(tag[i],"MN_COVAR"))continue;
      if(!strcmp(tag[i],"VBIAS_C"))continue;
      if(!strcmp(tag[i],"CVIR_FAC"))continue;
      if(!strcmp(tag[i],"ITRANS"))continue;
      if(!strcmp(tag[i],"IDUM_MCMC"))continue;
      if(!strcmp(tag[i],"FIX_PARAM"))continue;
      if(!strcmp(tag[i],"DEPROJECTED"))continue;
      if(!strcmp(tag[i],"OUTPUT"))continue;
      if(!strcmp(tag[i],"JENKINS_A"))continue;
      if(!strcmp(tag[i],"JENKINS_B"))continue;
      if(!strcmp(tag[i],"JENKINS_C"))continue;
      if(!strcmp(tag[i],"BEST_FIT"))continue;
      if(!strcmp(tag[i],"KAISER"))continue;
      if(!strcmp(tag[i],"SIGV"))continue;
      if(!strcmp(tag[i],"BETA"))continue;
      if(!strcmp(tag[i],"HUBBLE"))continue;
      if(!strcmp(tag[i],"OMEGA_B"))continue;
      
      //scale-dep bias
      if(!strcmp(tag[i],"HBIAS_C1"))continue;
      if(!strcmp(tag[i],"HBIAS_C2"))continue;
      if(!strcmp(tag[i],"HBIAS_D1"))continue;
      if(!strcmp(tag[i],"HBIAS_D2"))continue;
      
      if(!strcmp(tag[i],"MCMC"))continue;
      if(!strcmp(tag[i],"wp_minimize"))continue;
	  if(!strcmp(tag[i],"m2n_minimize"))continue;
      if(!strcmp(tag[i],"wp_format"))continue;
      if(!strcmp(tag[i],"n_wp"))continue;
      if(!strcmp(tag[i],"wp_npca"))continue;
      if(!strcmp(tag[i],"zspace_minimize"))continue;
      if(!strcmp(tag[i],"pi_max"))continue;
      if(!strcmp(tag[i],"esys"))continue;

      if(HOD.color)
	{
	  if(!strcmp(tag[i],"fblue0_cen") ||
	     !strcmp(tag[i],"sigma_fblue_cen") ||
	     !strcmp(tag[i],"fblue0_sat") ||
	     !strcmp(tag[i],"sigma_fblue_sat"))
	    {
	      fprintf(stderr,"Parameters for color HOD not specified.\n");
	      exit(0);
	    }
	  continue;
	}
	     
      if(!strcmp(tag[i],"color"))continue;
      if(!strcmp(tag[i],"fblue0_cen"))continue;
      if(!strcmp(tag[i],"sigma_fblue_cen"))continue;
      if(!strcmp(tag[i],"fblue0_sat"))continue;
      if(!strcmp(tag[i],"sigma_fblue_sat"))continue;


      if(!strcmp(tag[i],"free[0]"))continue;
      if(!strcmp(tag[i],"free[1]"))continue;
      if(!strcmp(tag[i],"free[2]"))continue;
      if(!strcmp(tag[i],"free[3]"))continue; 
      if(!strcmp(tag[i],"free[4]"))continue;
      if(!strcmp(tag[i],"free[5]"))continue;
      if(!strcmp(tag[i],"free[6]"))continue;
      if(!strcmp(tag[i],"free[7]"))continue;
      if(!strcmp(tag[i],"free[8]"))continue;
      if(!strcmp(tag[i],"free[9]"))continue;
      if(!strcmp(tag[i],"free[10]"))continue;
      if(!strcmp(tag[i],"free[11]"))continue;
      if(!strcmp(tag[i],"free[12]"))continue;
      if(!strcmp(tag[i],"free[13]"))continue;
      if(!strcmp(tag[i],"free[14]"))continue;
      if(!strcmp(tag[i],"free[15]"))continue;
      if(!strcmp(tag[i],"free[16]"))continue;
      if(!strcmp(tag[i],"free[17]"))continue;
      if(!strcmp(tag[i],"free[18]"))continue;
      if(!strcmp(tag[i],"free[19]"))continue;
      if(!strcmp(tag[i],"free[20]"))continue;
      if(!strcmp(tag[i],"free[21]"))continue;

      if(!strcmp(tag[i],"All"))continue;
      if(!strcmp(tag[i],"populate_sim"))continue;
      if(!strcmp(tag[i],"HaloFile"))continue;
	  if(!strcmp(tag[i],"CovarMCMCFile")) {
		//set to NOT use the file
		Files.i_CovarMCMC=0;
		continue;}
      if(!strcmp(tag[i],"HOD"))continue;
	  if(!strcmp(tag[i],"M2N"))continue;
	  if(!strcmp(tag[i],"M2N_type"))continue;
      if(!strcmp(tag[i],"PCA"))continue;
      if(!strcmp(tag[i],"PVD"))continue;
      if(!strcmp(tag[i],"matter_xi"))continue;
      if(!strcmp(tag[i],"matter_pk"))continue;
      if(!strcmp(tag[i],"sigma_r"))continue;
      if(!strcmp(tag[i],"kaiser_xi"))continue;
      if(!strcmp(tag[i],"cvir"))continue;
      if(!strcmp(tag[i],"massfunc"))continue;
      if(!strcmp(tag[i],"angular_xi"))continue;

      if(!strcmp(tag[i],"REDSHIFT_MIN"))continue;
      if(!strcmp(tag[i],"REDSHIFT_MAX"))continue;

      if(!strcmp(tag[i],"FIT_WTHETA"))continue;

      //If we want to do w(theta), but there's no fname_nz, give up
      if(!strcmp(tag[i],"fname_nz")) {
	if( (Task.angular_xi) || Task.FIT_WTHETA || Task.All ) {
	  fprintf(stderr,"No filename specified for n(z) data\n");
	  errorFlag=1;
	}
	continue;
      }

      if(!strcmp(tag[i],"TF_file"))
	{
	  if(ITRANS==11) {
	    sprintf(Files.TF_file,"CMBFAST_trans.dat");
	    fprintf(stderr,"No transfer function file, using [%s]\n",Files.TF_file);
	  }
	  continue;
	}

      if(!strcmp(tag[i],"fname_covar"))
	{
	  if(Task.wp_minimize) {
	    fprintf(stderr,"No filename specified for covariance matrix.\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"fname_wp"))
	{
	  if(Task.wp_minimize) {
	    fprintf(stderr,"No filename specified for wp data.\n");
	    errorFlag=1;
	  }
	  continue;
	}
	  if(!strcmp(tag[i],"m2n_filename"))
	{
	  if(Task.m2n_minimize) {
		fprintf(stderr,"No filename specified for m2n data.\n");
		errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"mn_fname_covar"))
    {
        if(Task.m2n_minimize) {
            fprintf(stderr,"No filename specified for m2n covariance matrix.\n");
            errorFlag=1;
        }
        continue;
    }
      if(!strcmp(tag[i],"M_cut"))
	{
	  if(HOD.pdfs==2 || HOD.pdfs==3){
	    fprintf(stderr,"No value for M_cut given for pdfs= 2/3\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"M_cen_max"))
	{
	  if(HOD.pdfc==7) {
	    fprintf(stderr,"No value for M_cen_max given for pdfc= 7\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"sigma_logM"))
	{
	  if(HOD.pdfc==2){
	    fprintf(stderr,"No value for sigma_logM given for pdfc=2\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"MaxCen"))
	{
	  if(HOD.pdfc==5){
	    fprintf(stderr,"No value for MaxCen given for pdfc=5\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(*tag[i])
	{
	  fprintf(stderr,"Error. I miss a value for tag '%s' in parameter file '%s'.\n",
		  tag[i],fname);
	  errorFlag=1;
	}
    }

  if(PCA==1 && COVAR==1)
    {
      fprintf(stderr,"ERROR: you have both PCA and COVAR set to 1.\n");
      errorFlag=1;
    }

  IDUM_MCMC=IDUM_MCMC_TEMP;
  MCMC = Task.MCMC;

  if(errorFlag)
    endrun("error in input_params ");

  /* Other initialization stuff.
   */
  MSTAR=mstar();

  //printf("NGALS xp %e\n",GALAXY_DENSITY);

  ctemp = HOD.color; 
  if(Task.wp_minimize)
    HOD.color = 0;
  wp.ngal=0;
  set_HOD_params();
  HOD.color = ctemp;
  if(XCORR)set_HOD2_params();

  //Set to run with cosmology correction if redshifts are available
  if(REDSHIFT_MIN >= 0 && REDSHIFT_MAX>REDSHIFT_MIN) GALDENS_CCORR=1;

  //printf("NGALS xx %e %e\n",GALAXY_DENSITY,wp.ngal);

#undef DOUBLE 
#undef STRING 
#undef INT 
#undef MAXTAGS
}
 

