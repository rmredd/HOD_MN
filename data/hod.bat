% Cosmological Parameters
%----------------------------------------------------------------------------

GAMMA		0.2
OMEGA_M		0.27
SIGMA_8		0.82
RHO_CRIT	2.78E11		% h^2 M_sol/Mpc^3
SPECTRAL_INDX	0.95
HUBBLE		0.7
OMEGA_B		0.04
DELTA_CRIT	1.686
ITRANS		5
TF_file		transfunc.WMAP3


% N-body simulation specifics
%----------------------------------------------------------------------------

REDSHIFT	0.05
DELTA_HALO	337.3
BOX_SIZE		0		% Mpc/h
RESOLUTION	0.375		% BOX_SIZE/Npart^(1/3)


% HOD Parameters
%----------------------------------------------------------------------------
% NB! -> If M_min is > 0, then GALAXY_DENSITY is calculated in the code.
%	 If M_min is <=0, then M_min is calculated for the given GALAXY_DENSITY
% 9.442036e+13 1.305428e+00
% 7.762375e+13 1.043098 1.470442e+13 0.123981
% 1.674134e+01 4.899611e+12 7.497712e+13 1.050227e+00 1.785872e+13 1.219268e-01  1.240018

M_min           0
M1              2.340e13
M_max           1.00E+16
M_cut           3.137e12
alpha           0.9500
GALAXY_DENSITY  0.003076

sigma_logM 	0.446


pdfs		12
pdfc		12

VBIAS		1
VBIAS_C		0.4
%CVIR_FAC	0.0664
CVIR_FAC	0.8187
M_cen_lin	3.0e14

EXCLUSION	4


% Free parameters for HOD fitting
% -----------------------------------------------------------------------------
% 1 is free

free[1]		0	% M_min
free[2]		1	% M1
free[3]		1	% alpha
free[4]		1	% M_cut
free[5]		1	% sigma_logM
free[6]		1	% CVIR_FAC
free[7]		0	% HOD.MaxCen
free[10]	1	% HOD.M_cen_lin -- pdfc==12
free[11]		1	% OMEGA_M
free[12]		1	% SIGMA_8
free[13]		0	% VBIAS
free[16]		0	% NSPEC
free[17]		0	% HUBBLE
% scale-dep halo bias parameters
free[18]		0	% HBIAS_C1
free[19]		0	% HBIAS_C2
free[20]		0	% HBIAS_D1
free[21]		0	% HBIAS_D2

% These are the different tasks that the program can perform
% If the All flag is hit, then the other flags are ignored and all the tasks
% are performed.
%-------------------------------------------------------------------------------

All		0
z_space_xi	0
real_space_xi	1
kaiser_xi	0
r_half		0
multipoles	0
wp_minimize	1
m2n_minimize	1	% Determines inclusion of M2N in fitting
COVAR		1
HOD		1
PVD		0
root_filename	b_vp_sc0.2_20.5_fittest
populate_sim	0
HaloFile		/home/tinker/LANL/WMAP/halo.WMAP.384
%HaloFile		/home/tinker/LANL/halo.768

fname_wp	/u/ki/rmredd/data/HOD_MN_dat/bolshoi_tests/bolshoi_vpeak_20.5_sc20/bolshoi_vpeak_20.5_s.20_c0_wp_jtcovar.dat
fname_covar	/u/ki/rmredd/data/HOD_MN_dat/bolshoi_tests/bolshoi_vpeak_20.5_sc20/bolshoi_vpeak_20.5_s.20_c0_covar_jtcovar.dat

% Files and flags for M2N data and minimization
%-------------------------------------------------------------------------------
m2n_filename	  /u/ki/rmredd/data/HOD_MN_dat/bolshoi_tests/bolshoi_vpeak_20.5_sc20/m2n_bolshoi_vpeak_20.5_sc20_rich.dat
mn_fname_covar    /u/ki/rmredd/data/HOD_MN_dat/bolshoi_tests/bolshoi_vpeak_20.5_sc20/m2n_bolshoi_vpeak_20.5_sc20_rich_covar.dat
M2N_type	  1	% 1=richness; 0=mass; 2=mass w/ Ntot, 3=richness w/ Ntot
MN_COVAR	  1 	% Whether to use M/N covariance matrix
M2N		  1	% Prints M2N model, richness or mass as specified

% Files and flags for zspace data and minimization
%-------------------------------------------------------------------------------

zspace_minimize	0
POWELL		0
MCMC		0

OUTPUT		0

