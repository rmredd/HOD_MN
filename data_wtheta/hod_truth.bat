% Cosmological Parameters
%----------------------------------------------------------------------------

GAMMA		0.2
OMEGA_M		0.286
SIGMA_8		0.82
RHO_CRIT	2.78E11		% h^2 M_sol/Mpc^3
SPECTRAL_INDX	0.96
HUBBLE		0.7
OMEGA_B		0.047
DELTA_CRIT	1.686
ITRANS		5
TF_file		transfunc.WMAP3


% N-body simulation specifics
%----------------------------------------------------------------------------

REDSHIFT	0.15
DELTA_HALO	392.7
BOX_SIZE		0		% Mpc/h
RESOLUTION	0.375		% BOX_SIZE/Npart^(1/3)


% HOD Parameters
%----------------------------------------------------------------------------
% NB! -> If M_min is > 0, then GALAXY_DENSITY is calculated in the code.
%	 If M_min is <=0, then M_min is calculated for the given GALAXY_DENSITY
% 9.442036e+13 1.305428e+00
% 7.762375e+13 1.043098 1.470442e+13 0.123981
% 1.674134e+01 4.899611e+12 7.497712e+13 1.050227e+00 1.785872e+13 1.219268e-01  1.240018

M_min           1.2868e12
M1              1.6710e13
M_max           1.00E+16
M_cut           3.9253e12
alpha           0.92208
GALAXY_DENSITY  0.0038588
GALDENS_ERR	0.000036

sigma_logM 	0.3970


pdfs		12
pdfc		12

VBIAS		1
VBIAS_C		0.4
%CVIR_FAC	0.0664
CVIR_FAC	0.1
M_cen_lin	1.0175e14
MaxCen		0.05

EXCLUSION	4


% Free parameters for HOD fitting
% -----------------------------------------------------------------------------
% 1 is free

free[0]		1	% GALAXY_DENSITY -- matches against data
free[1]		0	% M_min
free[2]		0	% M1
free[3]		0	% alpha
free[4]		0	% M_cut
free[5]		0	% sigma_logM
free[6]		0	% CVIR_FAC
free[7]		0	% HOD.MaxCen
free[10]	0	% HOD.M_cen_lin -- pdfc==12
free[11]		0	% OMEGA_M
free[12]		0	% SIGMA_8
free[13]		0	% VBIAS
free[16]		0	% NSPEC
free[17]		0	% HUBBLE

% These are the different tasks that the program can perform
% If the All flag is hit, then the other flags are ignored and all the tasks
% are performed.
%-------------------------------------------------------------------------------

FIT_WTHETA	0
All		0
z_space_xi	0
real_space_xi	1
kaiser_xi	0
angular_xi	1
r_half		0
multipoles	0
wp_minimize	1
m2n_minimize	1	% Determines inclusion of M2N in fitting
COVAR		0
HOD		1
PVD		0
root_filename	fit_truth_z_0.1_0.2
populate_sim	0
HaloFile		/home/tinker/LANL/WMAP/halo.WMAP.384

fname_wp	wp_z_0.1_0.2.dat
fname_covar	wp_covar_line_z_0.1_0.2.dat

fname_nz	nz_z_0.1_0.2_ztrue.dat % File with necessary n(z) data
%fname_nz	nz_test_z_0.1_0.2.dat % File with necessary n(z) data

% Files and flags for M2N data and minimization
%-------------------------------------------------------------------------------
m2n_filename	  mn_ch_m20.5_z_0.1_0.2.dat
mn_fname_covar	  mn_ch_covar_line_m20.5_z_0.1_0.2.dat
M2N_type	  2	% 1=rich, 0=mass; 2=mass w/ Ntot, 3=rich w/ Ntot
MN_COVAR	  0	% Whether to use M/N covariance matrix
M2N		  1	% Prints M2N model, richness or mass as specified

% Files and flags for zspace data and minimization
%-------------------------------------------------------------------------------

zspace_minimize	0
POWELL		0
MCMC		0

OUTPUT		0

