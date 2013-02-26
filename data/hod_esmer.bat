% Cosmological Parameters
%----------------------------------------------------------------------------

GAMMA		0.2
OMEGA_M		0.25
SIGMA_8		0.8
RHO_CRIT	2.78E11		% h^2 M_sol/Mpc^3
%RHO_CRIT	2.89E11		% h^2 M_sol/Mpc^3
SPECTRAL_INDX	1.0
HUBBLE		0.7
OMEGA_B		0.04
DELTA_CRIT	1.686
ITRANS		5
TF_file		transfunc.WMAP3

%HBIAS_C1	1.13
%HBIAS_C2	0.46
%HBIAS_D1	1.20
HBIAS_C1	1.17
HBIAS_C2	0.69
HBIAS_D1	1.49
HBIAS_D2	2.09

% N-body simulation specifics
%----------------------------------------------------------------------------

REDSHIFT	0.05
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

M_min           0
M1              1.65e14
M_max           1.00E+16
M_cut           1.76e12
alpha           0.7821
GALAXY_DENSITY  0.0003

sigma_logM 	0.440

pdfs		12
pdfc		2

VBIAS		1
VBIAS_C		0.4
%CVIR_FAC	0.45
CVIR_FAC	0.45
MaxCen		0.0726
%M_cen_lin	3.28e20

EXCLUSION	4


% Free parameters for HOD fitting
% -----------------------------------------------------------------------------
% 1 is free

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
massfunc	1
real_space_xi	1
kaiser_xi	0
matter_xi	0
r_half		0
multipoles	0
wp_minimize	0
m2n_minimize	1	% Determines inclusion of M2N in fitting
COVAR		1
HOD		1
PVD		0
root_filename	esmer_truth
populate_sim	0
HaloFile		/home/tinker/LANL/WMAP/halo.WMAP.384
%HaloFile		/home/tinker/LANL/halo.768

fname_wp	wp_esmer.dat
fname_covar	wp_esmer.covar

% Files and flags for M2N data and minimization
%-------------------------------------------------------------------------------
m2n_filename	  m2n_esmeralda_rich.dat
mn_fname_covar	  m2n_esmeralda_rich_covar.dat
M2N_type	  3	% 1=rich, 0=mass; 2=mass w/ Ntot, 3=rich w/ Ntot
MN_COVAR	  1	% Whether to use M/N covariance matrix
M2N		  1	% Prints M2N model, richness or mass as specified

% Files and flags for zspace data and minimization
%-------------------------------------------------------------------------------

zspace_minimize	0
POWELL		0
MCMC		0

OUTPUT		1

