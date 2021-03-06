
hd = $(HOME)/lib
LIB = -lm 

CC = gcc
CFLAGS = -O2
EXEC = HOD.x

parallel = 0
errors = 0

ifeq ($(parallel),1)
CC = mpicc
LIB = -lm -lmpich
CFLAGS = -O2 -DPARALLEL
EXEC = HOD.px
endif

ifeq ($(errors),1)
EXEC = HOD.x_err
endif

OBJ_DIR = ../obj

OBJ_HOD = $(addprefix $(OBJ_DIR)/, header.o main.o utility.o sigmac.o \
	transfnc.o transfunc_file.o \
	nonlinear_power_spectrum.o least_squares.o test.o \
	hod_functions.o \
	xi_matter.o one_halo_rspace.o \
	input_params.o dFdx.o mstar.o \
	halo_concentration.o growthfactor.o \
	halo_mass_conversion.o two_halo_rspace.o nfw_transform.o \
	pair_density.o tasks.o wp_minimization.o \
	output_params.o \
	one_halo_zspace.o galaxy_prob_vz.o \
	spherical_infall_velocity.o two_halo_zspace.o \
	xi_multipoles.o chi2_zspace.o integrated_bin.o \
	small_scale_measure.o \
	zspace_minimization.o linlin_bins.o xi2d_interp.o \
	kaiser_distortions.o \
	tf_eisenstein_hu.o density_dependence.o velocity_bias.o \
	vpf.o \
	neg_binomial.o populate_simulation.o color_samples.o \
	pvd.o pvz_temp.o \
	dark_matter_statistics.o hod_functions2.o ftread.o \
	ftwrite.o \
	velocity_moments.o cobenorm.o \
	covar_test.o \
	ml_ratios.o uber_comparison.o \
	nbody_xi.o meshlink2.o nbrsfind2.o i3tensor_2.o \
	chi2_grid.o jeans.o \
	fit_scale_bias.o aspen_breakout.o mcmc_color.o \
	populate_simulation_clf.o \
	mass2number.o m2n_mcmc_combined.o \
	wtheta.o drg_model.o mcmc_wtheta.o )


OBJ_ERR = $(addprefix $(OBJ_DIR)/, mcmc_with_errors.o \
	halo_mass_function_error.o \
	halo_bias_error.o)
OBJ_STD = $(addprefix $(OBJ_DIR)/, mcmc.o halo_mass_function.o halo_bias.o )

OBJ_ML  = $(addprefix $(OBJ_DIR)/, ml_minimization.o)

OBJ_NR = $(addprefix $(OBJ_DIR)/, nrutil.o qromo.o midpnt.o midinf.o \
	polint.o splint.o spline.o \
	zbrent.o qtrap.o trapzd.o cisi.o \
	complex.o amoeba.o amotry.o \
	gaussj.o powell.o linmin.o f1dim.o \
	mnbrak.o brent.o gasdev.o \
	ran1.o jacobi.o splin2.o splie2.o \
	ran2.o sort2.o )

ifeq ($(errors),1)
OBJS01 = $(OBJ_HOD) $(OBJ_NR) $(OBJ_ERR)
else
OBJS01 = $(OBJ_HOD) $(OBJ_NR) $(OBJ_STD) $(OBJ_ML)
endif


$(EXEC): $(OBJS01)
	$(CC) -o $@ $(OBJS01) $(LIB)
	mv *.x ..	
$(OBJ_DIR)/%.o: %.c header.h nrutil.h complex.h
	$(CC) -c -o $@ $< $(LIB)

OBJS11=	chi2_dist.o qromo.o midinf.o nrutil.o \
	polint.o midpnt.o
chi2_dist:	$(OBJS11)
	$(CC) -o $@ $(OBJS11) $(LIB)


clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f *.o