# The makefile compiles (type make) and makes an executable called XTANT.x
# This file was written by N.Medvedev
# in 2021-2024
#-----------------------------------------------------

EXEPATH = .

ifeq ($(db),y)
	EXE = $(EXEPATH)/XTANT_DEBUG.x
else
	EXE = $(EXEPATH)/XTANT.x
endif
LD = ld

# choose the compilers among gfortran, gfortran6.1.0, ifort2011, ifort2013
ifeq ($(c),gf)
	F90 = gfortran
	F90FLAGS = -ffree-line-length-0 -fdec-format-defaults -cpp -fdefault-real-8 -fdefault-double-8 -O2
	LDFLAGS = -lblas -llapack -fopenmp
	# need to link OpenMP,     ^^^^^^^ even for non OpenMP compilation. Can be fixed..
	ifeq (${OMP},y)
		F90FLAGS += -fopenmp
	endif
else ifeq ($(c),gf6)
	F90 = gfortran6.1.0
	F90FLAGS = -ffree-line-length-0 -fdec-format-defaults -cpp -fdefault-real-8 -fdefault-double-8 -O2
	LDFLAGS = -lblas -llapack -fopenmp
	ifeq (${OMP},y)
		F90FLAGS += -fopenmp
	endif
else ifeq ($(c),if13)
	F90 = ifort2013
	
	ifneq (${OMP},) # OMP=no means no omp
		F90FLAGS = -mkl
	else  # there is omp by default
		F90FLAGS = -qopenmp
		F90FLAGS += -mkl=parallel
	endif
	
	ifeq ($(db),y)
		# flags for debugging	
		#F90FLAGS += -g -fbounds-check
		F90FLAGS += -debug all -check all -fpe0 -fp-stack-check -O0 -g -fp-model precise -traceback -gen-interfaces -warn interfaces -fpp
	else
		# Flags for maximum performance
		F90FLAGS += -O5 -fpp -ipo -real-size 64
	endif
else
	F90 = ifort
	# F90 = ifort.exe
	
	ifneq (${OMP},) # OMP=no means no omp
		F90FLAGS = -mkl
	else  # there is omp by default
		F90FLAGS = -qopenmp
		F90FLAGS += -mkl=parallel
	endif
	
	ifeq ($(db),y)
		# flags for debugging	
		#F90FLAGS += -g -fbounds-check
		F90FLAGS += -debug all -check all -fpe0 -fp-stack-check -O0 -g -fp-model precise -traceback -gen-interfaces -warn interfaces -fpp
	else
		# Flags for maximum performance
		F90FLAGS += -O5 -fpp -ipo -real-size 64
	endif
endif

# list the files necessary for the compilation in the order of their dependencies
OBJS = Universal_constants.o Objects.o Variables.o MPI_subroutines.o BS_Cartesian_Gaussians.o BS_Spherical_Gaussians.o Periodic_table.o Algebra_tools.o BS_Basis_sets.o TB_Koster_Slater.o Dealing_with_files.o Gnuplotting.o Little_subroutines.o Dealing_with_EADL.o Dealing_with_DFTB.o Dealing_with_3TB.o Dealing_with_BOP.o Dealing_with_xTB.o Atomic_tools.o Atomic_thermodynamics.o MC_cross_sections.o  Electron_tools.o Nonadiabatic.o Dealing_with_POSCAR.o Dealing_with_mol2.o Dealing_with_CDF.o Read_input_data.o Van_der_Waals.o Coulomb.o ZBL_potential.o Exponential_wall.o Dealing_with_output_files.o Monte_Carlo.o TB_Fu.o TB_Pettifor.o TB_Molteni.o TB_NRL.o TB_DFTB.o TB_3TB.o TB_BOP.o TB_xTB.o TB.o Dealing_with_eXYZ.o Initial_configuration.o  Optical_parameters.o TB_complex.o Transport.o XTANT_MAIN_FILE.o

# print explanations
default:
	@echo "**************************************************"
ifeq ($(c),gf)
	@echo "Compiling XTANT-3 with gfortran"
else ifeq ($(c),gf6)
	@echo "Compiling XTANT-3 with gfortran 6.1.0"
else ifeq ($(c),if13)
	@echo "Compiling XTANT-3 with ifort 2013"
else
	@echo "Compiling XTANT-3 with ifort"
endif
#	@echo "Version 3"

# autoclean before the compilation
#	rm -f *.o
#	rm -f *.mod
	$(MAKE) $(EXE)
# autoclean after the compilation
#	rm *.o
#	rm *.mod

# compile and load
$(EXE):	$(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $(EXE)  $(OBJS)

.SUFFIXES: .f90 .o .mod
.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.x

# cleaning and removing all results - be careful!
# cleanall:
#	rm -f *.o
#	rm -f *.mod
#	rm -f *.x
#	rm -rf OUT*
