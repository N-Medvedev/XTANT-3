# The makefile compiles (type make) and makes executables for XTANT-results post-processing
# This file was written by N.Medvedev
# in 2023
#-----------------------------------------------------

default:
	@echo "**************************************************"
	@echo "Compiling XTANT_atomic_data_analysis code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_atomic_data_analysis.f90 -o XTANT_atomic_data_analysis.x
	@echo "Executable: XTANT_atomic_data_analysis.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_autocorrelators code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_autocorrelators.f90 -o XTANT_autocorrelators.x
	@echo "Executable: XTANT_autocorrelators.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_average_diffraction code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_average_diffraction.f90 -o XTANT_average_diffraction.x
	@echo "Executable: XTANT_average_diffraction.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_average_kappa code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_average_kappa.f90 -o XTANT_average_kappa.x
	@echo "Executable: XTANT_average_kappa.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_coupling_parameter code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_coupling_parameter.f90 -o XTANT_coupling_parameter.x
	@echo "Executable: XTANT_coupling_parameter.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_dielectric_function_analysis code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_dielectric_function_analysis.f90 -o XTANT_dielectric_function_analysis.x
	@echo "Executable: XTANT_dielectric_function_analysis.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_fragmentation code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_fragmentation.f90 -o XTANT_fragmentation.x
	@echo "Executable: XTANT_fragmentation.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_entropy code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_entropy.f90 -o XTANT_entropy.x
	@echo "Executable: XTANT_entropy.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_el_distribution_analysis code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_el_distribution_analysis.f90 -o XTANT_el_distribution_analysis.x
	@echo "Executable: XTANT_el_distribution_analysis.x"

	@echo "**************************************************"
	@echo "Compiling XTANT_average_identical code"
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics XTANT_average_identical.f90 -o XTANT_average_identical.x
	@echo "Executable: XTANT_average_identical.x"


clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.obj
	rm -f *.x
