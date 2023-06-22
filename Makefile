# The makefile calles for another makefile within Source_files directory
# which compiles and makes an executable called XTANT.x
# This file was written by N.Medvedev 
# in 2021-2023
#----------------------------------------------------- 

# To pass variables into the next make-file:
export

# Call makefile within the Source_files directory:
subsystem:
	cd Source_files && $(MAKE)

# Copy created executable into the parent directory:
	scp -r Source_files/XTANT.x XTANT.x

# Delete executable from the Source_files directory:
	rm -r Source_files/XTANT.x

# Compile post-processing subroutines:
	cd \!XTANT_ANALYSIS_SUBROUTINES && $(MAKE)


clean:
	rm -f Source_files/*.o
	rm -f Source_files/*.mod
	rm -f Source_files/*.obj
	rm -f Source_files/*.x

	rm -f \!XTANT_ANALYSIS_SUBROUTINES/*.o
	rm -f \!XTANT_ANALYSIS_SUBROUTINES/*.mod
	rm -f \!XTANT_ANALYSIS_SUBROUTINES/*.obj
	rm -f \!XTANT_ANALYSIS_SUBROUTINES/*.x
