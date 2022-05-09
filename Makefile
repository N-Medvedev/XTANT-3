# The makefile calles for another makefile within Source_files directory
# which compiles and makes an executable called XTANT.x
# This file was written by N.Medvedev 
# in 2021-2022
#----------------------------------------------------- 

# Call makefile within the Source_files directory:
subsystem:
	cd Source_files && $(MAKE)

# Copy created executable into the parent directory:
	scp -r Source_files/XTANT.x XTANT.x

# Delete executable from the Source_files directory:
	rm -r Source_files/XTANT.x


clean:
	cd Source_files &&$(MAKE clean)
	rm -f *.o
	rm -f *.mod
	rm -f XTANT.x
