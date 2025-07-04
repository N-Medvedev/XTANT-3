:: This file was written by N.Medvedev
:: in 2018-2024
:: -----------------------------------------------------
:: Brief notes on cmd programming :  https://ss64.com/nt/if.html
:: -----------------------------------------------------

@echo off
setlocal EnableDelayedExpansion

:: Go into the directory with source files:
cd Source_files

:: read argument from user
   SET arg1=%1

:: in case of empty argument, assume no debug
   IF "%1"=="" (
      SET arg1=NODEBUG
   )

:: shorthand expressions for debug options (convert them into DEBUGOMP option):
   IF /I %arg1%==DBG (
      SET arg1=DEBUGOMP
   )
   IF /I %arg1%==DB (
      SET arg1=DEBUGOMP
   )
:: shorthand expressions for no-optimization and no debug options (fast compiling):
   IF /I %arg1%==SLOW (
      SET arg1=FAST
   )


   SET "Starline=************************************************************************************"
   echo %Starline%
   echo Started compilation: %date% %time%

   :: Set the default compiler name:
   SET "Compiler=ifx"
   SET "Linking_options="

:: List of all program files to be compiled
   SET "List_of_files=Universal_constants.f90 Objects.f90 Variables.f90 MPI_subroutines.f90 BS_Cartesian_Gaussians.f90 BS_Spherical_Gaussians.f90 Periodic_table.f90 Algebra_tools.f90 BS_Basis_sets.f90 TB_Koster_Slater.f90 Dealing_with_files.f90 Gnuplotting.f90 Little_subroutines.f90 Dealing_with_EADL.f90 Dealing_with_DFTB.f90 Dealing_with_3TB.f90 Dealing_with_BOP.f90 Dealing_with_xTB.f90 Atomic_tools.f90 Atomic_thermodynamics.f90 MC_cross_sections.f90 Electron_tools.f90 Nonadiabatic.f90 Dealing_with_POSCAR.f90 Dealing_with_mol2.f90 Dealing_with_CDF.f90 ZBL_potential.f90 Read_input_data.f90 Van_der_Waals.f90 Coulomb.f90 Exponential_wall.f90  Dealing_with_output_files.f90 Monte_Carlo.f90 TB_Fu.f90 TB_Pettifor.f90 TB_Molteni.f90 TB_NRL.f90 TB_DFTB.f90 TB_3TB.f90 TB_BOP.f90 TB_xTB.f90 TB.f90 Dealing_with_eXYZ.f90 Initial_configuration.f90 Optical_parameters.f90 TB_complex.f90 Transport.f90 XTANT_MAIN_FILE.f90"

   :: List compiler options and the name of the executable:
   IF /I %arg1%==clean (
      echo %Starline%
      echo Cleaning up
      echo Started at: %date% %time%
      echo %Starline%

	  del *.mod
      del *.obj
      del *.optrpt
      del *.yaml
      del *.pdb
   )


   IF /I %arg1%==DEBUG (
      echo %Starline%
      echo Compiling with DEBUG option, no OpenMP or optimizations are included
      echo Started at: %date% %time%
      echo %Starline%

      :: List compiler options
      SET "Compile_options=/F9999999999 /QxHost /QaxAVX2 /fpp /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /fp:precise /standard-semantics /Qfp-stack-check"
      SET "Linking_options=/Qmkl=parallel"

      :: Set name of the executable:
      SET "Name_of_exe=XTANT_DEBUG.exe"
   )

   IF /I %arg1%==DEBUGOMP (
      echo %Starline%
      echo Compiling with DEBUGOMP option, OpenMP with /Qipo optimization are included
      echo Started at: %date% %time%
      echo %Starline%

      :: List compiler options
      ::SET "Compile_options=/F9999999999 /QxHost /QaxAVX2 /fpp /Qmkl=parallel /Qopenmp /real-size:64 /debug:all /O1 /Qipo /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /fp:precise /standard-semantics"
      SET "Compile_options=/F9999999999 /QxHost /QaxAVX2 /fpp /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /fp:precise /standard-semantics /Qfp-stack-check"
      SET "Linking_options=/Qopenmp /Qmkl=parallel"

      :: Set name of the executable:
      SET "Name_of_exe=XTANT_DEBUG_OMP.exe"
   )

   IF /I %arg1%==FAST (
      echo %Starline%
      echo Compiling with FAST option, OpenMP, no optimizations, no debug
      echo Started at: %date% %time%
      echo %Starline%

      :: List compiler options
      SET "Compile_options=/F9999999999 /fpp /real-size:64 /O1 /fpe:0 /fp:fast /Qipo /Qopt-report /standard-semantics"
      SET "Linking_options=/Qopenmp /Qmkl=parallel"

      :: Set name of the executable:
      SET "Name_of_exe=XTANT_OMP.exe"

      del *.pdb
   )

   IF /I %arg1%==NODEBUG (
      echo %Starline%
      echo Compiling for release, OpenMP and optimizations are included
      echo Started at: %date% %time%
      echo %Starline%

      :: List compiler options
      SET "Compile_options=/F9999999999 /fpp /real-size:64 /O3 /Qipo /standard-semantics /assume:nofpe_summary"
      SET "Linking_options=/Qopenmp /Qmkl=parallel"

      :: Set name of the executable:
      SET "Name_of_exe=XTANT.exe"
      del *.pdb
   )

   IF /I %arg1%==MPI (
      echo %Starline%
      echo Compiling with MPI parallelization
      echo Started at: %date% %time%
      echo %Starline%

      :: List compiler options (for release):
      SET "Compile_options=/F9999999999 /fpp /D MPI_USED /Qmkl=cluster /real-size:64 /O3 /Qip /fp:precise /standard-semantics /assume:nofpe_summary /gen-interfaces"'
      SET Linking_options= /NODEFAULTLIB:"mkl_intel_ilp64_dll.lib"

      :: List compiler options (for debug):
      ::SET Compile_options=/F9999999999 /QxHost /QaxAVX2 /fpp /D MPI_USED /Qmkl=cluster /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /fp:precise /standard-semantics /Qfp-stack-check

      :: Set name of the executable:
      SET "Name_of_exe=XTANT_MPI.exe"
      del *.pdb

      :: Set the compiler name:
      SET "Compiler=mpiifx"
   )

   :: compile modules one by one:
   IF /I not %arg1%==clean (
      FOR %%A IN (%List_of_files%) DO (
         :: Construct the command line for compilation of the current module:
         SET "Output=%Compiler% -c %Compile_options% %%A !Linking_options! 2>&1"
         echo * Compilation : !Output!
         IF ERRORLEVEL 1 (
            echo Error compiling %%A! See Captured Output above. Exiting...
            EXIT /B 1
         )
         :: Execute the compilation of the module:
         CALL !Output!
      )
      echo %Starline%

      echo %Starline%
      echo Assembling the files into executable: %Name_of_exe%
      echo Started at: %date% %time%
      echo %Starline%

:: Assemble the code from all created obj-files
::   ifx.exe %Compile_options% *.obj /exe:%Name_of_exe%
:: Construct the command line for creation fo executable:
      SET "Output=%Compiler% %Compile_options% *.obj /exe:%Name_of_exe% !Linking_options! 2>&1"
      IF ERRORLEVEL 1 (
         echo Error compiling %%A! See Captured Output above. Exiting...
         EXIT /B 1
      )
      echo * Executable : !Output!
      :: Combine modules into the executable:
      CALL !Output!

      echo %Starline%
::   echo Completed: %date% %time%
      echo The program %Name_of_exe% was created at %date% %time%
      echo %Starline%
      echo INSTRUCTIONS:
      IF /I %arg1%==MPI (
         echo Run the program as: mpiexec -np [n] ./%Name_of_exe%
         echo where [n] = number of processors to be used
      ) ELSE (
         echo Run the program as: %Name_of_exe%
      )
      echo %Starline%
   )

:: Go back into the parent directory from the source files:
cd ..\

:: Copy exe file from the source directory into the parent directory:
IF /I not %arg1%==clean (
   xcopy source_files\%Name_of_exe% %Name_of_exe%* /Y /Q

   :: Delete the exe file from the soure directory:
   del source_files\%Name_of_exe%
)
