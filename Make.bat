:: This file was written by N.Medvedev
:: in 2018-2022
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

   
   SET "Starline=************************************************************************************"
   echo %Starline%
   echo Started compilation: %date% %time%

:: List of all program files to be compiled
   SET "List_of_files=Universal_constants.f90 Objects.f90 Variables.f90 BS_Cartesian_Gaussians.f90 BS_Spherical_Gaussians.f90 Gnuplotting.f90 Periodic_table.f90 Algebra_tools.f90 BS_Basis_sets.f90  TB_Koster_Slater.f90 Dealing_with_files.f90 Little_subroutines.f90 Dealing_with_EADL.f90 Dealing_with_DFTB.f90 Dealing_with_3TB.f90 Dealing_with_BOP.f90 Dealing_with_xTB.f90 Atomic_tools.f90 Electron_tools.f90 Nonadiabatic.f90 Read_input_data.f90 Van_der_Waals.f90 Coulomb.f90 ZBL_potential.f90 Exponential_wall.f90 MC_cross_sections.f90 Dealing_with_output_files.f90 Monte_Carlo.f90 TB_Fu.f90 TB_Pettifor.f90 TB_Molteni.f90 TB_NRL.f90 TB_DFTB.f90 TB_3TB.f90 TB_BOP.f90 TB_xTB.f90 TB.f90 Initial_configuration.f90 Optical_parameters.f90 Transport.f90 XTANT_MAIN_FILE.f90"

:: List compiler options and the name of the executable:
   IF /I %arg1%==DEBUG (
      echo %Starline%
      echo Compiling with DEBUG option, no OpenMP or optimizations are included
      echo Started at: %date% %time%
      echo %Starline%

      :: List compiler options 
      SET "Compile_options=/F9999999999 /QxHost /QaxAVX  /fpp /Qmkl=parallel /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec /standard-semantics"

      :: Set name of the executable:
      SET "Name_of_exe=XTANT_DEBUG.exe"
   ) ELSE (
      IF /I %arg1%==DEBUGOMP (
         echo %Starline%
         echo Compiling with DEBUGOMP option, OpenMP but no optimizations are included
         echo Started at: %date% %time%
         echo %Starline%

         :: List compiler options 
         SET "Compile_options=/F9999999999 /QxHost /QaxAVX  /fpp /Qopenmp /D OMP_inside /Qmkl=parallel /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec /standard-semantics"

         :: Set name of the executable:
         SET "Name_of_exe=XTANT_DEBUG_OMP.exe"
     ) ELSE (
        IF /I %arg1%==FAST (
            echo %Starline%
            echo Compiling with FAST option, OpenMP, no optimizations, no debug
            echo Started at: %date% %time%
            echo %Starline%

            :: List compiler options
            SET "Compile_options=/F9999999999 /fpp /Qopenmp /D OMP_inside /Qmkl=parallel /real-size:64 /Od /fpe:0 /fp:precise /Qvec /standard-semantics"

            :: Set name of the executable:
            SET "Name_of_exe=XTANT_OMP.exe"
        ) ELSE (
            echo %Starline%
            echo Compiling for release, OpenMP and optimizations are included
            echo Started at: %date% %time%
            echo %Starline%

            :: List compiler options
            SET "Compile_options= /Qopenmp /D OMP_inside /Qmkl=parallel /O3 /fpp /Qvec /Qipo /real-size:64 /standard-semantics /F9999999999 "

            :: Set name of the executable:
            SET "Name_of_exe=XTANT.exe"

            del *.pdb
         )
       )
    )

:: Compile modules
   ifort.exe -c %Compile_options% %List_of_files%
   
   echo %Starline%
   echo Assembling the files into executable: %Name_of_exe%
   echo Started at: %date% %time%
   echo %Starline%

:: Assemble the code from all created obj-files
   ifort.exe %Compile_options% *.obj /exe:%Name_of_exe%

   echo %Starline%
::   echo Completed: %date% %time%
   echo The program %Name_of_exe% was created at %date% %time%
   echo %Starline%
   

:: Remove files that are no longer needed
del *.obj *.mod

:: Go back into the parent directory from the source files:
cd ..\

:: Copy exe file from the source directory into the parent directory:
xcopy source_files\%Name_of_exe% %Name_of_exe%* /Y /Q

:: Delete the exe file from the soure directory:
del source_files\%Name_of_exe%
