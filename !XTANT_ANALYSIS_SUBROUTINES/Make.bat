:: This file was written by N.Medvedev
:: in 2018-2022
:: -----------------------------------------------------
:: Brief notes on cmd programming :  https://ss64.com/nt/if.html
:: -----------------------------------------------------

@echo off
setlocal EnableDelayedExpansion

:: read argument from user
SET arg1=%1
    
:: in case of empty argument, assume compiling of all supporting subroutines
IF "%1"=="" (
   SET arg1=ALL
)

::del *.obj

SET "Starline=************************************************************************************"
echo %Starline%
echo Started compilation: %date% %time%

:: Program files to be compiled
SET "List_of_files= XTANT_atomic_data_analysis.f90"

:: List compiler options and the name of the executable:
IF /I %arg1%==DEBUG (
   echo %Starline%
   echo Compiling with DEBUG option, no OpenMP or optimizations are included
   echo Started at: %date% %time%
   echo %Starline%

   :: List compiler options
   SET "Compile_options=/F9999999999 /QxHost /QaxAVX  /fpp /Qmkl=parallel /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /standard-semantics"

   :: Set name of the executable:
   SET "Name_of_exe=XTANT_atomic_data_analysis.exe"
) ELSE (
   echo %Starline%
   echo Compiling for release, OpenMP and optimizations are included
   echo Started at: %date% %time%
   echo %Starline%

   :: List compiler options
   SET "Compile_options= /Qopenmp /D OMP_inside /Qmkl=parallel /O3 /fpp /Qipo /real-size:64 /standard-semantics /F9999999999 "

   :: Set name of the executable:
   SET "Name_of_exe=XTANT_atomic_data_analysis.exe"
)

:: Assemble the code from all created obj-files
ifort.exe %Compile_options% %List_of_files% /exe:%Name_of_exe%

echo %Starline%
::   echo Completed: %date% %time%
echo The program %Name_of_exe% was created at %date% %time%
echo %Starline%


:: Remove files that are no longer needed
del XTANT_atomic_data_analysis.obj


:: *********************************************************

echo Started compilation: %date% %time%

:: Program files to be compiled
SET "List_of_files= XTANT_autocorrelators.f90"

:: List compiler options and the name of the executable:
IF /I %arg1%==DEBUG (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_autocorrelators.exe"
) ELSE (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_autocorrelators.exe"
)

:: Assemble the code from all created obj-files
ifort.exe %Compile_options% %List_of_files% /exe:%Name_of_exe%

echo %Starline%
::   echo Completed: %date% %time%
echo The program %Name_of_exe% was created at %date% %time%
echo %Starline%


:: Remove files that are no longer needed
del XTANT_autocorrelators.obj


:: *********************************************************

echo Started compilation: %date% %time%

:: Program files to be compiled
SET "List_of_files= XTANT_coupling_parameter.f90"

:: List compiler options and the name of the executable:
IF /I %arg1%==DEBUG (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_coupling_parameter.exe"
) ELSE (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_coupling_parameter.exe"
)

:: Assemble the code from all created obj-files
ifort.exe %Compile_options% %List_of_files% /exe:%Name_of_exe%

echo %Starline%
::   echo Completed: %date% %time%
echo The program %Name_of_exe% was created at %date% %time%
echo %Starline%


:: Remove files that are no longer needed
del XTANT_coupling_parameter.obj


:: *********************************************************

echo Started compilation: %date% %time%

:: Program files to be compiled
SET "List_of_files= XTANT_dielectric_function_analysis.f90"

:: List compiler options and the name of the executable:
IF /I %arg1%==DEBUG (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_dielectric_function_analysis.exe"
) ELSE (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_dielectric_function_analysis.exe"

   del *.pdb
)

:: Assemble the code from all created obj-files
ifort.exe %Compile_options% %List_of_files% /exe:%Name_of_exe%

echo %Starline%
::   echo Completed: %date% %time%
echo The program %Name_of_exe% was created at %date% %time%
echo %Starline%

:: Remove files that are no longer needed
del XTANT_dielectric_function_analysis.obj

:: *********************************************************

echo Started compilation: %date% %time%

:: Program files to be compiled
SET "List_of_files= XTANT_fragmentation.f90"

:: List compiler options and the name of the executable:
IF /I %arg1%==DEBUG (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_fragmentation.exe"
) ELSE (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_fragmentation.exe"
)

:: Compiling:
ifort.exe -c %Compile_options% %List_of_files%

:: Assemble the code from all created obj-files
ifort.exe %Compile_options% *.obj /exe:%Name_of_exe%

echo %Starline%
::   echo Completed: %date% %time%
echo The program %Name_of_exe% was created at %date% %time%
echo %Starline%

:: Remove files that are no longer needed
del XTANT_fragmentation.obj

:: *********************************************************
echo Started compilation: %date% %time%

:: Program files to be compiled
SET "List_of_files= XTANT_entropy.f90"

:: List compiler options and the name of the executable:
IF /I %arg1%==DEBUG (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_entropy.exe"
) ELSE (
   :: Set name of the executable:
   SET "Name_of_exe=XTANT_entropy.exe"
)

:: Compiling:
ifort.exe -c %Compile_options% %List_of_files%

:: Assemble the code from all created obj-files
ifort.exe %Compile_options% *.obj /exe:%Name_of_exe%

echo %Starline%
::   echo Completed: %date% %time%
echo The program %Name_of_exe% was created at %date% %time%
echo %Starline%

:: Remove files that are no longer needed
del XTANT_entropy.obj


