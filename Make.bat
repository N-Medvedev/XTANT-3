:: This file was written by N.Medvedev
:: in 2018-2026
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

:: Read the ordered source list shared with the Unix Makefile
   SET "List_of_files="
   FOR /F "usebackq tokens=*" %%A IN ("sources.txt") DO (
      SET "List_of_files=!List_of_files! %%A"
   )

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
      SET "Compile_options=/F9999999999 /QxHost /QaxAVX2 /fpp /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /check:uninit /fpe:0 /fp:precise /standard-semantics /Qfp-stack-check"
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
      SET "Compile_options=/F9999999999 /QxHost /QaxAVX2 /fpp /real-size:64 /debug:all /debug:full /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /check:uninit /fpe:0 /fp:precise /standard-semantics /Qfp-stack-check"
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
      :: SET "Compile_options=/fpp /real-size:64 /O1 /fpe:0 /fp:fast /heap-arrays /standard-semantics /Qopenmp /Qmkl=parallel /F:1073741824"
      SET "Linking_options=/Qopenmp /Qmkl=parallel"

      :: Set name of the executable:
      SET "Name_of_exe=XTANT_OMP.exe"

      del *.pdb
   )

   IF /I %arg1%==HEAP (
      echo %Starline%
      echo Compiling with FAST option, OpenMP, no optimizations, no debug
      echo Started at: %date% %time%
      echo %Starline%

      :: List compiler options
      :: SET "Compile_options=/F9999999999 /fpp /real-size:64 /O1 /fpe:0 /fp:fast /Qipo /Qopt-report /standard-semantics"
      SET "Compile_options=/fpp /real-size:64 /O1 /fpe:0 /fp:fast /heap-arrays /standard-semantics /Qopenmp /Qmkl=parallel /F:1073741824"
      SET "Linking_options=/Qopenmp /Qmkl=parallel"

      :: Set name of the executable:
      SET "Name_of_exe=XTANT_HEAPED.exe"

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
         :: Execute the compilation of the module:
         CALL !Output!
         IF ERRORLEVEL 1 (
            echo Error compiling %%A! See Captured Output above. Exiting...
            EXIT /B 1
         )
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
      echo * Executable : !Output!
      :: Combine modules into the executable:
      CALL !Output!
      IF ERRORLEVEL 1 (
         echo Error linking %Name_of_exe%. Exiting...
         EXIT /B 1
      )

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
