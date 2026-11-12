@echo off
setlocal EnableDelayedExpansion

set "arg1=%~1"
if "%arg1%"=="" set "arg1=ALL"

set "Starline=************************************************************************************"
echo %Starline%
echo Started compilation: %date% %time%

if /I "%arg1%"=="DEBUG" (
   set "Compile_options=/F9999999999 /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /standard-semantics"
) else (
   set "Compile_options=/O3 /fpp /Qipo /real-size:64 /standard-semantics /F9999999999"
)

for /F "usebackq tokens=*" %%A in ("sources.txt") do (
   echo %Starline%
   echo Compiling %%A
   ifx.exe %Compile_options% "%%A" /exe:%%~nA.exe
   if errorlevel 1 (
      echo Error compiling %%A. Exiting...
      exit /b 1
   )
   echo Executable: %%~nA.exe
)

echo %Starline%
echo Finished compilation: %date% %time%
exit /b 0
