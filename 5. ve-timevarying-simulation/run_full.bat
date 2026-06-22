@echo off
REM ============================================================================
REM  run_full.bat - launch the COMPLETE time-varying VE(t) simulation pipeline
REM  (full tiered sweep over all methods + all summary tables and figures).
REM
REM  Just double-click this file, or run it from a terminal:  run_full.bat
REM
REM  Optional overrides (set before calling, e.g.  set VE_WORKERS=8 & run_full.bat):
REM    VE_WORKERS  number of parallel workers   (default: CPU cores - 1)
REM    VE_PARALLEL  1/0 turn parallelism on/off (default: 1)
REM    VE_METHODS  e.g. m1,m2,m3  to skip the slow M4/Stan model
REM    VE_FULL_R   replications per full-tier cell      (default 300)
REM    VE_HEAD_R   replications per headline-tier cell  (default 1000)
REM
REM  A full default run is LARGE (~58k replications; M4/Stan dominates) and can
REM  take many hours to days. For a quick check first, try:
REM      set VE_FULL_R=5 & set VE_HEAD_R=5 & run_full.bat
REM ============================================================================

setlocal enabledelayedexpansion

REM --- work from this script's own folder (the project root) ------------------
cd /d "%~dp0"

REM --- locate Rscript.exe -----------------------------------------------------
set "RSCRIPT="

REM 1) on PATH?
for /f "delims=" %%i in ('where Rscript.exe 2^>nul') do (
    if not defined RSCRIPT set "RSCRIPT=%%i"
)

REM 2) newest install under Program Files (64-bit then 32-bit)
if not defined RSCRIPT (
    for /f "delims=" %%d in ('dir /b /ad /o-n "%ProgramFiles%\R\R-*" 2^>nul') do (
        if not defined RSCRIPT if exist "%ProgramFiles%\R\%%d\bin\Rscript.exe" (
            set "RSCRIPT=%ProgramFiles%\R\%%d\bin\Rscript.exe"
        )
    )
)
if not defined RSCRIPT (
    for /f "delims=" %%d in ('dir /b /ad /o-n "%ProgramFiles(x86)%\R\R-*" 2^>nul') do (
        if not defined RSCRIPT if exist "%ProgramFiles(x86)%\R\%%d\bin\Rscript.exe" (
            set "RSCRIPT=%ProgramFiles(x86)%\R\%%d\bin\Rscript.exe"
        )
    )
)

if not defined RSCRIPT (
    echo [ERROR] Could not find Rscript.exe.
    echo         Install R from https://cran.r-project.org/ or add Rscript to PATH,
    echo         then run this file again.
    pause
    exit /b 1
)

echo Using R: "%RSCRIPT%"
echo Starting full simulation pipeline...
echo.

REM --- run the driver, tee a log file -----------------------------------------
"%RSCRIPT%" --vanilla run_full.R
set "RC=%ERRORLEVEL%"

echo.
if "%RC%"=="0" (
    echo ============================================================
    echo  DONE. Outputs in:  outputs\sweep  outputs\summary  outputs\figures
    echo ============================================================
) else (
    echo [ERROR] run_full.R exited with code %RC%. See messages above.
)

pause
exit /b %RC%
